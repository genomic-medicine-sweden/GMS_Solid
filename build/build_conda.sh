#!/usr/bin/env bash
# Build a self-contained bundle (pipeline + conda env + modules + configs + optional
# containers and reference files) for installation on a system without internet access.
#
# The build is resumable: every step writes a marker file into a state directory once it
# has completed successfully, and is skipped on a rerun. Use --force/--force-step to redo
# work, --status to see what is done, and --list-steps for the step names.
#
# Required environment variables:
#   PIPELINE_NAME          e.g. Twist_Solid
#   TAG_OR_BRANCH          pipeline tag or branch, e.g. v0.24.0
#   PIPELINE_GITHUB_REPO   clone url of the pipeline
#   PYTHON_VERSION         python version for the conda env, e.g. 3.9
#   CONFIG_NAME            e.g. GMS560_config
#   CONFIG_VERSION         config tag or branch, e.g. v1.18.0
#   CONFIG_GITHUB_REPO     clone url of the config repo
#
# Optional environment variables:
#   APPT_CACHE_STATUS         build | update | none   (default: none)
#                             build  = download containers and point the config at the cache
#                             update = only point the config at an existing remote cache
#   PATH_TO_apptainer_cache   required when APPT_CACHE_STATUS is build or update
#   GIT_CLONE_RETRIES         attempts per git clone (default: 3)
#   EXECUTOR_PLUGINS          comma-separated Snakemake executor plugins to add to the env:
#                             drmaa, slurm, sge, or "none". A plugin may carry a pip version
#                             spec, e.g. slurm==2.8.0. Needs PYTHON_VERSION >= 3.11. When unset,
#                             the script asks; without a terminal it keeps the previous choice
#                             (none on a first build).
#
# Any positional arguments are reference/design config files passed to
# `hydra-genetics references download`.

set -euo pipefail

# Steps in execution order. Forcing a step also invalidates every step after it.
STEPS=(
    clone_pipeline
    create_env
    stage_pipeline
    install_requirements
    install_executors
    pack_env
    clone_wrappers
    clone_modules
    clone_config
    update_config_paths
    containers
    pack_bundle
    references
    cleanup
)

HYDRA_MODULES=(
    alignment
    annotation
    biomarker
    cnv_sv
    compression
    filtering
    fusions
    misc
    mitochondrial
    parabricks
    prealignment
    qc
    reports
    snv_indels
    references
)

# Executor plugins EXECUTOR_PLUGINS may name; each maps to snakemake-executor-plugin-<name>.
KNOWN_EXECUTORS=(
    drmaa
    sge
    slurm
)

# Config files in the config repo that have ${TAG_OR_BRANCH} substituted into them.
CONFIG_SUBST_FILES=(
    profiles/Miarka/production/config.yaml
    profiles/Miarka/production_novaseqX/config.yaml
    profiles/Miarka/production_ctDNA/config.yaml
    profiles/Miarka/production_prio/config.yaml
    profiles/Miarka/production_prio_novaseqX/config.yaml
    config/Miarka/config_production_pipeline_miarka.yaml
)

FORCE_ALL=false
DRY_RUN=false
KEEP_WORKDIR=false
STATE_DIR=""
FORCED_STEPS=()
ONLY_STEPS=()
SKIP_STEPS=()
REFERENCE_CONFIGS=()

usage() {
    cat <<'USAGE'
Usage: bash build_conda.sh [options] [reference_config.yaml ...]

Options:
  -f, --force              Ignore all markers and rebuild everything from scratch.
      --force-step STEP    Redo STEP (and every step after it). Repeatable, comma-separated ok.
      --only STEP          Run only STEP, ignoring its marker. Repeatable, comma-separated ok.
      --skip STEP          Never run STEP. Repeatable, comma-separated ok.
      --state-dir DIR      Where markers are kept (default: <PIPELINE_NAME>/.build_state).
      --keep-workdir       Do not delete the conda env and staging directory at the end.
  -n, --dry-run            Print what would run, change nothing.
      --status             Print the state of every step and exit.
      --list-steps        Print the step names in execution order and exit.
  -h, --help               Show this help.

Reference config paths may be absolute, or relative to the pipeline clone
(<PIPELINE_NAME>/), the cloned config repo, the bundle root or the directory the script
was started from -- searched in that order, first match wins.

Steps are resumable. A completed step writes <state-dir>/<step>.done; the slow steps
(module clones, container downloads, reference downloads) additionally mark each repo,
container set and reference config separately, so an interrupted run resumes mid-step.
USAGE
}

log() { printf '[%s] %s\n' "$(date -u +%H:%M:%S)" "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

split_list() {
    # Split a comma-separated argument into words.
    printf '%s' "$1" | tr ',' ' '
}

step_exists() {
    local candidate=$1 step
    for step in "${STEPS[@]}"; do
        [[ $step == "$candidate" ]] && return 0
    done
    return 1
}

in_list() {
    local needle=$1 item
    shift
    for item in "$@"; do
        [[ $item == "$needle" ]] && return 0
    done
    return 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -f | --force)
            FORCE_ALL=true
            shift
            ;;
        --force-step)
            [[ $# -ge 2 ]] || die "--force-step requires a step name"
            for step in $(split_list "$2"); do FORCED_STEPS+=("$step"); done
            shift 2
            ;;
        --only)
            [[ $# -ge 2 ]] || die "--only requires a step name"
            for step in $(split_list "$2"); do ONLY_STEPS+=("$step"); done
            shift 2
            ;;
        --skip)
            [[ $# -ge 2 ]] || die "--skip requires a step name"
            for step in $(split_list "$2"); do SKIP_STEPS+=("$step"); done
            shift 2
            ;;
        --state-dir)
            [[ $# -ge 2 ]] || die "--state-dir requires a directory"
            STATE_DIR=$2
            shift 2
            ;;
        --keep-workdir)
            KEEP_WORKDIR=true
            shift
            ;;
        -n | --dry-run)
            DRY_RUN=true
            shift
            ;;
        --status)
            SHOW_STATUS=true
            shift
            ;;
        --list-steps)
            printf '%s\n' "${STEPS[@]}"
            exit 0
            ;;
        -h | --help)
            usage
            exit 0
            ;;
        --)
            shift
            REFERENCE_CONFIGS+=("$@")
            break
            ;;
        -*)
            die "unknown option: $1 (see --help)"
            ;;
        *)
            REFERENCE_CONFIGS+=("$1")
            shift
            ;;
    esac
done

SHOW_STATUS=${SHOW_STATUS:-false}

for step in ${FORCED_STEPS[@]+"${FORCED_STEPS[@]}"} ${ONLY_STEPS[@]+"${ONLY_STEPS[@]}"} \
    ${SKIP_STEPS[@]+"${SKIP_STEPS[@]}"}; do
    step_exists "$step" || die "unknown step: ${step} (see --list-steps)"
done

require_var() {
    local name=$1
    [[ -n ${!name:-} ]] || die "environment variable ${name} must be set (see --help)"
}

for var in PIPELINE_NAME TAG_OR_BRANCH PIPELINE_GITHUB_REPO PYTHON_VERSION \
    CONFIG_NAME CONFIG_VERSION CONFIG_GITHUB_REPO; do
    require_var "$var"
done

APPT_CACHE_STATUS=${APPT_CACHE_STATUS:-none}
case $APPT_CACHE_STATUS in
    build | update)
        require_var PATH_TO_apptainer_cache
        ;;
    none) ;;
    *)
        die "APPT_CACHE_STATUS must be build, update or none (got '${APPT_CACHE_STATUS}')"
        ;;
esac
GIT_CLONE_RETRIES=${GIT_CLONE_RETRIES:-3}

START_DIR=$PWD
WORK_DIR=${START_DIR}/${PIPELINE_NAME}
if [[ -n $STATE_DIR && $STATE_DIR != /* ]]; then
    STATE_DIR=${START_DIR}/${STATE_DIR}
fi
BUNDLE=${PIPELINE_NAME}_${TAG_OR_BRANCH}
ENV_DIR=${WORK_DIR}/${BUNDLE}_env
STAGE_DIR=${WORK_DIR}/${BUNDLE}
CONFIG_DIR=${STAGE_DIR}/${CONFIG_NAME}_${CONFIG_VERSION}
STAGED_CONFIG=${STAGE_DIR}/${PIPELINE_NAME}/config/config.yaml
APPT_CACHE_DIR=${WORK_DIR}/apptainer_cache
REF_DATA_DIR=${WORK_DIR}/design_and_ref_files
# Every step uses absolute paths, so a resumed run does not depend on the working directory.
: "${STATE_DIR:=${WORK_DIR}/.build_state}"

marker() { printf '%s/%s.done' "$STATE_DIR" "$1"; }
is_done() { [[ -f $(marker "$1") ]]; }

mark_done() {
    local path
    path=$(marker "$1")
    mkdir -p "$(dirname "$path")"
    date -u +"%Y-%m-%dT%H:%M:%SZ" >"$path"
}

clear_marker() {
    rm -f "$(marker "$1")"
    # Sub-markers of a step live in a directory named after it.
    rm -rf "${STATE_DIR:?}/$1"
}

# ---------------------------------------------------------------------------- status
if [[ $SHOW_STATUS == true ]]; then
    printf 'state directory: %s\n\n' "$STATE_DIR"
    for step in "${STEPS[@]}"; do
        if is_done "$step"; then
            printf '  %-20s done    %s\n' "$step" "$(cat "$(marker "$step")")"
        else
            printf '  %-20s pending\n' "$step"
        fi
    done
    exit 0
fi

# -------------------------------------------------------------------- executor plugins
PREVIOUS_EXECUTORS_FILE=${STATE_DIR}/install_executors/plugins
# The previous choice, as short names (slurm==2.8.0), for the default answer.
PREVIOUS_EXECUTORS=$(sed 's/snakemake-executor-plugin-//g' "$PREVIOUS_EXECUTORS_FILE" 2>/dev/null || true)

ask_executor_plugins() {
    local answer i default=${PREVIOUS_EXECUTORS:-none}
    {
        printf '\nSnakemake executor plugins to add to the env (needs python >= 3.11):\n'
        for i in "${!KNOWN_EXECUTORS[@]}"; do
            printf '  %d) %s\n' $((i + 1)) "${KNOWN_EXECUTORS[$i]}"
        done
        printf 'Numbers or names, comma/space separated, "none" for no plugins [%s]: ' "$default"
    } >&2
    read -r answer || answer=""
    answer=${answer:-$default}
    # Numbers become names; anything else (names, version specs) passes through.
    for i in $(split_list "$answer"); do
        if [[ $i =~ ^[0-9]+$ ]] && ((i >= 1 && i <= ${#KNOWN_EXECUTORS[@]})); then
            printf '%s ' "${KNOWN_EXECUTORS[$((i - 1))]}"
        else
            printf '%s ' "$i"
        fi
    done
}

if [[ -z ${EXECUTOR_PLUGINS+set} ]]; then
    if [[ -t 0 ]]; then
        EXECUTOR_PLUGINS=$(ask_executor_plugins)
    else
        EXECUTOR_PLUGINS=$PREVIOUS_EXECUTORS
    fi
fi

# Each entry is a pip requirement, e.g. snakemake-executor-plugin-slurm==2.8.0.
EXECUTOR_REQUIREMENTS=()
for spec in $(split_list "$EXECUTOR_PLUGINS"); do
    [[ $spec == none ]] && continue
    name=${spec%%[<>=!~]*}
    in_list "$name" "${KNOWN_EXECUTORS[@]}" ||
        die "unknown executor plugin '${name}' (known: ${KNOWN_EXECUTORS[*]})"
    EXECUTOR_REQUIREMENTS+=("snakemake-executor-plugin-${spec}")
done
if [[ ${#EXECUTOR_REQUIREMENTS[@]} -gt 0 ]]; then
    # Every executor plugin requires python >= 3.11; fail now rather than deep inside pip.
    IFS=. read -r py_major py_minor _ <<<"$PYTHON_VERSION"
    if ((py_major < 3 || (py_major == 3 && ${py_minor:-0} < 11))); then
        die "executor plugins need PYTHON_VERSION >= 3.11 (got ${PYTHON_VERSION})"
    fi
fi
# Sorted, so the same set in a different order counts as unchanged.
EXECUTOR_SET=$(printf '%s\n' ${EXECUTOR_REQUIREMENTS[@]+"${EXECUTOR_REQUIREMENTS[@]}"} | sort -u | paste -sd' ')
log "executor plugins: ${EXECUTOR_SET:-none}"

# ------------------------------------------------------------------- marker bookkeeping
# The build is a linear pipeline, so redoing a step invalidates everything after it.
invalidate_from() {
    local from=$1 step seen=false
    for step in "${STEPS[@]}"; do
        [[ $step == "$from" ]] && seen=true
        [[ $seen == true ]] && clear_marker "$step"
    done
}

if [[ $DRY_RUN == false ]]; then
    mkdir -p "$STATE_DIR"

    # A state directory only describes the build it was created for.
    FINGERPRINT="${PIPELINE_NAME} ${TAG_OR_BRANCH} ${CONFIG_NAME} ${CONFIG_VERSION} ${PYTHON_VERSION}"
    FINGERPRINT_FILE=${STATE_DIR}/build.fingerprint
    if [[ -f $FINGERPRINT_FILE && $(cat "$FINGERPRINT_FILE") != "$FINGERPRINT" ]]; then
        if [[ $FORCE_ALL == true ]]; then
            printf '%s\n' "$FINGERPRINT" >"$FINGERPRINT_FILE"
        else
            die "$(
                printf 'state directory %s belongs to a different build:\n' "$STATE_DIR"
                printf '  existing: %s\n' "$(cat "$FINGERPRINT_FILE")"
                printf '  current:  %s\n' "$FINGERPRINT"
                printf 'Rerun with --force, or pass --state-dir for this build.'
            )"
        fi
    else
        printf '%s\n' "$FINGERPRINT" >"$FINGERPRINT_FILE"
    fi

    if [[ $FORCE_ALL == true ]]; then
        log "--force given, clearing all markers"
        invalidate_from "${STEPS[0]}"
    fi

    for step in ${FORCED_STEPS[@]+"${FORCED_STEPS[@]}"}; do
        log "--force-step ${step}, clearing it and every later step"
        invalidate_from "$step"
    done

    for step in ${ONLY_STEPS[@]+"${ONLY_STEPS[@]}"}; do
        clear_marker "$step"
    done
fi

# The packed env must hold exactly the requested executor plugins, so a changed selection
# redoes install_executors and everything after it. A state directory from
# before install_executors existed has the env packed already, without any plugins.
if { is_done install_executors &&
    [[ $(cat "$PREVIOUS_EXECUTORS_FILE" 2>/dev/null) != "$EXECUTOR_SET" ]]; } ||
    { ! is_done install_executors && [[ -n $EXECUTOR_SET ]] && is_done pack_env; }; then
    if [[ $DRY_RUN == true ]]; then
        log "executor plugins changed, install_executors and every later step would be redone"
    else
        log "executor plugins changed, redoing install_executors and every later step"
        invalidate_from install_executors
    fi
fi

# ------------------------------------------------------------------------- environment
unset PYTHONPATH
export PYTHONNOUSERSITE=1
export TAG_OR_BRANCH

# conda's shell functions are not written for `set -u`.
set +u
eval "$(conda shell.bash hook)"
set -u

ENV_ACTIVE=false
ensure_env_active() {
    [[ $ENV_ACTIVE == true ]] && return 0
    [[ -d $ENV_DIR ]] || die "conda env ${ENV_DIR} is missing; rerun with --force-step create_env"
    set +u
    conda activate "$ENV_DIR"
    set -u
    ENV_ACTIVE=true
}

require_dir() {
    local path=$1 step=$2
    [[ -d $path ]] || die "${path} is missing; rerun with --force-step ${step}"
}

resolve_reference_config() {
    # A reference config can be given as an absolute path, or relative to any of the directories
    # the build works with: the pipeline clone the script runs from, the cloned config repo, the
    # bundle root (which holds the config repo and the staged pipeline clone), or the directory the
    # script was started from. First match wins, and the resolved path is logged before download.
    local given=$1 root candidate

    if [[ $given == /* ]]; then
        [[ -f $given ]] || die "reference config ${given} not found"
        printf '%s\n' "$given"
        return 0
    fi

    for root in "$PWD" "$CONFIG_DIR" "$STAGE_DIR" "$START_DIR"; do
        candidate=${root}/${given}
        if [[ -f $candidate ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done

    # Until 2025-11 the config repo was cloned unversioned into the pipeline clone, so paths of the
    # form GMS560_config/config/... (as in the README) predate the versioned CONFIG_DIR name.
    if [[ $given == "${CONFIG_NAME}/"* ]]; then
        candidate=${CONFIG_DIR}/${given#"${CONFIG_NAME}/"}
        if [[ -f $candidate ]]; then
            # stdout of this function is the resolved path, so the note goes to stderr.
            log "  ${given} -> ${CONFIG_NAME}_${CONFIG_VERSION} (the config repo is versioned now)" >&2
            printf '%s\n' "$candidate"
            return 0
        fi
    fi

    die "$(
        printf 'reference config %s not found. Looked in:\n' "$given"
        printf '  %s (the pipeline clone)\n' "$PWD"
        printf '  %s (the %s clone)\n' "$CONFIG_DIR" "$CONFIG_NAME"
        printf '  %s (the bundle)\n' "$STAGE_DIR"
        printf '  %s (where the script was started)' "$START_DIR"
    )"
}

git_clone() {
    # git_clone <url> <destination> [branch]
    local url=$1 dest=$2 branch=${3:-} attempt=1
    local -a args
    while :; do
        # A leftover directory means a previous attempt was interrupted.
        rm -rf "$dest"
        args=(clone --quiet)
        if [[ -n $branch ]]; then
            args+=(--branch "$branch")
        fi
        if git "${args[@]}" "$url" "$dest"; then
            return 0
        fi
        if ((attempt >= GIT_CLONE_RETRIES)); then
            die "failed to clone ${url} after ${attempt} attempt(s)"
        fi
        log "clone of ${url} failed, retrying ($((attempt + 1))/${GIT_CLONE_RETRIES})"
        attempt=$((attempt + 1))
        sleep 5
    done
}

# ------------------------------------------------------------------------------- steps

step_clone_pipeline() {
    # The build runs inside a clone of the pipeline; everything else is created below it.
    if [[ -d ${WORK_DIR}/.git ]]; then
        log "pipeline clone ${WORK_DIR} already present, keeping it"
    else
        git_clone "$PIPELINE_GITHUB_REPO" "$WORK_DIR" "$TAG_OR_BRANCH"
    fi
}

step_create_env() {
    rm -rf "$ENV_DIR"
    mamba create --prefix "$ENV_DIR" "python=${PYTHON_VERSION}" -y
}

step_stage_pipeline() {
    rm -rf "$STAGE_DIR"
    mkdir -p "$STAGE_DIR"
    git_clone "$PIPELINE_GITHUB_REPO" "${STAGE_DIR}/${PIPELINE_NAME}" "$TAG_OR_BRANCH"
}

step_install_requirements() {
    require_dir "${STAGE_DIR}/${PIPELINE_NAME}" stage_pipeline
    ensure_env_active
    # No --ignore-installed here: it makes pip reinstall the pip/setuptools that conda placed in
    # the prefix, and conda pack then refuses the env ("Files managed by conda were found to have
    # been deleted/overwritten"). PYTHONNOUSERSITE and the unset PYTHONPATH above already keep
    # packages from outside the prefix out of the install.
    "${ENV_DIR}/bin/pip3" install --no-cache-dir \
        -r "${STAGE_DIR}/${PIPELINE_NAME}/requirements.txt"
}

step_install_executors() {
    local package
    ensure_env_active
    # Drop known plugins that are no longer requested, so the env matches EXECUTOR_PLUGINS.
    for package in "${KNOWN_EXECUTORS[@]}"; do
        package=snakemake-executor-plugin-${package}
        if ! printf '%s\n' ${EXECUTOR_REQUIREMENTS[@]+"${EXECUTOR_REQUIREMENTS[@]}"} |
            grep -qE "^${package}([<>=!~]|$)" &&
            "${ENV_DIR}/bin/pip3" show --quiet "$package" >/dev/null 2>&1; then
            log "  removing ${package}"
            "${ENV_DIR}/bin/pip3" uninstall --yes "$package"
        fi
    done

    if [[ ${#EXECUTOR_REQUIREMENTS[@]} -eq 0 ]]; then
        log "  no executor plugins selected, nothing to install"
    else
        require_dir "${STAGE_DIR}/${PIPELINE_NAME}" stage_pipeline
        log "  installing ${EXECUTOR_REQUIREMENTS[*]}"
        # requirements.txt goes along so pip keeps its pins (snakemake>=9,<10) while resolving
        # the plugins' snakemake-interface-* dependencies.
        "${ENV_DIR}/bin/pip3" install --no-cache-dir \
            -r "${STAGE_DIR}/${PIPELINE_NAME}/requirements.txt" \
            "${EXECUTOR_REQUIREMENTS[@]}"
    fi

    mkdir -p "$(dirname "$PREVIOUS_EXECUTORS_FILE")"
    printf '%s\n' "$EXECUTOR_SET" >"$PREVIOUS_EXECUTORS_FILE"
}

step_pack_env() {
    ensure_env_active
    # requirements.txt pins packages that conda also ships in the prefix (setuptools), so pip
    # replaces some conda-managed files no matter what. The packed env is complete and only ever
    # used offline, so let conda pack proceed instead of aborting the build.
    conda pack --prefix "$ENV_DIR" -o "${STAGE_DIR}/env.tar.gz" --force --ignore-missing-files
}

step_clone_wrappers() {
    git_clone https://github.com/snakemake/snakemake-wrappers.git "${STAGE_DIR}/snakemake-wrappers"
}

step_clone_modules() {
    # One marker per module, so an interrupted run does not re-clone the whole set.
    local module dest
    mkdir -p "${STAGE_DIR}/hydra-genetics"
    for module in "${HYDRA_MODULES[@]}"; do
        dest=${STAGE_DIR}/hydra-genetics/${module}
        if is_done "clone_modules/${module}" && [[ -d ${dest}/.git ]]; then
            log "  ${module} already cloned"
            continue
        fi
        log "  cloning ${module}"
        git_clone "https://github.com/hydra-genetics/${module}.git" "$dest"
        mark_done "clone_modules/${module}"
    done
}

step_clone_config() {
    git_clone "$CONFIG_GITHUB_REPO" "$CONFIG_DIR" "$CONFIG_VERSION"
}

step_update_config_paths() {
    # Substitute ${TAG_OR_BRANCH} into the profiles and the production config.
    local relative path
    require_dir "$CONFIG_DIR" clone_config
    for relative in "${CONFIG_SUBST_FILES[@]}"; do
        path=${CONFIG_DIR}/${relative}
        [[ -f $path ]] || die "expected ${path} in ${CONFIG_NAME} ${CONFIG_VERSION}"
        log "  substituting TAG_OR_BRANCH in ${relative}"
        envsubst '$TAG_OR_BRANCH' <"$path" >"${path}.tmp"
        mv "${path}.tmp" "$path"
    done
}

step_containers() {
    if [[ $APPT_CACHE_STATUS == none ]]; then
        log "  APPT_CACHE_STATUS=none, nothing to do"
        return 0
    fi
    [[ -f $STAGED_CONFIG ]] || die "${STAGED_CONFIG} is missing; rerun with --force-step stage_pipeline"
    ensure_env_active

    if [[ $APPT_CACHE_STATUS == build ]]; then
        if is_done containers/download; then
            log "  containers already downloaded"
        else
            hydra-genetics prepare-environment create-singularity-files \
                -c "$STAGED_CONFIG" -o "$APPT_CACHE_DIR"
            mark_done containers/download
        fi
    fi

    # Point the staged config at the container cache as it will be seen on the target system.
    cp "$STAGED_CONFIG" "${STAGED_CONFIG}.copy"
    hydra-genetics prepare-environment container-path-update \
        -c "${STAGED_CONFIG}.copy" -n "$STAGED_CONFIG" -p "$PATH_TO_apptainer_cache"
}

step_pack_bundle() {
    require_dir "$STAGE_DIR" stage_pipeline
    rm -f "${WORK_DIR}/${BUNDLE}.tar.gz"
    tar -zcf "${WORK_DIR}/${BUNDLE}.tar.gz" -C "$WORK_DIR" "$BUNDLE"
}

step_references() {
    if [[ ${#REFERENCE_CONFIGS[@]} -eq 0 ]]; then
        log "  no reference configs given, nothing to do"
        return 0
    fi
    ensure_env_active

    # Resolved up front: hydra-genetics opens these paths itself, so a wrong one would otherwise
    # only surface as a FileNotFoundError traceback once the download is already running.
    local -a resolved=()
    local reference_config path key
    for reference_config in "${REFERENCE_CONFIGS[@]}"; do
        path=$(resolve_reference_config "$reference_config")
        resolved+=("$path")
    done

    for path in "${resolved[@]}"; do
        # Basename plus a checksum of the full path: readable, unique, and short enough to be a
        # filename whatever the path length.
        key=references/$(basename "$path")-$(printf '%s' "$path" | cksum | cut -d' ' -f1)
        if is_done "$key"; then
            log "  $(basename "$path") already downloaded"
            continue
        fi
        log "  downloading ${path}"
        hydra-genetics --debug references download -o "$REF_DATA_DIR" -v "$path"
        mark_done "$key"
    done

    rm -f "${WORK_DIR}/design_and_ref_files.tar.gz"
    tar -czf "${WORK_DIR}/design_and_ref_files.tar.gz" -C "$WORK_DIR" design_and_ref_files
}

step_cleanup() {
    if [[ $KEEP_WORKDIR == true ]]; then
        log "  --keep-workdir given, keeping ${ENV_DIR} and ${STAGE_DIR}"
        return 0
    fi
    if [[ $ENV_ACTIVE == true ]]; then
        set +u
        conda deactivate
        set -u
        ENV_ACTIVE=false
    fi
    rm -rf "$ENV_DIR" "$STAGE_DIR"
}

# --------------------------------------------------------------------------- run steps
ran_any=false

run_step() {
    local step=$1
    if in_list "$step" "${SKIP_STEPS[@]+"${SKIP_STEPS[@]}"}"; then
        log "skipping ${step} (--skip)"
        return 0
    fi
    if [[ ${#ONLY_STEPS[@]} -gt 0 ]] && ! in_list "$step" "${ONLY_STEPS[@]}"; then
        return 0
    fi
    if is_done "$step"; then
        log "skipping ${step} (already done)"
        return 0
    fi
    if [[ $DRY_RUN == true ]]; then
        if [[ $step == install_executors ]]; then
            log "would run ${step} (${EXECUTOR_SET:-no plugins})"
        else
            log "would run ${step}"
        fi
        return 0
    fi

    log "${step}"
    "step_${step}"
    mark_done "$step"
    ran_any=true
}

# The pipeline clone has to exist before anything else can run.
run_step "${STEPS[0]}"

# The rest of the build runs from inside the clone, so relative paths on the command line -- the
# reference configs, e.g. config/references/references.hg19.yaml -- resolve against it. Everything
# the script writes uses an absolute path, so this only affects paths the caller supplied.
if [[ $DRY_RUN == false ]]; then
    [[ -d $WORK_DIR ]] || die "${WORK_DIR} is missing; rerun with --force-step clone_pipeline"
    cd "$WORK_DIR"
fi

for step in "${STEPS[@]:1}"; do
    run_step "$step"
done

if [[ $DRY_RUN == true ]]; then
    log "dry run, nothing was changed"
elif [[ $ran_any == false ]]; then
    log "all steps already complete, nothing to do (see --status, --force)"
else
    log "build complete"
fi
