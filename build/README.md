# Building an offline bundle

`build_conda.sh` builds a self-contained bundle (pipeline, conda env, hydra-genetics modules,
snakemake-wrappers, an optional config repo, and optionally containers and reference files) to install
on a system without internet access.

Note! Before running the build script conda needs to be available. On HPC-clusters conda can usually be
accessed through the module system, e.g. with ```module load miniconda3```.

If you run on the compute node:
export SINGULARITY_CACHEDIR=/projects/wp4/nobackup/singularity_cache/{user}

## Example: building version v0.24.0

```bash
TAG_OR_BRANCH="v0.24.0" CONFIG_VERSION="v1.18.0" PIPELINE_NAME="Twist_Solid" PYTHON_VERSION="3.9" \
PIPELINE_GITHUB_REPO="https://github.com/genomic-medicine-sweden/Twist_Solid.git" \
CONFIG_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/GMS560_config.git" \
CONFIG_NAME="GMS560_config" \
PATH_TO_apptainer_cache="path/to/.sif/on/remote" \
APPT_CACHE_STATUS="build" \
EXECUTOR_PLUGINS="none" \
bash build_conda.sh config/references/design_files.hg19.yaml config/references/nextseq.hg19.pon.yaml \
config/references/references.hg19.yaml GMS560_config/config/reference_files.yaml
```

The script creates a working directory `<PIPELINE_NAME>/` (here `Twist_Solid/`) in the directory it is
started from; everything it produces ends up there.

Copy the following files and folders from `<PIPELINE_NAME>/` to the cluster (eg Miarka):
* Twist_Solid_{version}.tar.gz
* design_and_ref_files.tar.gz
* apptainer_cache

Next step: Extract the tar files and the env file inside the Twist_Solid_{version} folder:
```bash
cd Twist_Solid_{version}
mkdir Twist_Solid/venv/
tar -xzvf env.tar.gz -C Twist_Solid/venv/
```

## Configuration

All configuration is passed as environment variables; positional arguments are reference/design config
files handed to `hydra-genetics references download`.

| Variable | Required | Description |
| --- | --- | --- |
| `PIPELINE_NAME` | yes | e.g. `Twist_Solid` |
| `TAG_OR_BRANCH` | yes | pipeline tag or branch, e.g. `v0.24.0` |
| `PIPELINE_GITHUB_REPO` | yes | clone url of the pipeline |
| `PYTHON_VERSION` | yes | python version of the conda env, e.g. `3.9`. Snakemake 9 based releases need `>= 3.11`. |
| `CONFIG_NAME` | no | e.g. `GMS560_config` |
| `CONFIG_VERSION` | no | config repo tag or branch, e.g. `v1.18.0` |
| `CONFIG_GITHUB_REPO` | no | clone url of the config repo |
| `APPT_CACHE_STATUS` | no | `build` (download containers and point the config at the cache), `update` (only point the config at an existing remote cache) or `none` (default) |
| `PATH_TO_apptainer_cache` | with `build`/`update` | path of the container cache *on the target system* |
| `EXECUTOR_PLUGINS` | no | comma-separated Snakemake executor plugins to add to the env: `drmaa`, `slurm`, `sge` or `none`. A version spec is allowed, e.g. `slurm==2.8.0`. Needs `PYTHON_VERSION >= 3.11`. When unset the script asks interactively; without a terminal it reuses the previous choice (`none` on a first build). |
| `GIT_CLONE_RETRIES` | no | attempts per git clone (default `3`) |

### Building without a config repo

The config repo is all or nothing: set `CONFIG_NAME`, `CONFIG_VERSION` and `CONFIG_GITHUB_REPO`
together, or leave **all three unset** to build without one. Setting only some of them is an error.

Without a config repo the `clone_config` and `update_config_paths` steps still run but do nothing, and
reference config paths can only be absolute or relative to the pipeline clone (see below).

```bash
TAG_OR_BRANCH="v0.24.0" PIPELINE_NAME="Twist_Solid" PYTHON_VERSION="3.9" \
PIPELINE_GITHUB_REPO="https://github.com/genomic-medicine-sweden/Twist_Solid.git" \
EXECUTOR_PLUGINS="none" \
bash build_conda.sh config/references/references.hg19.yaml
```

(`--skip clone_config,update_config_paths` also leaves the config repo out of the bundle, but the build
still counts as a config-repo build, so reference configs prefixed with the config repo name will fail.
Unsetting the variables is the intended way.)

### Reference config paths

Reference config arguments are resolved against the repos the build cloned:

* `/absolute/path.yaml` — used as is
* `<CONFIG_NAME>_<CONFIG_VERSION>/rest.yaml` or `<CONFIG_NAME>/rest.yaml` — `rest.yaml` in the cloned
  config repo (a different `CONFIG_VERSION` in the prefix is an error)
* anything else — relative to the pipeline clone `<PIPELINE_NAME>/`, falling back to the directory the
  script was started from

All paths are checked before any download starts.

## Steps

The build is a linear list of steps, run in this order (`bash build_conda.sh --list-steps`):

| Step | What it does |
| --- | --- |
| `clone_pipeline` | clones the pipeline into `<PIPELINE_NAME>/` (kept if already present); the build runs from there |
| `create_env` | creates the conda env `<PIPELINE_NAME>/<PIPELINE_NAME>_<TAG>_env` with `mamba` |
| `stage_pipeline` | clones the pipeline again into the staging directory `<PIPELINE_NAME>/<PIPELINE_NAME>_<TAG>/` (the bundle content) |
| `install_requirements` | `pip install -r requirements.txt` into the env |
| `install_executors` | installs/removes executor plugins so the env matches `EXECUTOR_PLUGINS` |
| `pack_env` | `conda pack` the env into `<staging>/env.tar.gz` |
| `clone_wrappers` | clones `snakemake-wrappers` into the staging directory |
| `clone_modules` | clones the hydra-genetics modules into `<staging>/hydra-genetics/` |
| `clone_config` | clones the config repo into `<staging>/<CONFIG_NAME>_<CONFIG_VERSION>/` |
| `update_config_paths` | substitutes `${TAG_OR_BRANCH}` into the Miarka profiles and production config of the config repo |
| `containers` | `APPT_CACHE_STATUS=build`: downloads containers into `<PIPELINE_NAME>/apptainer_cache/`; `build`/`update`: rewrites container paths in the staged `config/config.yaml` to `PATH_TO_apptainer_cache` |
| `pack_bundle` | tars the staging directory into `<PIPELINE_NAME>/<PIPELINE_NAME>_<TAG>.tar.gz` |
| `references` | downloads the reference configs into `design_and_ref_files/` and packs `design_and_ref_files.tar.gz` |
| `cleanup` | removes the conda env and the staging directory (skip with `--keep-workdir`) |

## Resuming and re-running steps

Every completed step writes a marker `<state-dir>/<step>.done` (default state dir:
`<PIPELINE_NAME>/.build_state`) and is skipped on the next run. Rerunning the same command after a
failure therefore continues where it stopped. The slow steps also keep sub-markers, so an interrupted
step resumes mid-way: one per module (`clone_modules/<module>`), one for the container download
(`containers/download`) and one per reference config (`references/...`).

```bash
bash build_conda.sh --status        # show which steps are done / pending
bash build_conda.sh --list-steps    # step names in execution order
bash build_conda.sh -n              # dry run: print what would run, change nothing
```

`--status`, `--list-steps` and `-n` still need the required environment variables set, since the state
directory is derived from them.

| Option | Effect |
| --- | --- |
| `-f`, `--force` | ignore all markers and rebuild everything |
| `--force-step STEP` | redo `STEP` **and every step after it** |
| `--only STEP` | run only `STEP` (ignoring its marker); later steps are left alone |
| `--skip STEP` | never run `STEP` in this invocation |
| `--state-dir DIR` | keep markers in `DIR` instead of the default |
| `--keep-workdir` | do not delete the conda env and staging directory in `cleanup` |
| `-n`, `--dry-run` | print what would run, change nothing |
| `--status` | print the state of every step and exit |
| `--list-steps` | print the step names and exit |

`--force-step`, `--only` and `--skip` can be repeated or given a comma-separated list. Remember to pass
the same environment variables and reference configs as for the original build.

Examples:

```bash
# Re-download only the reference files
bash build_conda.sh --only references <reference configs...>

# Re-clone the modules and rebuild everything that depends on them (repack bundle, ...)
bash build_conda.sh --force-step clone_modules <reference configs...>

# Change a step in the middle and repack the bundle, without touching later steps
bash build_conda.sh --only update_config_paths,pack_bundle

# Rebuild the env from scratch but skip downloading the containers again
bash build_conda.sh --force-step create_env --skip containers <reference configs...>
```

Things to keep in mind:

* `--only` does not invalidate later steps, so a changed staging directory is **not** repacked unless
  `pack_bundle` is included (`--only X,pack_bundle`). Use `--force-step` when everything downstream
  should be redone.
* `cleanup` deletes the env and the staging directory (which also contains the config repo clone).
  Steps that need them fail afterwards with a hint to `--force-step` an earlier step. Pass
  `--keep-workdir` if you expect to rerun individual steps.
* The state directory is tied to `PIPELINE_NAME`, `TAG_OR_BRANCH`, `CONFIG_NAME`, `CONFIG_VERSION` and
  `PYTHON_VERSION`. Changing any of them aborts with an error; use `--force` to rebuild in place, or
  `--state-dir` to keep a separate state for the new build.
* Changing `EXECUTOR_PLUGINS` automatically redoes `install_executors` and every later step.
* `-n` can be combined with `--force`, `--force-step` and `--only` to preview exactly which steps
  would be redone; the markers on disk are left untouched.
