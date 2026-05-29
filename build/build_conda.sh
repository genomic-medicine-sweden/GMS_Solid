#!/usr/bin/env bash
set -e

eval "$(conda shell.bash hook)"

# Clone git
git clone --branch ${TAG_OR_BRANCH} ${PIPELINE_GITHUB_REPO}
cd ${PIPELINE_NAME}

# Create and activate conda envrionmnet in the current directory, then install pipeline requirements
mamba create --prefix ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env python=${PYTHON_VERSION} -y
conda activate ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env

if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH} ];
then
    rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}
fi

mkdir ${PIPELINE_NAME}_${TAG_OR_BRANCH}
# clone the required version of the pipeline
git clone --branch ${TAG_OR_BRANCH} ${PIPELINE_GITHUB_REPO} ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}
# install the requirements for the pipeline
unset PYTHONPATH
export PYTHONNOUSERSITE=1
./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env/bin/pip3 install --no-cache-dir -I -r ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/requirements.txt 
# pack the environment with the requriements installed
conda pack --prefix ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env -o ${PIPELINE_NAME}_${TAG_OR_BRANCH}/env.tar.gz


# # Clone snakemake-wrappers and hydra-genetics
mkdir -p ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics

git clone https://github.com/snakemake/snakemake-wrappers.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/snakemake-wrappers

git clone https://github.com/hydra-genetics/alignment.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/alignment
git clone https://github.com/hydra-genetics/annotation.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/annotation
git clone https://github.com/hydra-genetics/biomarker.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/biomarker
git clone https://github.com/hydra-genetics/cnv_sv.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/cnv_sv
git clone https://github.com/hydra-genetics/compression.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/compression
git clone https://github.com/hydra-genetics/filtering.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/filtering
git clone https://github.com/hydra-genetics/fusions.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/fusions
git clone https://github.com/hydra-genetics/misc.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/misc
git clone https://github.com/hydra-genetics/mitochondrial ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/mitochondrial
git clone https://github.com/hydra-genetics/parabricks ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/parabricks
git clone https://github.com/hydra-genetics/prealignment.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/prealignment
git clone https://github.com/hydra-genetics/qc.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/qc
git clone https://github.com/hydra-genetics/reports.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/reports
git clone https://github.com/hydra-genetics/snv_indels.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/snv_indels
git clone https://github.com/hydra-genetics/references.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/references


## Download the config files from the config repo
git clone --branch ${CONFIG_VERSION} ${CONFIG_GITHUB_REPO} ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}

# Set pipeline version in config paths
export TAG_OR_BRANCH=${TAG_OR_BRANCH}
envsubst '$TAG_OR_BRANCH' < ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production/config.yaml \
	> ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production/config.yaml.tmp \
	&& mv ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production/config.yaml.tmp \
	${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production/config.yaml

envsubst '$TAG_OR_BRANCH' < ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_novaseqX/config.yaml \
	> ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_novaseqX/config.yaml.tmp \
	&& mv ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_novaseqX/config.yaml.tmp \
	${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_novaseqX/config.yaml

envsubst '$TAG_OR_BRANCH' < ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_ctDNA/config.yaml \
	> ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_ctDNA/config.yaml.tmp \
	&& mv ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_ctDNA/config.yaml.tmp \
	${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_ctDNA/config.yaml

envsubst '$TAG_OR_BRANCH' < ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_prio/config.yaml \
	> ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_prio/config.yaml.tmp \
	&& mv ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_prio/config.yaml.tmp \
	${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_prio/config.yaml

envsubst '$TAG_OR_BRANCH' < ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_prio_novaseqX/config.yaml \
	> ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_prio_novaseqX/config.yaml.tmp \
	&& mv ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_prio_novaseqX/config.yaml.tmp \
	${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/profiles/Miarka/production_prio_novaseqX/config.yaml

envsubst '$TAG_OR_BRANCH' < ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/config/Miarka/config_production_pipeline_miarka.yaml \
	> ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/config/Miarka/config_production_pipeline_miarka.yaml.tmp \
	&& mv ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/config/Miarka/config_production_pipeline_miarka.yaml.tmp \
	${PIPELINE_NAME}_${TAG_OR_BRANCH}/${CONFIG_NAME}_${CONFIG_VERSION}/config/Miarka/config_production_pipeline_miarka.yaml

# # Download containers and update container path.
if [ ${APPT_CACHE_STATUS} == "build" ];
then
    hydra-genetics prepare-environment create-singularity-files -c ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/config.yaml -o apptainer_cache
    cp ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/config.yaml ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/config.yaml.copy
    hydra-genetics prepare-environment container-path-update -c ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/config.yaml.copy -n ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/config.yaml -p ${PATH_TO_apptainer_cache}
fi

# If apptainer already available remote and no new build is needed (specified by user). Only update path to cache in config.
if [ ${APPT_CACHE_STATUS} == "update" ];
then
    cp ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/config.yaml ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/config.yaml.copy
    hydra-genetics prepare-environment container-path-update -c ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/config.yaml.copy -n ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/config.yaml -p ${PATH_TO_apptainer_cache}
fi


# # Pack all cloned repositories
tar -zcvf ${PIPELINE_NAME}_${TAG_OR_BRANCH}.tar.gz ${PIPELINE_NAME}_${TAG_OR_BRANCH}

# Download references if given on command line
if [ "$#" != 0 ];
then
    for reference_config in "$@"
    do
        hydra-genetics --debug references download -o design_and_ref_files -v $reference_config
    done
    # Compress data
    tar -czvf design_and_ref_files.tar.gz design_and_ref_files
fi


conda deactivate

if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH}_env ];
then
    rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}_env
fi

if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH} ];
then
    rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}
fi
