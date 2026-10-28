#Building version v1.9.0

Note! Before running the build script conda needs to be available. On HPC-clusters conda can usually be accessed through the 
module system, e.g. with ```module load miniconda3```. 

If you run on the compute node:
export SINGULARITY_CACHEDIR=/projects/wp4/nobackup/singularity_cache/{user}

```bash
TAG_OR_BRANCH="v1.9.0" CONFIG_VERSION="v1.32.0" PIPELINE_NAME="Twist_Solid" PYTHON_VERSION="3.9" \
PIPELINE_GITHUB_REPO="https://github.com/genomic-medicine-sweden/Twist_Solid.git" \
CONFIG_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/GMS560_config.git" \
CONFIG_NAME="GMS560_config" \
PATH_TO_apptainer_cache="path/to/.sif/on/remote" \
APPT_CACHE_STATUS="build" \
PGX_VERSION="v0.4.0" PGX_GITHUB_REPO="https://github.com/genomic-medicine-sweden/pgx.git" \
bash build_conda.sh config/references/design_files.hg19.yaml config/references/nextseq.hg19.pon.yaml \
config/references/references.hg19.yaml GMS560_config/config/reference_files.yaml
```

The PGX add-on pipeline is included with `PGX_VERSION` and `PGX_GITHUB_REPO`. PGX gets its own environment, built with
`PGX_PYTHON_VERSION` (defaults to `PYTHON_VERSION`). Without `PGX_VERSION`, PGX is left out.


Copy the following files and folders to the cluster (eg Miarka):
* Twist_Solid_{version}.tar.gz
* design_and_ref_files.tar.gz
* apptainer_cache


Next step: Extract the tar files and the env file inside the Twist_Solid_{version} folder:
```bash
cd Twist_Solid_{version}
mkdir Twist_Solid/venv/
tar -xzvf env.tar.gz -C Twist_Solid/venv/
mkdir pgx_{pgx_version}/venv/
tar -xzvf pgx_env.tar.gz -C pgx_{pgx_version}/venv/
```
