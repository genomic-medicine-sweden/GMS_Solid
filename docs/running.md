# Running the pipeline

## Requirements
**Recommended hardware**

- CPU: >10 cores per sample
- Memory: 6GB per core
- Storage: >75GB per sample

**Note**: Running the pipeline with less resources may work, but has not been tested.

**Software**

- [python](https://www.python.org/), version 3.9 or newer
- [pip3](https://pypi.org/project/pip/)
- [virtualenv](https://docs.python.org/3/library/venv.html)
- [singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html)

**Nice to have**

- DRMAA compatible scheduler

## Installation
A list of releases of the Twist Solid pipeline can be found at: [Releases](https://github.com/genomic-medicine-sweden/Twist_Solid/releases).

### Clone the Twist Solid git repo
We recommend that the repository is cloned to your working directory. 
```bash
# Set up a working directory path
WORKING_DIRECTORY="/path_working_to_directory"
```

Fetch pipeline
```bash
# Set version
VERSION="v0.4.0"

# Clone selected version
git clone --branch ${VERSION} https://github.com/genomic-medicine-sweden/Twist_Solid.git ${WORKING_DIRECTORY}
```

### Create python environment
To run the Twist Solid pipeline a python virtual environment is needed.
```bash
# Enter working directory
cd ${WORKING_DIRECTORY}

# Create a new virtual environment
python3 -m venv ${WORKING_DIRECTORY}/virtual/environment
```

### Install pipeline requirements
Activate the virtual environment and install pipeline requirements specified in `requirements.txt`.
```bash
# Enter working directory
cd ${WORKING_DIRECTORY}

# Activate python environment
source virtual/environment/bin/activate

# Install requirements
pip install -r requirements.txt
```

### Setup required data and config

**Download data**
```bash
# make sure hydra-genetics is available
# make sure that TMPDIR points to a location with a lot of storage, it
# will be required to fetch reference data
# export TMPDIR=/PATH_TO_STORAGE
hydra-genetics --debug --verbose references download -o design_and_ref_files  -v config/references/references.hg19.yaml -v config/references/design_files.hg19.yaml -v config/references/nextseq.hg19.pon.yaml
```

**Update config**
```yaml 
# file config/config.data.hg19.yaml
# change rows:
PROJECT_DESIGN_DATA: "PATH_TO/design_and_ref_files" # parent folder for GMS560 design, ex GMS560/design
PROJECT_PON_DATA: "PATH_TO/design_and_ref_files" # artifact/background/PoN, ex GMS560/PoN
PROJECT_REF_DATA: "PATH_TO/design_and_ref_files" # parent folder for ref_data, ex ref_data/hg19

```

## Input sample files
The pipeline uses sample input files (`samples.tsv` and `units.tsv`) with information regarding sample information, sequencing meta information as well as the location of the fastq-files. Specification for the input files can be found at [Twist Solid schemas](https://github.com/genomic-medicine-sweden/Twist_Solid/blob/develop/workflow/schemas/). Using the python virtual environment created above it is possible to generate these files automatically using [hydra-genetics create-input-files](https://hydra-genetics.readthedocs.io/en/latest/run_pipeline/create_sample_files/):
```bash
hydra-genetics create-input-files -d path/to/fastq-files/
```

**Note**: Sample names cannot include "_" (underscore)!

## Profile and configuration
When running the pipeline, multiple configuration files should be used depending on the sample type. These can be specified in a Snakemake profile or directly in the command line using `--configfiles`.

### FFPE samples
For FFPE samples, the following configuration files are required:
- `config/config.yaml`
- `config/config.data.<genome_build>.yaml`

### ctDNA samples
For ctDNA samples, an additional configuration file is required to enable UMI processing:
- `config/config.yaml`
- `config/config.data.<genome_build>.yaml`
- `config/config.data.ctdna.<genome_build>.yaml`

### Specific sequencer configuration (optional)
For specific sequencers, an additional configuration file can be added:
- `config/config.data.<genome_build>.<sequencer>.yaml` (e.g., `config/config.data.hg19.novaseqX.yaml`)

## Run command
Using the activated python virtual environment created above, this is a basic command for running the pipeline:
```bash
snakemake --profile profiles/NAME_OF_PROFILE -s workflow/Snakefile
```  
<br />
There are many additional [snakemake running options](https://snakemake.readthedocs.io/en/stable/executing/cli.html#) some of which is listed below. However, options that are always used should be put in the [profile](https://hydra-genetics.readthedocs.io/en/latest/run_pipeline/profile/).

* --notemp - Saves all intermediate files. Good for development and testing different options.
* --until <rule> - Runs only rules dependent on the specified rule.

<br />
**Note:** Remember to have singularity and drmaa available on the system where the pipeline will be run.

<br />
