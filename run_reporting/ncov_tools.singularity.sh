#!/bin/bash 
#SBATCH --time=12:00:00
#SBATCH --mem=92G
#SBATCH -N 1 -n 20
#SBATCH -J ncov-tools_run

# Load modules
module purge 
module load singularity/3.6 python/3.7
echo "Finished loading modules..."

export ENVDIR="/genfs/projects/analyste_dev/python_venvs/snakemake/bin/activate"

# Activate python virtual env
source ${ENVDIR}

export NCOVTOOLS_SIF="/genfs/projects/analyste_dev/singularity/images/ncov-tools_v1.1.sif"

# Define directories
PROJ="REPLACE-PWD"
RUN_DIR="${PROJ}/ncov_tools"

# Run commands
cd ${RUN_DIR}

singularity exec \
    -B /genfs:/genfs \
    -B /lustre01:/lustre01 \
    -B /lustre02:/lustre02 \
    -B /lustre03:/lustre03 \
    -B /lustre04:/lustre04   \
    -B /cvmfs:/cvmfs \
    -B /project:/project \
    -B /scratch:/scratch \
    --cleanenv ${NCOVTOOLS_SIF} \
    bash ${RUN_DIR}/snakemake_run_all.sh


deactivate

