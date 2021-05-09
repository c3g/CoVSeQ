#!/bin/bash
set -eu -o pipefail

# FOR USE WITH ncov_tools
# A simple script that will use find to search for the appropriate files 
# In the output of an iVar illumina run and link then for ncov_tools
# Will only link datasets for which a complete set of files is found
# Assumes it is being run from the top of the genpipes project directory

export SAMPLEID=${1}

# Save the location of all files in a file with all the info
echo $SAMPLEID "," $(find alignment/ -name ${SAMPLEID}.sorted.filtered.primerTrim.bam -type f) \
    "," $(find consensus/ -name ${SAMPLEID}.consensus.*.fasta -type f) \
    "," $(find variant/ -name ${SAMPLEID}.*tsv -type f)

mkdir -p report/ncov_tools/data 
if [ "$(find alignment/ -name ${SAMPLEID}.sorted.filtered.primerTrim.bam -type f)" != "" ] && \
   [ "$(find consensus/ -name ${SAMPLEID}.consensus.*.fasta -type f)" != "" ] && \
   [ "$(find variant/ -name ${SAMPLEID}*tsv -type f)" != "" ]; then
        ln -fs $(pwd -P )/$(find alignment/ -name ${SAMPLEID}.sorted.filtered.primerTrim.bam -type f) \
            report/ncov_tools/data/${SAMPLEID}.mapped.primertrimmed.sorted.bam
        
        ln -fs $(pwd -P )/$(find consensus/ -name ${SAMPLEID}.consensus.*.fasta -type f) \
            report/ncov_tools/data/${SAMPLEID}.consensus.fasta
        
        ln -fs $(pwd -P )/$(find variant/ -name ${SAMPLEID}.*tsv -type f) \
            report/ncov_tools/data/${SAMPLEID}.variants.tsv
fi
