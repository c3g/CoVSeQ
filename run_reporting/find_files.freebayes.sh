#!/bin/bash
set -eu -o pipefail

# FOR USE WITH ncov_tools
# A simple script that will use find to search for the appropriate files 
# In the output of an freebayes illumina run and link then for ncov_tools
# Will only link datasets for which a complete set of files is found
# Assumes it is being run from the top of the genpipes project directory

export SAMPLEID=${1}

echo $SAMPLEID "," $(find alignment/ -name ${SAMPLEID}.sorted.filtered.primerTrim.bam -type f) \
    "," $(find consensus/ -name ${SAMPLEID}.consensus.*.fasta -type f) \
    "," $(find variant/ -name ${SAMPLEID}.freebayes_calling.consensus.vcf -type f)

mkdir -p report/ncov_tools/data 
if [ "$(find alignment/ -name ${SAMPLEID}.sorted.filtered.primerTrim.bam -type f)" != "" ] && \
   [ "$(find consensus/ -name ${SAMPLEID}.freebayes_calling.consensus.fasta -type f)" != "" ] && \
   [ "$(find variant/ -name ${SAMPLEID}.freebayes_calling.consensus.vcf -type f)" != "" ]; then
        ln -fs $(pwd -P )/$(find alignment/ -name ${SAMPLEID}.sorted.filtered.primerTrim.bam -type f) \
            report/ncov_tools/data/${SAMPLEID}.mapped.primertrimmed.sorted.bam
        
        ln -fs $(pwd -P )/$(find consensus/ -name ${SAMPLEID}.freebayes_calling.consensus.fasta -type f) \
            report/ncov_tools/data/${SAMPLEID}.consensus.fasta
        
        ln -fs $(pwd -P )/$(find variant/ -name ${SAMPLEID}.freebayes_calling.consensus.vcf -type f) \
            report/ncov_tools/data/${SAMPLEID}.variants.vcf
fi
