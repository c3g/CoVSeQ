#!/bin/bash
set -eu -o pipefail

# A simple script that will use find to search for the appropriate files in the output of a production run
# Assumes it is being run from the top of the genpipes project directory

export SAMPLEID=${1}

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
