#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l mem=32G
#PBS -l walltime=24:0:0
#PBS -N prepare_reporting
set -eu -o pipefail

# Load Modules
module purge 
module load mugqic/R_Bioconductor/3.5.3_3.8
PROJ="/lustre03/project/6007512/C3G/projects/Moreira_COVID19_Genotyping"
REPORT_TMPLTS="${PROJ}/hgalvez/sample_reporting"
GENPIPES_VERSION="covseq_v1.0beta"

# Define samples and names
RUN_NAME="test_run"
RUN_PATH="${PROJ}/hgalvez/sample_reporting/${RUN_NAME}"
mkdir -p report/sample_reports 
cd ${RUN_PATH}/report

cp ${REPORT_TMPLTS}/run_report.tmplt.Rmd run_report.Rmd
cp ${REPORT_TMPLTS}/generate_report_tables.tmplt.R generate_report_tables.R

echo "run_name , " ${RUN_NAME} > run_metadata.csv
echo "genpipes_version , " ${GENPIPES_VERSION} >> run_metadata.csv
grep "^cluster_server" ../CoVSeQ.config.trace.ini | sed s:=:,:g >> run_metadata.csv
grep "^assembly_" ../CoVSeQ.config.trace.ini | sed s:=:,:g >> run_metadata.csv
grep -m 1 "^sequencing_technology" ../CoVSeQ.config.trace.ini | sed s:=:,:g >> run_metadata.csv

grep "^module_" ../CoVSeQ.config.trace.ini | sed s:=:,:g > module_table.tmp.csv


#cat tracks.tmplt.R | sed s:REPLACE:$RUN_NAME:g > $RUN_NAME/tracks.${RUN_NAME}.R 
#cat metricsPlots.tmplt.R | sed s:REPLACE:$RUN_NAME:g > $RUN_NAME/metricsPlots.${RUN_NAME}.R 

Rscript generate_report_tables.R

Rscript -e "rmarkdown::render('run_report.Rmd')"

rm module_table.tmp.csv
