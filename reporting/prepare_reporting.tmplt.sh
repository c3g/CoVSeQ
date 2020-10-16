#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l mem=32G
#PBS -l walltime=24:0:0
#PBS -N prepare_reporting
set -eu -o pipefail

# Load Modules
module purge 
module load mugqic/R_Bioconductor/3.5.3_3.8
module load mugqic/bcftools/1.9
PROJ="/lustre03/project/6007512/C3G/projects/Moreira_COVID19_Genotyping"
REPORT_TMPLTS="${PROJ}/hgalvez/CoVSeQ/reporting"
COLLECT_METRICS="${PROJ}/hgalvez/CoVSeQ/metrics/covid_collect_metrics.sh"
GENPIPES_VERSION="genpipes/covid_release/1.0"

# Define samples and names
RUN_NAME=""
RUN_PATH="${PROJ}/hgalvez/sample_reporting/${RUN_NAME}"
cd ${RUN_PATH}
bash ${COLLECT_METRICS} ${RUN_PATH}/readset.txt

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

# Prepare problematic variants count

Rscript -e "rmarkdown::render('run_report.Rmd', output_format = 'all')"

rm module_table.tmp.csv
