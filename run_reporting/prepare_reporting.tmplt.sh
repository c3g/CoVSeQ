#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l mem=32G
#PBS -l walltime=24:0:0
#PBS -N prepare_reporting
set -eu -o pipefail

# Load Modules
echo "Loading modules..."
module purge 
module load mugqic/R_Bioconductor/3.5.3_3.8
module load mugqic/bcftools/1.9
PROJ="/lustre03/project/6007512/C3G/projects/Moreira_COVID19_Genotyping"
REPORT_TMPLTS="${PROJ}/hgalvez/CoVSeQ/reporting"
COLLECT_METRICS="${PROJ}/hgalvez/CoVSeQ/metrics/covid_collect_metrics.sh"
GENPIPES_VERSION="genpipes/covid_release/1.1_beta"

# Define samples and names
RUN_NAME=""
RUN_PATH="/genfs/projects/COVID_full_processing/illumina/${RUN_NAME}"
cd ${RUN_PATH}

# Create ncov_tools links
echo "Linking files for ncov_tools..."
echo "sample,bam.path,fasta.path,tsv.path" > output_file_paths.csv
for item in $(cut -f 1 readset.txt | grep -v Sample); do
    bash ${REPORT_TMPLTS}/find_files.sh ${item} >> output_file_paths.csv
done

# Run Collect Metrics Script
echo "Collecting metrics..."
cp ${REPORT_TMPLTS}/prepare_collect_metrics.tmplt.R prepare_collect_metrics.R 
Rscript prepare_collect_metrics.R
bash ${COLLECT_METRICS} ${RUN_PATH}/report.readset.txt

# Create reporting scripts
mkdir -p report/sample_reports 
cd ${RUN_PATH}/report

echo "Copying template scripts..."
cp ${REPORT_TMPLTS}/run_report.tmplt.Rmd run_report.Rmd
cp ${REPORT_TMPLTS}/generate_report_tables.tmplt.R generate_report_tables.R
cp ${REPORT_TMPLTS}/prepare_ncov_tools.tmplt.R prepare_ncov_tools.R 
echo "cd $(pwd -P)/ncov_tools" > ncov_tools/snakemake_run_all.sh 
cat ${REPORT_TMPLTS}/snakemake_run_all.tmplt.sh >> ncov_tools/snakemake_run_all.sh

# Prepare ncov_tools
echo "Preparing to run ncov_tools..."
Rscript prepare_ncov_tools.R
cat ${REPORT_TMPLTS}/ncov_tools.config.tmplt.yaml | sed s:REPLACE-RUN:${RUN_NAME}: |\
 sed s:REPLACE-NEG-CTLS:$(cat neg_controls.txt): > ncov_tools/config.yaml
cat ${REPORT_TMPLTS}/ncov_tools.singularity.tmplt.sh | sed s:REPLACE-PWD:$(pwd -P): > ncov_tools/ncov_tools.singluarity.sh 
sbatch ncov_tools/ncov_tools.singluarity.sh

# Prepare run metadata
echo "run_name , " ${RUN_NAME} > run_metadata.csv
echo "genpipes_version , " ${GENPIPES_VERSION} >> run_metadata.csv
grep "^cluster_server" ../CoVSeQ.config.trace.ini | sed s:=:,:g >> run_metadata.csv
grep "^assembly_" ../CoVSeQ.config.trace.ini | sed s:=:,:g >> run_metadata.csv
grep -m 1 "^sequencing_technology" ../CoVSeQ.config.trace.ini | sed s:=:,:g >> run_metadata.csv

grep "^module_" ../CoVSeQ.config.trace.ini | sed s:=:,:g > module_table.tmp.csv

# Generate report tables
Rscript generate_report_tables.R

# Prepare problematic variants count

# Render report
Rscript -e "rmarkdown::render('run_report.Rmd', output_format = 'all')"

# Cleanup
rm module_table.tmp.csv
