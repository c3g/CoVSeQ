#!/bin/bash
set -eu -o pipefail

##################################################
# Load modules
module purge
module load mugqic/R_Bioconductor/3.5.3_3.8 mugqic/bcftools/1.9

# Set file paths 
REPORT_PATH="/genfs/projects/COVID_full_processing/Report"
rm -r illumina_reports mgi_reports nanopore_reports
mkdir -p ${REPORT_PATH}/illumina_reports/{internal_metrics,ncov_tools}
mkdir -p ${REPORT_PATH}/mgi_reports/{internal_metrics,ncov_tools}
mkdir -p ${REPORT_PATH}/nanopore_reports/{internal_metrics,ncov_tools}

##################################################
# STEP 1: Prepare Illumina links and tables 
## Prepare full list of illumina samples
echo "Preparing Illumina reports..." 
echo "Preparing Illumina sample list" 
cd ${REPORT_PATH}/illumina_reports
find /genfs/projects/COVID_full_processing/illumina/ \
    -type f -wholename */alignment*/*/*.sorted.filtered.bam |\
    cut -d \/ -f 5,6,8 |\
    sed s:\/:,:g > illumina_samples.csv

## Prepare internal metrics reports
echo "Preparing Illumina internal metrics" 
cd ${REPORT_PATH}/illumina_reports/internal_metrics 
for item in $(cat ${REPORT_PATH}/illumina_reports/illumina_samples.csv | cut -f 2 -d , | sort | uniq); do 
  ln -sf $(find ../../../illumina/${item} -name report_metrics.csv -type f) \
    ${item}_report_metrics.csv
done   

## Prepare ncov_tools reports
echo "Preparing Illumina ncov_tools reports" 
cd ${REPORT_PATH}/illumina_reports/ncov_tools 
for item in $(find ../../../illumina/*/report -name *_summary_qc.tsv -type f); do ln -sf $item .; done
for item in $(find ../../../illumina/*/report -name *_negative_control_report.tsv -type f); do ln -sf $item .; done
for item in $(find ../../../illumina/*/report -name *_ambiguous_position_report.tsv -type f); do ln -sf $item .; done
for item in $(find ../../../illumina/*/report -name *_lineage_report.csv -type f); do ln -sf $item .; done

grep WARN *_negative_control_report.tsv > ${REPORT_PATH}/control_warnings/illumina_neg_control_warnings.txt 

##################################################
# STEP 2: Prepare MGI links and tables 
## Prepare full list of MGI runs
echo "Preparing MGI reports..." 
echo "Preparing MGI sample list" 
cd ${REPORT_PATH}/mgi_reports
find /genfs/projects/COVID_full_processing/mgi/ \
    -type f -wholename */alignment*/*/*.sorted.filtered.bam |\
    cut -d \/ -f 5,6,8 |\
    sed s:\/:,:g > mgi_samples.csv

## Prepare internal metrics reports
echo "Preparing MGI internal reports" 
cd ${REPORT_PATH}/mgi_reports/internal_metrics 
for item in $(cat ${REPORT_PATH}/mgi_reports/mgi_samples.csv | cut -f 2 -d , | sort | uniq); do 
  ln -sf $(find ../../../mgi/${item} -name report_metrics.csv -type f) \
    ${item}_report_metrics.csv
done   

## Prepare ncov_tools reports
# echo "Preparing MGI ncov_tools reports" 
#cd ${REPORT_PATH}/mgi_reports/ncov_tools 
#for item in $(find ../../../mgi/ -name *_summary_qc.tsv -type f); do ln -sf $item .; done
#for item in $(find ../../../mgi/ -name *_negative_control_report.tsv -type f); do ln -sf $item .; done
#for item in $(find ../../../mgi/ -name *_ambiguous_position_report.tsv -type f); do ln -sf $item .; done
#for item in $(find ../../../mgi/ -name *_lineage_report.csv -type f); do ln -sf $item .; done

##################################################
# STEP 3: Prepare Nanopore links and tables 
## Prepare full list of nanopore runs
echo "Preparing nanopore reports..." 
echo "Preparing nanopore sample lists"  
cd ${REPORT_PATH}/nanopore_reports
find /genfs/projects/COVID_full_processing/nanopore/ \
    -type f -wholename */analysis/*/*/*_commands.txt |\
    cut -d \/ -f 5,6,8,9 |\
    sed s:\/:,:g > nanopore_samples.csv

# Prepare internal reports links
echo "Preparing nanopore internal metrics"  
cd ${REPORT_PATH}/nanopore_reports/internal_metrics
cat ${REPORT_PATH}/nanopore.metrics.header.csv > full_nanopore_report.csv
for item in $(find ../../../nanopore/*/analysis/ -name *metrics.csv -type f); do
    echo $(tail -n 1 ${item})","${item} >> full_nanopore_report.csv 
done

## Prepare ncov_tools reports
echo "Preparing nanopore ncov_tools reports"  
cd ${REPORT_PATH}/nanopore_reports/ncov_tools 
for item in $(find ../../../nanopore/*/analysis/*/ncov_tools -name *_summary_qc.tsv -type f); do ln -sf $item .; done
for item in $(find ../../../nanopore/*/analysis/*/ncov_tools -name *_negative_control_report.tsv -type f); do ln -sf $item .; done
for item in $(find ../../../nanopore/*/analysis/*/ncov_tools -name *_ambiguous_position_report.tsv -type f); do ln -sf $item .; done
for item in $(find ../../../nanopore/*/analysis/*/ncov_tools -name *_lineage_report.csv -type f); do ln -sf $item .; done

##################################################
# Fetch latest freezeman report
echo "Fetch freezeman information" 
cd ${REPORT_PATH}/LIMS_reports
bash freezeman_fetch.sh $(date +"%Y-%m-%d") 
ln -sf $(date +"%Y-%m-%d")_freezeman_fetch.csv latest_freezeman_report.csv

##################################################
# Run R scripts to generate tables
echo "Parse LIMS reports"  
cd ${REPORT_PATH}/LIMS_reports
Rscript parse_LIMS_reports.R 


##################################################
# Generate release table 
echo "Generate release table" 
cd ${REPORT_PATH}/release_lists
grep "^>" latest_GISAID_fasta.fasta | sed s:\|:,:g | sed s:\>::g  > latest_release_list.csv
