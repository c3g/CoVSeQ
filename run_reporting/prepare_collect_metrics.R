# Load libraries
suppressPackageStartupMessages(library(tidyverse))
set.seed(123456789)

###############################################################################
# Create paths 
readset.file.path <- file.path("readset.txt")
output.files.path <- file.path("output_file_paths.csv") 

###############################################################################
# Read input data
readset.table <- readr::read_tsv(readset.file.path)
file.table <- readr::read_csv(output.files.path, na = c("", " ", "NULL"))

###############################################################################
# Produce ncov_tools tables 

complete.samples <- file.table %>% filter(!is.na(bam.path) & !is.na(fasta.path) & !is.na(tsv.path)) %>% pull(sample) 

tmp.readset <- readset.table %>% filter(Sample %in% complete.samples) %>% distinct(.,Sample,.keep_all = T)

# Prepare the ncov_tools metadata file
write_tsv(tmp.readset, path = file.path("report.readset.txt"))


