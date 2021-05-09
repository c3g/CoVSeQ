# Load libraries
suppressPackageStartupMessages(library(tidyverse))
set.seed(123456789)

###############################################################################
# Create paths 
readset.file.path <- file.path("..", "readset.txt")
output.files.path <- file.path("..", "output_file_paths.csv") 

###############################################################################
# Read input data
readset.table <- readr::read_tsv(readset.file.path)
file.table <- readr::read_csv(output.files.path, na = c("", " ", "NULL"))

###############################################################################
# Produce ncov_tools tables 

complete.samples <- file.table %>% filter(!is.na(bam.path) & !is.na(fasta.path) & !is.na(tsv.path)) %>% pull(sample) 

# For now, Ct and Date fields are empty
ncov.tools.metadata <- tibble(sample = complete.samples, ct = NA, date = NA)

# Write a string with the negative controls for later configuration setup
neg.controls <- ncov.tools.metadata %>% filter(str_detect(sample, regex("((negctrl|ext)|ntc)|ctrl_neg", ignore_case = TRUE))) %>% pull(sample)
neg.controls <- str_c(neg.controls, collapse = "\",\"")
neg.controls <- paste0("\"", neg.controls, "\"")
if (neg.controls == "\"\""){
    neg.controls <- ""
}

write_lines(neg.controls, 
    path = "neg_controls.txt")

# Prepare the ncov_tools metadata file
write_tsv(ncov.tools.metadata, path = file.path("ncov_tools", "metadata.tsv"))


