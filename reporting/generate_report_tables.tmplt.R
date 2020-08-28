# Load libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(tidyverse))
set.seed(123456789)

###############################################################################
# Define Parsers
CoV2.genome <- Seqinfo(seqnames = c("MN908947.3"), seqlengths = c(29903), isCircular=c(FALSE), genome="SARSCoV2")

read_vcf_plus <- function(sample.name){
  print(sample.name)
  sample.vcf.path <- file.path("..","variant", sample.name, paste0(sample.name,".sorted.filtered.primerTrim.vcf"))
  sample.vcf <- readVcf(sample.vcf.path, CoV2.genome)
  
  if (dim(sample.vcf)[1] == 0 ){
    
    variant.ids <- c(NaN)
    position <- c(NaN)
    ref.allele <- c(NaN)
    alt.allele <- c(NaN)
    variation <- c(NaN) 
    alt.FREQ <- c(NaN)
    
      } 
  
  else if (dim(sample.vcf)[1] >= 1){ 
    
    variant.ids <- geno(sample.vcf)$ALT_FREQ %>% rownames()
    alt.FREQ  <- geno(sample.vcf)$ALT_FREQ %>% as_tibble() %>% pull(1) %>% as.double()
    position <- start(sample.vcf) %>% as.double()
    ref.allele <- ref(sample.vcf) %>% as.character()
    alt.allele <- alt(sample.vcf) %>% unlist() %>% as.character() 
    variation <- paste0(ref.allele, position, alt.allele)
    
    }
  
  tibble(variant.ids, position, ref.allele, alt.allele, variation, alt.FREQ)
}


sample_status <- function(percent.N){
  if (is.na(percent.N)) {
    percent.N <- 100
  }  


  if (percent.N <= 1){
    "PASS"
  }
  else if (percent.N > 1 & percent.N <= 5) {
    "FLAG"
  }
  else if (percent.N > 5) { 
    "REJ"}
  else {
    "NA"
  }
}

###############################################################################
# Create paths 
readset.file.path <- file.path("..", "readset.txt")
metrics.table.path <- file.path("..", "metrics", "metrics.csv")
host.metrics.table.path <- file.path("..", "metrics", "host_contamination_metrics.tsv")
metadata.path <- file.path("run_metadata.csv")
module.path <- file.path("module_table.tmp.csv")


###############################################################################
# Read input data
readset.table <- readr::read_tsv(readset.file.path)
metrics.table <- readr::read_csv(metrics.table.path) 
host.metrics.table <- readr::read_tsv(host.metrics.table.path)
metadata.table <- readr::read_csv(metadata.path, col_names = c("category", "value"))
module.table <- readr::read_csv(module.path, col_names = c("category", "value"))

###############################################################################
# Produce software versions table 
module.table %>% 
  dplyr::rename("Software Versions" = value) %>% 
  dplyr::select("Software Versions") %>% 
  write_csv(path="software_versions.csv")

###############################################################################
# Produce metrics table

## Define final column names (like Sandrine Requested)
final.column.names <- c("Sample",
                        "Nb reads",
                        "Percent human reads",
                        "Nb clean reads",
                        "Mean coverage",
                        "Percent N",
                        "Length low cov region (<20X)",
                        "Percent consensus > 100X",
                        "Length consensus",
                        "Nb variants > 10 perc allele freq",
                        "Nb variants > 75 perc allele freq",
                        "PASS/FLAG/REJ")


## Calculate variant numbers 

variant.numbers <- tibble(
  Sample = character(), 
  var.num.10.more = numeric(), 
  var.num.75.more = numeric()
)

for (sample in readset.table$Sample) {
  vcf.table <- read_vcf_plus(sample)
  write_csv(vcf.table, path = file.path("sample_reports", paste0(sample, "_vcf_info.csv")))
  var.num.10 <- vcf.table %>% filter(alt.FREQ > 0.10) %>% tally() %>% pull(n)
  var.num.75 <- vcf.table %>% filter(alt.FREQ > 0.75) %>% tally() %>% pull(n)
  variant.numbers <- add_row(variant.numbers, Sample = sample, var.num.10.more = var.num.10, var.num.75.more = var.num.75) 
}


## Join all tables 
full.table <- left_join(readset.table, metrics.table, by = c("Sample" = "sample"))
full.table <- left_join(full.table, host.metrics.table, by = "Sample")
full.table <- left_join(full.table, variant.numbers, by = "Sample")

## Calculate last missing metrics 
full.table <- full.table %>% mutate(cons.per.N = as.numeric(cons.per.N)) 

sample.status <- full.table %>% pull(cons.per.N) %>% map(sample_status) %>% unlist()

full.table <- full.table %>% 
                  mutate(total.reads = Total_aligned + Unmapped_only) %>% 
                  mutate(status = sample.status) %>% 
                  mutate(Human_only_perc = round(Human_only_perc, 2)) %>% 
                  mutate(consensus.length = 29903 - (cons.per.N /100 * 29903)) %>% mutate(consensus.length = as.integer(consensus.length)) %>% 
                  mutate(length.lowcov = 29903 - (bam.perc.20x /100 * 29903)) %>% mutate(length.lowcov = as.integer(length.lowcov)) 
                  
                  
final.columns <- c("Sample",
                   "total.reads",
                   "Human_only_perc", 
                   "SARS_only",
                   "bam.mean.cov",
                   "cons.per.N",
                   "length.lowcov",
                   "bam.perc.100x",
                   "consensus.length",
                   "var.num.10.more",
                   "var.num.75.more",
                   "status")

final.table <- full.table %>% dplyr::select(final.columns) %>% rename_at(vars(final.columns), ~final.column.names)

write_csv(final.table, path = "report_metrics.csv")
write_tsv(final.table, path = "report_metrics.tsv")
###############################################################################
