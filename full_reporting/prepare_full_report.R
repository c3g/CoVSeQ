library(tidyverse)
library(lubridate)
library(ggplot2)

##################################################################
# FUNCTIONS 
find_site <- function(sample.id){
  prefix <- str_split_fixed(sample.id, "-", 2)[1]
  if (prefix %in% LSPQ.sites$Prefix){
    site <- prefix
  }else if(str_detect(sample.id, "^L00")){
    site <- "LSPQ"
  }else{
    site <- NA
  }
  site
}

remove_v_from_id <- function(sample.id){
  id_sans_v <- stringr::str_replace(sample.id, "_\\d$", "" )
  id_sans_v
}


classify_sample <- function(sample.id){
  sample.type <- "other"
  upper.id <- str_to_upper(sample.id)
  
  for (site.substring in site.prefixes){
    if (str_detect(upper.id, site.substring)){
      sample.type <- "clinical"
    }
  }
  
  for (td.substring in techdev.substrings){
    if (str_detect(upper.id, td.substring)){
      sample.type <- "techdev"
    }
  }
    
  for (ctrl.substring in control.substrings){
    if (str_detect(upper.id, ctrl.substring)){
        sample.type <- "control"
      }
  }
  sample.type   
}

extract_id_from_header <- function(gisaid.header){
  detected.id <- stringr::str_replace(gisaid.header, regex("/202\\d$", ignore_case = T), "")
  detected.id <- stringr::str_replace(detected.id, regex("^HCOV-19/CANADA/QC-LSPQ-", ignore_case = T), "")
  detected.id <- stringr::str_replace(detected.id, regex("^HCOV-19/CANADA/QC_", ignore_case = T), "")
  detected.id <- stringr::str_replace(detected.id, regex("^HCOV-19/CANADA/QC-", ignore_case = T), "")
  detected.id <- stringr::str_replace(detected.id, regex("^CANADA/QC-", ignore_case = T), "")
  detected.id 
}

assign_platform <- function(lims.run.name){
  platform <- NA
  if (is.na(lims.run.name)){
    platform <- NA
  }
  else if (stringr::str_detect(lims.run.name, regex("illumina", ignore_case = T))){
    platform <- "illumina"
  }
  else if (stringr::str_detect(lims.run.name, regex("nanopore", ignore_case = T))){
    platform <- "nanopore"
  }
  else if (stringr::str_detect(lims.run.name, regex("mgi", ignore_case = T))){
    platform <- "mgi"
  }
  platform
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
    "REJ"
  }
  else {
    "NA"
  }
}


neg_ctrl_check <- function(plate.id){
  ncov_tools.neg.ctrl.warn = "PASS" 
  if (plate.id %in% problem.negs){
    ncov_tools.neg.ctrl.warn = "FLAG"  
  }
  ncov_tools.neg.ctrl.warn
}


flag_neg_ctrl <- function(pass.flag.rej, neg.ctrl.flag){
  if (is.na(neg.ctrl.flag) | is.na(pass.flag.rej)){
    pass.flag.rej = pass.flag.rej
  }
  else if (neg.ctrl.flag == "FLAG" & pass.flag.rej == "PASS"){
    pass.flag.rej = "FLAG"
  }
  pass.flag.rej
}

flag_frameshift <- function(pass.flag.rej, frameshift.detect){
  if (is.na(frameshift.detect) | is.na(pass.flag.rej)){
    pass.flag.rej = pass.flag.rej
  }
  else if (frameshift.detect == TRUE & pass.flag.rej == "PASS"){
    pass.flag.rej = "FLAG"
  }
  pass.flag.rej
}

flag_low_cov <- function(pass.flag.rej, length.low.cov){
  if (is.na(length.low.cov) | is.na(pass.flag.rej)){
    pass.flag.rej = pass.flag.rej
  }
  else if (pass.flag.rej == "PASS" & length.low.cov > 300){
    pass.flag.rej = "FLAG"
  }
  pass.flag.rej
}


##################################################################
# PROLOGUE
## Parse list of sites
LSPQ.sites <- readr::read_csv(file.path("LSPQ_metadata", "ListeCH_prefix.csv"))

## Prepare sample classifier lists
techdev.substrings <- c("POOL", "SPIKEIN", "SPIKE", "CANCOV", "5ULRT")
control.substrings <- c("CTRL", "CTL", "EXT", "ACCUKIT", "POS", "NEG", "NTC", "PTC")
site.prefixes <- LSPQ.sites %>% pull(Prefix) 
site.prefixes <- c(site.prefixes, "L00")

## Plates with negative control issues (manually curated based on reports)
problem.negs <- c("GQ_013", "GQ_027", "GQ_034", "GQ_042", "GQ_062",
                  "GQ_078", "GQ_085", "GQ_087", "GQ_095", "GQ_096", 
                  "GQ_104", "GQ_105")

##################################################################
# BUILD INDIVIDUAL REPORT TABLES
## By decreasing order, from most, to least samples... 
## Logic is that the first table will contain *all* samples, 
## then progressively each additional table will contain subsets of that totality

##################################################################
# TABLE 1: Create a list of all samples and all runs for all platforms
## Import all nanopore run info
nanopore.runs <- readr::read_csv(file.path("nanopore_reports", "nanopore_samples.csv"), 
                                 col_names = c("platform", "run", "analysis", "sample"))
nanopore.runs <- nanopore.runs %>% mutate(alias.id = sample, sample = remove_v_from_id(sample))
nanopore.runs <- nanopore.runs %>% distinct() # Remove any duplicates
## Import all illumina run info
illumina.runs <- readr::read_csv(file.path("illumina_reports", "illumina_samples.csv"), 
                                 col_names = c("platform", "run", "sample")) %>% 
  mutate(analysis = run)
illumina.runs <- illumina.runs %>% distinct() # Remove any duplicates
### Import all MGI run info
mgi.runs <- readr::read_csv(file.path("mgi_reports", "mgi_samples.csv"), 
                            col_names = c("platform", "run", "sample")) %>% 
  mutate(analysis = run)
mgi.runs <- mgi.runs %>% distinct() # Remove any duplicates

### Merge all run info
run.list <- dplyr::bind_rows(nanopore.runs, illumina.runs)
run.list <- dplyr::bind_rows(run.list, mgi.runs)

### Parse run type (From a manually curated list)
run.type <- readr::read_csv(file.path("run_lists","run_type.dates.csv"), 
                            col_names = c("run", "run.type", "run.start.date"), 
                            col_types = cols(
                              run = col_character(),
                              run.type = col_character(),
                              run.start.date = col_date(format = "%Y%m%d")))

run.type <- run.type%>% mutate(run.start.month = lubridate::month(run.start.date))

### Add run type information
run.list <- dplyr::left_join(run.list, run.type)

### Add sample type information
run.list <- run.list %>% mutate(sample.type = map_chr(run.list %>% pull(sample), classify_sample))
run.list <- run.list %>% mutate(site = map_chr(run.list %>% pull(sample), find_site)) %>% 
  select(sample, platform, sample.type, 
         site, run, run.type, run.start.date, run.start.month, 
         analysis, alias.id)

### Segregate by Production/Non-production
production.list <- run.list %>% filter(run.type == "production")
non.production.list <- run.list %>% filter(run.type != "production")

##################################################################
# TABLE 2: LIMS table (should include failed samples that don't get sequenced)
# Import latest LIMS report (sent by Haig) and parsed by my scripts
LIMS.latest.report.path <- file.path("LIMS_reports", "latest_LIMS_full_report.csv")
LIMS.report <- readr::read_csv(LIMS.latest.report.path) 
LIMS.report <- LIMS.report

LIMS.production.samples <- LIMS.report %>% 
  filter(seq.stage == "Illumina Sequencing (NovaSeq) 1.0 McGill 1.0")
LIMS.non.production.samples <- LIMS.report %>% 
  filter(seq.stage != "Illumina Sequencing (NovaSeq) 1.0 McGill 1.0" | is.na(seq.stage)) 

##################################################################
# TABLE 3: LSPQ metadata tables (shared by Eric)
# # Import latest LSPQ metadata table (with Ct values and collection dates)
# LSPQ.latest.full.metadata.path <- file.path("LSPQ_metadata", "latest_LSPQmetadata.tsv")
# LSPQ.full.metadata <- readr::read_tsv(LSPQ.latest.full.metadata.path, col_types = cols(
#     DTNAISS = col_date(format = ""), SEXE = col_factor(), BIOBANK_ID = col_character(),
#     COVBANK_ID = col_character(), SAMPLE_DATE = col_date(format = ""), CT = col_double() )) %>% 
#   rename(birth.date = DTNAISS, sex = SEXE, biobank.id = BIOBANK_ID, 
#     covbank.id = COVBANK_ID, sampling.date = SAMPLE_DATE) %>% 
#   mutate(sample = biobank.id) %>%
#   select(sample, biobank.id, covbank.id, sampling.date, CT, sex) 
# 
# # Import latest MINIMAL LSPQ metadata table (only includes sampling date)
# LSPQ.latest.minimal.metadata.path <- file.path("LSPQ_metadata", "latest_minimal_LSPQmetadata.tsv")
# LSPQ.minimal.metadata <- readr::read_tsv(LSPQ.latest.minimal.metadata.path) %>% 
#   rename(biobank.id = SAMPLE, sampling.date = DT_PRELEV) %>% 
#   mutate(sample = biobank.id)
# LSPQ.minimal.metadata <- LSPQ.minimal.metadata %>% distinct()

# Import latest MINIMAL LSPQ metadata "DUMP" table (includes 3 sample IDs and sampling date)
LSPQ.latest.dump.metadata.path <- file.path("LSPQ_metadata", "latest_minimal_LSPQmetadata_dump.tsv")
LSPQ.minimal.metadata.dump <- readr::read_tsv(LSPQ.latest.dump.metadata.path, 
                                              na = c("", "NA", "not available"),
                                              col_types = cols(
                                                SEXE = col_factor(),
                                                SAMPLING_DATE = col_date(),
                                                DTNAISS = col_date())) %>% 
  rename(sample.id.1 = SAMPLE_ID_1, 
         sample.id.2 = SAMPLE_ID_2, 
         sample.id.3 = SAMPLE_ID_3, 
         sampling.date = SAMPLING_DATE, 
         sex = SEXE, 
         birth.date = DTNAISS, 
         lineage = LIGNEE_PANGOLIN) %>% 
  mutate(sample = sample.id.1) %>% 
  mutate(age = lubridate::dseconds(lubridate::int_length(birth.date %--% sampling.date))) %>%
  select(sample, sampling.date, sex, age)
LSPQ.minimal.metadata.dump <- LSPQ.minimal.metadata.dump %>% distinct()


# Import latest Pangolin lineages provided by LSPQ
LSPQ.latest.pangolin.path <- file.path("LSPQ_metadata", "latest_Pangolin_Report.csv")
LSPQ.latest.pangolin <- readr::read_csv(LSPQ.latest.pangolin.path, na = c("", "NA"),
                                        col_types = cols(
                                          note = col_character())
                                        ) %>% 
  rename(fasta.header = taxon, 
         pangoLEARN.version = pangoLEARN_version) %>% 
  mutate(sample = extract_id_from_header(fasta.header)) %>% 
  select(sample, fasta.header, lineage, probability, pangoLEARN.version, status)


# Join LSPQ metadata tables
LSPQ.metadata <- dplyr::left_join(LSPQ.minimal.metadata.dump, LSPQ.latest.pangolin)
LSPQ.metadata <- LSPQ.metadata %>% distinct() # Remove duplicates

##################################################################
# TABLE 4: Illumina internal metrics and ncov_tools metrics
# Import all illumina internal metrics report files
illumina.metrics.path.list <- list.files(path = "illumina_reports/internal_metrics",
                                       pattern = "*report_metrics.csv",
                                       full.names = T)
illumina.metrics <- illumina.metrics.path.list %>%
  map_df(~ readr::read_csv(.) %>% 
           mutate(filename = .x )) %>% 
  mutate(platform = "illumina") 
illumina.metrics <- illumina.metrics %>% 
  mutate(metrics.run = str_split(filename, "/", simplify = T)[,3]) %>%
  mutate(metrics.run = str_replace(metrics.run, "_report_metrics.csv", "")) %>% 
           select(-filename)
illumina.metrics <- illumina.metrics %>% distinct() # Remove duplicates

# Import all ncov_tools metrics files
illumina.ncov.tools.path.list <- list.files(path = "illumina_reports/ncov_tools",
                                         pattern = "*_summary_qc.tsv",
                                         full.names = T, 
                                         recursive = T)
illumina.ncov.tools <- illumina.ncov.tools.path.list %>%
  map_df(~ readr::read_tsv(., col_types = cols(run_name = col_character())) %>% 
           mutate(filename = .x)) %>% 
  mutate(ncov_tools.pass = str_detect(qc_pass, "PASS"), 
         ncov_tools.partial = str_detect(qc_pass, "PARTIAL_GENOME"), 
         ncov_tools.incomplete = str_detect(qc_pass, "INCOMPLETE_GENOME"),
         ncov_tools.frameshift = str_detect(qc_pass, "POSSIBLE_FRAMESHIFT_INDELS"), 
         ncov_tools.ambiguity = str_detect(qc_pass, "EXCESS_AMBIGUITY"),
         platform = "illumina")
illumina.ncov.tools <- illumina.ncov.tools %>% 
  mutate(ncov.tools.run = str_split(filename, "/", simplify = T)[,3]) %>%
  mutate(ncov.tools.run = str_replace(ncov.tools.run, "_summary_qc.tsv", "")) %>% 
  select(-filename)
illumina.ncov.tools <- illumina.ncov.tools %>% distinct() # Remove duplicates

# # Import all ncov_tools lineage files
# illumina.ncov.lineage.path.list <- list.files(path = "illumina_reports/ncov_tools",
#                                             pattern = "*_lineage_report.csv",
#                                             full.names = T, 
#                                             recursive = T)
# illumina.ncov.lineage <- illumina.ncov.lineage.path.list %>%
#   map_df(~ readr::read_csv(.) %>% 
#            mutate(filename = .x)) %>% 
#   rename(fasta.header = taxon, preelim.lineage = lineage, lineage.prob = probability, version.pangoLEARN = pangoLEARN_version) %>%
#   mutate(platform = "illumina", sample = extract_id_from_header(fasta.header)) %>% 
#   select(sample, fasta.header, preelim.lineage, lineage.prob, version.pangoLEARN, platform, filename)
# illumina.ncov.lineage <- illumina.ncov.lineage %>% 
#   mutate(lineage.run = str_split(filename, "/", simplify = T)[,3]) %>%
#   mutate(lineage.run = str_replace(lineage.run, "_lineage_report.csv", "")) %>% 
#   select(-filename)

illumina.full.metrics <- dplyr::full_join(illumina.metrics, illumina.ncov.tools, 
                                          by = c("Sample" = "sample", 
                                                 "platform" = "platform",
                                                 "metrics.run" = "ncov.tools.run"))
# illumina.full.metrics <- dplyr::left_join(illumina.full.metrics, illumina.ncov.lineage, 
#                                           by = c("Sample" = "sample", 
#                                                  "platform" = "platform",
#                                                  "metrics.run" = "lineage.run"))


##################################################################
# TABLE 5: Nanopore internal metrics and ncov_tools metrics
# Import all nanopore metrics reports 
nanopore.metrics <- readr::read_csv(file.path("nanopore_reports", "internal_metrics", "full_nanopore_report.csv")) 
nanopore.metrics <- nanopore.metrics %>% 
  mutate(alias.id = sample, 
         sample = remove_v_from_id(sample), 
         platform = "nanopore", 
         `Nb reads` = NA, 
         `Percent human reads` = NA, 
         `Nb clean reads` = NA, 
         `Length low cov region (<20X)` = NA, 
         `Nb variants > 10 perc allele freq` = NA, 
         `Nb variants > 75 perc allele freq` = NA, 
         `PASS/FLAG/REJ` = map_chr(nanopore.metrics %>% pull(cons.perc.N), sample_status)) %>% 
  rename(Sample = sample, 
         `Mean coverage` = bam.mean.cov, 
         `Percent N` = cons.perc.N, 
         `Percent consensus > 100X` = bam.perc.100x, 
         `Length consensus` = cons.len) %>% 
  select(Sample, `Nb reads`, `Percent human reads`, `Nb clean reads`, `Mean coverage`, `Percent N`,
         `Length low cov region (<20X)`, `Percent consensus > 100X`, `Length consensus`,
         `Nb variants > 10 perc allele freq`, `Nb variants > 75 perc allele freq`, `PASS/FLAG/REJ`, filename, platform
         )
nanopore.metrics <- nanopore.metrics %>% 
  mutate(metrics.run = str_split(filename, "/", simplify = T)[,5]) %>%
  select(-filename)
nanopore.metrics <- nanopore.metrics %>% distinct() # Remove duplicates

# Import all ncov_tools metrics files
nanopore.ncov.tools.path.list <- list.files(path = "nanopore_reports/ncov_tools", 
                                       pattern = "*_summary_qc.tsv", 
                                       full.names = T, 
                                       recursive = T)
nanopore.ncov.tools <- nanopore.ncov.tools.path.list %>% 
  map_df(~ readr::read_tsv(., col_types = cols(run_name = col_character())) %>%
           mutate(filename = .x)) %>%
  mutate(ncov_tools.pass = str_detect(qc_pass, "PASS"), 
         ncov_tools.partial = str_detect(qc_pass, "PARTIAL_GENOME"), 
         ncov_tools.incomplete = str_detect(qc_pass, "INCOMPLETE_GENOME"),
         ncov_tools.frameshift = str_detect(qc_pass, "POSSIBLE_FRAMESHIFT_INDELS"), 
         ncov_tools.ambiguity = str_detect(qc_pass, "EXCESS_AMBIGUITY"),
         platform = "nanopore")
nanopore.ncov.tools <- nanopore.ncov.tools %>% 
  mutate(ncov.tools.run = str_split(filename, "/", simplify = T)[,3]) %>%
  mutate(ncov.tools.run = str_replace(ncov.tools.run, "_summary_qc.tsv", "")) %>% 
  select(-filename)
nanopore.ncov.tools <- nanopore.ncov.tools %>% distinct() # Remove duplicates

# # Import all ncov_tools lineage files
# nanopore.ncov.lineage.path.list <- list.files(path = "nanopore_reports/ncov_tools",
#                                               pattern = "*_lineage_report.csv",
#                                               full.names = T, 
#                                               recursive = T)
# nanopore.ncov.lineage <- nanopore.ncov.lineage.path.list %>%
#   map_df(~ readr::read_csv(.) %>% 
#            mutate(filename = .x )) %>% 
#   rename(fasta.header = taxon, preelim.lineage = lineage, lineage.prob = probability, version.pangoLEARN = pangoLEARN_version) %>%
#   mutate(platform = "nanopore", sample = extract_id_from_header(fasta.header)) %>% 
#   select(sample, fasta.header, preelim.lineage, lineage.prob, version.pangoLEARN, platform, filename)
# nanopore.ncov.lineage <- nanopore.ncov.lineage %>% 
#   mutate(lineage.run = str_split(filename, "/", simplify = T)[,3]) %>%
#   mutate(lineage.run = str_replace(lineage.run, "_lineage_report.csv", "")) %>% 
#   select(-filename)

nanopore.full.metrics <- dplyr::full_join(nanopore.metrics, nanopore.ncov.tools, 
                                          by = c("Sample" = "sample", 
                                                 "platform" = "platform"))
# nanopore.full.metrics <- dplyr::left_join(nanopore.full.metrics, nanopore.ncov.lineage, 
#                                           by = c("Sample" = "sample", 
#                                                  "platform" = "platform",
#                                                  "ncov.tools.run" = "lineage.run"))

##################################################################
# TABLE 6: MGI internal metrics 
# Import all MGI metrics reports 
mgi.metrics.path.list <- list.files(path = "mgi_reports/internal_metrics",
                                         pattern = "*report_metrics.csv",
                                         full.names = T)
mgi.metrics <- mgi.metrics.path.list %>%
  map_df(~ readr::read_csv(.) %>% 
           mutate(filename = .x)) %>% 
  mutate(platform = "mgi")
mgi.metrics <- mgi.metrics %>% 
  mutate(metrics.run = str_split(filename, "/", simplify = T)[,3]) %>%
  mutate(metrics.run = str_replace(metrics.run, "_report_metrics.csv", "")) %>% 
  select(-filename)
mgi.metrics <- mgi.metrics %>% distinct() # Remove duplicates

##################################################################
# TABLE 7: All joined metrics
full.metrics <- bind_rows(illumina.full.metrics, mgi.metrics)
full.metrics <- bind_rows(full.metrics, nanopore.full.metrics)

##################################################################
# TABLE 9: GISAID release information
GISAID.release.info <- readr::read_csv(file.path("release_lists", "latest_release_list.csv"), 
                                       col_names = c("gisaid.header", "gisaid.accession", "gisaid.sampling.date"))
GISAID.release.info <- GISAID.release.info %>% 
                        mutate(release = T, sample = extract_id_from_header(gisaid.header)) %>% 
                        select(sample, release, gisaid.header, gisaid.accession, gisaid.sampling.date)


##################################################################
# JOIN TABLES IN ORDER

# Join Production portion of report
final.production.report <- dplyr::left_join(production.list, LIMS.production.samples, 
                                          by = c("sample" = "sample"))

# Join non-production portion of report 
final.non.production.report <- dplyr::left_join(non.production.list, LIMS.non.production.samples, 
                                          by = c("sample" = "sample"))

# Join both
final.LSPQ.report <- dplyr::bind_rows(final.production.report, final.non.production.report)

# Add samples with no run information
final.no.run.report <- anti_join(LIMS.non.production.samples, final.LSPQ.report, 
                                 by = c("sample" = "sample")) %>% 
                        filter(!is.na(PCR.date)) # Remove duplicates by only keeping samples with PCR dates
final.LSPQ.report <- dplyr::bind_rows(final.LSPQ.report, final.no.run.report)

## Now, we progressively add all the extra information, starting with all metrics 
final.LSPQ.report <- dplyr::left_join(final.LSPQ.report, full.metrics, 
                                      by = c("sample" ="Sample", 
                                             "run" = "metrics.run"))

# Add LSPQ metadata information 
final.LSPQ.report <- dplyr::left_join(final.LSPQ.report, LSPQ.metadata, 
                                      by = c("sample" = "sample"))

# Finally add GISAID Release information
final.LSPQ.report <- dplyr::left_join(final.LSPQ.report, GISAID.release.info, 
                                      by = c("sample" = "sample"))

final.LSPQ.report <- final.LSPQ.report %>% 
  rename(`Sample Name` = sample, `Plate ID` = lspq.plate.id, preelim.QC = `PASS/FLAG/REJ`) %>% 
  mutate(ncov_tools.neg_ctrl = map_chr(final.LSPQ.report %>% pull(lspq.plate.id), neg_ctrl_check), 
         `PASS/FLAG/REJ` = map2_chr(preelim.QC, ncov_tools.frameshift, flag_frameshift),
         `PASS/FLAG/REJ` = map2_chr(`PASS/FLAG/REJ`, ncov_tools.neg_ctrl, flag_neg_ctrl), 
         `PASS/FLAG/REJ` = map2_chr(`PASS/FLAG/REJ`, `Length low cov region (<20X)`, flag_low_cov), 
         `PASS/FLAG/REJ` = factor(`PASS/FLAG/REJ`, levels = c("PASS", "FLAG", "REJ"),
                                  labels = c("Pass", "Flag", "Reject"), ordered = T)) 

final.LSPQ.report <- final.LSPQ.report %>% distinct()

##################################################################
# EXPORT TABLES

# Export full report
write_csv(x = final.LSPQ.report, 
          file = file.path("report_archive", paste0(format(Sys.Date(), "%Y-%m-%d"), "_LSPQReport.csv")), 
          na = "NA", 
          append = F)
write_tsv(x = final.LSPQ.report, 
          file = file.path("report_archive", paste0(format(Sys.Date(), "%Y-%m-%d"), "_LSPQReport.tsv")), 
          na = "NA", 
          append = F)

# ## Export Pilot/TechDev portion
# report.pilot.techdev <- final.LSPQ.report %>% 
#   filter(run.type %in% c("pilot", "techdev"))
# write_csv(x = report.pilot.techdev, 
#           file = file.path("report_archive", paste0(format(Sys.Date(), "%Y-%m-%d"), "_LSPQ_Pilot_Report.csv")), 
#           na = "NA", 
#           append = F)
# write_tsv(x = report.pilot.techdev, 
#           file = file.path("report_archive", paste0(format(Sys.Date(), "%Y-%m-%d"), "_LSPQ_Pilot_Report.tsv")), 
#           na = "NA", 
#           append = F)
# 
# ## Export Production/Outbreak portion
# report.prod.outbreak <- final.LSPQ.report %>% 
#   filter(run.type %in% c("production", "outbreak"))
# write_csv(x = report.prod.outbreak, 
#           file = file.path("report_archive", paste0(format(Sys.Date(), "%Y-%m-%d"), "_LSPQ_Production_Report.csv")), 
#           na = "NA", 
#           append = F)
# write_tsv(x = report.prod.outbreak, 
#           file = file.path("report_archive", paste0(format(Sys.Date(), "%Y-%m-%d"), "_LSPQ_Production_Report.tsv")), 
#           na = "NA", 
#           append = F)
# 
# ## Not-sequenced portion
# report.na.run.type <- final.LSPQ.report %>% 
#   filter(is.na(run.type))
# write_csv(x = report.na.run.type, 
#           file = file.path("report_archive", paste0(format(Sys.Date(), "%Y-%m-%d"), "_LSPQ_Unsequenced_Report.csv")), 
#           na = "NA", 
#           append = F)
# write_tsv(x = report.na.run.type, 
#           file = file.path("report_archive", paste0(format(Sys.Date(), "%Y-%m-%d"), "_LSPQ_Unsequenced_Report.tsv")), 
#           na = "NA", 
#           append = F)

