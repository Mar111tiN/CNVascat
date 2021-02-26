library(tidyverse)
library(glue)
library(ASCAT)
library(gtools)
########### INSTALL REFPHASE ############################################
# load refphase code locally to path
# home <- "/Users/mahtin"
# home <- "/Users/martinscience"
library(devtools)
# refphase.path <- glue(home, "Dropbox/Icke/Work/somVar/CNV/ASCATrefphase/refphase", .sep="/")
# devtools::install(refphase.path)
devtools::install_bitbucket("schwarzlab/refphase", auth_user = "martin37szyska", password = "hJ9WmjSn6Fpzygw7MSSu")
library(refphase, warn.conflicts = FALSE)

### PATHS and HELPER FUNCTIONS

{
  home <- "/Users/mahtin"
  # home <- "/Users/martinscience"
  base.dir <- glue(home, "Dropbox/Icke/Work/somVar/CNV/ASCATrefphase", .sep="/")
  data.dir <- glue(base.dir, "data", .sep="/")
  input.dir <- glue(data.dir, "pos", .sep="/")
  R.dir <- glue(base.dir, "R", .sep="/")
  output.dir <- glue(data.dir, "results", .sep="/")
  list.files(path=input.dir)
  source(glue(R.dir, "cnv_helpers.R", .sep="/"))
}

### SNAKE INPUT #########
{
  # from wildcards
  sample = "67"
  tumor = "A"
  normal = "B"
  # ------------------------------------------------------------------------------
  # GLOBAL CONFIG
  # ------------------------------------------------------------------------------
  cutoff = 0.4
  min.reads = 50
  snp.window.size = 100
  max.snps.in.window = 2
  gender = "XX"
  ascat.penalty = 25
  gamma = 1.
  calculate.logR.4.normal = FALSE
}

{
  # create the sample directory
  sample.dir <- get.path(output.dir, sample)
  dir.create(sample.dir, showWarnings = FALSE)
  
  # get the sample names from wildcards
  t <- glue(sample, tumor, .sep="_")
  n <- glue(sample, normal, .sep="_")
  tn <- glue(t, normal, .sep="-")
}


merged.df <- cnv.preprocess(
  tumor=t,
  normal=n, 
  input.dir=input.dir,
  output.dir=sample.dir,
  min.reads=min.reads,
  window.size=snp.window.size,
  max.snp = max.snps.in.window,
  calculate.normal=calculate.logR.4.normal
)

# ------------------------------------------------------------------------------
# RUN ASCAT
# ------------------------------------------------------------------------------
ascat.output <- run.ascat(tn, sample.dir=sample.dir)


# ------------------------------------------------------------------------------
# RUN REFPHASE
# ------------------------------------------------------------------------------  

# write the segment file from ascat output
ascat.output$segments %>% 
  dplyr::rename(chrom=chr, start=startpos, end=endpos, cn_major=nMajor, cn_minor=nMinor) %>% 
  drop_na() %>% 
  select(-c(sample)) %>% 
  write_tsv(get.path(sample.dir, glue(tn, "_ascat_segs.tsv")))

(file.path <- get.path(sample.dir, tn))

rphase_sample_data <- data.frame(sample_id = tn, 
                                 segmentation = toString(glue(file.path, "_ascat_segs.tsv")), 
                                 snps = toString(glue(file.path, "_snps.tsv")), 
                                 purity = ascat.output$aberrantcellfraction[[1]],
                                 ploidy = ascat.output$ploidy[[1]]
                                 ) %>% 
  write_tsv(get.path(sample.dir, "rphase_sample_data.tsv"))

# this does not work on a tibble so this has to be re-read as a data.frame
sample_data <- read.delim(get.path(sample.dir, "rphase_sample_data.tsv"), stringsAsFactors = FALSE)

sample_data
# Pre-load all experimental data and sample metadata
exper <- refphase_load(sample_data, segs_format = "refphase", snps_format = "refphase")

exper <- center_baf(exper)

exper <- fit_logr_to_ascat(exper)

# Run rphase on the experimental data
results <- refphase(
  exper, cn_method = "ascat",
  merge_min_width = 1e+065,
  max_merge_width = 1e+06,
  gap_fill_width = 2e+06,
  verbosity = "INFO",
  phase_all_segs = TRUE
)

refphase:::write_segs(results$phased_segs, file = glue(file.path, "_refphase_segs.tsv"))

# (optional) output the SNPs, including phasing information
refphase:::write_snps(results$phased_snps, file = glue(file.path, "_refphase_phased_snps.tsv"))

# Ploidy might have changed, if we have updated any copy numbers
write.table(results$sample_data, file = glue(file.path, "_refphase-sample-data.tsv"), sep = "\t", row.names = FALSE)

# Plot the results as PDFs (these files can be huge, over 100 MB in some cases)
pdf(glue(file.path, "_refphase_genome.pdf"), width = 10, height = 4 * nrow(sample_data))
plot_genome(results$sample_data, results$phased_snps, results$phased_segs)
dev.off()

pdf(glue(file.path, "_refphase_chromosomes.pdf"), width = 10, height = 2 + 4 * nrow(sample_data))
plot_all_chromosomes(results$sample_data, results$phased_snps, list(refphase = results$phased_segs, orig = exper$segs))
dev.off()


