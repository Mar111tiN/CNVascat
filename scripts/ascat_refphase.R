# General refphase workflow for detecting SCNA in tumor samples.
# See https://bitbucket.org/schwarzlab/refphase/src/master/ for further details.


############# IMPORTS ######################
suppressMessages({
    library(glue)
    library(tidyverse, warn.conflicts = FALSE)
    library(gtools, warn.conflicts = FALSE)
    library(ASCAT, warn.conflicts = FALSE)
})

## Install refphase from R
# R
# > library(devtools)
# > devtools::install("./scripts/refphase")
# use tryCatch to check for installed library refphase
tryCatch(library(refphase, warn.conflicts = FALSE), 
    error = function(e){
        message("\"",e,"Refphase not installed --> installing now")
        # load necessary library devtools
        library(devtools)
        # install refphase from scripts/refphase
        devtools::install(glue("{snakemake@scriptdir}/refphase"))
        suppressMessages({
          library(refphase, warn.conflicts = FALSE)
        })  
    })


# load the helper functions
source(glue("{snakemake@scriptdir}/cnv_helpers.R"))


########## SNAKEMAKE GENERICS #####################
{
    i <- snakemake@input
    o <- snakemake@output
    w <- snakemake@wildcards
    p <- snakemake@params
    c.a <- snakemake@config[['ascat']]
    c.r <- snakemake@config[['refphase']]

    tumor <- w[['tumor']]
    normal <- w[['normal']]
    sample <- w[["sample"]]

    tumor.pos <- i[['tumor']]
    normal.pos <- i[['normal']]
}


########## PARAMS/CONFIG #########################
{
    cutoff <- c.a[['cutoff']]
    # gender; can be the same for first try (XY)
    min.reads <- c.a[['min_reads']]
    snp.window.size <- c.a[['snp_window_size']]
    max.snps.in.window <- c.a[['max_snps_in_window']]
    gender <- c.a[['gender']]
    ascat.penalty <- c.a[['ascat_penalty']]
    gamma <- c.a[['gamma']]
    calculate.normal <- c.a[['calculate_normal']]
    merge_min_width <- c.r[['merge_min_width']]
    max_merge_width <- c.r[['max_merge_width']]
    gap_fill_width <- c.r[['gap_fill_width']]
}


## get the file paths and names
{
  # create the sample directory
  sample.dir <- dirname(toString(o))

  # snakemake should take care of folder creation
  # dir.create(sample.dir, showWarnings = FALSE)
  
  # get the sample names from wildcards
  sample.name <- glue(sample, "_", tumor, "-", normal)
  # get the full path to sample base name
  # e.g. "cnv/67A/67_A-B"
  sample.file <- get.path(sample.dir, sample.name)
}

# ------------------------------------------------------------------------------
# RUN PREPROCESSING
# ------------------------------------------------------------------------------
merged.df <- cnv.preprocess(
  tumor.pos=tumor.pos,
  normal.pos=normal.pos,
  sample.name=sample.name,
  output.dir=sample.dir,
  min.reads=min.reads,
  window.size=snp.window.size,
  max.snp = max.snps.in.window,
  calculate.normal=calculate.normal
)

# ------------------------------------------------------------------------------
# RUN ASCAT
# ------------------------------------------------------------------------------
ascat.output <- run.ascat(sample.name, sample.dir=sample.dir)


# ------------------------------------------------------------------------------
# RUN REFPHASE
# ------------------------------------------------------------------------------  

# write the segment file from ascat output
ascat.output$segments %>% 
  dplyr::rename(chrom=chr, start=startpos, end=endpos, cn_major=nMajor, cn_minor=nMinor) %>% 
  drop_na() %>% 
  select(-c(sample)) %>% 
  write_tsv(glue(sample.file, "_ascat_segs.tsv.gz"))


rphase_sample_data <- data.frame(sample_id = sample.name, 
                                 segmentation = toString(glue(sample.file, "_ascat_segs.tsv.gz")), 
                                 snps = toString(glue(sample.file, "_snps.tsv.gz")), 
                                 purity = ascat.output$aberrantcellfraction[[1]],
                                 ploidy = ascat.output$ploidy[[1]]
                                 ) %>% 
  write_tsv(get.path(sample.dir, "rphase_sample_data.tsv"))

# this does not work on a tibble so this has to be re-read as a data.frame
sample_data <- read.delim(get.path(sample.dir, "rphase_sample_data.tsv"), stringsAsFactors = FALSE)

# Pre-load all experimental data and sample metadata
exper <- refphase_load(sample_data, segs_format = "refphase", snps_format = "refphase")

exper <- center_baf(exper)

exper <- fit_logr_to_ascat(exper)

# Run rphase on the experimental data
results <- refphase(
  exper, cn_method = "ascat",
  merge_min_width = merge_min_width,
  max_merge_width = max_merge_width,
  gap_fill_width = gap_fill_width,
  verbosity = "INFO",
  phase_all_segs = TRUE
)

refphase:::write_segs(results$phased_segs, file = glue(sample.file, "_refphase_segs.tsv.gz"))

# (optional) output the SNPs, including phasing information
refphase:::write_snps(results$phased_snps, file = glue(sample.file, "_refphase_phased_snps.tsv.gz"))

# Ploidy might have changed, if we have updated any copy numbers
write.table(results$sample_data, file = glue(sample.file, "_refphase-sample-data.tsv"), sep = "\t", row.names = FALSE)

# Plot the results as PDFs (these files can be huge, over 100 MB in some cases)
pdf(glue(sample.file, "_refphase_genome.pdf"), width = 10, height = 4 * nrow(sample_data))
plot_genome(results$sample_data, results$phased_snps, results$phased_segs)
dev.off()

pdf(glue(sample.file, "_refphase_chromosomes.pdf"), width = 10, height = 2 + 4 * nrow(sample_data))
plot_all_chromosomes(results$sample_data, results$phased_snps, list(refphase = results$phased_segs, orig = exper$segs))
dev.off()
