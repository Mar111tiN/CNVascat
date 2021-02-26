# library(tidyverse)
# library(glue)

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS FOR ASCAT WORKFLOW
# ------------------------------------------------------------------------------
# 
# -----------------------------UTILS--------------------------------------------
get.path <- function(...) {
  # convenience function
  return(glue(..., .sep="/"))
}

# -----------------------------LOAD SAMPLE----------------------------------------
load.pos.file <- function(pos.file, min.reads=50) {
  # function for loading and adjusting the position files
  # uses tibble as dataframe and uses factor for chromosome
  
  # chrom_names for sorting
  chrom_names <- paste0("chr", c(1:22, "X", "Y"))
  
  # load into tibble
  sample.df <- read_delim(pos.file, delim = "\t") %>% 
    # turn chrom into factor
    mutate(chrom = factor(chrom, levels=chrom_names))
  
  filtered.df <- sample.df %>% 
    # Removal of insertion 
    distinct(chrom, pos, .keep_all=TRUE) %>%
    # Remove low quality SNPs
    filter(ref + alt > min.reads)
  # verbose
  print(glue("Loading ", pos.file, ": ", nrow(sample.df), " ->[total reads > ", min.reads, "]-> ", nrow(filtered.df)))
  return(sample.df)
}

# -----------------------------FILTER SNP-DENSE REGIONS----------------------------------------
# ------implemented function with flexible max.snp
filter.snp.window <- function(df, window.size=100, max.snp=2) {
  if (window.size==0) {
    return(df)
  }
  # create the indicator row as start of bad window
  df.filtered <- df %>% 
    # good =TRUE if distance to <max.snp>-ahead ( lead() ) position is less than window size
    # abs() is required for chrom gaps (creating large negative pos differences)
    mutate(good = abs(lead(pos, n=max.snp) - pos) > window.size, keep=good)
  # expand the keep from the <max.snp>-behind ( lag() ) positions
  for (i in seq(max.snp)) {
    df.filtered  <- df.filtered %>% mutate(keep=lag(good, n=i) & keep)
  }
  df.filtered <- df.filtered  %>% 
    # TRUE replace the first rows that have no lag
    replace_na(list(keep=TRUE)) %>% filter(keep) %>% 
    # remove out and bad columns
    select(-keep, -good) 
  print(glue("Filtered SNP-dense regions: ", nrow(df), "->[SNPs per ", window.size, " bases< ",max.snp," ]-> ", nrow(df.filtered)))
  return(df.filtered)
}

# ---------------------MERGE TUMOR AND NORMAL AND DERIVE COLUMNS---------
merge.TN <- function(tumor.df="", normal.df="", calculate.normal=FALSE) {
  
  # merge the two tibbles
  merged <- merge(tumor.df, normal.df, by = c("chrom", "pos"), 
                  suffixes = c("_tumor", "_normal")) %>% 
    # sort by chromosome
    arrange(chrom, pos) %>% 
    # add the derived columns
    mutate(
      total_t = alt_tumor + ref_tumor,
      total_n = alt_normal + ref_normal,
      baf_tumor = alt_tumor / total_t,
      # merged$baf_normal <- merged$alt_normal / (merged$alt_normal + merged$ref_normal)
      baf_normal = alt_normal / total_n
    ) 
  
  # options to include the normal logR by using the calculate.normal function argument
  if (calculate.normal) {
    merged <- merged %>%
      mutate(
        # THIS DOES NOT WORK!!!!!
        logr_tumor = log2(total_t/mean(total_t)),
        logr_normal = log2(total_n/mean(total_n))
      )
  } else {
    merged <- merged %>%
      mutate(
        # get the aggregates
        sum_t = sum(alt_tumor, na.rm = TRUE) + sum(ref_tumor, na.rm = TRUE),
        sum_n = sum(alt_normal, na.rm = TRUE) + sum(ref_normal, na.rm = TRUE),
        #  merged$logr_tumor <- log2(((merged$alt_tumor + merged$ref_tumor) / total_reads_tumor)
        # / ((merged$alt_normal + merged$ref_normal) / total_reads_normal))
        logr_tumor = log2(total_t/total_n) - log2(sum_t / sum_n),
        # logR of normal should be 0 in this case as normalization includes normal coverage
        logr_normal = 0
      ) %>%
      # remove the total reads columns
      select(-c(sum_t, sum_n))
  } 
  
  # go on with remaining calculations
  merged <- merged %>% 
    ## Remove NAs and +/-Inf
    filter_all(all_vars(is.finite(.))) %>% 
    # ASCAT only takes integer chromosome names
    mutate(chrom = as.numeric(str_replace_all(chrom, c(
      "chr" = "",
      "X" = "23",
      "Y" = "24"))))
  return(merged)
}


######### LOAD T & N ################
load.sample.pair <- function(
    tumor.pos=".",
    normal.pos=".", 
    min.reads=50,
    window.size=100, 
    max.snp=2,
    calculate.normal=FALSE # whether logR in normal should be used
) {
  # Load tumor 
  tumor <- load.pos.file(tumor.pos, min.reads=min.reads)
  # Load normal 
  normal <- load.pos.file(normal.pos) %>% 
    # Checks if in a window of 100bp are more than 2 SNPs and remove them
    filter.snp.window(window.size=window.size, max.snp=max.snp)
  # merge tumor and normal
  merged.df <- merge.TN(tumor.df=tumor, normal.df=normal, calculate.normal=calculate.normal)
  return(merged.df)
}


# ---------------------REFPHASE OUTPUT ---------
write.refphase <- function(df, sample.name, out.dir=".", cutoff=0.4) {
  out.file <- glue(out.dir, "/", sample.name, "_snps.tsv.gz")
  df %>% 
    # Set positions with germline BAF < 0.1 or BAF > 0.9 as germline hom
    mutate(germline_zygosity = ifelse(abs(df$baf_normal - 0.5) > cutoff, "hom", "het")) %>% 
    dplyr::rename(baf = baf_tumor, logr=logr_tumor) %>% 
    select(c(chrom,pos,baf,logr,germline_zygosity)) %>% 
    write_tsv(out.file)
  # return original file
  return(df)
}

# ---------------------ASCAT OUTPUT ---------
# Write the strange input format files that ASCAT expects
write.ascat.type <- function(df, sample.name, out.dir=".", type =".") {
  out.file <- get.path(out.dir, glue(sample.name, "_", type, ".tsv"))
  df %>%
    dplyr::rename(chrs=chrom, !!sample.name:= type) %>% 
    select(c(chrs, pos, !!sample.name)) %>% 
    # write to table
    write.table(
      file = out.file, 
      sep = "\t", 
      quote = FALSE, 
      row.names = TRUE, 
      col.names = NA
    )
}

write.ascat <- function(df, sample.name, out.dir=".") {
  # turn index into SNPindexer
  ascat.df <- df %>% 
    mutate(row.names = paste0("SNP",row_number())) %>% 
    column_to_rownames(var= "row.names")
  
  # cycle through the types and write the output file
  for (type in c("baf_tumor", "logr_tumor", "baf_normal", "logr_normal")) {
    write.ascat.type(ascat.df, sample.name=sample.name, out.dir=out.dir, type=type)
  }
  return(df)
}

# ---------------------PREPROCESSING AND PREPRO-OUTPUT---------
# load the pos files, do the computations and 
cnv.preprocess <- function(
  tumor.pos="",
  normal.pos="",
  sample.name="", 
  output.dir=".",
  min.reads=50,
  window.size=100,
  max.snp = 2,
  calculate.normal = FALSE
) {
  # loading the merged tumor normal file with all preprocessings
  merged.df <- load.sample.pair(
    tumor.pos=tumor.pos,
    normal.pos=normal.pos, 
    min.reads=min.reads,
    window.size=window.size,
    max.snp = max.snp,
    calculate.normal=calculate.normal
  ) %>% 
    write.refphase(sample.name, out.dir=output.dir) %>% 
    write.ascat(sample.name, out.dir=output.dir)
  return(merged.df)
}


# ---------------------LOAD FILES FOR ASCAT OBJECT ---------
ascat.load <- function(sample.name, input.dir=".", gender="XX") {
  # helper for loading the ascat file
  # make the filename base
  sample.file <- get.path(input.dir, sample.name)
  ascat.bc <- ascat.loadData(Tumor_LogR_file = glue(sample.file, "_logr_tumor.tsv"),
                             Tumor_BAF_file = glue(sample.file, "_baf_tumor.tsv"),
                             Germline_LogR_file = glue(sample.file, "_logr_normal.tsv"),
                             Germline_BAF_file = glue(sample.file, "_baf_normal.tsv"),
                             gender = gender) 
  return(ascat.bc)
}

# ------------------------------------------------------------------------------
# RUN ASCAT
# ------------------------------------------------------------------------------

run.ascat <- function(sample.name, sample.dir=".") {
  ascat.bc <- ascat.load(sample.name, input.dir=sample.dir)
  
  ascat.plotRawData(ascat.bc, img.dir=sample.dir)
  
  ascat.bc <- ascat.aspcf(ascat.bc, penalty = ascat.penalty, out.dir=sample.dir)
  
  ascat.plotSegmentedData(ascat.bc, img.dir=sample.dir)
  
  ascat.output <- ascat.runAscat(
    ascat.bc, 
    gamma = 1.0,    # actual log-diff between n=1 and n=2 --> 1 for NGS
    img.dir=sample.dir
  )
  return(ascat.output)
}

