# ==============================================================================
# rPhase_GeneralWorkflow.R
#
# Author: Raphael Hablesreiter (raphael.hablesreiter@charite.de)
#
# Adapted from Matt Huska (AG Schwarz - MDC Berlin)
#
# Description:
# General rPhase workflow for detecting SCNA in tumor samples.
# See https://bitbucket.org/schwarzlab/rphase/src/master/ for further details.
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------

# library(devtools)
# devtools::install_bitbucket("schwarzlab/rphase")

library(readr)
library(gtools)
library(ASCAT)
library(refphase)
library(dplyr)

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

chrom2integer <- function(x) {
    print(unique(x))
    as.integer(gsub("M", "25", 
                    gsub("X", "23", 
                         gsub("Y", "24", gsub("chr", "", x)))))
}


'%!in%' <- function(x,y)!('%in%'(x,y))

# ------------------------------------------------------------------------------
# Global parameters
# ------------------------------------------------------------------------------

work_dir = "/Users/raphaelhablesreiter/Projects/2020_CBF/20200916/result"
file_dir = "/Users/raphaelhablesreiter/Projects/2020_CBF/20200916/input"

# bed_file = "/Users/raphaelhablesreiter/Projects/2020_Test_rPhase/SS_HAEv6r2_hg38_Covered_nochr.bed"
# bed_file <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE) 
# bed_file <- bed_file[,c(1:3)]
# colnames(bed_file) <- c("chr", "start", "end")

patients <- unlist(strsplit(list.files(path = file_dir), "-"))
patients <- unique(patients[grep("AML", patients)])

cutoff = 0.4
minreads_normal = 50
minreads_tumor = 50
gender = "XX"
patients = c("0604_135")
pat_tmp = patients
for (pat_tmp in patients)
{
    # ------------------------------------------------------------------------------
    # Global settings
    # ------------------------------------------------------------------------------
    patient = pat_tmp
    normal_sample = paste0(pat_tmp, "-CR1")
    tumor_samples <- c(paste0(pat_tmp, "-D"), paste0(pat_tmp, "-Rel1"))
    
    
    # tumor_samples <- c(paste0(pat_tmp, "_Rel1"))
        
    dir.create(paste0(work_dir, "/", patient))
    setwd(paste0(work_dir, "/", patient))
    
    # ------------------------------------------------------------------------------
    # Preprocessing
    # ------------------------------------------------------------------------------
    
    column_names <- c("chrom", "pos", "ref", "alt")
    chrom_names <- c(1:22, "X", "Y")
    normal <- as.data.frame(read.csv(paste0(file_dir, "/", normal_sample, ".pos.gz"), 
                                     sep = "\t", header = FALSE, 
                                     stringsAsFactors = FALSE))
    colnames(normal) <- column_names
    
    # Removal of insertion 
    no_snp <- paste0(normal$chrom, "_", normal$pos)
    normal <- normal[no_snp %!in% no_snp[duplicated(no_snp)],]

    # Remove low quality SNPs
    normal <- normal[normal$ref + normal$alt > minreads_normal,]
    
    # Checks if in a window of 100bp are more than 2 SNPs and removes them
    filter <- as.data.frame(matrix(ncol = 2))
    colnames(filter) <- c("chrom", "pos")
    for (chrom in unique(normal$chrom))
    {
        print(chrom)
        tmp <- normal[normal$chrom %in% chrom,2]
        call_start <- NA
        call_end <- NA
        
        for (i in 1:length(tmp)){
            call_start <- c(call_start, seq(-99,0) + tmp[i])
            call_end <- c(call_end, seq(0,99) + tmp[i])
        }
        
        windows <- as.data.frame(unique(list(call_start, call_end)))
        windows <- windows[-1,]
        tmp_out <- ""
        
        for (i in 1:nrow(windows)){
            tmp_snps <- tmp[tmp > windows[i,1] & tmp < windows[i,2]]
            if (length(tmp_snps) > 2){
                tmp_out <- c(tmp_out, tmp_snps)
            }
        }
       
        tmp_out <- unique(tmp_out)
        tmp_chrom <- rep(chrom, length(tmp_out))
        
        tmp_filter <- as.data.frame(cbind(tmp_chrom, tmp_out))
        colnames(tmp_filter) <-c("chrom", "pos")
        
        filter <- rbind(filter, tmp_filter)
    }
    
    normal <- normal[paste0(normal$chrom,normal$pos) %!in% paste0(filter$chrom,filter$pos),]
    

    
    for (tumor_sample in tumor_samples) {
        
        tumor <- as.data.frame(read.csv(paste0(file_dir, "/", tumor_sample, ".pos.gz"), 
                                        sep = "\t", header = FALSE,
                                        stringsAsFactors = FALSE))
        colnames(tumor) <- column_names
        
        no_snp <- paste0(tumor$chrom, "_", tumor$pos)
        tumor <- tumor[no_snp %!in% no_snp[duplicated(no_snp)],]
        
        # Filter low quality SNPs
        tumor <- tumor[tumor$ref + tumor$alt > minreads_tumor,]
        
        
        merged <- merge(normal, tumor, by = c("chrom", "pos"), 
                        suffixes = c("_normal", "_tumor"))
        
        merged$pos <- as.integer(merged$pos)
        merged <- merged[order(chrom2integer(merged$chrom), merged$pos), ]
        
        total_reads_tumor <- sum(merged$alt_tumor, na.rm = TRUE) + 
            sum(merged$ref_tumor, na.rm = TRUE)
        total_reads_normal <- sum(merged$alt_normal, na.rm = TRUE) + 
            sum(merged$ref_normal, na.rm = TRUE)
        
        merged$baf_tumor <- merged$alt_tumor / (merged$alt_tumor + merged$ref_tumor)
        merged$logr_tumor <- log2(((merged$alt_tumor + merged$ref_tumor) / total_reads_tumor) / ((merged$alt_normal + merged$ref_normal) / total_reads_normal))
        
        merged$baf_normal <- merged$alt_normal / 
            (merged$alt_normal + merged$ref_normal)
        
        # logR of normal should always be 0
        merged$logr_normal <- merged$logr_tumor
        merged$logr_normal[] <- 0
        
        ## Remove NAs and +/-Inf
        keep_numeric <- apply(merged[, c("logr_tumor", "logr_normal", 
                                         "baf_normal", "baf_tumor")], 1, 
                              function(x) all(unlist(Map(is.finite, x))))
        merged <- merged[keep_numeric, ]
        
        # ASCAT only takes integer chromosome names
        merged$chrom <- gsub("Y", "24",
                             gsub("X", "23",
                                  gsub("chr", "", merged$chrom)))        
        
        # Write the specific tsv format that we can later use with rphase
        # Set positions with germline BAF < 0.1 or BAF > 0.9 as germline homozygous,
        # and the rest as heterozygous
        write.table(data.frame(chrom = merged$chrom, pos = merged$pos, 
                               baf = merged$baf_tumor, logr = merged$logr_tumor, 
                               germline_zygosity = 
                                   ifelse(abs(merged$baf_normal - 0.5) > cutoff, 
                                          "hom", "het")), 
                    file = paste0(tumor_sample, "_rphase_snps.tsv"), sep = "\t", 
                    quote = FALSE, row.names = FALSE)
        

        
        # Write the strange input format files that ASCAT expects
        rownames(merged) <- paste0("SNP", seq_len(nrow(merged)))
        write.table(data.frame(chrs = merged$chrom, pos = merged$pos, 
                               sample = merged$baf_tumor, 
                               row.names = rownames(merged)), 
                    file = paste0(tumor_sample, "_baf_tumor.tsv"), sep = "\t", 
                    quote = FALSE, row.names = TRUE, col.names = NA)
        
        # Tumor
        dat <- data.frame(chrs = merged$chrom, pos = merged$pos, 
                          sample = merged$baf_tumor, row.names = rownames(merged))
        colnames(dat)[[3]] <- tumor_sample
        
        write.table(dat, file = paste0(tumor_sample, "_baf_tumor.tsv"), sep = "\t", 
                    quote = FALSE, row.names = TRUE, col.names = NA)
        
        dat[[tumor_sample]] <- merged$logr_tumor
        write.table(dat, file = paste0(tumor_sample, "_logr_tumor.tsv"), sep = "\t", 
                    quote = FALSE, row.names = TRUE, col.names = NA)
        
        # Normal
        dat[[tumor_sample]] <- merged$baf_normal
        write.table(dat, file = paste0(tumor_sample, "_baf_normal.tsv"), sep = "\t", 
                    quote = FALSE, row.names = TRUE, col.names = NA)
        
        dat[[tumor_sample]] <- merged$logr_normal
        write.table(dat, file = paste0(tumor_sample, "_logr_normal.tsv"), sep = "\t", 
                    quote = FALSE, row.names = TRUE, col.names = NA)
        
    }
    
    rphase_sample_data <- data.frame(sample_id = tumor_samples, 
                                     segmentation = paste0(tumor_samples, "_rphase_segs.tsv"), 
                                     snps = paste0(tumor_samples, "_rphase_snps.tsv"), 
                                     purity = NA, ploidy = NA, row.names = tumor_samples)
    ascat_penalty <- c(25,25)
    for (tumor_sample in tumor_samples) {
        ascat_bc <- ascat.loadData(Tumor_LogR_file = paste0(tumor_sample, "_logr_tumor.tsv"),
                                   Tumor_BAF_file = paste0(tumor_sample, "_baf_tumor.tsv"),
                                   Germline_LogR_file = paste0(tumor_sample, "_logr_normal.tsv"),
                                   Germline_BAF_file = paste0(tumor_sample, "_baf_normal.tsv"),
                                   gender = gender) 
        # ascat_bc <- ascat.GCcorrect(ascat_bc, GCcontentfile =
        #                                 paste0(file_dir, "/", patient, ".GCcontent.txt"))
        ascat.plotRawData(ascat_bc)
        if (tumor_sample == tumor_samples[1])
        {
            ascat_bc <- ascat.aspcf(ascat_bc, penalty = ascat_penalty[1])
        } else {
            ascat_bc <- ascat.aspcf(ascat_bc, penalty = ascat_penalty[2])
        }
        
        ascat.plotSegmentedData(ascat_bc)
        ascat_output <- ascat.runAscat(ascat_bc, gamma = 1.0)
        segs <- ascat_output$segments
        
        # Save the segmentation in the default ASCAT format
        write.table(segs, file = paste0(tumor_sample, "_rphase_segs.tsv"), 
                    sep = "\t", quote = FALSE, row.names = FALSE)
        
        rphase_sample_data[tumor_sample, "purity"] <- ascat_output$aberrantcellfraction[[1]]
        rphase_sample_data[tumor_sample, "ploidy"] <- ascat_output$ploidy[[1]]
    }
    
    write.table(rphase_sample_data, file = "rphase_sample_data.tsv", sep = "\t", 
                quote = FALSE, row.names = FALSE)
    
    # rename header of "_rphase_segs.tsv"
    sample_data <- read.delim("rphase_sample_data.tsv", stringsAsFactors = FALSE)
    
    for (i in 1:nrow(sample_data))
    {
        for (j in 2:3)
        {
            tmp <- read.delim(sample_data[i,j], stringsAsFactors = FALSE)
            if (j == 2)
            {
                colnames(tmp) <- c("sample", "chrom", "start",
                                   "end", "cn_major", "cn_minor")
            } else {
                colnames(tmp) <- c("chrom", "pos", "baf",
                                   "logr", "germline_zygosity")
            }
            tmp <- tmp[complete.cases(tmp[,"chrom"]),]
            write.table(tmp, file = sample_data[i,j], sep = "\t", quote = FALSE,
                        row.names = FALSE, col.names = TRUE)
        }
    }
    
    # Pre-load all experimental data and sample metadata
    exper <- refphase_load(sample_data, segs_format = "refphase", snps_format = "refphase")
    
    # (Optional) If your data shows reference bias, which presents itself as BAFs
    # in regions with a balanced copy number that are systematically shifted away
    # from 0.5 (usually something like 0.47), this can try to correct for that.
    exper <- center_baf(exper)
    
    # (Optional) Fit SNP logr data to improve copy number re-estimation in refphase,
    # when using the default ASCAT formula-based method fo re-estimating copy
    # numbers
    exper <- fit_logr_to_ascat(exper)
    
    # Run rphase on the experimental data
    results <- refphase(exper, cn_method = "ascat",
                        merge_min_width = 1e+06,
                        max_merge_width = 1e+06,
                        gap_fill_width = 2e+06, verbosity = "INFO")
    
    
    refphase:::write_segs(results$phased_segs, file = paste0(patient, "-refphase-segmentation.tsv"))
    
    # (optional) output the SNPs, including phasing information
    refphase:::write_snps(results$phased_snps, file = paste0(patient, "-refphase-phased-snps.tsv.gz"))
    
    # Ploidy might have changed, if we have updated any copy numbers
    write.table(results$sample_data, file = paste0(patient, "-refphase-sample-data-updated.tsv"), sep = "\t", row.names = FALSE)
    
    # Plot the results as PDFs (these files can be huge, over 100 MB in some cases)
    pdf(paste0(patient, "-refphase-genome.pdf"), width = 10, height = 4 * nrow(sample_data))
    plot_genome(results$sample_data, results$phased_snps, results$phased_segs)
    dev.off()
    
    pdf(paste0(patient, "-refphase-chromosomes.pdf"), width = 10, height = 2 + 4 * nrow(sample_data))
    plot_all_chromosomes(results$sample_data, results$phased_snps, list(refphase = results$phased_segs, orig = exper$segs))
    dev.off()
    
    # # Plot the results as individual PNG files
    # png(paste0(patient, "-refphase-genome.png"), width = 600, height = 200 + 300 * nrow(sample_data))
    # plot_genome(results$sample_data, results$phased_snps, results$phased_segs)
    # dev.off()
    # 
    # png(paste0(patient, "-refphase-chromosomes-%02d.png"), width = 700, height = 300 * nrow(sample_data))
    # plot_all_chromosomes(results$sample_data, results$phased_snps, list(refphase = results$phased_segs, orig = exper$segs))
    # dev.off()
}




