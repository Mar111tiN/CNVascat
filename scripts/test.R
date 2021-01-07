
#
# Description:
# General refphase workflow for detecting SCNA in tumor samples.
# See https://bitbucket.org/schwarzlab/refphase/src/master/ for further details.
# Adapted from Matt Huska (AG Schwarz - MDC Berlin)

## create conda environment using the env file R-env.yaml
# conda env create -n R-env -f ../env/R-env.yaml
## Install refphase from R
# cd scripts; R
# > library(devtools)
# > devtools::install("refphase")

############# IMPORTS ######################
library(tidyverse)
library(gtools)
library(ASCAT)
library(refphase)

############# CUSTOM FUNCTIONS #################

chrom2integer <- function(x) {
    print(unique(x))
    as.integer(gsub("M", "25", 
                    gsub("X", "23", 
                         gsub("Y", "24", gsub("chr", "", x)))))
}

tumor <- snakemake@input[["tumor"]]
params <- snakemake@params
cutoff <- params[['cutoff']]

print(tumor)
print(cutoff)


