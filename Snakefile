import yaml
import os
import math
import pandas as pd

# ############ SETUP ##############################
configfile: "config/config_AMLM7.yaml"
# configfile: "configs/config.json"
workdir: config['workdir']
snakedir = os.path.dirname(workflow.snakefile)
# include helper functions
# include: "includes/io.snk"
include: "includes/utils.snk"

chrom_list = get_chrom_list(config)

# get the samples
sample_sheet = os.path.join(snakedir, config['samples']['sheet'])
sample_df = get_bam_files(config['bam_folder'], sample_sheet)
print(sample_df)


# ############ INCLUDES ##############################  
include: "includes/preprocessing.snk"
include: "includes/ascat.snk"

# specified wildcards have to match the regex
wildcard_constraints:
    # eg sample cannot contain _ or / to prevent ambiguous wildcards
    sample = "[^/._]+",
    tumor = "[AR][1-3]?",
    normal = "B",
    type = "[^/._]+"

# ############## MASTER RULE ##############################################

rule all:
    input:
        expand("pos/{sample}.pos.gz", sample=sample_df.index),
        # quick fix for tumor-normal pairs
        expand("cnv/{sample}/{sample}-B_refphase_genome.pdf", sample=[s for s in sample_df.index if s.split("_")[1] not in config['samples']['normal']])

###########################################################################

# print out of installed tools
onstart:
    print("    CNV ASCAT PIPELINE STARTING.......")
    # write config to the results directory
    path_to_config = os.path.join(config['workdir'], "config.yaml")
    with open(path_to_config, 'w+') as stream:
        yaml.dump(config, stream, default_flow_style=False)
    # create logs folder


onsuccess:
    # shell("export PATH=$ORG_PATH; unset ORG_PATH")
    print("Workflow finished - everything ran smoothly")

    # cleanup
    if config['cleanup']:
        print('I would like to cleanup but I don\'t know how!')
