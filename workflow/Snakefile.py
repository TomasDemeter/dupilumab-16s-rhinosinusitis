####################################################
# Snakemake pipeline for 16S metagenomics analysis #
####################################################

####################
# Python pacakages #
####################
import pandas as pd

# Configuration
configfile: "../config/config.yaml"
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]
samplefile = config["refs"]["metadata"]

##########################
# Samples and conditions #
##########################
samples = pd.read_csv(samplefile)
SAMPLES = samples["Sample_ID"].tolist()

#########
# rules #
#########
include: "rules/dada2_preprocessing.smk"
include: "rules/dada2_filter.smk"
include: "rules/fastqc.smk"
include: "rules/multiqc.smk"
include: "rules/dada2_taxonomy.smk"
include: "rules/dada2_phyloseq.smk"

###################
# Desired outputs #
###################
MULTIQC_RAW         = rules.multiqc_raw.output.report
MULTIQC_FILTERED    = rules.multiqc_filtered.output.report
PHYLOSEQ            = rules.Dada2_to_phyloseq.output.phyloseq
TRACKING_PLOT       = rules.Dada2_loss_plotting.output.tracking_reads_summary_plot
PRIMER_DETECTION    = expand(rules.detect_primers.output.txt_file, sample = SAMPLES)

############
# Pipeline #
############
rule all:
    input:
        MULTIQC_RAW,
        PRIMER_DETECTION,
        MULTIQC_FILTERED,
        PHYLOSEQ,
        TRACKING_PLOT
    message:
        "16S metagenomic pipeline run complete!"