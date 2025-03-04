#!/usr/bin/env Rscript

# Rscript 06_dada2_chimeras_taxonomy.R -i ../../results/dada2/03_tabletax/seqtab.rds -d ../../inputs/database/silva_nr99_v138.1_train_set.fa.gz -l 400 -u 435

###########################
# Load required libraries #
###########################
library(dada2)
library(tidyverse)
library(tibble)
library(phyloseq)
library(ShortRead)
library(optparse)

#############
# Functions #
#############
remove_chimeras_run_taxonomy <- function(input_path, database, lower_cutoff, upper_cutoff, minFoldParentOverAbundance, multithread){
    # Read in RDS 
    seqtab <- readRDS(input_path)

    # Remove sequences that are not within the expected length range
    seqtab <- seqtab[,nchar(colnames(seqtab)) %in% lower_cutoff:upper_cutoff]

    input_path <- dirname(input_path)

    # Remove chimeras
    seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread = multithread, minFoldParentOverAbundance = minFoldParentOverAbundance)

    # Assign taxonomy
    tax <- assignTaxonomy(seqtab_nochim, database, tryRC = TRUE,
                        multithread = multithread)

    saveRDS(seqtab_nochim, paste0(input_path, "/seqtab_final.rds"))
    saveRDS(tax, paste0(input_path, "/tax_final.rds"))

}

##################################################
# Main function to handle command-line arguments #
##################################################
main <- function() {
    # Define options
    option_list <- list(
        make_option(c("-i", "--input_path"), type = "character", help = "Path to input files"),
        make_option(c("-d", "--database"), type = "character", help = "Path to database"),
        make_option(c("-l", "--lower_cutoff"), type = "integer", default = 0, help = "Lower cutoff for sequence length"),
        make_option(c("-u", "--upper_cutoff"), type = "integer", default = 1000, help = "Upper cutoff for sequence length"),
        make_option(c("-m", "--multithread"), type = "logical", default = TRUE, help = "Use multiple threads"),
        make_option(c("-p", "--minFoldParentOverAbundance"), type = "integer", default = 2, help = " Only sequences greater than this-fold more abundant than a sequence can be its parent")
    )

    # Parse options
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)

    # Check if all required arguments are provided
    required_args <- c("input_path", "database")
    missing_args <- setdiff(required_args, names(opt[!sapply(opt, is.null)]))

    if (length(missing_args) > 0) {
        stop("Missing arguments: ", paste(missing_args, collapse = ", "))
    }

    # Run the main function
    remove_chimeras_run_taxonomy(opt$input_path, opt$database, opt$lower_cutoff, opt$upper_cutoff, opt$minFoldParentOverAbundance, opt$multithread)
}

# Execute the main function if the script is run directly
if (!interactive()) {
    main()
}