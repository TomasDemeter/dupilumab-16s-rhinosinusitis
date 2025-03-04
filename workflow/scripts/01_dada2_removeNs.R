#!/usr/bin/env Rscript

###########################
# Load required libraries #
###########################
library(dada2)
library(tidyverse)
library(tibble)
library(phyloseq)
library(ShortRead)
library(optparse)

remove_Ns <- function(input_fp, results_fp,  pattern_FWD_read, pattern_REV_read, maxN, multithread){

    # Create output directories if they do not exist
    if (!dir.exists(results_fp)) {
        dir.create(results_fp, recursive = TRUE)
    }

    preprocess_fp <- file.path(results_fp, "01_preprocess")
    if (!dir.exists(preprocess_fp)) {
        dir.create(preprocess_fp, recursive = TRUE)
    }

    filtN_fp <- file.path(preprocess_fp, "filtN")
    if (!dir.exists(filtN_fp)) {
        dir.create(filtN_fp, recursive = TRUE)
    }

    fnFs <- sort(list.files(input_fp, pattern = pattern_FWD_read, full.names = TRUE))
    fnRs <- sort(list.files(input_fp, pattern = pattern_REV_read, full.names = TRUE))

    # Name the N-filtered files to put them in filtN/ subdirectory
    fnFs_filtN <- file.path(filtN_fp, basename(fnFs))
    fnRs_filtN <- file.path(filtN_fp, basename(fnRs))

    # filter Ns
    filterAndTrim(fnFs, fnFs_filtN, fnRs, fnRs_filtN, maxN = maxN, multithread = multithread)
}

##################################################
# Main function to handle command-line arguments #
##################################################
main <- function() {
    # Define options
    option_list <- list(
        make_option(c("-i", "--input_path"), type = "character", help = "Path to input files"),
        make_option(c("-o", "--output_path"), type = "character", help = "Path for output files"),
        make_option(c("-f", "--pattern_FWD_read"), type = "character", default = "_1.fq.gz", help = "Pattern for forward read"),
        make_option(c("-r", "--pattern_REV_read"), type = "character", default = "_2.fq.gz", help = "Pattern for reverse read"),
        make_option(c("-m", "--multithread"), type = "logical", default = TRUE, help = "Use multiple threads"),
        make_option(c("-n", "--maxN"), type = "numeric", default = 0, help = "Maximum number of Ns allowed in a read")
    )

    # Parse options
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)

    # Check if all required arguments are provided
    required_args <- c("input_path", "output_path")
    missing_args <- setdiff(required_args, names(opt[!sapply(opt, is.null)]))

    if (length(missing_args) > 0) {
        stop(sprintf("Missing required arguments: %s", paste(missing_args, collapse = ", ")), call. = FALSE)
    }

    # Main workflow
    remove_Ns(opt$input_path, opt$output_path, opt$pattern_FWD_read, opt$pattern_REV_read, opt$maxN, opt$multithread)
}

# Execute the main function if the script is run directly
if (!interactive()) {
    main()
}