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


#############
# Functions #
#############
# Create output directories and symlink input files
create_output_directories <- function(input_path, output_path, pattern_FWD_read, pattern_REV_read) {
    # Create output directories
    filter_fp <- normalizePath(file.path(output_path, "02_filter"), mustWork = FALSE)
    dir.create(filter_fp, recursive = TRUE, showWarnings = FALSE)

    ## Get the forward and reverse read files
    fnFs_cut <- sort(list.files(input_path, pattern = pattern_FWD_read, full.names = TRUE))
    fnRs_cut <- sort(list.files(input_path, pattern = pattern_REV_read, full.names = TRUE))

    # Create subdirectories for preprocessed files
    subF_fp <- normalizePath(file.path(filter_fp, "preprocessed_F"), mustWork = FALSE)
    subR_fp <- normalizePath(file.path(filter_fp, "preprocessed_R"), mustWork = FALSE)
    dir.create(subF_fp, recursive = TRUE, showWarnings = FALSE)
    dir.create(subR_fp, recursive = TRUE, showWarnings = FALSE)

    # Symlink input files to output directories, only if they do not already exist
    for (fnF in fnFs_cut) {
        fnF_Q <- normalizePath(file.path(subF_fp, basename(fnF)), mustWork = FALSE)
        if (!file.exists(fnF_Q)) {
            file.symlink(from = normalizePath(fnF, mustWork = FALSE), to = fnF_Q)
        }
    }
    for (fnR in fnRs_cut) {
        fnR_Q <- normalizePath(file.path(subR_fp, basename(fnR)), mustWork = FALSE)
        if (!file.exists(fnR_Q)) {
            file.symlink(from = normalizePath(fnR, mustWork = FALSE), to = fnR_Q)
        }
    }

    return(list(subF_fp = subF_fp, subR_fp = subR_fp))
}

# Filter and trim reads
filter_and_trim_reads <- function(subF_fp, subR_fp, pattern_FWD_read, pattern_REV_read, trunc_length_FWD, trunc_length_REV, maxEE_FWD, maxEE_REV, truncQ, multithread) {
    # Create output directories
    filtpathF <- file.path(subF_fp, "filtered")
    filtpathR <- file.path(subR_fp, "filtered")
    dir.create(filtpathF, recursive = TRUE, showWarnings = FALSE)
    dir.create(filtpathR, recursive = TRUE, showWarnings = FALSE)

    # create dir for saving stats about lost reads
    lost_reads_fp <- file.path(dirname(subF_fp))

    # Get forward and reverse read files
    fastqFs <- sort(list.files(subF_fp, pattern = pattern_FWD_read))
    fastqRs <- sort(list.files(subR_fp, pattern = pattern_REV_read))
    if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

    # Filter and trim reads
    filt_out <- filterAndTrim(fwd = file.path(subF_fp, fastqFs), filt = file.path(filtpathF, fastqFs),
                                rev = file.path(subR_fp, fastqRs), filt.rev = file.path(filtpathR, fastqRs),
                                truncLen = c(trunc_length_FWD, trunc_length_REV), maxEE = c(maxEE_FWD, maxEE_REV), truncQ = truncQ,
                                rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = multithread)

    print(filt_out)
    saveRDS(filt_out, file.path(lost_reads_fp, "filt_out.rds"))
    return(list(filt_out = filt_out, filtpathF = filtpathF, filtpathR = filtpathR))
}




generate_quality_plots <- function(subF_fp, subR_fp, filtpathF, filtpathR, pattern_FWD_read, pattern_REV_read) {
    # Get forward and reverse read files before trimming
    fastqFs <- sort(list.files(subF_fp, pattern = pattern_FWD_read, full.names = TRUE))
    fastqRs <- sort(list.files(subR_fp, pattern = pattern_REV_read, full.names = TRUE))
    num_samples <- length(fastqFs)
    sample_indices <- if(num_samples <= 20) 1:num_samples else sample(1:num_samples, 20)

    # Generate quality plots for files before trimming
    fwd_qual_plots_before <- plotQualityProfile(fastqFs[sample_indices])
    rev_qual_plots_before <- plotQualityProfile(fastqRs[sample_indices])

    # Get forward and reverse read files after trimming
    fastqFs <- sort(list.files(filtpathF, full.names = TRUE))
    fastqRs <- sort(list.files(filtpathR, full.names = TRUE))
    num_samples <- length(fastqFs)
    sample_indices <- if(num_samples <= 20) 1:num_samples else sample(1:num_samples, 20)

    # Generate quality plots for files after trimming
    fwd_qual_plots_after <- plotQualityProfile(fastqFs[sample_indices])
    rev_qual_plots_after <- plotQualityProfile(fastqRs[sample_indices])

    # Save plots
    saveRDS(list(before = fwd_qual_plots_before, after = fwd_qual_plots_after), file.path(subF_fp, "fwd_qual_plots.rds"))
    saveRDS(list(before = rev_qual_plots_before, after = rev_qual_plots_after), file.path(subR_fp, "rev_qual_plots.rds"))

    # Save plots as images
    ggsave(plot = fwd_qual_plots_before, filename = file.path(subF_fp, "fwd_qual_plots_before.png"), width = 12, height = 10, dpi = "retina")
    ggsave(plot = rev_qual_plots_before, filename = file.path(subR_fp, "rev_qual_plots_before.png"), width = 12, height = 10, dpi = "retina")
    ggsave(plot = fwd_qual_plots_after, filename = file.path(subF_fp, "fwd_qual_plots_after.png"), width = 12, height = 10, dpi = "retina")
    ggsave(plot = rev_qual_plots_after, filename = file.path(subR_fp, "rev_qual_plots_after.png"), width = 12, height = 10, dpi = "retina")
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
        make_option(c("-F", "--trunc_length_FWD"), type = "integer", help = "Length to truncate forward reads"),
        make_option(c("-R", "--trunc_length_REV"), type = "integer", help = "Length to truncate reverse reads"),
        make_option(c("-e", "--maxEE_FWD"), type = "numeric", default = 2, help = "Maximum expected errors for forward reads"),
        make_option(c("-E", "--maxEE_REV"), type = "numeric", default = 2, help = "Maximum expected errors for reverse reads"),
        make_option(c("-q", "--truncQ"), type = "numeric", default = 2, help = "Quality threshold for truncation"),
        make_option(c("-m", "--multithread"), type = "logical", default = TRUE, help = "Use multiple threads for filtering and trimming")
    )

    # Parse options
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)

    # Check if all required arguments are provided
    required_args <- c("input_path", "output_path", "trunc_length_FWD", "trunc_length_REV")
    missing_args <- setdiff(required_args, names(opt[!sapply(opt, is.null)]))

    if (length(missing_args) > 0) {
        stop(sprintf("Missing required arguments: %s", paste(missing_args, collapse = ", ")), call. = FALSE)
    }

    # Create output directories and symlink input files
    dirs <- create_output_directories(opt$input_path, opt$output_path, opt$pattern_FWD_read, opt$pattern_REV_read)


    # Main workflow
    filter_trim_results <- filter_and_trim_reads(dirs$subF_fp, dirs$subR_fp, opt$pattern_FWD_read, opt$pattern_REV_read, opt$trunc_length_FWD, opt$trunc_length_REV, opt$maxEE_FWD, opt$maxEE_REV, opt$truncQ, opt$multithread)

    # Generate quality plots
    generate_quality_plots(dirs$subF_fp, dirs$subR_fp, filter_trim_results$filtpathF, filter_trim_results$filtpathR, opt$pattern_FWD_read, opt$pattern_REV_read)
}

# Execute the main function if the script is run directly
if (!interactive()) {
    main()
}