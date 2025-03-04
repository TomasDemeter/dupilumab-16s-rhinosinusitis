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
dada2_sequence_table <- function(filtpathF, filtpathR, errF, errR, pattern_FWD_read, pattern_REV_read, multithread, pool) {

    parent_dir <- dirname(dirname(dirname(errF)))
    table_fp <- file.path(parent_dir, "03_tabletax")
    if(!dir.exists(table_fp)) dir.create(table_fp)

    # read RDS files
    errF <- readRDS(errF)
    errR <- readRDS(errR)

    filtFs <- list.files(filtpathF, pattern = pattern_FWD_read, full.names = TRUE)
    filtRs <- list.files(filtpathR, pattern = pattern_REV_read, full.names = TRUE)

    # Get the sample names
    sample_names <- basename(filtFs) # doesn't drop fastq.gz
    sample_names <- gsub(pattern_FWD_read, "", sample_names)
    sample_namesR <- basename(filtRs) # doesn't drop fastq.gz 
    sample_namesR <- gsub(pattern_REV_read, "", sample_namesR)

    # Double check
    if(!identical(sample_names, sample_namesR)) stop("Forward and reverse files do not match.")
    names(filtFs) <- sample_names
    names(filtRs) <- sample_names

    # make lists to hold the loop output
    mergers <- vector("list", length(sample_names))
    names(mergers) <- sample_names
    ddF <- vector("list", length(sample_names))
    names(ddF) <- sample_names
    ddR <- vector("list", length(sample_names))
    names(ddR) <- sample_names

    # For each sample, get a list of merged and denoised sequences
    for(sam in sample_names) {
        cat("Processing:", sam, "\n")
        # Dereplicate forward reads
        derepF <- derepFastq(filtFs[[sam]])
        # Infer sequences for forward reads
        dadaF <- dada(derepF, err = errF, multithread = multithread, pool = pool)
        ddF[[sam]] <- dadaF
        # Dereplicate reverse reads
        derepR <- derepFastq(filtRs[[sam]])
        # Infer sequences for reverse reads
        dadaR <- dada(derepR, err = errR, multithread = multithread, pool = pool)
        ddR[[sam]] <- dadaR
        # Merge reads together
        merger <- mergePairs(ddF[[sam]], derepF, ddR[[sam]], derepR)
        mergers[[sam]] <- merger
    }

    # Merge sequences
    seqtab <- makeSequenceTable(mergers)

    # Save the output
    saveRDS(ddF, file.path(table_fp, "ddF.rds"))
    saveRDS(ddR, file.path(table_fp, "ddR.rds"))
    saveRDS(mergers, file.path(table_fp, "mergers.rds"))
    saveRDS(seqtab, file.path(table_fp, "seqtab.rds"))

    return(list(seqtab = seqtab, ddF = ddF, ddR = ddR, table_fp = table_fp))
}


# Plot the read length distribution
plot_length_distribution <- function(seqtab_data) {

    # Get the sequence table
    seqtab <- seqtab_data$seqtab

    # Get the table file path
    table_fp <- seqtab_data$table_fp
    
    # Get the read length distribution
    length_distribution <- as.data.frame(table(nchar(getSequences(seqtab))))

    # Plot the read length distribution
    length_distribution_plot <- ggplot(length_distribution, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    labs(title = "Read Lengths Distribution",
            x = "Read Length",
            y = "Count") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
            panel.background = element_rect(fill = "white"))

    # Save the plot
    ggsave(plot = length_distribution_plot, filename = paste0(table_fp, "/length_distribution_plot.png"), width = 15, height = 15, dpi = "retina", bg = "white")
}

##################################################
# Main function to handle command-line arguments #
##################################################
main <- function() {
    # Define options
    option_list <- list(
        make_option(c("-F", "--filtpathF"), type = "character", help = "Forward filtered reads"),
        make_option(c("-R", "--filtpathR"), type = "character", help = "Reverse filtered reads"),
        make_option(c("-e", "--errF"), type = "character", help = "Forward error rates"),
        make_option(c("-E", "--errR"), type = "character", help = "Reverse error rates"),
        make_option(c("-f", "--pattern_FWD_read"), type = "character", default = "_1.fq.gz", help = "Pattern for forward read"),
        make_option(c("-r", "--pattern_REV_read"), type = "character", default = "_2.fq.gz", help = "Pattern for reverse read"),
        make_option(c("-m", "--multithread"), type = "logical", default = TRUE, help = "Use multiple threads"),
        make_option(c("-p", "--pool"), type = "logical", default = FALSE, help = "Pool reads")
    )

    # Parse options
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
    # Check if all required arguments are provided
    required_args <- c("filtpathF", "filtpathR", "errF", "errR", "pattern_FWD_read", "pattern_REV_read")
    missing_args <- setdiff(required_args, names(opt[!sapply(opt, is.null)]))

    if (length(missing_args) > 0) {
        stop("Missing arguments: ", paste(missing_args, collapse = ", "))
    }

    # Main workflow
    seqtab_data <- dada2_sequence_table(opt$filtpathF, opt$filtpathR, opt$errF, opt$errR, opt$pattern_FWD_read, opt$pattern_REV_read, opt$multithread, opt$pool)
    plot_length_distribution(seqtab_data)
}

# Execute the main function if the script is run directly
if (!interactive()) {
    main()
}
