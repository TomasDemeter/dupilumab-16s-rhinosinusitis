#!/usr/bin/env Rscript

###########################
# Load required libraries #
###########################
library(dada2)
library(dplyr)
library(tibble)
library(phyloseq)
library(Biostrings)
library(optparse)

############
# Function #
############
create_phyloseq <- function(counts, taxa, metadata, sample_names, output_path, experiment) {
    # Read the data
    counts_raw <- readRDS(counts)
    taxa_raw <- readRDS(taxa)
    
    # Create short ASV IDs first
    asv_ids <- paste0("ASV", seq(ncol(counts_raw)))
    
    # Store the sequence-to-ASV mapping
    seq_to_asv <- setNames(asv_ids, colnames(counts_raw))
    
    # Process counts table with ASV IDs
    counts_obj <- counts_raw
    colnames(counts_obj) <- seq_to_asv[colnames(counts_obj)]
    counts_obj <- otu_table(counts_obj, taxa_are_rows = FALSE)
    
    # Process taxa table with ASV IDs
    taxa_obj <- taxa_raw
    rownames(taxa_obj) <- seq_to_asv[rownames(taxa_obj)]
    taxa_obj <- tax_table(as.matrix(taxa_obj))
    
    # Read metadata
    metadata <- read.csv(metadata, row.names = sample_names) %>%
        sample_data()
    
    # Create initial phyloseq object
    ps <- phyloseq(counts_obj, taxa_obj, metadata)
    
    # Add reference sequences
    # Store original sequences as refseq
    dna <- Biostrings::DNAStringSet(names(seq_to_asv))
    names(dna) <- seq_to_asv
    ps@refseq <- dna
    
    # Save phyloseq object
    output_filename <- paste0(experiment, "_phyloseq.rds")
    saveRDS(ps, file.path(output_path, output_filename))
    
    # Print some validation info
    print(paste("Created phyloseq object with", ntaxa(ps), "taxa and", nsamples(ps), "samples"))
    print("Example of first few taxa:")
    print(head(taxa_names(ps)))
}

##################################################
# Main function to handle command-line arguments #
##################################################

main <- function() {
    # Define options
    option_list <- list(
        make_option(c("-c", "--counts"), type="character", default=NULL, help="Path to the counts file"),
        make_option(c("-t", "--taxa"), type="character", default=NULL, help="Path to the taxa file"),
        make_option(c("-m", "--metadata"), type="character", default=NULL, help="Path to the metadata file"),
        make_option(c("-o", "--output_path"), type="character", default=NULL, help="Path to the output directory"),
        make_option(c("-s", "--sample_names"), type="character", default = "Sample_ID", help="Name of the column containing sample names (default: Sample_ID)"),
        make_option(c("-e", "--experiment"), type="character", default="phyloseq", help="Experiment name for output file (default: phyloseq)")
)

    # Parse command line arguments
    opt_parser <- OptionParser(option_list=option_list)
    opt <- parse_args(opt_parser)

    required_args <- c("counts", "taxa", "metadata", "output_path")
    missing_args <- required_args[!required_args %in% names(opt)]
    
    if (length(missing_args) > 0) {
        stop(paste("Missing required arguments:", paste(missing_args, collapse=", ")))
    }
    
    # Call the create_phyloseq function
    create_phyloseq(opt$counts, opt$taxa, opt$metadata, opt$sample_names, opt$output_path, opt$experiment)
}

# Execute the main function if the script is run directly
if (!interactive()) {
    main()
}