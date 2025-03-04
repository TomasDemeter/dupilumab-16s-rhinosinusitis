#!/usr/bin/env Rscript

###########################
# Load required libraries #
###########################
suppressPackageStartupMessages({
    library(dada2)
    library(tidyverse)
    library(tibble)
    library(phyloseq)
    library(ShortRead)
    library(optparse)
})

#############
# Functions #
#############
all_orients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}

remove_primers <- function(input_path_fwd, input_path_rev, output_path_fwd, output_path_rev, forward_primer_seq, reverse_primer_seq, minimum_length, nextseq_trim, cores) {
    # Get all orientations of the forward and reverse primers
    FWD_orients <- all_orients(forward_primer_seq)
    REV_orients <- all_orients(reverse_primer_seq)

    # Create directory to hold the output from cutadapt if it doesn't exist
    output_dir <- dirname(output_path_fwd)
    if (!dir.exists(output_dir)) dir.create(output_dir)

    # Save the reverse complements of the primers to variables
    FWD_RC <- dada2:::rc(FWD_orients)
    REV_RC <- dada2:::rc(REV_orients)

    ##  Create the cutadapt flags ##
    # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
    R1_flags <- paste("-g", FWD_orients, "-a", REV_RC, "--minimum-length", minimum_length)

    # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
    R2_flags <- paste("-G", REV_orients, "-A", FWD_RC, "--minimum-length", minimum_length)

    # Run Cutadapt on the single pair of files
    system2("cutadapt", args = c("--cores", cores, "--nextseq-trim", as.character(nextseq_trim), R1_flags, R2_flags, "-n", "2", # -n 2 required to remove FWD and REV from reads
                                "-o", output_path_fwd, "-p", output_path_rev,
                                input_path_fwd, input_path_rev))
}

##################################################
# Main function to handle command-line arguments #
##################################################
main <- function() {
    # Define options
    option_list <- list(
        make_option(c("-i", "--input_path_fwd"), type = "character", help = "Path to forward input files"),
        make_option(c("-I", "--input_path_rev"), type = "character", help = "Path to reverse input files"),
        make_option(c("-o", "--output_path_fwd"), type = "character", help = "Path for forward output files"),
        make_option(c("-O", "--output_path_rev"), type = "character", help = "Path for reverse output files"),
        make_option(c("-F", "--forward_primer_seq"), type = "character", help = "Forward primer sequence"),
        make_option(c("-R", "--reverse_primer_seq"), type = "character", help = "Reverse primer sequence"),
        make_option(c("-m", "--minimum_length"), type = "integer", default = 50, help = "Minimum length of reads after primer removal"),
        make_option(c("-n", "--nextseq_trim"), type = "integer", default = 20, help = "Number of bases to trim from the end of NextSeq reads"),
        make_option(c("-j", "--cores"), type = "integer", default = 0, help = "Number of cores to use")
        )

    # Parse options
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)

    # Check if all required arguments are provided
    required_args <- c("input_path_fwd", "input_path_rev", "output_path_fwd", "output_path_rev", "forward_primer_seq", "reverse_primer_seq")
    missing_args <- setdiff(required_args, names(opt[!sapply(opt, is.null)]))

    if (length(missing_args) > 0) {
        stop(sprintf("Missing required arguments: %s", paste(missing_args, collapse = ", ")), call. = FALSE)
    }

    # Main workflow
remove_primers(opt$input_path_fwd, opt$input_path_rev, opt$output_path_fwd, opt$output_path_rev, opt$forward_primer_seq, opt$reverse_primer_seq, opt$minimum_length, opt$nextseq_trim, opt$cores)
}

# Execute the main function if the script is run directly
if (!interactive()) {
    main()
}