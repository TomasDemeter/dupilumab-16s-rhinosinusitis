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
    library(Biostrings)
    library(utils)
})


all_orients <- function(primer) {
    # Create all orientations of the input sequence
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}

primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}



detect_primers <- function(FWD_before_removal, REV_before_removal, FWD_after_removal, REV_after_removal, forward_primer_seq, reverse_primer_seq, output_path) {
    # Get all orientations of the forward and reverse primers
    FWD_orient <- all_orients(forward_primer_seq)
    REV_orient <- all_orients(reverse_primer_seq)
    
    before_primer_removal <- rbind(
        FWD.ForwardReads = sapply(FWD_orient, primerHits, fn = FWD_before_removal),
        FWD.ReverseReads = sapply(FWD_orient, primerHits, fn = REV_before_removal),
        REV.ForwardReads = sapply(REV_orient, primerHits, fn = FWD_before_removal),
        REV.ReverseReads = sapply(REV_orient, primerHits, fn = REV_before_removal)
    )

    after_primer_removal <- rbind(
        FWD.ForwardReads = sapply(FWD_orient, primerHits, fn = FWD_after_removal),
        FWD.ReverseReads = sapply(FWD_orient, primerHits, fn = REV_after_removal),
        REV.ForwardReads = sapply(REV_orient, primerHits, fn = FWD_after_removal),
        REV.ReverseReads = sapply(REV_orient, primerHits, fn = REV_after_removal)
    )

    # Combine results into a list
    results <- list(
    Before_Primer_Removal = before_primer_removal,
    After_Primer_Removal = after_primer_removal
    )

    # Write the results to a text file
    writeLines(c("Before Primer Removal:", 
                capture.output(print(results$Before_Primer_Removal)),
                "",
                "After Primer Removal:",
                capture.output(print(results$After_Primer_Removal))),
                con = output_path)

    # Return the results invisibly
    invisible(results)
}


##################################################
# Main function to handle command-line arguments #
##################################################
main <- function() {
    # Define options
    option_list <- list(
        make_option(c("-i", "--FWD_before_removal"), type = "character", help = "Path to forward input files"),
        make_option(c("-I", "--REV_before_removal"), type = "character", help = "Path to reverse input files"),
        make_option(c("-o", "--FWD_after_removal"), type = "character", help = "Path to forward output files"),
        make_option(c("-O", "--REV_after_removal"), type = "character", help = "Path to reverse output files"),
        make_option(c("-F", "--forward_primer_seq"), type = "character", help = "Forward primer sequence"),
        make_option(c("-R", "--reverse_primer_seq"), type = "character", help = "Reverse primer sequence"),
        make_option(c("-p", "--output_path"), type = "character", help = "Path to output txt file")
        )

    # Parse options
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)

    # Check if all required arguments are provided
    required_args <- c("FWD_before_removal", "REV_before_removal", "FWD_after_removal", "REV_after_removal", "forward_primer_seq", "reverse_primer_seq", "output_path")
    missing_args <- setdiff(required_args, names(opt[!sapply(opt, is.null)]))

    if (length(missing_args) > 0) {
        stop(sprintf("Missing required arguments: %s", paste(missing_args, collapse = ", ")), call. = FALSE)
    }

    # Main workflow
detect_primers(opt$FWD_before_removal, opt$REV_before_removal, opt$FWD_after_removal, opt$REV_after_removal, opt$forward_primer_seq, opt$reverse_primer_seq, opt$output_path)
}

# Execute the main function if the script is run directly
if (!interactive()) {
    main()
}