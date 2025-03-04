#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(phyloseq)
  library(parallel)
  library(doParallel)
  library(optparse)
})

# Set up command line arguments
option_list <- list(
  make_option(c("-i", "--input"), type="character", 
              help="Input RDS file containing phyloseq object"),
  make_option(c("-o", "--output"), type="character", 
              help="Output RDS file path"),
  make_option(c("-n", "--iterations"), type="integer", default=100,
              help="Number of iterations [default=%default]"),
  make_option(c("-s", "--sample_size"), type="integer", default=NULL,
              help="Sample size for rarefaction [default=minimum sample sum]"),
  make_option(c("-c", "--cores"), type="integer", default=NULL,
              help="Number of cores to use [default=detected cores - 1]"),
  make_option(c("-r", "--replace"), action="store_true", default=FALSE,
              help="Sample with replacement [default=%default]"),
  make_option(c("-t", "--trim"), action="store_true", default=TRUE,
              help="Trim OTUs [default=%default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Verbose output [default=%default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Function definition
parallel_rarefy_phyloseq <- function(physeq_obj, n_iterations = 100, sample_size = NULL, 
                                   rngseed = 42, replace = FALSE, trimOTUs = TRUE, 
                                   verbose = FALSE, n_cores = NULL) {
  
  # Set up parallel processing
  if (is.null(n_cores)) {
    n_cores <- detectCores() - 1
  }
  
  if(verbose) cat(sprintf("Using %d cores\n", n_cores))
  
  # Initialize cluster
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Export required packages and objects to all cores
  clusterEvalQ(cl, {
    library(phyloseq)
  })
  
  # Export the phyloseq object to all cores
  clusterExport(cl, c("physeq_obj", "sample_size", "replace"), environment())
  
  # Set sample size if not provided
  if (is.null(sample_size)) {
    sample_size <- min(sample_sums(physeq_obj))
  }
  
  if(verbose) cat(sprintf("Sample size: %d\n", sample_size))
  
  # Generate seeds
  set.seed(rngseed)
  seeds <- sample.int(1e6, n_iterations)
  
  # Parallel rarefaction
  if(verbose) cat("Starting rarefaction...\n")
  rarefied_list <- parLapply(cl, seeds, function(seed) {
    rarefy_even_depth(physeq = physeq_obj,
                     sample.size = sample_size,
                     rngseed = seed,
                     replace = replace,
                     trimOTUs = FALSE,
                     verbose = FALSE)
  })
  
  # Stop cluster
  stopCluster(cl)
  
  if(verbose) cat("Merging results...\n")
  
  # Merge the rarefied OTU tables
  merged_otu <- Reduce('+', lapply(rarefied_list, function(x) otu_table(x)))
  
  # Average and round the merged OTU table
  avg_otu <- round(merged_otu / n_iterations)
  
  # Create new phyloseq object with averaged OTU table
  avg_physeq <- phyloseq(otu_table(avg_otu, taxa_are_rows = taxa_are_rows(physeq_obj)),
                        sample_data(physeq_obj),
                        tax_table(physeq_obj),
                        refseq(physeq_obj))
  
  # Optionally trim OTUs
  if (trimOTUs) {
    if(verbose) cat("Trimming OTUs...\n")
    avg_physeq <- prune_taxa(taxa_sums(avg_physeq) > 0, avg_physeq)
  }
  
  return(avg_physeq)
}

# Main execution
if(opt$verbose) cat("Loading phyloseq object...\n")
physeq <- readRDS(opt$input)

if(opt$verbose) cat("Starting parallel rarefaction process...\n")
result <- parallel_rarefy_phyloseq(
  physeq_obj = physeq,
  n_iterations = opt$iterations,
  sample_size = opt$sample_size,
  replace = opt$replace,
  trimOTUs = opt$trim,
  verbose = opt$verbose,
  n_cores = opt$cores
)

if(opt$verbose) cat("Saving results...\n")
saveRDS(result, file = opt$output)

if(opt$verbose) cat("Done!\n")
