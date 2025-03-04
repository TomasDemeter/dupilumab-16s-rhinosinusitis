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
# learning error rates with NovaSeq data
loessErrfun <- function(trans) {

    # this function was found at https://github.com/ErnakovichLab/dada2_ernakovichlab?tab=readme-ov-file
    # from JacobRPrice alter loess arguments (weights and span and enforce monotonicity)
    # https://github.com/benjjneb/dada2/issues/1307
    qq <- as.numeric(colnames(trans))
    est <- matrix(0, nrow=0, ncol=length(qq))
    for(nti in c("A","C","G","T")) {
        for(ntj in c("A","C","G","T")) {
        if(nti != ntj) {
            errs <- trans[paste0(nti,"2",ntj),]
            tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
            rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
            rlogp[is.infinite(rlogp)] <- NA
            df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)

            # Gulliem Salazar's solution
            # https://github.com/benjjneb/dada2/issues/938
            mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)

            pred <- predict(mod.lo, qq)
            maxrli <- max(which(!is.na(pred)))
            minrli <- min(which(!is.na(pred)))
            pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
            pred[seq_along(pred)<minrli] <- pred[[minrli]]
            est <- rbind(est, 10^pred)
            } # if(nti != ntj)
        } # for(ntj in c("A","C","G","T"))
    } # for(nti in c("A","C","G","T"))
    
    # HACKY
    MAX_ERROR_RATE <- 0.25
    MIN_ERROR_RATE <- 1e-7
    est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
    est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE

    # enforce monotonicity
    # https://github.com/benjjneb/dada2/issues/791
    estorig <- est
    est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                                . >= X40 ~ .))) %>% as.matrix()
    rownames(est) <- rownames(estorig)
    colnames(est) <- colnames(estorig)

    # Expand the err matrix with the self-transition probs
    err <- rbind(1-colSums(est[1:3,]), est[1:3,],
                est[4,], 1-colSums(est[4:6,]), est[5:6,],
                est[7:8,], 1-colSums(est[7:9,]), est[9,],
                est[10:12,], 1-colSums(est[10:12,]))
    rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
    colnames(err) <- colnames(trans)
    # Return
    return(err)
}


learn_error_rates <- function(filtpathF,filtpathR, pattern_FWD_read, pattern_REV_read, nbases, multithread){
    
    # Create output directory if it does not exist
    parent_dir <- dirname(dirname(filtpathF))
    error_rates_dir <- file.path(parent_dir, "error_rates")
    dir.create(error_rates_dir)

    # File parsing
    filtFs <- list.files(filtpathF, pattern = pattern_FWD_read, full.names = TRUE)
    filtRs <- list.files(filtpathR, pattern = pattern_REV_read, full.names = TRUE)

    # Sample names in order
    sample_names <- basename(filtFs) # doesn't drop fq.gz
    sample_names <- gsub(pattern_FWD_read, "", sample_names)
    sample_namesR <- basename(filtRs) # doesn't drop fq.gz 
    sample_namesR <- gsub(pattern_REV_read, "", sample_namesR)

    # Double check
    if(!identical(sample_names, sample_namesR)) stop("Forward and reverse files do not match.")
    names(filtFs) <- sample_names
    names(filtRs) <- sample_names

    # Learn errors for forward and reverse reads
    errF <- learnErrors(
        filtFs,
        multithread = multithread,
        nbases = nbases,
        errorEstimationFunction = loessErrfun,
        verbose = TRUE
    )

    errR <- learnErrors(
        filtRs,
        multithread = multithread,
        nbases = nbases,
        errorEstimationFunction = loessErrfun,
        verbose = TRUE
    )


    # Save error rate objects for downstream analysis
    saveRDS(errF, paste0(error_rates_dir, "/errF.rds"))
    saveRDS(errR, paste0(error_rates_dir, "/errR.rds"))

    # Save error rate plots
    errF_plot <- plotErrors(errF, nominalQ = TRUE)
    errR_plot <- plotErrors(errR, nominalQ = TRUE)

    saveRDS(errF_plot, paste0(error_rates_dir, "/errF_plot.rds"))
    saveRDS(errR_plot, paste0(error_rates_dir, "/errR_plot.rds"))

    ggsave(plot = errF_plot, filename = paste0(error_rates_dir, "/errF_plot.png"), 
            width = 10, height = 10, dpi = "retina")
    ggsave(plot = errR_plot, filename = paste0(error_rates_dir, "/errR_plot.png"), 
            width = 10, height = 10, dpi = "retina")
}

##################################################
# Main function to handle command-line arguments #
##################################################
main <- function() {
    # Define options
    option_list <- list(
        make_option(c("-F", "--filtpathF"), type = "character", help = "Path to forward filtered files"),
        make_option(c("-R", "--filtpathR"), type = "character", help = "Path to reverse filtered files"),
        make_option(c("-f", "--pattern_FWD_read"), type = "character", default = "_1.fq.gz", help = "Pattern for forward read"),
        make_option(c("-r", "--pattern_REV_read"), type = "character", default = "_2.fq.gz", help = "Pattern for reverse read"),
        make_option(c("-n", "--nbases"), type = "numeric", default = 1e8, help = "Number of bases to use for error learning"),
        make_option(c("-m", "--multithread"), type = "logical", default = TRUE, help = "Use multiple threads")
    )

    # Parse options
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)

    # Check if all required arguments are provided
    required_args <- c("filtpathF", "filtpathR", "pattern_FWD_read", "pattern_REV_read")
    missing_args <- setdiff(required_args, names(opt[!sapply(opt, is.null)]))

    if (length(missing_args) > 0) {
        stop(sprintf("Missing required arguments: %s", paste(missing_args, collapse = ", ")), call. = FALSE)
    }

    # Main workflow
    learn_error_rates(opt$filtpathF, opt$filtpathR, opt$pattern_FWD_read, opt$pattern_REV_read, opt$nbases, opt$multithread)
}

# Execute the main function if the script is run directly
if (!interactive()) {
    main()
}