#!/usr/bin/env Rscript

# Rscript 08_dada2_loss_plotting.R -F ../../results/dada2/02_filter/filt_out.rds -D ../../results/dada2/03_tabletax/ddF.rds -R ../../results/dada2/03_tabletax/ddR.rds -M ../../results/dada2/03_tabletax/mergers.rds -S ../../results/dada2/03_tabletax/seqtab_final.rds -o ../../results/dada2/

###########################
# Load required libraries #
###########################
library(dada2)
library(tidyverse)
library(tibble)
library(phyloseq)
library(ShortRead)
library(optparse)
library(Hmisc)


#############
# Functions # 
#############
getN <- function(x) sum(getUniques(x)) # function to grab sequence counts from output objects


generate_track_df <- function(filtering_output, ddF, ddR, mergers, seqtab_nochim, output_path, pattern_FWD_read, pattern_REV_read){

    # Load RDS files
    filtering_output <- readRDS(filtering_output)
    ddF <- readRDS(ddF)
    ddR <- readRDS(ddR)
    mergers <- readRDS(mergers)
    seqtab_nochim <- readRDS(seqtab_nochim)

    sample_names <- filtering_output %>%
        data.frame() %>%
        rownames(.) %>%
        gsub(pattern_FWD_read, "", .)

    filt_out_track <- filtering_output %>%
        data.frame() %>%
        dplyr::mutate(Sample = gsub(pattern_FWD_read, "",rownames(.))) %>%
        dplyr::rename(input = reads.in, filtered = reads.out)

    rownames(filt_out_track) <- filt_out_track$Sample


    ddF_track <- data.frame(denoisedF = sapply(ddF[sample_names], getN)) %>%
        mutate(Sample = row.names(.))
    ddR_track <- data.frame(denoisedR = sapply(ddR[sample_names], getN)) %>%
        mutate(Sample = row.names(.))
    merge_track <- data.frame(merged = sapply(mergers, getN)) %>%
        mutate(Sample = row.names(.))
    chim_track <- data.frame(nonchim = rowSums(seqtab_nochim)) %>%
        mutate(Sample = row.names(.))


    track <- left_join(filt_out_track, ddF_track, by = "Sample") %>%
        left_join(ddR_track, by = "Sample") %>%
        left_join(merge_track, by = "Sample") %>%
        left_join(chim_track, by = "Sample") %>%
        replace(., is.na(.), 0) %>%
        select(Sample, everything())
    
    row.names(track) <- track$Sample

    track_pct <- track %>% 
    data.frame() %>%
    mutate(Sample = rownames(.),
         filtered_pct = ifelse(filtered == 0, 0, 100 * (filtered/input)),
         denoisedF_pct = ifelse(denoisedF == 0, 0, 100 * (denoisedF/filtered)),
         denoisedR_pct = ifelse(denoisedR == 0, 0, 100 * (denoisedR/filtered)),
         merged_pct = ifelse(merged == 0, 0, 100 * merged/((denoisedF + denoisedR)/2)),
         nonchim_pct = ifelse(nonchim == 0, 0, 100 * (nonchim/merged)),
         total_pct = ifelse(nonchim == 0, 0, 100 * nonchim/input)) %>%
    select(Sample, ends_with("_pct"))

    # summary stats of tracked reads averaged across samples
    track_pct_avg <- track_pct %>% summarize_at(vars(ends_with("_pct")), 
                                                list(avg = mean))

    track_pct_med <- track_pct %>% summarize_at(vars(ends_with("_pct")), 
                                                list(avg = stats::median))

    
    saveRDS(track, paste0(output_path, "/tracking_reads.rds"))
    saveRDS(track_pct, paste0(output_path, "/tracking_reads_percentage.rds"))

    return(list(track = track, track_pct_avg = track_pct_avg, track_pct_med = track_pct_med))
}



generate_track_plot <- function(track_list, output_path) {

    # Extracting the data from the list
    track <- track_list$track
    track_pct_avg <- track_list$track_pct_avg
    track_pct_med <- track_list$track_pct_med

    # Plotting each sample's reads through the pipeline
    track_plot <- track %>% 
    data.frame() %>%
    mutate(Sample = rownames(.)) %>%
    gather(key = "Step", value = "Reads", -Sample) %>%
    mutate(Step = factor(Step, 
                            levels = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim"))) %>%
    ggplot(aes(x = Step, y = Reads)) +
    geom_line(aes(group = Sample), alpha = 0.2) +
    geom_point(alpha = 0.5, position = position_jitter(width = 0)) + 
    stat_summary(fun.y = median, geom = "line", group = 1, color = "steelblue", size = 1, alpha = 0.5) +
    stat_summary(fun.y = median, geom = "point", group = 1, color = "steelblue", size = 2, alpha = 0.5) +
    stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5), 
                    geom = "ribbon", group = 1, fill = "steelblue", alpha = 0.2) +
    geom_label(data = t(track_pct_avg[1:5]) %>% data.frame() %>% 
                    dplyr::rename(Percent = 1) %>%
                    mutate(Step = c("filtered", "denoisedF", "denoisedR", "merged", "nonchim"),
                        Percent = paste(round(Percent, 2), "%")),
             aes(label = Percent), y = 1.1 * max(track[,2])) +
    geom_label(data = track_pct_avg[6] %>% data.frame() %>%
                dplyr::rename(total = 1),
                aes(label = paste("Total\nRemaining:\n", round(track_pct_avg[1,6], 2), "%")), 
                y = mean(track[,6]), x = 6.5) +
    expand_limits(y = 1.1 * max(track[,2]), x = 7) +
    theme_classic()

    saveRDS(track_plot, paste0(output_path, "/tracking_reads_summary_plot.rds"))
    ggsave(plot = track_plot, filename = paste0(output_path, "/tracking_reads_summary_plot.png"), width = 15, height = 15, dpi = "retina", bg = "white")
}


##################################################
# Main function to handle command-line arguments #
##################################################
main <- function() {
    # Define options
    option_list <- list(
        make_option(c("-F", "--filtering_output"), type = "character", help = "filt_out.rds from filter_and_trim.R"),
        make_option(c("-D", "--ddF"), type = "character", help = "ddF.rds from dada2_sequence_table.R"), 
        make_option(c("-R", "--ddR"), type = "character", help = "ddR.rds from dada2_sequence_table.R"), 
        make_option(c("-M", "--mergers"), type = "character", help = "mergers.rds from dada2_sequence_table.R"),
        make_option(c("-S", "--seqtab_final"), type = "character", help = "seqtab_final.rds from dada2_sequence_table.R"),
        make_option(c("-o", "--output_path"), type = "character", help = "Path for output files"),
        make_option(c("-f", "--pattern_FWD_read"), type = "character", default = "_1.fq.gz", help = "Pattern for forward read"),
        make_option(c("-r", "--pattern_REV_read"), type = "character", default = "_2.fq.gz", help = "Pattern for reverse read")
    )

    # Parse options
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)

    # Check if all required arguments are provided
    required_args <- c("filtering_output", "ddF", "ddR", "mergers", "seqtab_final", "output_path")
    missing_args <- setdiff(required_args, names(opt[!sapply(opt, is.null)]))

    if (length(missing_args) > 0) {
        stop(sprintf("Missing required arguments: %s", paste(missing_args, collapse = ", ")), call. = FALSE)
    }

    # Main workflow
    track_list <- generate_track_df(opt$filtering_output, opt$ddF, opt$ddR, opt$mergers, opt$seqtab_final, opt$output_path, opt$pattern_FWD_read, opt$pattern_REV_read)
    generate_track_plot(track_list, opt$output_path)
}

# Execute the main function if the script is run directly
if (!interactive()) {
    main()
}