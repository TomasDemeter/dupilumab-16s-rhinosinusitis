rule Dada2_to_phyloseq:
    input:
        seqtab_final    = rules.Dada2_chimeras_taxonomy.output.seqtab_final,
        taxonomy        = rules.Dada2_chimeras_taxonomy.output.taxonomy,
        metadata        = config["refs"]["metadata"]
    output:
        phyloseq        = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/04_phyloseq/" + config["refs"]["experiment"] + "_phyloseq.rds"
    params:
        output_path     = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/04_phyloseq/",
        sample_names    = config["Dada2_to_phyloseq"]["sample_names"],
        experiment      = config["refs"]["experiment"]
    conda:
        "dada2_env"
    message:
        "Creating phyloseq object"
    shell:
        "Rscript scripts/07_dada2_to_phyloseq.R "
        "--counts {input.seqtab_final} "
        "--taxa {input.taxonomy} "
        "--metadata {input.metadata} "
        "--sample_names {params.sample_names} "
        "--experiment {params.experiment} "
        "--output_path {params.output_path}"

rule Dada2_loss_plotting:
    input:
        filtered_out_stats          = rules.Dada2_filter_and_trim.output.filtered_out_stats,
        ddF                         = rules.Dada2_sequence_table.output.ddF,
        ddR                         = rules.Dada2_sequence_table.output.ddR,
        mergers                     = rules.Dada2_sequence_table.output.mergers,
        seqtab_final                = rules.Dada2_chimeras_taxonomy.output.seqtab_final
    output:
        tracking_reads_summary_plot = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/tracking_reads_summary_plot.png"
    params:
        pattern_FWD_read            = config["Dada2_global"]["pattern_FWD_read"],
        pattern_REV_read            = config["Dada2_global"]["pattern_REV_read"],
        output_path                 = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/"
    conda:
        "dada2_env"
    message:
        "Creating tracking reads summary plot"
    shell:
        "Rscript scripts/08_dada2_loss_plotting.R "
        "--filtering_output {input.filtered_out_stats} "
        "--ddF {input.ddF} "
        "--ddR {input.ddR} "
        "--mergers {input.mergers} "
        "--seqtab_final {input.seqtab_final} "
        "--pattern_FWD_read {params.pattern_FWD_read} "
        "--pattern_REV_read {params.pattern_REV_read} "
        "--output_path {params.output_path}"
    