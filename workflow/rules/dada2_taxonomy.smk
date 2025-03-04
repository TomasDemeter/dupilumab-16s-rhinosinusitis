rule Dada2_sequence_table:
    input:
        errF                = rules.Dada2_learn_error_rates.output.errF,
        errR                = rules.Dada2_learn_error_rates.output.errR
    output:
        seqtab              = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/03_tabletax/seqtab.rds",
        ddF                 = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/03_tabletax/ddF.rds",
        ddR                 = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/03_tabletax/ddR.rds",
        mergers             = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/03_tabletax/mergers.rds"
    params:
        filtpathF           = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/02_filter/preprocessed_F/filtered",
        filtpathR           = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/02_filter/preprocessed_R/filtered",
        pattern_FWD_read    = config["Dada2_global"]["pattern_FWD_read"],
        pattern_REV_read    = config["Dada2_global"]["pattern_REV_read"],
        multithread         = config["Dada2_global"]["multithread"],
        pooling             = config["Dada2_sequence_table"]["pooling"]
    conda:
        "dada2_env"
    message:
        "Creating sequence table"
    shell:
        "Rscript scripts/05_dada2_sequence_table.R "
        "--filtpathF {params.filtpathF} "
        "--filtpathR {params.filtpathR} "
        "--pattern_FWD_read {params.pattern_FWD_read} "
        "--pattern_REV_read {params.pattern_REV_read} "
        "--pool {params.pooling} "
        "--errF {input.errF} "
        "--errR {input.errR} "
        "--multithread {params.multithread}"

rule Download_silva_training_set:
    output:
        silva_training_set      = config["Download_silva_training_set"]["training_set"]
    params:
        silva_training_set_link = config["Download_silva_training_set"]["training_set_link"]
    message:
        "Downloading Silva training set"
    shell:
        "wget -O {output.silva_training_set} {params.silva_training_set_link}"

    
rule Dada2_chimeras_taxonomy:
    input:
        seqtab              = rules.Dada2_sequence_table.output.seqtab,
        database            = rules.Download_silva_training_set.output.silva_training_set
    output:
        seqtab_final        = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/03_tabletax/seqtab_final.rds",
        taxonomy            = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/03_tabletax/tax_final.rds"
    params:
        multithread         = config["Dada2_global"]["multithread"],
        lower_cutoff        = config["Dada2_chimeras_taxonomy"]["lower_cutoff"],
        upper_cutoff        = config["Dada2_chimeras_taxonomy"]["upper_cutoff"],
        min_fold            = config["Dada2_chimeras_taxonomy"]["minFoldParentOverAbundance"]
    conda:
        "dada2_env"
    message:
        "Identifying chimeras and assigning taxonomy"
    shell:
        "Rscript scripts/06_dada2_chimeras_taxonomy.R "
        "--input_path {input.seqtab} "
        "--database {input.database} "
        "--lower_cutoff {params.lower_cutoff} "
        "--upper_cutoff {params.upper_cutoff} "
        "--minFoldParentOverAbundance {params.min_fold} "
        "--multithread {params.multithread}"
