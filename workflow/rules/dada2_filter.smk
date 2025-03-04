rule Dada2_filter_and_trim:
    input:
        output_read_1           = expand(rules.Dada2_remove_primers.output.output_read_1, sample = SAMPLES),
        output_read_2           = expand(rules.Dada2_remove_primers.output.output_read_2, sample = SAMPLES)
    output:
        filtered_out_stats      = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/02_filter/filt_out.rds"
    params:
        input_path              = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/01_preprocess/trimmed/",
        output_path             = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/" ,
        pattern_FWD_read        = config["Dada2_global"]["pattern_FWD_read"],
        pattern_REV_read        = config["Dada2_global"]["pattern_REV_read"],
        multithread             = config["Dada2_global"]["multithread"],
        trunc_length_FWD        = config["Dada2_filter_and_trim"]["trunc_length_FWD"],
        trunc_length_REV        = config["Dada2_filter_and_trim"]["trunc_length_REV"],
        maxEE_FWD               = config["Dada2_filter_and_trim"]["maxEE_FWD"],
        maxEE_REV               = config["Dada2_filter_and_trim"]["maxEE_REV"],
        truncQ                  = config["Dada2_filter_and_trim"]["truncQ"]
    conda:
        "dada2_env"
    message:
        "Filtering and trimming reads"
    shell:
        "Rscript scripts/03_dada2_filter_and_trim.R "
        "--input_path {params.input_path} "
        "--output_path {params.output_path} "
        "--pattern_FWD_read {params.pattern_FWD_read} "
        "--pattern_REV_read {params.pattern_REV_read} "
        "--trunc_length_FWD {params.trunc_length_FWD} "
        "--trunc_length_REV {params.trunc_length_REV} "
        "--maxEE_FWD {params.maxEE_FWD} "
        "--maxEE_REV {params.maxEE_REV} "
        "--truncQ {params.truncQ} "
        "--multithread {params.multithread}"

rule Dada2_learn_error_rates:
    input:
        Dada2_filter_and_trim   = rules.Dada2_filter_and_trim.output.filtered_out_stats
    output:
        errF                    = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/02_filter/error_rates/errF.rds",
        errR                    = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/02_filter/error_rates/errR.rds"
    params:
        filtpathF               = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/02_filter/preprocessed_F/filtered",
        filtpathR               = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/02_filter/preprocessed_R/filtered",
        pattern_FWD_read        = config["Dada2_global"]["pattern_FWD_read"],
        pattern_REV_read        = config["Dada2_global"]["pattern_REV_read"],
        multithread             = config["Dada2_global"]["multithread"],
        nbases                  = config["Dada2_learn_error_rates"]["nbases"]
    conda:
        "dada2_env"
    message:
        "Learning error rates"
    shell:
        "Rscript scripts/04_dada2_learn_error_rates.R "
        "--filtpathF {params.filtpathF} "
        "--filtpathR {params.filtpathR} "
        "--pattern_FWD_read {params.pattern_FWD_read} "
        "--pattern_REV_read {params.pattern_REV_read} "
        "--nbases {params.nbases} "
        "--multithread {params.multithread}"
