########################
# removes N from reads #
########################
rule Dada2_removeNs:
    input:
        raw_read_1          = expand(WORKING_DIR + "raw_reads/" + config["refs"]["experiment"] + "/{sample}" + config["Dada2_global"]["pattern_FWD_read"], sample = SAMPLES),
        raw_read_2          = expand(WORKING_DIR + "raw_reads/" + config["refs"]["experiment"] + "/{sample}" + config["Dada2_global"]["pattern_REV_read"], sample = SAMPLES)
    output:
        output_read_1       = expand(RESULT_DIR + config["refs"]["experiment"] + "/Dada2/01_preprocess/filtN/{sample}" + config["Dada2_global"]["pattern_FWD_read"], sample = SAMPLES),
        output_read_2       = expand(RESULT_DIR + config["refs"]["experiment"] + "/Dada2/01_preprocess/filtN/{sample}" + config["Dada2_global"]["pattern_REV_read"], sample = SAMPLES)
    params:
        input_path          = WORKING_DIR + "raw_reads/" + config["refs"]["experiment"],
        output_path         = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/",
        pattern_FWD_read    = config["Dada2_global"]["pattern_FWD_read"],
        pattern_REV_read    = config["Dada2_global"]["pattern_REV_read"],
        multithread         = config["Dada2_global"]["multithread"],
        maxN                = config["Dada2_removeNs"]["maxN"]
    conda:
        "dada2_env"
    message:
        "Removing Ns from reads for all files"
    shell:
        "Rscript scripts/01_dada2_removeNs.R "
        "--input_path {params.input_path} "
        "--output_path {params.output_path} "
        "--pattern_FWD_read {params.pattern_FWD_read} "
        "--pattern_REV_read {params.pattern_REV_read} "
        "--maxN {params.maxN} "
        "--multithread {params.multithread}"

##############################
# removes primers from reads #
##############################

rule Dada2_remove_primers:
    input:
        read_1              = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/01_preprocess/filtN/{sample}" + config["Dada2_global"]["pattern_FWD_read"],
        read_2              = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/01_preprocess/filtN/{sample}" + config["Dada2_global"]["pattern_REV_read"]
    output:
        output_read_1       = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/01_preprocess/trimmed/{sample}" + config["Dada2_global"]["pattern_FWD_read"],
        output_read_2       = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/01_preprocess/trimmed/{sample}" + config["Dada2_global"]["pattern_REV_read"]
    params:
        forward_primer_seq  = config["Dada2_remove_primers"]["forward_primer_seq"],
        reverse_primer_seq  = config["Dada2_remove_primers"]["reverse_primer_seq"],
        minimum_length      = config["Dada2_remove_primers"]["minimum_length"],
        nextseq_trim        = config["Dada2_remove_primers"]["nextseq_trim"]
    conda:
        "dada2_env"
    message:
        "Removing primers from reads"
    shell:
        "Rscript scripts/02_dada2_remove_primers.R "
        "--input_path_fwd {input.read_1} "
        "--input_path_rev {input.read_2} "
        "--output_path_fwd {output.output_read_1} "
        "--output_path_rev {output.output_read_2} "
        "--forward_primer_seq {params.forward_primer_seq} "
        "--reverse_primer_seq {params.reverse_primer_seq} "
        "--minimum_length {params.minimum_length} "
        "--nextseq_trim {params.nextseq_trim}"

rule detect_primers:
    input:
        before_removal_FWD  = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/01_preprocess/filtN/{sample}" + config["Dada2_global"]["pattern_FWD_read"],
        before_removal_REV  = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/01_preprocess/filtN/{sample}" + config["Dada2_global"]["pattern_REV_read"],
        after_removal_FWD   = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/01_preprocess/trimmed/{sample}" + config["Dada2_global"]["pattern_FWD_read"],
        after_removal_REV   = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/01_preprocess/trimmed/{sample}" + config["Dada2_global"]["pattern_REV_read"]
    output:
        txt_file            = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/01_preprocess/stats/{sample}.txt"
    params:
        output_dir          = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/01_preprocess/stats/",
        forward_primer_seq  = config["Dada2_remove_primers"]["forward_primer_seq"],
        reverse_primer_seq  = config["Dada2_remove_primers"]["reverse_primer_seq"],
    conda:
        "dada2_env"
    message:
        "Checking primer contamination before and after primer removal step"
    shell:
        "mkdir -p {params.output_dir}; "
        "Rscript scripts/primer_detection.R "
        "--FWD_before_removal {input.before_removal_FWD} "
        "--REV_before_removal {input.before_removal_REV} "
        "--FWD_after_removal {input.after_removal_FWD} "
        "--REV_after_removal {input.after_removal_REV} "
        "--forward_primer_seq {params.forward_primer_seq} "
        "--reverse_primer_seq {params.reverse_primer_seq} "
        "--output_path {output.txt_file}"