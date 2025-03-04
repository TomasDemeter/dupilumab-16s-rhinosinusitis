############################################
# Running FastQC on raw and filtered reads #
############################################
rule FastQC_raw:
    input:
        raw_read_1          = WORKING_DIR + "raw_reads/" + config["refs"]["experiment"] + "/{sample}" + config["Dada2_global"]["pattern_FWD_read"],
        raw_read_2          = WORKING_DIR + "raw_reads/" + config["refs"]["experiment"] + "/{sample}" + config["Dada2_global"]["pattern_REV_read"]
    output:
        html_1              = RESULT_DIR + config["refs"]["experiment"] + "/FastQC/raw/{sample}_1_fastqc.html",
        zip_file_1          = RESULT_DIR + config["refs"]["experiment"] + "/FastQC/raw/{sample}_1_fastqc.zip",
        html_2              = RESULT_DIR + config["refs"]["experiment"] + "/FastQC/raw/{sample}_2_fastqc.html",
        zip_file_2          = RESULT_DIR + config["refs"]["experiment"] + "/FastQC/raw/{sample}_2_fastqc.zip"        
    params:
        output_dir          = directory(RESULT_DIR + config["refs"]["experiment"] + "/FastQC/raw/")
    message:
        "Running FastQC on raw reads "
    conda: 
        "multiqc_env"
    shell:
        "mkdir -p {params.output_dir}; "
        "fastqc "
        "{input.raw_read_1} {input.raw_read_2} "
        "--outdir {params.output_dir} "
        "--threads {resources.cpus_per_task}"

rule FastQC_filtered:
    input:
        dada2_filtered_done = rules.Dada2_filter_and_trim.output.filtered_out_stats
    output:
        html_1              = RESULT_DIR + config["refs"]["experiment"] + "/FastQC/filtered/{sample}_1_fastqc.html",
        zip_file_1          = RESULT_DIR + config["refs"]["experiment"] + "/FastQC/filtered/{sample}_1_fastqc.zip",
        html_2              = RESULT_DIR + config["refs"]["experiment"] + "/FastQC/filtered/{sample}_2_fastqc.html",
        zip_file_2          = RESULT_DIR + config["refs"]["experiment"] + "/FastQC/filtered/{sample}_2_fastqc.zip"        
    params:
        filtered_read_1     = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/01_preprocess/trimmed/{sample}" + config["Dada2_global"]["pattern_FWD_read"],
        filtered_read_2     = RESULT_DIR + config["refs"]["experiment"] + "/Dada2/01_preprocess/trimmed/{sample}" + config["Dada2_global"]["pattern_REV_read"],
        output_dir          = directory(RESULT_DIR + config["refs"]["experiment"] + "/FastQC/filtered/")
    message:
        "Running FastQC on filtered reads "
    conda: 
        "fastqc_env"
    shell:
        "mkdir -p {params.output_dir}; "
        "fastqc "
        "{params.filtered_read_1} {params.filtered_read_2} "
        "--outdir {params.output_dir} "
        "--threads {resources.cpus_per_task}"
