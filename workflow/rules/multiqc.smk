##################
# MultiQC report #
##################
rule multiqc_raw:
    input:
        fastqc_read1    = expand(rules.FastQC_raw.output.zip_file_1, sample = SAMPLES),
        fastqc_read2    = expand(rules.FastQC_raw.output.zip_file_2, sample = SAMPLES)
    output:
        report          = RESULT_DIR + config["refs"]["experiment"] + "/MultiQC/raw/multiqc_report.html"
    params:
        fastqc_results  = RESULT_DIR + config["refs"]["experiment"] + "/FastQC/raw/",
        outdir          = directory(RESULT_DIR + config["refs"]["experiment"] + "/MultiQC/raw/")
    conda: 
        "multiqc_env"
    message: "Summarising reports with multiqc"
    shell:
        "mkdir -p {params.outdir}; "
        "multiqc --force "
        "--outdir {params.outdir} "
        "{params.fastqc_results}"

rule multiqc_filtered:
    input:
        fastqc_read1    = expand(rules.FastQC_filtered.output.zip_file_1, sample = SAMPLES),
        fastqc_read2    = expand(rules.FastQC_filtered.output.zip_file_2, sample = SAMPLES)
    output:
        report          = RESULT_DIR + config["refs"]["experiment"] + "/MultiQC/filtered/multiqc_report.html"
    params:
        fastqc_results  = RESULT_DIR + config["refs"]["experiment"] + "/FastQC/filtered/",
        outdir          = directory(RESULT_DIR + config["refs"]["experiment"] + "/MultiQC/filtered/")
    conda: 
        "multiqc_env"
    message: "Summarising reports with multiqc"
    shell:
        "mkdir -p {params.outdir}; "
        "multiqc --force "
        "--outdir {params.outdir} "
        "--interactive "
        "{params.fastqc_results}"