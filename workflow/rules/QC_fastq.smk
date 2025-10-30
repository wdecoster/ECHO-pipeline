# rules/fastq_QC.smk

rule QC_fastq_pre_filtering:
    input:
        fastq=raw_fastq
    output:
        summary_stats=f"{OUTPUT_DIR}/qc/fastq/pre_filtering/{{sample}}/{{sample}}_fastq_NanoStats.txt",
        html_report=f"{OUTPUT_DIR}/qc/fastq/pre_filtering/{{sample}}/{{sample}}_fastq_NanoPlot-report.html"
    params:
        results_dir=f"{OUTPUT_DIR}/qc/fastq/pre_filtering/{{sample}}/",
        sample_name=f"{{sample}}_fastq_"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/QC_fastq_pre_filtering/{{sample}}.log"
    threads: 48
    singularity:
       "docker://quay.io/biocontainers/nanoplot:1.46.1--pyhdfd78af_0"
    shell:
        """
        NanoPlot -t {threads} \
            -o {params.results_dir} \
            --dpi 250 -c blue -f png --N50 \
            --title {params.sample_name} \
            --prefix {params.sample_name} \
            --fastq {input.fastq} \
            > {log} 2>&1
        """


# filtering reads based on length and quality in fastq files
if DO_FILTER:
    rule filter_fastq:
        input:
            fastq=raw_fastq
        output:
            filtered_fastq=f"{OUTPUT_DIR}/00_raw_data/basecalled/fastq/{{sample}}/{{sample}}_filtered.fastq"
        log:
            f"{OUTPUT_DIR}/logs/snakemake_rules/filter_fastq/{{sample}}.log"
        threads: 24
        singularity:
            "docker://quay.io/biocontainers/chopper:0.11.0--hcdda2d0_0"
        shell:
            """
            chopper -q {MIN_READ_QUAL} -l {MIN_READ_LENGTH} \
                    -i {input.fastq} -t {threads} > {output.filtered_fastq} 2> {log}
            """

    rule QC_fastq_post_filtering:
        input:
            fastq=f"{OUTPUT_DIR}/00_raw_data/basecalled/fastq/{{sample}}/{{sample}}_filtered.fastq"
        output:
            summary_stats=f"{OUTPUT_DIR}/qc/fastq/post_filtering/{{sample}}/{{sample}}_filtered_fastq_NanoStats.txt",
            html_report=f"{OUTPUT_DIR}/qc/fastq/post_filtering/{{sample}}/{{sample}}_filtered_fastq_NanoPlot-report.html"
        params:
            results_dir=f"{OUTPUT_DIR}/qc/fastq/post_filtering/{{sample}}/",
            sample_name=f"{{sample}}_filtered_fastq_"
        log:
            f"{OUTPUT_DIR}/logs/snakemake_rules/QC_fastq_post_filtering/{{sample}}.log"
        threads: 48
        singularity:
            "docker://quay.io/biocontainers/nanoplot:1.46.1--pyhdfd78af_0"
        shell:
            """
            NanoPlot -t {threads} \
                -o {params.results_dir} \
                --dpi 250 -c blue -f png --N50 \
                --title {params.sample_name} \
                --prefix {params.sample_name} \
                --fastq {input.fastq} \
                > {log} 2>&1
            """



