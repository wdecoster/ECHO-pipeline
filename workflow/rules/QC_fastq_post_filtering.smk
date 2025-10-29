# rules/QC_fastq_post_filtering.smk

if DO_FILTER:
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

