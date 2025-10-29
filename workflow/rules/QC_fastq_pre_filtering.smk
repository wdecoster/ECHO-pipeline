# rules/QC_fastq_pre_filtering.smk

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

