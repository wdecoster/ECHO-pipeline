# multiQC

if DO_FILTER:
    rule multiqc:
        input:
            fastq_pre_report=f"{OUTPUT_DIR}/qc/fastq/pre_filtering/{{sample}}/{{sample}}_fastq_NanoPlot-report.html",
            fastq_post_report=f"{OUTPUT_DIR}/qc/fastq/post_filtering/{{sample}}/{{sample}}_filtered_fastq_NanoPlot-report.html", 
            phased_bam_report=f"{OUTPUT_DIR}/qc/phased_bam/{{sample}}/{REFERENCE_NAME}/nanoplot/{{sample}}_{REFERENCE_NAME}_phased_bam_NanoPlot-report.html"
        output:
            multiqc_out=f"{OUTPUT_DIR}/qc/multiqc/{{sample}}/{REFERENCE_NAME}/multiqc_report.html"
        singularity:
            "docker://quay.io/biocontainers/multiqc:1.31--pyhdfd78af_0"
        log:
            f"{OUTPUT_DIR}/logs/snakemake_rules/multiqc/{{sample}}.log"
        shell:
            """
            outdir=$(dirname {output.multiqc_out})
            mkdir -p "$outdir"

            multiqc \
            --force \
            -n multiqc_report.html \
            -o "$outdir" \
            {OUTPUT_DIR}/qc/fastq/pre_filtering/{wildcards.sample} \
            {OUTPUT_DIR}/qc/fastq/post_filtering/{wildcards.sample} \
            {OUTPUT_DIR}/qc/phased_bam/{wildcards.sample}/{REFERENCE_NAME}/nanoplot \
            > {log} 2>&1
            """
else:
    rule multiqc:
        input:
            fastq_pre_report=f"{OUTPUT_DIR}/qc/fastq/pre_filtering/{{sample}}/{{sample}}_fastq_NanoPlot-report.html",
            phased_bam_report=f"{OUTPUT_DIR}/qc/phased_bam/{{sample}}/{REFERENCE_NAME}/nanoplot/{{sample}}_{REFERENCE_NAME}_phased_bam_NanoPlot-report.html"
        output:
            multiqc_out=f"{OUTPUT_DIR}/qc/multiqc/{{sample}}/{REFERENCE_NAME}/multiqc_report.html"
        singularity:
            "docker://quay.io/biocontainers/multiqc:1.31--pyhdfd78af_0"
        log:
            f"{OUTPUT_DIR}/logs/snakemake_rules/multiqc/{{sample}}.log"
        shell:
            """
            outdir=$(dirname {output.multiqc_out})
            mkdir -p "$outdir"

            multiqc \
            --force \
            -n multiqc_report.html \
            -o "$outdir" \
            {OUTPUT_DIR}/qc/fastq/pre_filtering/{wildcards.sample} \
            {OUTPUT_DIR}/qc/phased_bam/{wildcards.sample}/{REFERENCE_NAME}/nanoplot \
            > {log} 2>&1
            """

