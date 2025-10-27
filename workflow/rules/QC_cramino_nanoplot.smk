# rules/QC_cramino_nanoplot.smk

# QC of phased bam file from cramino output
rule QC_cramino_nanoplot:
    input:
        arrow=f"{OUTPUT_DIR}/qc/phased_bam/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_out.arrow"
    output:
        nanoplot_output=f"{OUTPUT_DIR}/qc/phased_bam/{{sample}}/{REFERENCE_NAME}/nanoplot/{{sample}}_{REFERENCE_NAME}_phased_bam_NanoPlot-report.html"
    params:
        output_dir=f"{OUTPUT_DIR}/qc/phased_bam/{{sample}}/{REFERENCE_NAME}",
        sample_name=f"{{sample}}_{REFERENCE_NAME}_phased_bam_"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/QC_phased_bam/{{sample}}_nanoplot.log"
    threads: 16
    singularity:
        "docker://quay.io/biocontainers/nanoplot:1.46.1--pyhdfd78af_0"
    shell:
        """
        (
            NanoPlot \
                --threads {threads} \
                --title {params.sample_name} \
                --prefix {params.sample_name} \
                --outdir {params.output_dir}/nanoplot \
                --arrow {input.arrow} \
                --no_static \
        )  > {log} 2>&1
        """
