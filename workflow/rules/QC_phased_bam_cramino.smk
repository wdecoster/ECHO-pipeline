# rules/QC_phased_bam_cramino.smk

# QC of phased bam file
rule QC_phased_bam_cramino:
    input:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam"
    output:
        cramino_output=f"{OUTPUT_DIR}/qc/phased_bam/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_cramino_output.txt",
        arrow_output=f"{OUTPUT_DIR}/qc/phased_bam/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_out.arrow"
    params:
        output_dir=f"{OUTPUT_DIR}/qc/phased_bam/{{sample}}/{REFERENCE_NAME}",
        hist_output=f"{OUTPUT_DIR}/qc/phased_bam/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_histogram.out",
        sample_name=f"{{sample}}_{REFERENCE_NAME}_phased_bam_"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/QC_phased_bam/{{sample}}_cramino.log"
    threads: 16
    singularity:
        "docker://quay.io/biocontainers/cramino:1.1.0--h3dc2dae_0"
    shell:
        """
        cramino \
            --hist {params.hist_output} \
            --arrow {output.arrow_output} \
            --karyotype \
            --phased \
            --format text \
            {input.phased_bam} > {output.cramino_output} 2> {log}
        """
