# rules/phased_bam_index.smk

# phasing: bgzip and index files
rule phased_bam_index:
    input:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam"
    output:
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam.bai"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/phasing/{{sample}}_index_bam.log"
    threads: 16
    singularity:
        "docker://staphb/samtools:1.22"
    shell:
        """
        (
            samtools index {input.phased_bam}
        ) > {log} 2>&1
        """
