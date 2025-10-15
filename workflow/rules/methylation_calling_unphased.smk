# rules/methylation_calling_unphased.smk

# methylation calling
rule methylation_calling_unphased:
    input:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam.bai",
        reference=ancient(REFERENCE)
    params:
        unphased_bed=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/unphased/{{sample}}_{REFERENCE_NAME}_unphased.bed"
    output:
        unphased_bed=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/unphased/{{sample}}_{REFERENCE_NAME}_unphased.bed"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/methylation_calling/{{sample}}_pileup_unphased.log"
    threads: 16
    singularity:
        "docker://quay.io/biocontainers/ont-modkit:0.5.0--hcdda2d0_2"
    shell:
        """
        (
            modkit pileup {input.phased_bam} {params.unphased_bed} \
                --threads {threads} \
                --ref {input.reference} \
                --combine-strands \
                --cpg
        ) > {log} 2>&1
        """
