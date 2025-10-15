# rules/methylation_calling_phased.smk

# methylation calling
rule methylation_calling_phased:
    input:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam.bai",
        reference=ancient(REFERENCE)
    params:
        out_dir=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/",
        hap1_bed=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_1.bed",
        hap2_bed=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_2.bed",
        ungrouped_bed=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_ungrouped.bed",
        sample_name=f"{{sample}}_{REFERENCE_NAME}_haplotype"
    threads: 16
    output:
        hap1_bed=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_1.bed",
        hap2_bed=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_2.bed",
        ungrouped_bed=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_ungrouped.bed"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/methylation_calling/{{sample}}_pileup_phased.log"
    singularity:
        "docker://quay.io/biocontainers/ont-modkit:0.5.0--hcdda2d0_2"
    shell:
        """
        (
            modkit pileup {input.phased_bam} {params.out_dir}/phased/ \
                --threads {threads} \
                --ref {input.reference} \
                --combine-strands \
                --cpg \
                --partition-tag HP \
                --prefix {params.sample_name}
        ) > {log} 2>&1
        """
