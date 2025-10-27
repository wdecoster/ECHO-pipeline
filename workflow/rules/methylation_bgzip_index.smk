# rules/methylation_bgzip_index.smk

# bgzip and index the bedmethyl files
rule methylation_bgzip_index:
    input:
        hap1_bed=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_1.bed",
        hap2_bed=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_2.bed",
        ungrouped_bed=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_ungrouped.bed",
        unphased_bed=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/unphased/{{sample}}_{REFERENCE_NAME}_unphased.bed"
    output:
        hap1_bed_gz=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_1.bed.gz",
        hap1_bed_tbi=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_1.bed.gz.tbi",
        hap2_bed_gz=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_2.bed.gz",
        hap2_bed_tbi=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_2.bed.gz.tbi",
        ungrouped_bed_gz=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_ungrouped.bed.gz",
        ungrouped_bed_tbi=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_ungrouped.bed.gz.tbi",
        unphased_bed_gz=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/unphased/{{sample}}_{REFERENCE_NAME}_unphased.bed.gz",
        unphased_bed_tbi=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/unphased/{{sample}}_{REFERENCE_NAME}_unphased.bed.gz.tbi"
    threads: 8
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/methylation_calling/{{sample}}_bgzip_tabix.log"
    singularity:
        "docker://quay.io/biocontainers/htslib:1.22.1--h566b1c6_0"
    shell:
        """
        (
            bgzip -@ {threads} {input.hap1_bed}
            bgzip -@ {threads} {input.hap2_bed}
            bgzip -@ {threads} {input.ungrouped_bed}
            bgzip -@ {threads} {input.unphased_bed}
            tabix {output.hap1_bed_gz}
            tabix {output.hap2_bed_gz}
            tabix {output.ungrouped_bed_gz}
            tabix {output.unphased_bed_gz}
        ) > {log} 2>&1
        """
