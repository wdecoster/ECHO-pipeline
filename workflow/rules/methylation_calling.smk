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
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/methylation_calling_phased/{{sample}}_methylation_calling_phased.tsv"
    shell:
        """
        (
            modkit pileup {input.phased_bam} {params.out_dir}/phased/ \
                --threads {threads} \
                --ref {input.reference} \
                --ignore h \
                --combine-strands \
                --cpg \
                --partition-tag HP \
                --prefix {params.sample_name}
        ) > {log} 2>&1
        """

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
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/methylation_calling_phased/{{sample}}_methylation_calling_unphased.tsv"
    shell:
        """
        (
            modkit pileup {input.phased_bam} {params.unphased_bed} \
                --threads {threads} \
                --ref {input.reference} \
                --ignore h \
                --combine-strands \
                --cpg
        ) > {log} 2>&1
        """

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
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/methylation_calling_phased/{{sample}}_methylation_bgzip_index.tsv"
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



