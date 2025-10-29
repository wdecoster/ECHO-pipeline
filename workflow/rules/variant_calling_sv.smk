# rules/variant_calling_sv.smk

# variant calling: structural variants
rule variant_calling_sv:
    input:
        aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam",
        bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam.bai",
        reference=ancient(REFERENCE)
    output:
        sv_vcf=temp(f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_SV_unphased.vcf")
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/variant_calling_sv/{{sample}}.log"
    threads: 16
    singularity:
        "docker://quay.io/biocontainers/sniffles:2.6.3--pyhdfd78af_0"
    shell:
        """
        sniffles \
            --input {input.aligned_bam} \
            --reference {input.reference} \
            --vcf {output.sv_vcf} \
            --threads {threads} \
            --output-rnames \
        > {log} 2>&1
        """

