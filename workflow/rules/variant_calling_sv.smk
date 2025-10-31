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
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/variant_calling_sv/{{sample}}_variant_calling_sv.tsv"
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

# SV vcf compress and index
rule sv_vcf_compress_index:
    input:
        sv_vcf=f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_SV_unphased.vcf"
    output:
        sv_vcf_gz=f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_SV_unphased.vcf.gz",
        sv_vcf_index=f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_SV_unphased.vcf.gz.tbi"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/SV_vcf_compress_index/{{sample}}.log"
    singularity:
        "docker://quay.io/biocontainers/htslib:1.22.1--h566b1c6_0"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/variant_calling_sv/{{sample}}_sv_vcf_compress_index.tsv"
    shell:
        """
        (
            bgzip -@ 8 -c {input.sv_vcf} > {output.sv_vcf_gz}
            tabix -f -p vcf {output.sv_vcf_gz}
        ) > {log} 2>&1
        """

