# rules/sv_vcf_compress_index

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
    shell:
        """
        (
            bgzip -@ 8 -c {input.sv_vcf} > {output.sv_vcf_gz}
            tabix -f -p vcf {output.sv_vcf_gz}
        ) > {log} 2>&1
        """
