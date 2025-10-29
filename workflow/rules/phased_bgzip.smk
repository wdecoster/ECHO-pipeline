# rules/phased_bgzip.smk

# gzip and index phased VCF files
rule phased_bgzip:
    input:
        SNV_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased.vcf",
        SV_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_SV.vcf"
    params:
        out_prefix=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}"
    output:
        SNV_vcf_gz=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased.vcf.gz",
        SNV_vcf_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased.vcf.gz.tbi",
        SV_vcf_gz=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_SV.vcf.gz",
        SV_vcf_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_SV.vcf.gz.tbi",
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/phasing/{{sample}}_bgzip_phased_vcf.log"
    threads: 16
    singularity:
        "docker://quay.io/biocontainers/htslib:1.22.1--h566b1c6_0"
    shell:
        """
        (
            bgzip -@ 8 {params.out_prefix}_phased.vcf
            bgzip -@ 8 {params.out_prefix}_phased_SV.vcf
            tabix {params.out_prefix}_phased.vcf.gz
            tabix {params.out_prefix}_phased_SV.vcf.gz
        ) > {log} 2>&1
        """
