# rules/VCF_filtering.smk

# filter SNPs based on PASS flag and at least 5 reads used for variant calling
rule filter_snps:
    input:
        SNV_vcf_gz=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased.vcf.gz"
    output:
        SNV_filt_vcf_gz=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/filt/{{sample}}_{REFERENCE_NAME}_SNP_filt.vcf.gz",
        SNV_filt_vcf_gz_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/filt/{{sample}}_{REFERENCE_NAME}_SNP_filt.vcf.gz.tbi"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/VCF_filtering_snps/{{sample}}.log"
    threads: 4
    singularity:
        "docker://leenaputzeys/te_methylation:v1.0"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/VCF_filtering_snps/{{sample}}.tsv"
    shell:
        r"""
        mkdir -p $(dirname {output.SNV_filt_vcf_gz})
        (
            bcftools view {input.SNV_vcf_gz} \
                -i 'FILTER="PASS" & FORMAT/DP>5' \
                --threads {threads} \
                -O z -o {output.SNV_filt_vcf_gz}
            tabix {output.SNV_filt_vcf_gz}
        ) > {log} 2>&1
        """

rule filter_SVs:
    input:
        SV_vcf_gz=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_SV.vcf.gz"
    output:
        SV_filt_vcf_gz=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/filt/{{sample}}_{REFERENCE_NAME}_SV_filt.vcf.gz",
        SV_filt_vcf_gz_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/filt/{{sample}}_{REFERENCE_NAME}_SV_filt.vcf.gz.tbi"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/VCF_filtering_svs/{{sample}}.log"
    threads: 4
    singularity:
        "docker://leenaputzeys/te_methylation:v1.0"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/VCF_filtering_svs/{{sample}}.tsv"
    shell:
        r"""
        mkdir -p $(dirname {output.SV_filt_vcf_gz})
        (
            bcftools view {input.SV_vcf_gz} \
                -i 'FILTER="PASS" & INFO/SUPPORT>5' \
                --threads {threads} \
                -O z -o {output.SV_filt_vcf_gz}
            tabix {output.SV_filt_vcf_gz}
        ) > {log} 2>&1
        """

