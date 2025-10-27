# rules/longphase_phase.smk

# phasing: co-phasing SNP, SV and methylation
rule longphase_phase:
    input:
        aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam",
        bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam.bai",
        snp_vcf_gz=f"{OUTPUT_DIR}/02_variant_calling/SNVs_Indels/{{sample}}/{REFERENCE_NAME}/phased_merge_output.vcf.gz",
        snp_vcf_index=f"{OUTPUT_DIR}/02_variant_calling/SNVs_Indels/{{sample}}/{REFERENCE_NAME}/phased_merge_output.vcf.gz.tbi",
        sv_vcf_gz=f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_SV_unphased.vcf.gz",
        sv_vcf_index=f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_SV_unphased.vcf.gz.tbi",
        modcall_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_modcall.vcf",
        reference=ancient(REFERENCE)
    params:
        out_prefix=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}"
    output:
        SNV_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased.vcf",
        SV_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_SV.vcf",
        modcall_vcf_phased=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_mod.vcf"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/phasing/{{sample}}_phase.log"
    threads: 32
    singularity:
        "docker://quay.io/biocontainers/longphase:1.7.3--hf5e1c6e_0"
    shell:
        """
        (
            longphase phase \
                --snp-file {input.snp_vcf_gz} \
                --mod-file {params.out_prefix}_modcall.vcf \
                --sv-file {input.sv_vcf_gz} \
                -b {input.aligned_bam} \
                -r {input.reference} \
                -o {params.out_prefix}_phased \
                -t {threads} \
                --ont

        ) > {log} 2>&1
        """
