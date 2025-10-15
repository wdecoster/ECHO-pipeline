# rules/longphase_haplotag.smk

# phasing: haplotagging
rule longphase_haplotag:
    input:
        aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam",
        bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam.bai",
        SNV_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased.vcf.gz",
        SV_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_SV.vcf.gz",
        modcall_vcf_phased=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_mod.vcf",
        reference=ancient(REFERENCE)
    params:
        out_prefix=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}"
    output:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/phasing/{{sample}}_haplotagging.log"
    threads: 16
    singularity:
        "docker://quay.io/biocontainers/longphase:2.0--h13024bc_0"
    shell:
        """
        (
            longphase haplotag \
                -b {input.aligned_bam} \
                -r {input.reference} \
                -s {input.SNV_vcf} \
                --sv-file {input.SV_vcf} \
                --mod-file {input.modcall_vcf_phased} \
                -o {params.out_prefix}_phased_alignment
        ) > {log} 2>&1
        """
