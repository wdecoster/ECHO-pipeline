# rules/TR_calling.smk

# tandem repeats calling using LongTR
rule TR_calling:
    input:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam.bai",
        TR_catalog=ancient(TR_CATALOG),
        reference=ancient(REFERENCE)
    output:
        TR_vcf=f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}.vcf.gz"
    params:
        haploid_chrs=HAPLOID_CHRS,
        type_of_tr=TYPE_OF_TR,
        sample_name=f"{{sample}}_{REFERENCE_NAME}"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/TR_calling/{{sample}}.log"
    singularity:
        "docker://quay.io/biocontainers/longtr:1.2--h077b44d_1"
    shell:
        """
        LongTR \
            --bams {input.phased_bam} \
            --regions {input.TR_catalog} \
            --fasta {input.reference} \
            --tr-vcf  {output.TR_vcf}\
            --bam-samps {params.sample_name} \
            --bam-libs {params.sample_name} \
            --haploid {params.haploid_chrs} \
            --phased-bam \
            > {log} 2>&1
        """


