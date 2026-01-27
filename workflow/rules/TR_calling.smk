# rules/TR_calling.smk

# tandem repeats calling using LongTR
rule TR_calling:
    input:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam.bai",
        TR_catalog=ancient(TR_CATALOG),
        sex = f"{OUTPUT_DIR}/qc/phased_bam/{{sample}}/{REFERENCE_NAME}/{{sample}}_sex_inference.csv",
        reference=ancient(REFERENCE)
    output:
        TR_vcf=f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}.vcf.gz"
    params:
        haploid_chrs= lambda wildcards, input: get_haploid_chromosomes(input.sex),
        type_of_tr=TYPE_OF_TR,
        sample_name=f"{{sample}}_{REFERENCE_NAME}"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/TR_calling/{{sample}}.log"
    singularity:
        "docker://quay.io/biocontainers/longtr:1.2--h077b44d_1"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/TR_calling/{{sample}}_TR_calling.tsv"
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


rule TR_calling_tabix:
    input:
        vcf=f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}.vcf.gz"
    output:
        tbi=f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}.vcf.gz.tbi"
    threads: 1
    singularity:
        "docker://leenaputzeys/tr_methylation:v1.0"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/TR_calling_tabix/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        tabix -f {input.vcf} > {log} 2>&1
        """
