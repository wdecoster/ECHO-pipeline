# rules/TR_methylation_calling.smk

# tandem repeats methylation calling
rule TR_methylation_calling:
    input:
        TR_vcf=f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}.vcf.gz",
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam.bai",
        reference=ancient(REFERENCE)
    output:
        out_file=f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}/{{sample}}_{REFERENCE_NAME}_TR_methylation_summary.tsv"
    params:
        out_dir=f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}",
        flanking_length_bp=FLANKING_LENGTH_BP,
        haploid_chrs=HAPLOID_CHRS,
        extension=CONSENSUS_EXTENSION,
        type_of_TR=TYPE_OF_TR,
        sample_name=f"{{sample}}_{REFERENCE_NAME}"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/TR_methylation_calling/{{sample}}.log"
    threads: 32
    singularity:
         "docker://leenaputzeys/tr_methylation:v1.0"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/TR_methylation_calling/{{sample}}_TR_methylation_calling.tsv"
    shell:
        """
        bash {TR_LONGTR_METH_SCRIPT_PATH} \
            -v {input.TR_vcf} \
            -r {input.reference} \
            -i {input.phased_bam} \
            -o {params.out_dir} \
            -s {params.sample_name} \
            -e {params.extension} \
            -t {threads} \
            -f {params.flanking_length_bp} \
            -h {params.haploid_chrs} \
            > {log} 2>&1
        """

