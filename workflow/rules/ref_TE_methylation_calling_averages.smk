# rules/ref_TE_methylation_calling_averages.smk

# reference TE average methylation calculations
rule ref_TE_methylation_calling_averages:
    input:
        phased_meth_bed_gz=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_1.bed.gz",
        SNV_vcf_gz=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased.vcf.gz"
    output:
        output_file=f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/{{sample}}_{REFERENCE_NAME}_methylation_summary.bed"
    params:
        phased_dir=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased",
        unphased_dir=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/unphased",
        out_dir=f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}",
        phased_variation_dir=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}",
        TE_catalog=TE_CATALOG,
        type_of_TE=TYPE_OF_TE,
        flanking_length_bp=FLANKING_LENGTH_BP,
        sample_name=f"{{sample}}_{REFERENCE_NAME}"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/ref_TE_methylation_calling_averages/{{sample}}.log"
    threads: 64
    singularity:
        "docker://leenaputzeys/te_methylation:v1.0"
    shell:
        """
        bash {REF_TE_METH_AVERAGES_SCRIPT_PATH} \
            -p {params.phased_dir} \
            -u {params.unphased_dir} \
            -v {params.phased_variation_dir} \
            -c {params.TE_catalog} \
            -t {params.type_of_TE} \
            -s {params.sample_name} \
            -x {threads} \
            -o {params.out_dir} \
            -f {params.flanking_length_bp} \
            > {log} 2>&1
        """

