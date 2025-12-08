# rules/ref_TE_methylation_calling_cpg_res.smk

# rules/refactor_ref_TE_methylation_calling_averages.smk

rule ref_TE_meth_call_cpg_res_phased:
    input:
        pileup=lambda w: (
            f"{OUTPUT_DIR}/04_methylation_calling/"
            f"{w.sample}/{REFERENCE_NAME}/phased/"
            f"{w.sample}_{REFERENCE_NAME}_haplotype_{w.pileup}.bed.gz"
        ),
        catalog=lambda w: (
            f"{OUTPUT_DIR}/08_TE_methylation_calling/"
            f"{w.sample}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/"
            f"{TYPE_OF_TE}_{w.catalog}.bed"
        )
    output:
        f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/"
        f"cpg_resolution/mod_phased/{{sample}}_{REFERENCE_NAME}_{TYPE_OF_TE}_{{catalog}}_pileup_{{pileup}}.bed"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/ref_TE_methylation_calling_cpg_res_phased/"
        f"{{sample}}_phased_haplo_{{pileup}}_{{catalog}}.log"
    threads: 32
    singularity:
        "docker://leenaputzeys/te_methylation:v1.0"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/ref_TE_meth_call_cpg_res_phased/"
        f"{{sample}}_ref_TE_meth_call_cpg_res_phased_{{catalog}}_{{pileup}}.tsv"
    shell:
        """
        zcat {input.pileup} | bedtools intersect -a - -b {input.catalog} -wa -wb > {output}
        """


rule ref_TE_meth_call_cpg_res_unphased:
    input:
        pileup=lambda w: (
            f"{OUTPUT_DIR}/04_methylation_calling/"
            f"{w.sample}/{REFERENCE_NAME}/unphased/"
            f"{w.sample}_{REFERENCE_NAME}_unphased.bed.gz"
        ),
        catalog=lambda w: (
            f"{OUTPUT_DIR}/08_TE_methylation_calling/"
            f"{w.sample}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/"
            f"{TYPE_OF_TE}_{w.catalog}.bed"
        )
    output:
        f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/"
        f"cpg_resolution/mod_unphased/{{sample}}_{REFERENCE_NAME}_{TYPE_OF_TE}_{{catalog}}_pileup_unphased.bed"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/ref_TE_methylation_calling_cpg_res_unphased/"
        f"{{sample}}_unphased_{{catalog}}.log"
    threads: 32
    singularity:
        "docker://leenaputzeys/te_methylation:v1.0"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/ref_TE_meth_call_cpg_res_unphased/"
        f"{{sample}}_ref_TE_meth_call_cpg_res_unphased_{{catalog}}.tsv"
    shell:
        """
        zcat {input.pileup} | bedtools intersect -a - -b {input.catalog} -wa -wb > {output}
        """


