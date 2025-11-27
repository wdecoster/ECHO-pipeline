# rules/ref_TE_methylation_calling_cpg_res.smk

# rules/refactor_ref_TE_methylation_calling_averages.smk

# reference TE average methylation calculations
rule prep_ref_TE_methylation_calling_averages:
    input:
        phased_meth_bed_gz_1=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_1.bed",
        phased_meth_bed_gz_2=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_2.bed",
        unphased_meth_bed_gz=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/unphased/{{sample}}_{REFERENCE_NAME}_unphased.bed",
        snp_vcf_gz=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased.vcf.gz",
        sv_vcf_gz=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_SV.vcf.gz"
    output:
        upstream_te_catalog=f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/{TYPE_OF_TE}_upstream.bed",
        downstream_te_catalog=f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/{TYPE_OF_TE}_downstream.bed",
        snp_intersect=f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/variants/{{sample}}_{REFERENCE_NAME}_SNPs_intersect.bed",
        sv_intersect=f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/variants/{{sample}}_{REFERENCE_NAME}_SV_intersect.bed",
        snp_intersect_count=f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/variants/{{sample}}_{REFERENCE_NAME}_SNPs_intersect_count.bed",
        sv_intersect_count=f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/variants/{{sample}}_{REFERENCE_NAME}_SV_intersect_count.bed"
    params:
        phased_dir=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased",
        unphased_dir=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/unphased",
        out_dir=f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}",
        phased_variation_dir=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}",
        sample_name=f"{{sample}}_{REFERENCE_NAME}"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/ref_TE_methylation_calling_averages/prep_{{sample}}.log"
    threads: 4
    singularity:
        "docker://leenaputzeys/te_methylation:v1.0"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/ref_TE_methylation_calling/prep_{{sample}}_ref_TE_methylation_calling_averages.tsv"
    shell:
        """
        bash {WORKFLOW_ROOT}/scripts/prep_ref_TE_avg_meth.sh \
            -p {params.phased_dir} \
            -u {params.unphased_dir} \
            -v {params.phased_variation_dir} \
            -c {TE_CATALOG} \
            -t {TYPE_OF_TE} \
            -s {params.sample_name} \
            -x {threads} \
            -o {params.out_dir} \
            -f {FLANKING_LENGTH_BP} \
            > {log} 2>&1
        """

# a simple rule to make the full catalogue follow the naming converntions of the small catalogue
rule cp_catalogue:
    input:
        f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/{TYPE_OF_TE}_upstream.bed",  # this is not used but means it runs after prev
        full_cat=TE_CATALOG
    output:
        f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/{TYPE_OF_TE}_full.bed"
    threads: 
        1
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/ref_TE_methylation_calling/prep_{{sample}}_cp_catalogue.tsv"
    shell:
        """
        cp {input.full_cat} {output}
        """

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


