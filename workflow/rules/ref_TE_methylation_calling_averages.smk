# rules/refactor_ref_TE_methylation_calling_averages.smk

# reference TE average methylation calculations
rule prep_ref_TE_methylation_calling_averages:
    input:
        phased_meth_bed_gz_1=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_1.bed",
        phased_meth_bed_gz_2=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_2.bed",
        unphased_meth_bed_gz=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/unphased/{{sample}}_{REFERENCE_NAME}_unphased.bed",
        snp_vcf_gz=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased.vcf.gz",
        snp_vcf_gz=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_SV.vcf.gz"
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

rule modkit_ref_TE_meth_call_phased:
    input:
        pileup = f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_{{pileup}}.bed.gz",
        catalog = f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/{TYPE_OF_TE}_{{catalog}}.bed"
    output:
        tsv = f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/"
        f"mod_phased/{{sample}}_{REFERENCE_NAME}_{{catalog}}_stats_{{pileup}}.tsv"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/ref_TE_methylation_calling_averages/modkit_{{sample}}_phased_haplo{{pileup}}_{{catalog}}.log"
    threads: 32
    singularity:
        "docker://quay.io/biocontainers/ont-modkit:0.5.0--hcdda2d0_2"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/ref_TE_methylation_calling/modkit_{{sample}}_phased_haplo{{pileup}}_{{catalog}}_ref_TE_methylation_calling_averages.tsv"
    shell:
        """
        modkit stats \
            --force \
            -t {threads} \
            --regions {input.catalog} \
            -c m \
            -o {output.tsv} \
            {input.pileup} \
            > {log} 2>&1
        """

rule modkit_ref_TE_meth_call_unphased:
    input:
        pileup = f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/unphased/{{sample}}_{REFERENCE_NAME}_unphased.bed.gz",
        catalog = f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/{TYPE_OF_TE}_{{catalog}}.bed"
    output:
        tsv = f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/"
        f"mod_unphased/{{sample}}_{REFERENCE_NAME}_{{catalog}}_stats_unphased.tsv"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/ref_TE_methylation_calling_averages/modkit_{{sample}}_unphased__{{catalog}}.log"
    threads: 32
    singularity:
        "docker://quay.io/biocontainers/ont-modkit:0.5.0--hcdda2d0_2"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/ref_TE_methylation_calling/modkit_{{sample}}_unphased_{{catalog}}_ref_TE_methylation_calling_averages.tsv"
    shell:
        """
        modkit stats \
            --force \
            -t {threads} \
            --regions {input.catalog} \
            -c m \
            -o {output.tsv} \
            {input.pileup} \
            > {log} 2>&1
        """

# all catalogues used
CATALOGS=["upstream", "downstream", "full"]

rule summarise_ref_TE_methylation_calling_averages:
    input:
        expand(f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/mod_phased/{{sample}}_{REFERENCE_NAME}_{{catalog}}_stats_{{pileup}}.tsv", 
               sample=SAMPLES, pileup=['1', '2'], catalog=CATALOGS),
        expand(f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/mod_unphased/{{sample}}_{REFERENCE_NAME}_{{catalog}}_stats_unphased.tsv", 
               sample=SAMPLES, catalog=CATALOGS)
    output:
        f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/{{sample}}_{REFERENCE_NAME}_methylation_summary.bed"
    params:
        out_dir=f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}",
        phased_variation_dir=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}",
        TE_catalog=TE_CATALOG,
        sample_name=f"{{sample}}_{REFERENCE_NAME}"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/ref_TE_methylation_calling_averages/summarise_{{sample}}.log"
    threads: 1
    singularity:
        "docker://leenaputzeys/te_methylation:v1.0"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/ref_TE_methylation_calling/summarise_{{sample}}_ref_TE_methylation_calling_averages.tsv"
    shell:
        """
        bash {WORKFLOW_ROOT}/scripts/summarise_ref_TE_avg_meth.sh \
            -c {params.TE_catalog} \
            -s {params.sample_name} \
            -o {params.out_dir} \
            > {log} 2>&1
        """
