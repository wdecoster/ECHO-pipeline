# rules/non_ref_TE_methylation_calling.smk

# Non reference TE methylation calling based on TLDR output
rule non_ref_TE_methylation_calling:
    input:
        TLDR_txt_output=f"{OUTPUT_DIR}/05_non_ref_TE_calling/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}.table.txt"
    output:
        output_file=f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/non_ref_TE/{{sample}}_{REFERENCE_NAME}.table.pass.summary.meth.phased.txt"
    params:
        out_dir=f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/non_ref_TE",
        flanking_length_bp=FLANKING_LENGTH_BP,
        sample_name=f"{{sample}}_{REFERENCE_NAME}"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/non_ref_TE_methylation_calling/{{sample}}.log"
    threads: 32
#    WARNING: PARALELISATION DOES NOT WORK WITH CONTAINER - FIX  
    singularity:
        "/mnt/ngs/projects/Nanopore_pipeline/users/lfw156/container_creation/non_ref_te_meth_2.0.sif"
#    conda:
#        "../envs/te_methylation.yaml"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/non_ref_TE_methylation_calling/{{sample}}_non_ref_TE_methylation_calling.tsv"
    shell:
        """
        bash {TLDR_METH_SCRIPT_PATH} \
            -i {input.TLDR_txt_output} \
            -o {params.out_dir} \
            -s {params.sample_name} \
            -f {params.flanking_length_bp} \
            -t {threads} \
            > {log} 2>&1
        """
