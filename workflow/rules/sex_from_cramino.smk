# rules/sex_from_cramino.smk

# infer sex from cramino karyotyping 
rule sex_from_cramino:
    input:
        cramino_out=f"{OUTPUT_DIR}/qc/phased_bam/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_cramino_output.txt"
    output:
        out_file=f"{OUTPUT_DIR}/qc/phased_bam/{{sample}}/{REFERENCE_NAME}/{{sample}}_sex_inference.csv"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/sex_from_cramino/{{sample}}.log"
    singularity:
         "docker://biocontainers/pandas:1.5.1_cv1"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/sex_from_cramino/{{sample}}_sex.tsv"
    shell:
        """
        python3 {SEX_SCRIPT_PATH} {input.cramino_out} -o {output.out_file} > {log} 2>&1
        """

