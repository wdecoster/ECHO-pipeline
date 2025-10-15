# rules/alignment.smk

# alignment
rule alignment:
    input:
        fastq=fastq_for_pipeline,
        reference=ancient(REFERENCE)
    output:
        sam=temp(f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}.sam")
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/alignment/{{sample}}.log"
    threads: 80
    singularity:
        "docker://quay.io/biocontainers/minimap2:2.30--h577a1d6_0"
    shell:
        """
         minimap2 -ax map-ont -y -2 -t {threads} {input.reference} {input.fastq} > {output.sam} 2> {log}
        """


