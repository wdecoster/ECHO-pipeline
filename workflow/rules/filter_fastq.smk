# rules/filter_fastq.smk

# filtering reads based on length and quality in fastq files
if DO_FILTER:
    rule filter_fastq:
        input:
            fastq=raw_fastq
        output:
            filtered_fastq=f"{OUTPUT_DIR}/00_raw_data/basecalled/fastq/{{sample}}/{{sample}}_filtered.fastq"
        log:
            f"{OUTPUT_DIR}/logs/snakemake_rules/filter_fastq/{{sample}}.log"
        threads: 24
        singularity:
            "docker://quay.io/biocontainers/chopper:0.11.0--hcdda2d0_0"
        shell:
            """
            chopper -q {MIN_READ_QUAL} -l {MIN_READ_LENGTH} \
                    -i {input.fastq} -t {threads} > {output.filtered_fastq} 2> {log}
            """


