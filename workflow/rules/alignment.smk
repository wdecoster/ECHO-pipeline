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
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/alignment/{{sample}}_alignment.tsv"
    shell:
        """
         minimap2 -ax map-ont -y -2 -t {threads} {input.reference} {input.fastq} > {output.sam} 2> {log}
        """

# SAM convert, sort and index
rule bam_sort_index:
    input:
        sam=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}.sam"
    output:
        bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam",
        bai=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam.bai"
    threads: 24
    singularity:
        "docker://staphb/samtools:1.22"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/sort_index/{{sample}}.log"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/alignment/{{sample}}_bam_sort_index.tsv"
    shell:
        """
        (
            samtools sort -@ {threads} {input.sam} > {output.bam}
            samtools index {output.bam}
        ) > {log} 2>&1
        """


