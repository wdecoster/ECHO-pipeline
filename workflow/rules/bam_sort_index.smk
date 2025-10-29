# rules/bam_sort_index.smk

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
    shell:
        """
        (
            samtools sort -@ {threads} {input.sam} > {output.bam}
            samtools index {output.bam}
        ) > {log} 2>&1
        """

