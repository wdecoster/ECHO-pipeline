# rules/TE_calling.smk

# TE calling
rule TE_calling:
    input:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam.bai",
        reference_TE=ancient(REFERENCE_TE),
        reference=ancient(REFERENCE)
    output:
        out_file=f"{OUTPUT_DIR}/05_non_ref_TE_calling/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}.table.txt"
    params:
        chr_file=f"{OUTPUT_DIR}/05_non_ref_TE_calling/{{sample}}/{REFERENCE_NAME}/{{sample}}_chr.txt",
        out_dir=directory(f"{OUTPUT_DIR}/05_non_ref_TE_calling/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}"),
        consensus_extension=CONSENSUS_EXTENSION
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/TE_calling/{{sample}}.log"
    threads: 48
    singularity:
        "docker://leenaputzeys/tldr:v1.2.2-7c6dfda"
    shell:
        """
        (
        for chr in {{1..22}} M X Y; do
            echo "chr${{chr}}" >> {params.chr_file}
        done

        tldr \
            -b {input.phased_bam} \
            -e {input.reference_TE} \
            -r {input.reference} \
            -p {threads} \
            -c {params.chr_file} \
            --outbase {params.out_dir} \
            --detail_output \
            --embed_minreads 3 \
            --methylartist \
            --max_cluster_size 500 \
            --extend_consensus {params.consensus_extension}
        ) > {log} 2>&1
        rm {params.chr_file}
        """

