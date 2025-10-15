# rules/longphase_modcall.smk

# phasing: longphase modcall
rule longphase_modcall:
    input:
        aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam",
        reference=ancient(REFERENCE)
    output:
        modcall_vcf=temp(f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_modcall.vcf")
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/phasing/{{sample}}_modcall.log"
    params:
        out_prefix=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}"
    threads: 16
    singularity:
        "docker://quay.io/biocontainers/longphase:2.0--h13024bc_0"
    shell:
        """
        (
            longphase modcall \
                -b {input.aligned_bam} \
                -r {input.reference} \
                -o {params.out_prefix}_modcall \
                -t {threads}
        ) > {log} 2>&1
        """
