# rules/variant_calling_snps_indels.smk

# variant calling: SNVs and Indels
rule variant_calling_snps_indels:
    input:
        aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam",
        bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam.bai",
        reference=ancient(REFERENCE)
    params:
        out_dir=f"{OUTPUT_DIR}/02_variant_calling/SNVs_Indels/{{sample}}/{REFERENCE_NAME}",
        model_path="/usr/local/bin/models/r1041_e82_400bps_sup_v500"
    output:
        snp_vcf_gz=f"{OUTPUT_DIR}/02_variant_calling/SNVs_Indels/{{sample}}/{REFERENCE_NAME}/merge_output.vcf.gz",
        snp_vcf_index=f"{OUTPUT_DIR}/02_variant_calling/SNVs_Indels/{{sample}}/{REFERENCE_NAME}/merge_output.vcf.gz.tbi"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/variant_calling_snps_indels/{{sample}}.log"
    threads: 48
    singularity:
        "docker://quay.io/biocontainers/clair3:1.2.0--py310h779eee5_0"
    shell:
        """
        export OMP_NUM_THREADS={threads}
        run_clair3.sh \
            --bam_fn {input.aligned_bam} \
            --ref_fn {input.reference} \
            --output {params.out_dir} \
            --threads {threads} \
            --platform="ont" \
            --model_path {params.model_path} \
            > {log} 2>&1
        """
