# rules/variant_calling_snps_indels.smk

# variant calling: SNVs and Indels
rule variant_calling_snps_indels:
    input:
        aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam",
        bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam.bai",
        reference=ancient(REFERENCE)
    params:
        out_dir=f"{OUTPUT_DIR}/02_variant_calling/SNVs_Indels/{{sample}}/{REFERENCE_NAME}",
        model_path="/usr/local/bin/models/r1041_e82_400bps_sup_v500",
        tmp_root = CLAIR3_TMP_ROOT
    output:
        snp_vcf_gz=f"{OUTPUT_DIR}/02_variant_calling/SNVs_Indels/{{sample}}/{REFERENCE_NAME}/phased_merge_output.vcf.gz",
        snp_vcf_index=f"{OUTPUT_DIR}/02_variant_calling/SNVs_Indels/{{sample}}/{REFERENCE_NAME}/phased_merge_output.vcf.gz.tbi"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/variant_calling_snps_indels/{{sample}}.log"
    threads: 16
    singularity:
        "docker://quay.io/biocontainers/clair3:1.2.0--py310h779eee5_0"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/variant_calling_snps_indels/{{sample}}_variant_calling_snps_indels.tsv"
    shell:
        """
        set -euo pipefail
        # Job-local Clair3 temporary directory
        TMPDIR="{params.tmp_root}/{wildcards.sample}"
        export TMPDIR
        export TMP="$TMPDIR"
        export TEMP="$TMPDIR"

        mkdir -p "$TMPDIR"
        cleanup() {{
             rm -rf "$TMPDIR"
         }}
        trap cleanup EXIT 
        
        echo "[INFO] Clair3 TMPDIR=$TMPDIR" >> {log}

        export OMP_NUM_THREADS={threads}
        run_clair3.sh \
            --bam_fn {input.aligned_bam} \
            --ref_fn {input.reference} \
            --output {params.out_dir} \
            --threads {threads} \
            --platform="ont" \
            --enable_phasing \
            --longphase_for_phasing \
            --use_longphase_for_final_output_phasing \
            --model_path {params.model_path} \
            --remove_intermediate_dir \
            > {log} 2>&1
        """

