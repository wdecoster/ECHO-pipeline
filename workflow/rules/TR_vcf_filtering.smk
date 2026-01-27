# rules/TR_vcf_filtering.smk
# Produces a canonical VCF used by TR methylation:
# - If CpG filtering is OFF -> copy-through
# - If CpG filtering is ON  -> bcftools view -R with CpG STR BED

rule TR_vcf_for_methylation:
    input:
        tr_vcf=f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}.vcf.gz",
        tr_vcf_tbi=f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}.vcf.gz.tbi"
    params:
        cpg_bed=lambda wc: CPG_STR_BED  # may be None if filtering is off
    output:
        out_vcf=f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}_for_methylation.vcf.gz"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/TR_vcf_for_methylation/{{sample}}.log"
    threads: 1
    singularity:
        "docker://leenaputzeys/tr_methylation:v1.0"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/TR_vcf_for_methylation/{{sample}}_TR_vcf_for_methylation.tsv"
    run:
        if not FILTER_TO_CPG_STR:
            shell(r"""
                set -euo pipefail
                echo "[TR_vcf_for_methylation] CpG filter OFF -> copy-through" > {log}
                cp {input.tr_vcf} {output.out_vcf} >> {log} 2>&1
                tabix -f {output.out_vcf} >> {log} 2>&1 || true
            """)
        else:
            if not params.cpg_bed:
                raise ValueError("CpG STR filtering enabled but cpg_str_bed is missing.")
            shell(r"""
                set -euo pipefail
                echo "[TR_vcf_for_methylation] CpG filter ON -> bcftools view -R" > {log}

                N_VARS=$(bcftools view -H {input.tr_vcf} | wc -l) >> {log} 2>&1 || true
                echo "N_VARS=$N_VARS" >> {log}

                if (( N_VARS == 0 )); then
                    echo "No variants in input VCF; writing empty filtered VCF" >> {log}
                    bcftools view -h {input.tr_vcf} | bgzip -c > {output.out_vcf} 2>> {log}
                    tabix -f {output.out_vcf} >> {log} 2>&1
                    exit 0
                fi

                bcftools view -R {params.cpg_bed} -Oz -o {output.out_vcf} {input.tr_vcf} >> {log} 2>&1
                tabix -f {output.out_vcf} >> {log} 2>&1
            """)

