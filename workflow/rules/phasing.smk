# rules/phasing.smk

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
        "docker://quay.io/biocontainers/longphase:1.7.3--hf5e1c6e_0"
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

# phasing: longphase phase
rule longphase_phase:
    input:
        aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam",
        bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam.bai",
        snp_vcf_gz=f"{OUTPUT_DIR}/02_variant_calling/SNVs_Indels/{{sample}}/{REFERENCE_NAME}/phased_merge_output.vcf.gz",
        snp_vcf_index=f"{OUTPUT_DIR}/02_variant_calling/SNVs_Indels/{{sample}}/{REFERENCE_NAME}/phased_merge_output.vcf.gz.tbi",
        sv_vcf_gz=f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_SV_unphased.vcf.gz",
        sv_vcf_index=f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_SV_unphased.vcf.gz.tbi",
        modcall_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_modcall.vcf",
        reference=ancient(REFERENCE)
    params:
        out_prefix=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}"
    output:
        SNV_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased.vcf",
        SV_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_SV.vcf",
        modcall_vcf_phased=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_mod.vcf"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/phasing/{{sample}}_phase.log"
    threads: 32
    singularity:
        "docker://quay.io/biocontainers/longphase:1.7.3--hf5e1c6e_0"
    shell:
        """
        (
            longphase phase \
                --snp-file {input.snp_vcf_gz} \
                --mod-file {params.out_prefix}_modcall.vcf \
                --sv-file {input.sv_vcf_gz} \
                -b {input.aligned_bam} \
                -r {input.reference} \
                -o {params.out_prefix}_phased \
                -t {threads} \
                --ont

        ) > {log} 2>&1
        """

# phasing: haplotagging
rule longphase_haplotag:
    input:
        aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam",
        bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam.bai",
        SNV_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased.vcf.gz",
        SV_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_SV.vcf.gz",
        modcall_vcf_phased=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_mod.vcf",
        reference=ancient(REFERENCE)
    params:
        out_prefix=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}"
    output:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/phasing/{{sample}}_haplotagging.log"
    threads: 48
    singularity:
        "docker://quay.io/biocontainers/longphase:1.7.3--hf5e1c6e_0"
    shell:
        """
        (
            longphase haplotag \
                -b {input.aligned_bam} \
                -r {input.reference} \
                -s {input.SNV_vcf} \
                --sv-file {input.SV_vcf} \
                --mod-file {input.modcall_vcf_phased} \
                -o {params.out_prefix}_phased_alignment
        ) > {log} 2>&1
        """

# phasing: bgzip and index files
rule phased_bam_index:
    input:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam"
    output:
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam.bai"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/phasing/{{sample}}_index_bam.log"
    threads: 16
    singularity:
        "docker://staphb/samtools:1.22"
    shell:
        """
        (
            samtools index {input.phased_bam}
        ) > {log} 2>&1
        """

# gzip and index phased VCF files
rule phased_bgzip:
    input:
        SNV_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased.vcf",
        SV_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_SV.vcf"
    params:
        out_prefix=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}"
    output:
        SNV_vcf_gz=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased.vcf.gz",
        SNV_vcf_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased.vcf.gz.tbi",
        SV_vcf_gz=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_SV.vcf.gz",
        SV_vcf_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_SV.vcf.gz.tbi",
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/phasing/{{sample}}_bgzip_phased_vcf.log"
    threads: 16
    singularity:
        "docker://quay.io/biocontainers/htslib:1.22.1--h566b1c6_0"
    shell:
        """
        (
            bgzip -@ 8 {params.out_prefix}_phased.vcf

            bgzip -@ 8 {params.out_prefix}_phased_SV.vcf
            tabix {params.out_prefix}_phased.vcf.gz
            tabix {params.out_prefix}_phased_SV.vcf.gz
        ) > {log} 2>&1
        """




