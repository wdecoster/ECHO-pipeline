# rules/TR_methylation_calling.smk

# TR meth pep - fixes header
rule TR_methyl_call_vcf_prep:
    input:
        tr_vcf=f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}.vcf.gz"
    output:
        out_vcf=f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}_input_sorted.vcf.gz"
    log: 
        f"{OUTPUT_DIR}/logs/snakemake_rules/TR_methyl_call_vcf_prep/{{sample}}.log"
    threads: 1
    singularity:
        "docker://staphb/bcftools:1.22"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/TR_methyl_call_vcf_prep/{{sample}}_TR_methyl_call_vcf_prep.tsv"
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1

        BODY_N=$(bcftools view -H {input.tr_vcf} | wc -l)
        echo "VCF body rows: $BODY_N"
        bcftools annotate \
            --header-lines <(echo '##FORMAT=<ID=DFLANKINDEL,Number=1,Type=Integer,Description="Total number of reads with an indel in the regions flanking the STR">') \
            -Ov {input.tr_vcf} | \
            bcftools sort - -Oz -o {output.out_vcf}
        """


# splits vcf into N vcfs for parallel rules - chunk size set by var in setup
checkpoint TR_methyl_call_vcf_chunking:
    input:
        fixed_vcf=f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}_input_sorted.vcf.gz"
    params:
        out_pre=f"{{sample}}_chunk"
    output:
        chunk_dir=directory(f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}/chunked_vcfs"),
        done=f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}/chunked_vcfs/.done"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/TR_methyl_call_vcf_chunking/{{sample}}.log"
    threads: 1
    singularity:
        "docker://staphb/bcftools:1.22"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/TR_methyl_call_vcf_chunking/{{sample}}_TR_methyl_call_vcf_chunking.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > {log} 2>&1

        echo "[INFO] Starting TR_methyl_call_vcf_chunking for sample {wildcards.sample}"
        echo "[INFO] Input VCF: {input.fixed_vcf}"

        mkdir -p {output.chunk_dir}

        N_VARS=$(bcftools view -H {input.fixed_vcf} | wc -l)
        echo "$N_VARS variants in input vcf"

        CHUNK_VARS=$(( ( N_VARS / {N_VCF_CHUNKS} ) + ( N_VARS % {N_VCF_CHUNKS} > 0 ) ))

        if (( N_VARS == 0 )); then
            echo "No variants were found in the TR vcf file"
            for chunk_n in $(seq 0 $(({N_VCF_CHUNKS} - 1))); do
                cp {input.fixed_vcf} {output.chunk_dir}/{params.out_pre}${{chunk_n}}.vcf.gz
            done
            touch {output.done}
            exit 0
        fi

        echo "putting $CHUNK_VARS variants per chunked vcf"
        bcftools +scatter {input.fixed_vcf} -Oz -n $CHUNK_VARS -o {output.chunk_dir} -p {params.out_pre}
        touch {output.done}
        echo "Rule completed"
        """


def get_chunks(wildcards):
    ckpt = checkpoints.TR_methyl_call_vcf_chunking.get(sample=wildcards.sample)
    chunk_dir = ckpt.output["chunk_dir"]

    pattern = os.path.join(chunk_dir, f"{wildcards.sample}_chunk{{chunk}}.vcf.gz")
    chunks = glob_wildcards(pattern).chunk
    return sorted(chunks, key=lambda x: int(x))


# tandem repeats methylation calling
rule TR_methylation_calling:
    input:
        chunk_vcf=lambda wc: os.path.join(
            checkpoints.TR_methyl_call_vcf_chunking.get(sample=wc.sample).output["chunk_dir"],
            f"{wc.sample}_chunk{wc.chunk}.vcf.gz"
        ),
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam.bai",
        sex=f"{OUTPUT_DIR}/qc/phased_bam/{{sample}}/{REFERENCE_NAME}/{{sample}}_sex_inference.csv",
        reference=ancient(REFERENCE)
    output:
        unsorted_vcf=f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}/chunked_vcfs/{{sample}}_chunk{{chunk}}_methylated.vcf"
    params:
        tmp_out_dir=f"{SYS_TMP}/{REFERENCE_NAME}_{{sample}}_{TYPE_OF_TR}_methyl_call_chunk{{chunk}}", # puts out dir on tmp
        flanking_length_bp=FLANKING_LENGTH_BP,
        haploid_chrs= lambda wildcards, input: get_haploid_chromosomes(input.sex),
        extension=CONSENSUS_EXTENSION,
        type_of_TR=TYPE_OF_TR,
        sample_name=f"{{sample}}_{REFERENCE_NAME}",
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/TR_methylation_calling/{{sample}}_chunk{{chunk}}.log"
    threads: 8
    singularity:
         "docker://leenaputzeys/tr_methylation:v1.0"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/TR_methylation_calling/{{sample}}_chunk{{chunk}}_TR_methylation_calling.tsv"
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1

        bash {WORKFLOW_ROOT}/scripts/TR-longTR-methylation.sh \
            -v {input.chunk_vcf} \
            -r {input.reference} \
            -i {input.phased_bam} \
            -o {params.tmp_out_dir} \
            -s {params.sample_name} \
            -e {params.extension} \
            -t {threads} \
            -f {params.flanking_length_bp} \
            -h {params.haploid_chrs} \
            > {log} 2>&1
        """

# sorting the output annotated vcf
rule TR_methyl_sort_vcf_chunk:
    input:
        fixed_vcf=f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}_input_sorted.vcf.gz",
        chunk_vcf=f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}/chunked_vcfs/{{sample}}_chunk{{chunk}}_methylated.vcf"
    output:
        out_vcf=f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}/chunked_vcfs/{{sample}}_chunk{{chunk}}_methylated.vcf.gz"
    log:
        f"{OUTPUT_DIR}/logs/snakemake_rules/TR_methyl_sort_vcf_chunk/{{sample}}_chunk{{chunk}}.log"
    threads: 1
    singularity:
         "docker://leenaputzeys/tr_methylation:v1.0"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/TR_methyl_sort_vcf_chunk/{{sample}}_chunk{{chunk}}_TR_methyl_sort_vcf_chunk.tsv"
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1

        # If TR vcf (output of LongTR) does not contain anz variant skip the rule and output empty files 
        N_VARS=$(bcftools view -H {input.chunk_vcf} | wc -l)
        if (( N_VARS == 0 )); then
            echo "No variants were found in the TR vcf file: skip bcftools sort"
            cp {input.chunk_vcf} {output.out_vcf}
            exit 0
        fi

        bcftools sort -Oz -o {output.out_vcf} {input.chunk_vcf}
        tabix {output.out_vcf}
        """

# merges the chunked vcfs and creates the summary tsv
rule TR_methyl_cat_vcfs:
    input: 
        fixed_vcf=f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}_input_sorted.vcf.gz",
        vcf_list=lambda wc: expand(
            f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}/chunked_vcfs/{{sample}}_chunk{{chunk}}_methylated.vcf.gz",
            sample=wc.sample,
            chunk=get_chunks(wc))        
    output:
        out_vcf=f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}_TR_methylation.vcf.gz",
        out_tsv=f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}_TR_methylation_summary.tsv"
    params:
        header_vcf=f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}_input_sorted.vcf.gz"
    log: 
        f"{OUTPUT_DIR}/logs/snakemake_rules/TR_methyl_cat_vcfs/{{sample}}.log"
    threads: 1
    singularity:
        "docker://staphb/bcftools:1.22"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/TR_methyl_cat_vcfs/{{sample}}_TR_methyl_cat_vcfs.tsv"
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1

        # init file with header
        echo -e "CHROM\tPOS\tID\tREF_ALLELE\tALT_ALLELES\tREF_MOTIF\tGT\tTR_LEN\tTR_N_CPG\tTR_PATTERN\tTR_AM\tTR_N_METH_VALID\tUPSTREAM_TR_AM\tUPSTREAM_TR_N_METH_VALID\tDOWNSTREAM_TR_AM\tDOWNSTREAM_TR_N_METH_VALID\tTR_CPG_METH_HP1\tTR_CPG_DEPTH_HP1\tTR_CPG_METH_HP2\tTR_CPG_DEPTH_HP2" > {output.out_tsv}
        #if TR vcf (output of LongTR) does not contain anz variant skip the rule and output empty files
        N_VARS=$(bcftools view -H {input.fixed_vcf} | wc -l)
        if (( N_VARS == 0 )); then
            echo "No variants were found in the TR vcf file: skip generation summary file"
            bcftools view -h {params.header_vcf} > {output.out_vcf}
            exit 0
        fi
        
        bcftools concat -a {input.vcf_list} | \
            bcftools annotate --remove 'INFO/NSKIP,INFO/NFILT,INFO/INEXACT_ALLELE,INFO/BPDIFFS,INFO/DP,INFO/DSNP,INFO/DFLANKINDEL,INFO/REFAC,INFO/AC' -Ov | \
            bcftools sort - -Oz -o {output.out_vcf}
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%MOTIF\t[%GT\t%TR_LEN\t%TR_N_CPG\t%TR_PATTERN\t%TR_AM\t%TR_N_METH_VALID\t%UPSTREAM_TR_AM\t%UPSTREAM_TR_N_METH_VALID\t%DOWNSTREAM_TR_AM\t%DOWNSTREAM_TR_N_METH_VALID\t%TR_CPG_METH_HP1\t%TR_CPG_DEPTH_HP1\t%TR_CPG_METH_HP2\t%TR_CPG_DEPTH_HP2\n]' {output.out_vcf} >> {output.out_tsv}
        """
