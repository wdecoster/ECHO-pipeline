SAMPLES = ["HG001_subset"]
OUTPUT_DIR = "/ifs/data/research/unique/projects/snakemake_test"
REFERENCE = "/ifs/data/research/unique/leena/references/T2T-CHM13/T2T-CHM13v2-renamed.fasta"



rule all:
    input:
        expand(f"{OUTPUT_DIR}/00_raw_data/basecalled/bam/{{sample}}_unaligned.bam", sample=SAMPLES),
	expand(f"{OUTPUT_DIR}/qc/qc_basecalling/{{sample}}_summary.tsv", sample=SAMPLES),
	expand(f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam", sample=SAMPLES), 
	expand(f"{OUTPUT_DIR}/02_variant_calling/SNPs_Indels/{{sample}}/phased_merge_output.vcf.gz", sample=SAMPLES), 
	expand(f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}_SV_unphased.vcf.gz", sample=SAMPLES), 
	expand(f"{OUTPUT_DIR}/03_phasing/{{sample}}_phased_alignment.bam", sample=SAMPLES), 
	expand(f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}_haplotype_1.bed", sample=SAMPLES), 
	expand(f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}_haplotype_2.bed", sample=SAMPLES), 
	expand(f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}_haplotype_ungrouped.bed", sample=SAMPLES), 
	expand(f"{OUTPUT_DIR}/05_TE_calling/{{sample}}", sample=SAMPLES), 
	expand(f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{{sample}}_TRs.vcf", sample=SAMPLES), 
	REFERENCE


# Basecalling runs on GPU
rule basecalling:
    input:
        pod5_dir=f"{OUTPUT_DIR}/00_raw_data/pod5/",
	model_dir="/ifs/software/research/unique/brando/dorado_models/"
    output:
        unaligned_bam=f"{OUTPUT_DIR}/00_raw_data/basecalled/bam/{{sample}}_unaligned.bam",
	summary_file=f"{OUTPUT_DIR}/qc/qc_basecalling/{{sample}}_summary.tsv"
    shell:
        """
	module load bioinf/dorado
	dorado basecaller sup,5mCG_5hmCG {input.pod5_dir} --trim --models-directory {input.model_dir} > {output.unaligned_bam}
	dorado summary {output.unaligned_bam} > {output.summary_file}
	"""


# Alignment runs on CPU
rule alignment:
    input:
        unaligned_bam=f"{OUTPUT_DIR}/00_raw_data/basecalled/bam/{{sample}}_unaligned.bam",
	reference=REFERENCE
    output:
        fastq=f"{OUTPUT_DIR}/00_raw_data/fastq/{{sample}}.fastq",
	aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam",
	bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam.bai"
    conda:
        "/ifs/software/research/unique/brando/pipeline/conda_env_yaml/minimap2.yaml"
    shell:
        """
	samtools fastq -T 'MM,ML' {input.unaligned_bam} > {output.fastq}
	minimap2 -y -t 64 -ax map-ont {input.reference} {output.fastq} | samtools sort -o {output.aligned_bam} -
	samtools index -o {output.bam_index} {output.aligned_bam}
	"""


# Variant calling: SNPs and Indels
rule variant_calling_snps_indels:
    input:
        aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam",
        bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam.bai",
	reference=REFERENCE,
	mode_path="/ifs/software/research/unique/leena/conda-envs/clair3/bin/models/r1041_e82_400bps_sup_v500"
    params:
        out_dir=f"{OUTPUT_DIR}/02_variant_calling/SNPs_Indels/{{sample}}"    
    output:
        snp_vcf_gz=f"{OUTPUT_DIR}/02_variant_calling/SNPs_Indels/{{sample}}/phased_merge_output.vcf.gz",
        snp_vcf_index=f"{OUTPUT_DIR}/02_variant_calling/SNPs_Indels/{{sample}}/phased_merge_output.vcf.gz.tbi"
    resources:
        cpus=32
    conda:
        "/ifs/software/research/unique/brando/pipeline/conda_env_yaml/clair3.yaml"
    shell:
        """
	export OMP_NUM_THREADS={resources.cpus}
	run_clair3.sh \
        --bam_fn {input.aligned_bam} \
        --ref_fn {input.reference} \
        --output {params.out_dir} \
        --threads {resources.cpus} \
        --platform="ont" \
        --model_path {input.mode_path} \
        --enable_phasing \
        --longphase_for_phasing
	"""


# Variant calling: STRUCTURAL VARIANTS
rule variant_calling_sv:
    input:
        aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam",
        bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam.bai",
        reference=REFERENCE
    output:
        sv_vcf_gz=f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}_SV_unphased.vcf.gz",
        sv_vcf_index=f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}_SV_unphased.vcf.gz.tbi"
    resources:
        cpus=8
    conda:
        "/ifs/software/research/unique/brando/pipeline/conda_env_yaml/sniffles2.yaml"
    shell:
        """
        sniffles \
        --input {input.aligned_bam} \
        --reference {input.reference} \
        --vcf {output.sv_vcf_gz}.tmp \
        --threads {resources.cpus} \
        --output-rnames

        bgzip -c {output.sv_vcf_gz}.tmp > {output.sv_vcf_gz} 
        tabix -f -p vcf {output.sv_vcf_gz} 
        rm -f {output.sv_vcf_gz}.tmp
        """


# Phasing
rule phasing:
    input:
        aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam",
        bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam.bai",
	snp_vcf_gz=f"{OUTPUT_DIR}/02_variant_calling/SNPs_Indels/{{sample}}/phased_merge_output.vcf.gz",
	snp_vcf_index=f"{OUTPUT_DIR}/02_variant_calling/SNPs_Indels/{{sample}}/phased_merge_output.vcf.gz.tbi",
	sv_vcf_gz=f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}_SV_unphased.vcf.gz",
        sv_vcf_index=f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}_SV_unphased.vcf.gz.tbi",
	reference=REFERENCE
    params:
        out_prefix=f"{OUTPUT_DIR}/03_phasing/{{sample}}"    
    output:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}_phased_alignment.bam.bai",
    resources:
        cpus=16
    conda:
        "/ifs/software/research/unique/brando/pipeline/conda_env_yaml/longphase.yaml"
    shell:
        """
	/ifs/software/research/unique/brando/pipeline/tools/longphase/longphase modcall \
	-b {input.aligned_bam} \
	-r {input.reference} \
	-o {params.out_prefix}_modcall \
	-t {resources.cpus}

	/ifs/software/research/unique/brando/pipeline/tools/longphase/longphase phase \
	--snp-file {input.snp_vcf_gz} \
	--mod-file {params.out_prefix}_modcall.vcf \
	--sv-file {input.sv_vcf_gz} \
	-b {input.aligned_bam} \
	-r {input.reference} \
	-o {params.out_prefix}_phased \
	-t {resources.cpus} \
	--ont

	/ifs/software/research/unique/brando/pipeline/tools/longphase/longphase haplotag \
	-b {input.aligned_bam} \
	-r {input.reference} \
	-s {params.out_prefix}_phased.vcf \
	--sv-file {params.out_prefix}_phased_SV.vcf \
	--mod-file {params.out_prefix}_phased_mod.vcf \
	-o {params.out_prefix}_phased_alignment 
        
	samtools index -o {output.phased_bam_index} {output.phased_bam} 
	"""


# Methylation calling
rule methylation_calling:
    input:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam.bai",
        reference=REFERENCE
    params:
        out_dir=f"{OUTPUT_DIR}/04_methylation_calling/"
    output:
        hap1_bed=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}_haplotype_1.bed",
        hap2_bed=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}_haplotype_2.bed",
        ungrouped_bed=f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}_haplotype_ungrouped.bed"
    conda:
        "/ifs/software/research/unique/brando/pipeline/conda_env_yaml/modkit.yaml"
    shell:
        """
        modkit pileup {input.phased_bam} {params.out_dir}/ \
        --ref {input.reference} \
        --combine-strands \
        --cpg \
        --partition-tag HP \
        --prefix {wildcards.sample}_haplotype
	"""


# TE calling
rule TE_calling:
    input:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam.bai",
        te_library="/ifs/data/research/unique/repeat-catalogs/TEs/teref.ont.human.fa",
	reference=REFERENCE
    output:
        out_dir=directory(f"{OUTPUT_DIR}/05_TE_calling/{{sample}}")
    conda:
        "/ifs/software/research/unique/brando/pipeline/conda_env_yaml/tldr.yaml"
    shell:
        """
        /ifs/software/research/unique/brando/pipeline/tools/tldr/tldr \
	-b {input.phased_bam} \
        -e {input.te_library} \
        -r {input.reference} \
	-p 32 \
        --outbase {output.out_dir} \
        --detail_output \
        --methylartist \
        --max_cluster_size 500 \
        --extend_consensus 200
	"""



# Tandem Repeats calling
rule TR_calling:
    input:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam.bai",
        STR_catalog="/ifs/data/research/unique/repeat-catalogs/TRs/T2T-CHM13/genome-wide/repeat_catalog_v1.hg38_liftOver_T2T.adjusted.1_to_1000bp_motifs.longTR.bed",
	reference=REFERENCE
    output:
        TR_vcf=f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{{sample}}_TRs.vcf"
    conda:
        "/ifs/software/research/unique/brando/pipeline/conda_env_yaml/longtr.yaml"
    shell:
        """
        /ifs/software/research/unique/leena/longTR/LongTR \
        --bams {input.phased_bam} \
        --regions {input.STR_catalog} \
        --fasta {input.reference} \
        --tr-vcf  {output.TR_vcf}\
        --bam-samps {wildcards.sample} \
        --bam-libs {wildcards.sample} \
        --haploid chrX,chrY \
        --phased-bam
        """



