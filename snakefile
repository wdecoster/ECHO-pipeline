configfile: "config.yaml"


SAMPLES = config["samples"]
START_FROM = config.get("start_from", "pod5")  # default to "pod5"
OUTPUT_DIR = config["output_dir"]
INPUT_DIR = config["input_dir"]
REFERENCE = config["reference"]
REFRENCE_TE = config["reference_TE"]
TR_CATALOG = config["tr_catalog"]
TE_CATALOG = config["te_catalog"]
FLANKING_LENGTH_BP = config["flanking_length_bp"]
HAPLOID_CHRS = config["haploid_chrs"]
TYPE_OF_TE = config["type_of_te"]


all_inputs = []

if START_FROM == "pod5":
    all_inputs.extend([
        expand(f"{OUTPUT_DIR}/00_raw_data/basecalled/ubam/{{sample}}_unaligned.bam", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/qc/qc_basecalling/{{sample}}_summary.tsv", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam", sample=SAMPLES),
    ])
elif START_FROM == "ubam":
    all_inputs.extend([
        expand(f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam", sample=SAMPLES),
    ])
elif START_FROM == "bam":
    # no basecalling or alignment outputs expected
    pass
else:
    print(f"ERROR: Invalid 'start_from' value in config.yaml: {START_FROM}")
    exit(1)


all_inputs.extend([
    expand(f"{OUTPUT_DIR}/02_variant_calling/SNPs_Indels/{{sample}}/phased_merge_output.vcf.gz", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}_SV_unphased.vcf.gz", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/03_phasing/{{sample}}/{{sample}}_phased_alignment.bam", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/04_methylation_calling/phased/{{sample}}_haplotype_1.bed.gz", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/04_methylation_calling/phased/{{sample}}_haplotype_2.bed.gz", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/04_methylation_calling/unphased/{{sample}}_unphased.bed.gz", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/05_TE_calling/{{sample}}/non_ref_TE/{{sample}}.table.txt", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{{sample}}_TRs.vcf.gz", sample=SAMPLES),
    directory(expand(f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}", sample=SAMPLES)),
    expand(f"{OUTPUT_DIR}/05_TE_calling/{{sample}}/non_ref_TE/{{sample}}/{{sample}}.table.pass.summary.meth.phased.txt", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/05_TE_calling/{{sample}}/ref_TE/{{sample}}_summary_per_ref_{TYPE_OF_TE}.txt", sample=SAMPLES),
    REFERENCE
])

rule all:
    input: all_inputs


#Basecalling
if START_FROM == "pod5":
    rule basecalling:
        input:
            pod5_dir=f"{INPUT_DIR}/00_raw_data/pod5/{{sample}}/",
	    model_dir="/ifs/software/research/unique/brando/dorado_models/"
        output:
            unaligned_bam=f"{OUTPUT_DIR}/00_raw_data/basecalled/ubam/{{sample}}_unaligned.bam",
	    summary_file=f"{OUTPUT_DIR}/qc/qc_basecalling/{{sample}}_summary.tsv" 
        shell:
            """
	    module load bioinf/dorado
	    dorado basecaller sup,5mCG_5hmCG {input.pod5_dir} --trim --models-directory {input.model_dir} > {output.unaligned_bam}
	    dorado summary {output.unaligned_bam} > {output.summary_file}
	    """


# Alignment
if START_FROM in ["pod5", "ubam"]:
    rule alignment:
        input:
            unaligned_bam=f"{OUTPUT_DIR}/00_raw_data/basecalled/ubam/{{sample}}_unaligned.bam",
    	    reference=REFERENCE
        output:
            fastq=f"{OUTPUT_DIR}/00_raw_data/basecalled/fastq/{{sample}}.fastq",
    	    aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam",
    	    bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam.bai"
        conda:
            "conda_env_yaml/minimap2.yaml"
        shell:
            """
    	    samtools fastq -T 'MM,ML' {input.unaligned_bam} > {output.fastq}
    	    minimap2 -y -t 64 -ax map-ont {input.reference} {output.fastq} | samtools sort -o {output.aligned_bam} -
    	    samtools index -o {output.bam_index} {output.aligned_bam}
    	    """


# Variant calling: SNPs and Indels
if START_FROM in ["pod5", "ubam", "bam"]:
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
            "conda_env_yaml/clair3.yaml"
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
if START_FROM in ["pod5", "ubam", "bam"]:
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
            "conda_env_yaml/sniffles2.yaml"
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
        out_prefix=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{{sample}}"    
    output:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{{sample}}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{{sample}}_phased_alignment.bam.bai",
	snp_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{{sample}}_phased.vcf.gz"
    resources:
        cpus=16
    conda:
        "conda_env_yaml/longphase.yaml"
    shell:
        """
	/ifs/software/research/unique/pipeline_tools/longphase/longphase modcall \
	-b {input.aligned_bam} \
	-r {input.reference} \
	-o {params.out_prefix}_modcall \
	-t {resources.cpus}

	/ifs/software/research/unique/pipeline_tools/longphase/longphase phase \
	--snp-file {input.snp_vcf_gz} \
	--mod-file {params.out_prefix}_modcall.vcf \
	--sv-file {input.sv_vcf_gz} \
	-b {input.aligned_bam} \
	-r {input.reference} \
	-o {params.out_prefix}_phased \
	-t {resources.cpus} \
	--ont

	/ifs/software/research/unique/pipeline_tools/longphase/longphase haplotag \
	-b {input.aligned_bam} \
	-r {input.reference} \
	-s {params.out_prefix}_phased.vcf \
	--sv-file {params.out_prefix}_phased_SV.vcf \
	--mod-file {params.out_prefix}_phased_mod.vcf \
	-o {params.out_prefix}_phased_alignment 
        
	samtools index -o {output.phased_bam_index} {output.phased_bam} 
	bgzip {params.out_prefix}_phased.vcf
	bgzip {params.out_prefix}_phased_SV.vcf 
	tabix {params.out_prefix}_phased.vcf.gz
	tabix {params.out_prefix}_phased_SV.vcf.gz
	"""



# Methylation calling
rule methylation_calling:
    input:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{{sample}}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{{sample}}_phased_alignment.bam.bai",
        reference=REFERENCE
    params:
        out_dir=f"{OUTPUT_DIR}/04_methylation_calling/",
        hap1_bed=f"{OUTPUT_DIR}/04_methylation_calling/phased/{{sample}}_haplotype_1.bed",
        hap2_bed=f"{OUTPUT_DIR}/04_methylation_calling/phased/{{sample}}_haplotype_2.bed",
        ungrouped_bed=f"{OUTPUT_DIR}/04_methylation_calling/phased/{{sample}}_haplotype_ungrouped.bed",
        unphased_bed=f"{OUTPUT_DIR}/04_methylation_calling/unphased/{{sample}}_unphased.bed"
    output:
        hap1_bed_gz=f"{OUTPUT_DIR}/04_methylation_calling/phased/{{sample}}_haplotype_1.bed.gz",
        hap2_bed_gz=f"{OUTPUT_DIR}/04_methylation_calling/phased/{{sample}}_haplotype_2.bed.gz",
        unphased_bed_gz=f"{OUTPUT_DIR}/04_methylation_calling/unphased/{{sample}}_unphased.bed.gz"
    conda:
        "conda_env_yaml/modkit.yaml"
    shell:
        """
        modkit pileup {input.phased_bam} {params.out_dir}/phased/ \
        --ref {input.reference} \
        --combine-strands \
        --cpg \
        --partition-tag HP \
        --prefix {wildcards.sample}_haplotype
	
	modkit pileup {input.phased_bam} {params.unphased_bed} \
        --ref {input.reference} \
        --combine-strands \
        --cpg 
        
        echo start	
	bgzip {params.hap1_bed}
	bgzip {params.hap2_bed}
	bgzip {params.ungrouped_bed}
	bgzip {params.unphased_bed}
	echo finished
	tabix {params.hap1_bed}.gz
        tabix {params.hap2_bed}.gz
       	tabix {params.ungrouped_bed}.gz
       	tabix {params.unphased_bed}.gz
        """


# TE calling
rule TE_calling:
    input:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{{sample}}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{{sample}}_phased_alignment.bam.bai",
        reference_TE=REFRENCE_TE,
	reference=REFERENCE
    output:
        chr_file=f"{OUTPUT_DIR}/05_TE_calling/{{sample}}/non_ref_TE/{{sample}}/chr.txt",
	out_file=f"{OUTPUT_DIR}/05_TE_calling/{{sample}}/non_ref_TE/{{sample}}.table.txt"
    params:
        out_dir=directory(f"{OUTPUT_DIR}/05_TE_calling/{{sample}}/non_ref_TE/{{sample}}"),
        flanking_length_bp=FLANKING_LENGTH_BP
    conda:
        "conda_env_yaml/tldr.yaml"
    shell:
        """
        # Generate chromosome list file
        for chr in {{1..22}} M X Y; do
            echo "chr${{chr}}" >> {output.chr_file}
        done

        /ifs/software/research/unique/pipeline_tools/tldr/tldr \
	-b {input.phased_bam} \
        -e {input.reference_TE} \
        -r {input.reference} \
	-p 32 \
        -c {output.chr_file} \
        --outbase {params.out_dir} \
        --detail_output \
        --methylartist \
        --max_cluster_size 500 \
        --extend_consensus {params.flanking_length_bp}
	
	"""


# Tandem Repeats calling
rule TR_calling:
    input:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{{sample}}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{{sample}}_phased_alignment.bam.bai",
        TR_catalog=TR_CATALOG,
	reference=REFERENCE
    output:
        TR_vcf=f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{{sample}}_TRs.vcf.gz"
    params:
        haploid_chrs=HAPLOID_CHRS
    container:
        "/ifs/software/research/unique/containers/longtr_2025_11_03.sif"
    shell:
        """
        /bin/LongTR \
        --bams {input.phased_bam} \
        --regions {input.TR_catalog} \
        --fasta {input.reference} \
        --tr-vcf  {output.TR_vcf}\
        --bam-samps {wildcards.sample} \
        --bam-libs {wildcards.sample} \
        --haploid {params.haploid_chrs} \
        --phased-bam
        """


# Tandem Repeats methylation calling
rule TR_methylation_calling:
    input:
        TR_vcf=f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{{sample}}_TRs.vcf.gz",
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{{sample}}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{{sample}}_phased_alignment.bam.bai",
        reference=REFERENCE
    output:
        out_dir=directory(f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}")
    params:
        flanking_length_bp=FLANKING_LENGTH_BP,
	haploid_chrs=HAPLOID_CHRS
    conda:
        "conda_env_yaml/TR_longTR_methylation.yaml"
    shell:
        """
        bash scripts/TR-longTR-methylation_v3.sh \
	-v {input.TR_vcf} \
	-r {input.reference} \
	-i {input.phased_bam} \
	-o {output.out_dir} \
	-s {wildcards.sample} \
	-f {params.flanking_length_bp} \
	-h {params.haploid_chrs}
        """


# Non reference TE methylation calling
rule non_ref_TE_methylation_calling:
    input:
        TLDR_txt_output=f"{OUTPUT_DIR}/05_TE_calling/{{sample}}/non_ref_TE/{{sample}}.table.txt"
    output:
        output_file=f"{OUTPUT_DIR}/05_TE_calling/{{sample}}/non_ref_TE/{{sample}}/{{sample}}.table.pass.summary.meth.phased.txt"
    params:
        out_dir=f"{OUTPUT_DIR}/05_TE_calling/{{sample}}/non_ref_TE/{{sample}}",
        flanking_length_bp=FLANKING_LENGTH_BP
    conda:
        "conda_env_yaml/modkit.yaml"
    shell:
        """
        bash scripts/TLDR-methylation_v1.sh \
        -i {input.TLDR_txt_output} \
        -o {params.out_dir} \
        -s {wildcards.sample} \
        -f {params.flanking_length_bp} \
        """

# Non reference TE methylation calling
rule ref_TE_methylation_calling:
    input:
        phased_meth_bed_gz=f"{OUTPUT_DIR}/04_methylation_calling/phased/{{sample}}_haplotype_1.bed.gz",
        snp_vcf=f"{OUTPUT_DIR}/03_phasing/{{sample}}/{{sample}}_phased.vcf.gz"
    output:
        output_file=f"{OUTPUT_DIR}/05_TE_calling/{{sample}}/ref_TE/{{sample}}_summary_per_ref_{TYPE_OF_TE}.txt"
    params:
        phased_dir=f"{OUTPUT_DIR}/04_methylation_calling/phased",
        unphased_dir=f"{OUTPUT_DIR}/04_methylation_calling/unphased",
        phased_variation_dir=f"{OUTPUT_DIR}/03_phasing/{{sample}}",
        out_dir=f"{OUTPUT_DIR}/05_TE_calling/{{sample}}/ref_TE",
        TE_catalog=TE_CATALOG, 
	type_of_TE=TYPE_OF_TE,
        flanking_length_bp=FLANKING_LENGTH_BP
    conda:
        "conda_env_yaml/modkit.yaml"
    shell:
        """
        bash scripts/refTE-meth_v1.sh \
        -p {params.phased_dir} \
        -u {params.unphased_dir} \
        -v {params.phased_variation_dir} \
        -c {params.TE_catalog} \
	-t {params.type_of_TE} \
	-s {wildcards.sample} \
	-o {params.out_dir} \
	-f {params.flanking_length_bp}
        """

