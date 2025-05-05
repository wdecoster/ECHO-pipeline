configfile: "config.yaml"


SAMPLES = config["samples"]
START_FROM = config.get("start_from", "pod5")
OUTPUT_DIR = config["output_dir"]
INPUT_DIR = config["input_dir"]
REFERENCE = config["reference"]
TR_CATALOG = config["tr_catalog"]
TE_CATALOG = config["te_catalog"]
FLANKING_LENGTH_BP = config["flanking_length_bp"]
HAPLOID_CHRS = config["haploid_chrs"]

# debugging in case --config is not parsed correctly
print(f"START_FROM = {START_FROM}")

all_inputs = []

if START_FROM == "pod5":
    all_inputs.extend([
        expand(f"{OUTPUT_DIR}/00_raw_data/basecalled/ubam/{{sample}}_unaligned.bam", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/qc/qc_basecalling/{{sample}}_summary.tsv", sample=SAMPLES),
    ])
elif START_FROM == "ubam":
    all_inputs.extend([
        expand(f"{OUTPUT_DIR}/qc/unaligned/{{sample}}_unaligned_bam_NanoPlot-report.html", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/qc/unaligned/QC_LRS_{{sample}}_fastq.html", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/qc/aligned/{{sample}}_LRS_aligned_bam.html", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/qc/aligned/{{sample}}_aligned_bam_NanoPlot-report.html", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/qc/aligned/{{sample}}_coverage.chr.stat.gz", sample=SAMPLES),
    ])
elif START_FROM == "bam":
    # no basecalling or alignment outputs expected
    pass
else:
    print(f"ERROR: Invalid 'start_from' value in config.yaml: {START_FROM}")
    exit(1)
    # This is redundant and will never be executed since you define the fallback value in config.get("start_from", "pod5")

# /ifs/data/research/unique/projects/snakemake_test/QC_reports
all_inputs.extend([
    expand(f"{OUTPUT_DIR}/02_variant_calling/SNPs_Indels/{{sample}}/phased_merge_output.vcf.gz", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}_SV_unphased.vcf.gz", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/03_phasing/{{sample}}_phased_alignment.bam", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}_haplotype_1.bed", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}_haplotype_2.bed", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}_haplotype_ungrouped.bed", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/05_TE_calling/{{sample}}", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{{sample}}_TRs.vcf.gz", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}", sample=SAMPLES),
    REFERENCE
])

rule all:
    input:
        all_inputs


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


# QC for unaligned and Alignment 
if START_FROM in ["pod5", "ubam"]:
    
    # QC for unaligned bam with -- NanoPlot --
    rule QC_NanoP_unaligned:
        input:
            unaligned_bam=f"{OUTPUT_DIR}/00_raw_data/basecalled/ubam/{{sample}}_unaligned.bam",
        params:
            results_dir=f"{OUTPUT_DIR}/qc/unaligned/",
            name_id=f"unaligned_bam"
        output:
            summary_stats=f"{OUTPUT_DIR}/qc/unaligned/{{sample}}_unaligned_bam_NanoStats.txt",
            html_report=f"{OUTPUT_DIR}/qc/unaligned/{{sample}}_unaligned_bam_NanoPlot-report.html"
        threads:24
        log:
            f"{OUTPUT_DIR}/qc/log/NanoP_{{sample}}_unaligned.log"
        conda:
            "conda_env_yaml/Nanoplot_env.yaml"
        shell:
            """
            NanoPlot -t {threads} \
                -o {params.results_dir} \
                --dpi 250 -c green -f png --N50 \
                --title {wildcards.sample}_{params.name_id} \
                --prefix {wildcards.sample}_{params.name_id}_ \
                --ubam {input.unaligned_bam} \
                >& {log}
            """
        
    # ALIGNMENT 
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
# QC for fastq file -- LongReadSum --
    rule QC_LRS_fastq:
        input:
            sum_txt=f"{OUTPUT_DIR}/00_raw_data/basecalled/fastq/{{sample}}.fastq"
        output:
            sum_report=f"{OUTPUT_DIR}/qc/unaligned/QC_LRS_{{sample}}_fastq.html",
            sum_txt=f"{OUTPUT_DIR}/qc/unaligned/QC_LRS_{{sample}}_summary.txt"
        params:
            name_id=f"QC_LRS_{{sample}}_",
            out_dir=f"{OUTPUT_DIR}/qc/unaligned"
        log:
            f"{OUTPUT_DIR}/qc/log/LRS_{{sample}}_unaligned.log"
        threads: 24
        conda:
            "conda_env_yaml/Longreadsum_env.yaml"
        shell:
            """
            longreadsum fq \
                -i {input} \
                -o {params.out_dir} \
                -t {threads} \
                --outprefix {params.name_id} \
                -u 33 
            mv {params.out_dir}/FASTQ_summary.txt {params.out_dir}/{params.name_id}summary.txt
            """
            # -u quality offset for bases in fastq, default 33
            # 8 min for HG001_subset, 
    
    
    
    # QC for aligned bam with -- LongReadSum --
    rule QC_LRS_aligned:
        input:
            aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam"
        params:
            results_dir=f"{OUTPUT_DIR}/qc/aligned", # not aligned/ that will fail the commandline separation in shell
            name_id=f"{{sample}}_LRS_aligned_"
        output:
            summary_bam=f"{OUTPUT_DIR}/qc/aligned/{{sample}}_LRS_aligned_bam_summary.txt", 
            report=f"{OUTPUT_DIR}/qc/aligned/{{sample}}_LRS_aligned_bam.html"
        threads: 24
        log:
            f"{OUTPUT_DIR}/qc/log/LRS_{{sample}}_aligned.log"
        conda:
            "conda_env_yaml/Longreadsum_env.yaml"
        shell:
            """
            longreadsum bam -i {input.aligned_bam} -o {params.results_dir} \
                --log {log} \
                -t {threads} \
                --outprefix {params.name_id}
            mv {params.results_dir}/bam_summary.txt {params.results_dir}/{params.name_id}bam_summary.txt
            """

    # QC for aligned bam with -- NanoPlot --
    rule QC_NanoP_aligned:
        input:
            aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam"
        params:
            results_dir=f"{OUTPUT_DIR}/qc/aligned/",
            name_id=f"aligned_bam"
        output:
            summary_stats=f"{OUTPUT_DIR}/qc/aligned/{{sample}}_aligned_bam_NanoStats.txt",
            html_report=f"{OUTPUT_DIR}/qc/aligned/{{sample}}_aligned_bam_NanoPlot-report.html" 
        threads: 24
        log:
            f"{OUTPUT_DIR}/qc/log/NanoP_{{sample}}_aligned.log"
        conda:
            "conda_env_yaml/Nanoplot_env.yaml"
        shell:
            """
            NanoPlot -t {threads} \
                -o {params.results_dir} \
                --dpi 250 -c blue -f png --N50 \
                --title {wildcards.sample}_{params.name_id} \
                --prefix {wildcards.sample}_{params.name_id}_ \
                --bam {input.aligned_bam} \
                >& {log}
            """
            
    # QC coverage calculation with -- PanDepth --
    rule coverage:
        input:
            aligned_bam=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam"
        output:
            coverage_zip=f"{OUTPUT_DIR}/qc/aligned/{{sample}}_coverage.chr.stat.gz" # .chr.stat.gz als default suffix
        params:
            prefix=f"{OUTPUT_DIR}/qc/aligned/{{sample}}_coverage"
        conda:
            "conda_env_yaml/PanDepth_env.yaml"
        threads: 24
        shell:
            """
            pandepth \
                -i {input} \
                -t {threads} \
                -o {params.prefix}
            """
            # cd /ifs/software/research/unique/leonard/tools/PanDepth/bin --> only if included in PATH variable of bashrc



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
        out_prefix=f"{OUTPUT_DIR}/03_phasing/{{sample}}"
    output:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/03_phasing/{{sample}}_phased_alignment.bam.bai",
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
        "conda_env_yaml/modkit.yaml"
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
        TE_catalog=TE_CATALOG,
        reference=REFERENCE
    output:
        out_dir=directory(f"{OUTPUT_DIR}/05_TE_calling/{{sample}}")
    params:
        flanking_length_bp=FLANKING_LENGTH_BP
    conda:
        "conda_env_yaml/tldr.yaml"
    shell:
        """
        /ifs/software/research/unique/pipeline_tools/tldr/tldr \
        -b {input.phased_bam} \
            -e {input.TE_catalog} \
            -r {input.reference} \
        -p 32 \
            --outbase {output.out_dir} \
            --detail_output \
            --methylartist \
            --max_cluster_size 500 \
            --extend_consensus {params.flanking_length_bp}
        """


# Tandem Repeats calling
rule TR_calling:
    input:
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam.bai",
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
        vcf_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam.bai",
        phased_bam=f"{OUTPUT_DIR}/03_phasing/{{sample}}_phased_alignment.bam",
        phased_bam_index=f"{OUTPUT_DIR}/01_alignment/{{sample}}_sorted.bam.bai",
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
