# ECHO: a nanopore sequencing-based workflow for (epi)genetic profiling of the human repeatome

## Introduction 

Repetitive DNA elements make up more than half of the human genome and include both tandem repeats (TRs) and transposable elements (TEs). These elements are highly polymorphic and tightly regulated by epigenetic mechanisms such as DNA methylation. They play key roles in genome regulation, evolution, and disease. However, their repetitive nature has made them difficult to analyze with short-read sequencing approaches.
Oxford Nanopore long-read sequencing provides the unique ability to span full-length repeat regions while simultaneously detecting native DNA methylation. This opens the door to comprehensive analyses of both the genetic and epigenetic landscape of the human repeatome in a single experiment.

Here we introduce **ECHO**, a comprehensive [Snakemake](https://snakemake.readthedocs.io/en/stable/)-based pipeline for the (**E**pi)genomic **C**haracterisation of **H**uman Repetitive Elements using **O**xford Nanopore Sequencing. It integrates state-of-the-art tools for QC, mapping, variant detection, phasing and methylation calling into a single reproducible workflow. With dedicated modules for both TR and TE analysis, ECHO enables joint profiling of sequence variation and CpG methylation across the full spectrum of repetitive elements. 

## Pipeline overview
### **Schematic overview**
![Pipeline schematic](https://github.com/leenput/repeatome_pipeline/blob/readme/docs/figures/DAG-pipeline.jpg)

### **Overview of tools**
**Raw read manipulation**  
- Basecalling with [dorado](https://github.com/nanoporetech/dorado) (*optional, GPU availability needed*)   
- Conversion from ubam to fastq format using [samtools](https://www.htslib.org/) (*optional*)     
- Read filtering using [chopper](https://github.com/wdecoster/chopper) (*optional*)   

**QC**  
- pre- and post-filtered read QC, evaluation of mapping metrics and phasing with [cramino](https://github.com/wdecoster/cramino) and [nanoplot](https://github.com/wdecoster/NanoPlot)  

**Alignment**  
- Read alignment to human genome reference ([GRCh38](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/)/[T2T-CHM13v2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/)) with [minimap2](https://github.com/lh3/minimap2)    

**Variant calling**  
- Small variant calling (SNVs and Indels) using [Clair3](https://github.com/HKU-BAL/Clair3)  
- Structial variant (SV) calling using [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)  
- Filter variants with [BCFtools](https://samtools.github.io/bcftools/bcftools.html)  

**Phasing**    
- Phase and haplotag reads with [LongPhase](https://github.com/twolinin/longphase)  

**Methylation**  
- Generate methylation pileups at CpG sites with [modkit](https://github.com/nanoporetech/modkit)  

**TR characterization**  
- TR genotyping using [LongTR](https://github.com/gymrek-lab/LongTR)    
- Using a custom [script](workflow/scripts/TR-longTR-methylation_v4.sh), methylation information is incorporated for each TR allele called by LongTR

**TE characterization**  
- non-ref TE analysis: identification of TE insertions not present in the reference genome using [TLDR](https://github.com/adamewing/tldr), followed by methylation information extraction using [script](workflow/scripts/TLDR-methylation_v2.sh)  
- ref TE analysis: [script](workflow/scripts/ref_TE_avg_meth.sh) to analyse and summarize sequence variants and methylation information across annotated TEs in the reference genome   


## Manual to Launch Snakemake Pipeline

## 🧪 Environment Setup

To run the Snakemake pipeline, first load the required Conda and Singularity environments:

```bash
module load bioinf/conda
. /cm/shared/apps/bioinf/conda/23.10.0/etc/profile.d/conda.sh
conda activate /ifs/software/research/unique/leena/conda-envs/snakemake-env
module load bioinf/singularity
```

---

## 📥 Input File Requirements

To run the pipeline, input files must be in one of the following formats:

- `.pod5`
- `.ubam`
- `.bam`

---

## 🗂 Project Directory Structure

If you start the pipeline from `pod5`, ensure that the `pod5` files follow this directory structure:

```
/ifs/data/research/unique/projects/{project_name}/00_raw_data/pod5/{sample_name}
```

If you start the pipeline from `ubam`, the `ubam` file should be in this directory:

```
/ifs/data/research/unique/projects/{project_name}/00_raw_data/basecalled/ubam/{sample_name}
```

If you start the pipeline from `bam`, the `bam` and the index `bai` files should be in this directory:

```
/ifs/data/research/unique/projects/{project_name}/01_alignment/{sample_name}
```



Replace `{project_name}` and `{sample_name}` with your actual project and sample identifiers.

---

## ⚙️ Configuration File

Before running the pipeline, you need to create a `config.yaml` file that includes the following:

- Sample ID
- Input format (`.pod5`, `.ubam`, or `.bam`)
- Input directory
- Output directory (preferably the same of input directory)
- Reference genome path
- TE catalog path
- TR catalog path
- Length of flanking regions for TE and TR analysis

> An example `config.yaml` file is provided in profiles/slurm_profile/. You can create your own slurm profile directory in profiles/ and copz the config.yaml there an customize it for your own analysis.

---

## 🚀 Running the Pipeline

To run the pipeline, use the following command:

```bash
snakemake -s workflow/snakefile --profile profiles/slurm_profile 
```

---

