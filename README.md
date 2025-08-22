# ECHO: a nanopore sequencing-based workflow for (epi)genetic profiling of the human repeatome

---
## TABLE OF CONTENTS
- [Introduction](#introduction)
- [Pipeline overview](#pipeline-overview)
- [Repeat catalogs](#repeat-catalogs)
- [Setting up](#setting-up-the-pipeline)
  - [Installation](#installation)
  - [Input Files](#prepare-input-files)
  - [Directory Structure](#project-directory-structure)
  - [Configuration](#set-up-configuration-file)
- [Running the pipeline](#running-the-pipeline)
- [Outputs](#outputs)
- [Questions](#questions)

## INTRODUCTION

Repetitive DNA elements make up more than half of the human genome and include both tandem repeats (TRs) and transposable elements (TEs). These elements are highly polymorphic and tightly regulated by epigenetic mechanisms such as DNA methylation. They play key roles in genome regulation, evolution, and disease. However, their repetitive nature has made them difficult to analyze with short-read sequencing approaches.
Oxford Nanopore long-read sequencing provides the unique ability to span full-length repeat regions while simultaneously detecting native DNA methylation. This opens the door to comprehensive analyses of both the genetic and epigenetic landscape of the human repeatome in a single experiment.

Here we introduce **ECHO**, a comprehensive [Snakemake](https://snakemake.readthedocs.io/en/stable/)-based pipeline for the (**E**pi)genomic **C**haracterisation of **H**uman Repetitive Elements using **O**xford Nanopore Sequencing. It integrates state-of-the-art tools for QC, mapping, variant detection, phasing and methylation calling into a single reproducible workflow. With dedicated modules for both TR and TE analysis, ECHO enables joint profiling of sequence variation and CpG methylation across the full spectrum of repetitive elements. 

---

## PIPELINE OVERVIEW
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

---

## REPEAT CATALOGS
The ECHO pipeline includes multiple repeat catalogs specifically designed to capture the repetitive elements of interest, both tandem repeats (TRs) and transposable elements (TEs), for the selected reference genome (GRCh38 or T2T CHM13v2). These catalogs, which are hosted within this [zenodo repository](https://zenodo.org/records/16925640), are compiled from published resources and adapted to ensure compatibility with the pipeline. They define the genomic loci where genotyping and/or methylation profiling is performed, enabling analysis of the human repeatome. 

### Tandem repeats
For tandem repeat (TR) analysis, catalogs specify the 1-based genomic coordinates of repeat loci to be profiled in the following format.
```bash
chr  start  end  TRmotif  TR_ID
```

The following catalogs are provided with the pipeline:
- **STRs associated with human diseases**:
    - 73 loci
    - available for GRCh38 and T2T-CHM13v2
    - based on [STRchive](https://strchive.org/loci/) (v2.2.1)
- **forensically relevant STRs**:
    - 77 loci
    - available for GRCh38
    - based on [STRbase](https://strbase.nist.gov/Loci) (v2.0)
- **genome-wide TR panel**:
    - 1,177,430 loci
    - available for GRCh38
    - derived from ~1.8M TRs in the [project adotto catalog v1.2.1](https://zenodo.org/records/13987414), released as part of the GIAB tandem repeat benchmark variant set ([English et al. 2025](https://www.nature.com/articles/s41587-024-02225-z). This panel was filtered to remove homopolymers, repeats with motif lengths >100 bp (STRs and VNTRs), and regions containing multiple overlapping TRs, which cannot be robustly evaluated with Truvari (inspired by approach of [Loughleed et al., 2025](https://www.biorxiv.org/content/10.1101/2025.03.25.645269v1.full).
- **STRs with CpG sites in their repeat motif (1-6 bp)**:
    - 7,394 loci
    - available for GRCh38
    - derived from the genome-wide panel above with additional filtering for TRs with motif lengths < 7bp and the occurance of CpG sites in their motifs 

### Transposable elements
TE catalogues are derived from RepeatMasker annotations, which were obtained from the UCSC Genome Browser:  

- **GRCh38**: [hg38.fa.out.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz)  
- **T2T-CHM13v2**: [hs1.repeatMasker.out.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz)  

The RepeatMasker (RM) outputs were filtered to retain only bona fide TEs and uncertain classifications (entries containing “?”) were removed from the panel.

ECHO provides:  
- a **genome-wide TE catalog** covering all annotated TEs,  
- class-specific catalogs for **DNA transposons, RC/Helitrons, LINEs, SINEs, LTR retrotransposons, and SVAs**, based on the RM classification.  

In addition, we provide the file `reref.ont.human.fa`, a FASTA reference used by the TLDR tool to annotate the most relevant TE families in the human genome.

---

## SETTING UP THE PIPELINE

### Installation
  
To obtain the ECHO pipline, use:
```bash
git clone https://github.com/leenput/repeatome_pipeline.git # clone the repository
cd repeatome_pipeline
bash workflow/scripts/download_catalogs.sh # download the repeat catalogs
```
  
### **Prepare input files**  
To run the pipeline, ONT input files must be in one of the following formats:  

| Format  | Description                                                                       |
| ------- | --------------------------------------------------------------------------------- |
| `.pod5` | Raw signal-level data (GPU required for basecalling)                              |
| `.ubam` | Pre-basecalled, unaligned data (dorado, methylation-aware model: sup,5mCG\_5hmCG) |
| `.bam`  | Pre-basecalled (dorado) + aligned ONT data (aligned to GRCh38 or T2T-CHM13v2)              |


### Project directory structure

Depending on start point of the pipeline, ensure that your input files are stored in the following data structures (paths are tailored to usage on our Abacus HPC system):

- `pod5`:   

```
/ifs/data/research/unique/projects/{project_name}/00_raw_data/pod5/{sample_name}/<your-file.pod5>
```

- `ubam`:  

```
/ifs/data/research/unique/projects/{project_name}/00_raw_data/basecalled/ubam/{sample_name}/<your-file.bam>
```

- `bam` (with index `.bai`):  

```
/ifs/data/research/unique/projects/{project_name}/01_alignment/{sample_name}/<your-file.bam>
/ifs/data/research/unique/projects/{project_name}/01_alignment/{sample_name}/<your-file.bam.bai>
```

Replace `{project_name}` and `{sample_name}` with your actual project and sample identifiers.
  

### Set up configuration file
Before running the pipeline, you need to create a `config.yaml` file that includes the following:

- Sample ID
- Input format (`.pod5`, `.ubam`, or `.bam`)
- project input directory
- project output directory (preferably the same of input directory)
- reference genome path
- TE catalog path
- TR catalog path
- length of flanking regions for TE and TR analysis

*Note: an example `config.yaml` file is provided in profiles/slurm_profile/. You can create your own slurm profile directory in profiles/ and copz the config.yaml there an customize it for your own analysis.*

---

## RUNNING THE PIPELINE
  
### Abacus HPC environment

In abacus, first load the required Conda and Singularity environments:

```bash
module load bioinf/conda
. /cm/shared/apps/bioinf/conda/23.10.0/etc/profile.d/conda.sh
conda activate /ifs/software/research/unique/leena/conda-envs/snakemake-env
module load bioinf/singularity
```
  
### Launch the pipeline
  
Finally, to launch the pipeline, use the following command:

```bash
snakemake -s workflow/snakefile --profile profiles/slurm_profile 
```

---

## QUESTIONS?
Please leave any feedback, issue or question on the [Issues section](https://github.com/leenput/repeatome_pipeline/issues).   
