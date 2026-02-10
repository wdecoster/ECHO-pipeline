# ECHO: a nanopore sequencing-based workflow for (epi)genetic profiling of the human repeatome

---
## TABLE OF CONTENTS
- [Introduction](#introduction)
- [Pipeline overview](#pipeline-overview)
- [Setting up](#setting-up-the-pipeline)
  - [Installation](#installation)
  - [Input Files](#prepare-input-files)
  - [Directory Structure](#project-directory-structure)
  - [Configuration](#set-up-configuration-file)
- [Running the pipeline](#running-the-pipeline)
- [Citation](#citation)
- [Questions](#questions)

## INTRODUCTION

Repetitive DNA elements make up more than half of the human genome and include both tandem repeats (TRs) and transposable elements (TEs). These elements are highly polymorphic and tightly regulated by epigenetic mechanisms such as DNA methylation. They play key roles in genome regulation, evolution, and disease. However, their repetitive nature has made them difficult to analyze with short-read sequencing approaches.
Oxford Nanopore long-read sequencing provides the unique ability to span full-length repeat regions while simultaneously detecting native DNA methylation. This opens the door to comprehensive analyses of both the genetic and epigenetic landscape of the human repeatome in a single experiment.

Here we introduce **ECHO**, a comprehensive [Snakemake](https://snakemake.readthedocs.io/en/stable/)-based pipeline for the (**E**pi)genomic **C**haracterisation of **H**uman Repetitive Elements using **O**xford Nanopore Sequencing. It integrates state-of-the-art tools for QC, mapping, variant detection, phasing and methylation calling into a single reproducible workflow. With dedicated modules for both TR and TE analysis, ECHO enables joint profiling of sequence variation and CpG methylation across the full spectrum of repetitive elements. 

---

## PIPELINE OVERVIEW
### **Schematic overview**
![Pipeline schematic](https://github.com/leenput/repeatome_pipeline/)

For more information on all the tools used, see [`docs/tools.md`](docs/tools.md). 

---

## SETTING UP THE PIPELINE

### Installation
  
To install the ECHO pipline, use:
```bash
git clone https://github.com/leenput/repeatome_pipeline.git # clone the repository
cd repeatome_pipeline
bash scripts/download_repeat_catalogs.sh # download the repeat catalogs in current directory
```
  
### **Prepare input files**  
To run the pipeline, ONT input files must be in one of the following formats:  

| Format  | Description                                                                       |
| ------- | --------------------------------------------------------------------------------- |
| `.pod5` | Raw signal-level data (GPU required for basecalling)                              |
| `.ubam` | Pre-basecalled, unaligned data (dorado, methylation-aware model: sup,5mCG\_5hmCG) |
| `.bam`  | Pre-basecalled (dorado) + aligned ONT data (aligned to GRCh38 or T2T-CHM13v2)              |


### Project directory structure

The project directory will be organized in the following way:
```
projectID/                        
├── 00_raw_data/                  
├── 01_alignment/             
├── 02_variant_calling/      
├── 03_phasing/              
├── 04_methylation_calling/
├── 05_non_ref_TE_calling/
├── 06_TR_calling/
├── 07_TR_methylation_calling/
├── 08_TE_methylation_calling/
├── qc/
├── logs/
```

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
- project input directory (*format: /ifs/data/research/unique/projects/projectID*)  
- project output directory (preferably the same of input directory)
- reference genome path (path to GRCh38 or T2T-CHM13v2 fasta file) 
- TE catalog path 
- TR catalog path
- length of flanking regions for TE and TR analysis

*Note: an example `config.yaml` file is provided in profiles/slurm_profile/. You can create your own slurm profile directory in profiles/ and copy the config.yaml there an customize it for your own analysis.*

---

## RUNNING THE PIPELINE
  
### *Abacus* HPC environment

In abacus, first load the required Conda and Singularity environments:

```bash
module load bioinf/conda
. /cm/shared/apps/bioinf/conda/23.10.0/etc/profile.d/conda.sh
conda activate /ifs/software/research/unique/leena/conda-envs/snakemake_env_v9  # contains snakemake v9.13.4
module load bioinf/singularity
```

### *Other* HPC environment
For other HPC systems, ensure that the following are installed and available in your environment:

- Conda (version ≥23.3, tested on 23.10.0)
- Singularity (version ≥3.7 and <4.0, tested on 3.7.0)
- Snakemake (version ≥7.0 and <9.0, tested on 7.32.4)
 
### Launch the pipeline
  
Finally, to launch the pipeline, use the following command:

```bash
snakemake -s workflow/snakefile --profile profiles/slurm_profile 
```

---
## REPEAT CATALOGS
For a detailed description of the repeat catalogs bundled with ECHO, see  
📄 [`docs/repeat_catalogs.md`](docs/repeat_catalogs.md)


---
## OUTPUT

In your project folder, numerous output files are provided, with the most important ones explained [here](docs/output.md)

---
## CITATION

TBD

---

## QUESTIONS?
Please leave any feedback, issue or question on the [Issues section](https://github.com/leenput/repeatome_pipeline/issues).   
