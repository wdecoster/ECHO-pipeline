# ECHO: a nanopore sequencing-based workflow for (epi)genetic profiling of the human repeatome

## Introduction 

Repetitive DNA elements make up more than half of the human genome and include both tandem repeats (TRs) and transposable elements (TEs). These elements are highly polymorphic and tightly regulated by epigenetic mechanisms such as DNA methylation. They play key roles in genome regulation, evolution, and disease. However, their repetitive nature has made them difficult to analyze with short-read sequencing approaches.
Oxford Nanopore long-read sequencing provides the unique ability to span full-length repeat regions while simultaneously detecting native DNA methylation. This opens the door to comprehensive analyses of both the genetic and epigenetic landscape of the human repeatome in a single experiment.

Here we introduce **ECHO**, a comprehensive Snakemake-based pipeline for the (**E**pi)genomic **C**haracterisation of **H**uman Repetitive Elements using **O**xford Nanopore Sequencing. It integrates state-of-the-art tools for QC, mapping, variant detection, phasing and methylation calling into a single reproducible workflow. With dedicated modules for both TR and TE analysis, ECHO enables joint profiling of sequence variation and CpG methylation across the full spectrum of repetitive elements. 

## Pipeline overview


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

