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
![Pipeline schematic](docs/DAG-pipeline.jpg)

For more information on all the tools used, see [`docs/tools.md`](docs/tools.md). 

---

## SETTING UP THE PIPELINE

--

### Installation
  
To install the ECHO pipline, use:
```bash
git clone https://github.com/leenput/repeatome_pipeline.git # clone the repository
cd repeatome_pipeline
bash scripts/download_repeat_catalogs.sh # download the repeat catalogs in current directory
```

---
  
### **Prepare input files**  
To run the pipeline, ONT input files must be in one of the following formats:  

| Format  | Description                                                                       |
| ------- | --------------------------------------------------------------------------------- |
| `.pod5` | Raw signal-level data (GPU required for basecalling)                              |
| `.ubam` | Pre-basecalled, unaligned data (dorado, methylation-aware model: sup,5mCG\_5hmCG) |
| `.bam`  | Pre-basecalled (dorado) + aligned ONT data (aligned to GRCh38 or T2T-CHM13v2)              |

---

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

Depending on start point of the pipeline, ensure that your input files are stored in the following data structures:

- `pod5`:   

```
/path/to/projects/{project_name}/00_raw_data/pod5/{sample_name}/<your-file.pod5>
```

- `ubam`:  

```
/path/to/projects/{project_name}/00_raw_data/basecalled/ubam/{sample_name}/<your-file.bam>
```

- `bam` (with index `.bai`):  

```
/path/to/projects/{project_name}/01_alignment/{sample_name}/<your-file.bam>
/path/to/projects/{project_name}/01_alignment/{sample_name}/<your-file.bam.bai>
```

Replace `{project_name}` and `{sample_name}` with your actual project and sample identifiers.
  
---

### Set up configuration file

Before running the pipeline, you must generate a `config.yaml` file.
This is done using the provided helper script `make_config_tiny.py`,
which creates and validates a Snakemake configuration for ECHO.

The configuration defines:
- Sample IDs
- Input format (`.pod5`, `.ubam`, or `.bam`)
- Project input and output directories
- Reference genome (GRCh38 or T2T-CHM13v2)
- TR and TE catalogs
- Key analysis parameters (e.g. flanking length, read filters)  

A minimal configuration using the bundled ECHO repeat catalogs **defaults** can be generated automatically:

```bash
python scripts/make_config_tiny.py init \
  --output configs/<config-name>.yaml \
  --samples SAMPLE1 SAMPLE2 \
  --start-from ubam \
  --input-dir /path/to/project \
  --output-dir /path/to/project \
  --reference /path/to/GRCh38.fa \
  --reference-name GRCh38 \
  --use-bundled-db
```

📄 For full configuration details and advanced usage, see  
[`docs/configuration.md`](docs/configuration.md)


#### Default behaviour (bundled database mode)

When using the configuration script with `--use-bundled-db`, ECHO applies the
following defaults unless explicitly overridden:

| Category | Setting | Default value | Notes |
|--------|--------|---------------|-------|
| Reference | Reference build | As specified by `--reference-name` | `GRCh38` or `T2T-CHM13v2` |
| TR analysis | TR catalog | Adotto genome-wide longTR | Sensitive, genome-wide TR catalog |
| TR analysis | TR type | `genome-wide` | Used for output folder naming |
| TR methylation | CpG filtering | Enabled | Restricts analysis to canonical STRs containing CpGs |
| TE analysis | TE catalog | Genome-wide (`all`) | All TE classes from UCSC RepeatMasker |
| Read filtering | Minimum read quality | 7 | Applied to FASTQ/UBAM inputs |
| Read filtering | Minimum read length | 500 bp | Shorter reads are discarded |
| Analysis | Flanking region length | 250 bp | Used for TR and TE analyses |
| Analysis | Repeat consensus extension | 1000 bp | Extension for repeat consensus building |

These defaults are chosen to provide a **sensitive, genome-wide analysis**
while keeping computational requirements manageable.
All defaults can be modified via the configuration script, see [here](docs/configuration.md).

---

### Configure the execution profile

ECHO is executed using a Snakemake *profile*, which defines how jobs are submitted
to the compute environment (e.g. SLURM settings, partitions, resources).

After generating `<config-name>.yaml`, you must update the profile configuration so that:  
1. the profile points to the generated `config.yaml`  
2. cluster-specific settings match your local compute infrastructure  

#### Point the profile to your config.yaml

In the profile directory (e.g. `profiles/slurm_profile/`), edit `config.yaml`
so that it references the configuration file you generated:

```yaml
configfile: /full/path/to/your/<config-name>.yaml
```

#### Adjust cluster-specific settings
In the same profile config, you should adapt cluster-specific settings to match your available  
computational infrastructure, including:  
- singularity bind mounts to the project folder 
- partition or queue names  
- default memory or runtime limits  

---

## RUNNING THE PIPELINE
  
### *Abacus* HPC environment

In Abacus, first load the required Conda and Singularity environments:

```bash
module load bioinf/conda
. /cm/shared/apps/bioinf/conda/23.10.0/etc/profile.d/conda.sh
conda activate /ifs/software/research/unique/leena/conda-envs/snakemake_env_v9  # contains snakemake v9.13.4
module load bioinf/singularity
```
And then launch it:
```bash
snakemake -s workflow/snakefile --profile profiles/slurm_profile
```

### *Other* HPC environment
Ensure the following are installed and available in your environment:

- Conda (version ≥23.3, tested on 23.10.0)
- Singularity (version ≥3.7 and <4.0, tested on 3.7.0)
- Snakemake (version ≥7.0 and <9.0, tested on 7.32.4)
   
To launch the pipeline, use the following command:

```bash
snakemake -s workflow/snakefile --profile profiles/slurm_profile 
```

### *Local* execution

ADD

---

## REPEAT CATALOGS
 
For a detailed description of the repeat catalogs bundled with ECHO, or how to use custom catalogs, see  
📄 [`docs/repeat_catalogs.md`](docs/repeat_catalogs.md)

---

## OUTPUT FILES 

In your project folder, numerous output files are provided, with the most important ones explained [here](docs/output.md)

---

## TEST DATA

ADD

---
## CITATION

ADD

---

## QUESTIONS?

Please leave any feedback, issue or question on the [Issues section](https://github.com/leenput/repeatome_pipeline/issues).   
