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
![Pipeline schematic](https://github.com/leenput/repeatome_pipeline/blob/a607e743ae1e6a5f04cc4b91df79611345c5e504/DAG-pipeline.jpg)

### **Overview of tools**
**Raw read manipulation**  
- Basecalling with [dorado](https://github.com/nanoporetech/dorado) (*optional, GPU availability needed*)   
- Conversion from ubam to fastq format using [samtools](https://www.htslib.org/) (*optional*)     
- Read filtering using [chopper](https://github.com/wdecoster/chopper) (*optional*)   

**QC**  
- Pre- and post-filtered read QC and evaluation of mapping and phasing metrics with [cramino](https://github.com/wdecoster/cramino) and [nanoplot](https://github.com/wdecoster/NanoPlot)  

**Alignment**  
- Read alignment to human reference genome ([GRCh38](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/)/[T2T-CHM13v2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/)) with [minimap2](https://github.com/lh3/minimap2)    

**Variant calling**  
- Small variant calling (SNVs and Indels) using [Clair3](https://github.com/HKU-BAL/Clair3)  
- Structural variant (SV) calling using [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)  
- Filter variants with [BCFtools](https://samtools.github.io/bcftools/bcftools.html)  

**Phasing**    
- Phase and haplotag reads with [LongPhase](https://github.com/twolinin/longphase)  

**Methylation**  
- Generate methylation pileups at CpG sites with [modkit](https://github.com/nanoporetech/modkit)  

**TR characterization**  
- TR genotyping using [LongTR](https://github.com/gymrek-lab/LongTR)    
- Using a custom [script](workflow/scripts/TR-longTR-methylation_v4.sh), methylation information is incorporated for each TR allele called by LongTR
- Motif decomposition using [uTR](https://github.com/morisUtokyo/uTR).

**TE characterization**  
- non-ref TE analysis: identification of TE insertions not present in the reference genome using [TLDR](https://github.com/adamewing/tldr), followed by methylation information extraction using [script](workflow/scripts/TLDR-methylation_v2.sh)  
- ref TE analysis: two scripts ([script1](workflow/scripts/ref_TE_avg_meth.sh) and [script2](workflow/scripts/ref_TE_cpg_res.sh)) to analyse and summarize sequence variants and methylation information across annotated TEs in the reference genome

---

## REPEAT CATALOGS
The ECHO pipeline includes multiple repeat catalogs specifically designed to capture the repetitive elements of interest, both tandem repeats (TRs) and transposable elements (TEs), for the selected reference genome (GRCh38 or T2T CHM13v2). These catalogs, which are hosted within this [zenodo repository](https://zenodo.org/records/16925640), are compiled from published resources and adapted to ensure compatibility with the pipeline. They define the genomic loci where genotyping and/or methylation profiling is performed, enabling analysis of the human repeatome. 

### Tandem repeats
For tandem repeat (TR) analysis, catalogs specify the 1-based genomic coordinates of repeat loci to be profiled in the following format:
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
    - derived from ~1.8M TRs in the [project adotto catalog v1.2.1](https://zenodo.org/records/13987414), released as part of the GIAB tandem repeat benchmark variant set ([English et al. 2025](https://www.nature.com/articles/s41587-024-02225-z)). This panel was filtered to remove homopolymers, repeats with motif lengths >100 bp (STRs and VNTRs), and regions containing multiple overlapping TRs, which cannot be robustly evaluated with Truvari (inspired by approach of [Loughleed et al., 2025](https://www.biorxiv.org/content/10.1101/2025.03.25.645269v1.full).
- **STRs with CpG sites in their repeat motif (1-6 bp)**:
    - 7,394 loci
    - available for GRCh38
    - derived from the genome-wide panel above with additional filtering for TRs with motif lengths < 7bp and the occurence of CpG sites in their motifs 

### Transposable elements
TE catalogues are derived from RepeatMasker annotations, which were obtained from the UCSC Genome Browser:  

- **GRCh38**: [hg38.fa.out.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz)  
- **T2T-CHM13v2**: [hs1.repeatMasker.out.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz)  

The RepeatMasker (RM) outputs were filtered to retain only bona fide TEs and uncertain classifications (entries containing “?”) were removed from the panel.

ECHO can be configured to use different TE catalogs depending on the user's choice:  
- class-specific catalogs for **DNA transposons, RC/Helitrons, LINEs, SINEs, LTR retrotransposons, and SVAs**, based on the RM classification.
- a **genome-wide TE catalog** covering all annotated TEs,  

In addition, we provide the file `reref.ont.human.fa`, a FASTA reference used by the TLDR tool to annotate the most relevant TE families in the human genome.

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
## OUTPUT

In your project folder, numerous output files are provided, with the most important ones explained below.

**0. raw data files**
```
projectID/                
├── 00_raw_data/
    ├── pod5/
        ├── sampleID/
            ├── sampleID.pod5 
    ├── basecalled/
        ├── ubam
            ├── sampleID/
                ├── sampleID_unaligned.bam
        ├── fastq
            ├── sampleID/
                ├── sampleID.fastq
                ├── sampleID_filtered.fastq # optional
```

**1. alignment data**
```
projectID/
├── 01_alignment/
    ├── sampleID/
        ├── RefGenome/
            ├── sampleID_RefGenome_sorted.bam
            ├── sampleID_RefGenome_sorted.bam.bai
```

**2. variant calling**
```
projectID/                
├── 02_variant_calling/
    ├── SNVs_Indels/ 
        ├── sampleID/
            ├── RefGenome/
                ├── *.vcf.gz     # Clair3 output and intermediate files
    ├── SVs
        ├── sampleID/
            ├── RefGenome/
                ├── *.vcf.gz     # Sniffles2 ouput files 
                        
```

**3. phasing**

```
projectID/      
├── 03_phasing/
        ├── sampleID/
            ├── RefGenome/
                ├── sampleID_RefGenome_phased_alignment.bam         # phased and haplotagged alignment file 
                ├── sampleID_RefGenome_phased_alignment.bam.bai     # index of alignment file
                ├── sampleID_RefGenome_SNP_filt.vcf.gz(.tbi)              # phased VCF and index of Clair3-identified SNVs and Indels (filtered for PASS variants)
                ├── sampleID_RefGenome_SV_filt.vcf.gz(.tbi)               # phased VCF of index Sniffles2-identified SVs (filtered for PASS variants)

```

**4. whole-genome methylation pile-ups**

```
projectID/
├── 04_methylation_calling/  
        ├── sampleID/
            ├── RefGenome/
                ├── phased/
                    ├── sampleID_RefGenome_haplotype_1.bed.gz(.tbi)           # CpG methylation pileup in bedmethyl format (and index) across the genome for HP:1
                    ├── sampleID_RefGenome_haplotype_2.bed.gz(.tbi)           # CpG methylation pileup in bedmethyl format (and index) across the genome for HP:2
                    ├── sampleID_RefGenome_haplotype_ungrouped.bed.gz(.tbi)   # CpG methylation pileup in bedmethyl format (and index) across the genome for unphased data
                ├── unphased/
                    ├── sampleID_RefGenome_unphased.bed.gz(.tbi)              # haplotype-unaware CpG methylation calls in bedmethylformat (and index)

```

More information on the bedmethyl format can be found [here](https://nanoporetech.github.io/modkit/intro_pileup.html#bedmethyl-column-descriptions).  

**5. calling non-reference TE insertions**

```
projectID/
├── 05_non_ref_TE_calling/
        ├── sampleID/
            ├── RefGenome/
                ├── SampleID_RefGenome.table.txt                                   # summary table of TDLR
                ├── SampleID_RefGenome/                                            # detailed output folder of TLDR, with seperate files for each TE insertion (with own UUID)
                    ├── SampleID_RefGenome_phased_alignment.uuid.te.bam(.bai)      # local alignment file (and index) of the TE insertion
                    ├── uuid.cons.ref.fa(.fai)                                     # consensus sequence (and index) of the TE insertion                   

```
More information on the output files generated by TLDR can be obtained [here](https://github.com/adamewing/tldr).

**6. TR calling**

```
projectID/
├── 06_TR_calling/
        ├── sampleID/
            ├── RefGenome/
                ├── sampleID_RefGenome_TRs_TRCatalog.vcf.gz # Phased VCF generated by LongTR 
```

Specifics of the VCF file generated by LongTR can be found [here](https://github.com/gymrek-lab/LongTR?tab=readme-ov-file).

**7. Methylation analysis of TR regions**  
Based on the LongTR VCF output, methylation information is added and stored in the following output files:

```
projectID/
├── 07_TR_methylation_calling/
        ├── sampleID/
            ├── RefGenome/
                ├── TR_catalog/
                      ├─ sampleID_refGenome_TR_methylation.vcf          # LongTR VCF output, extended with TR characteristics and methylation information in the FORMAT fields (see below)
                      ├─ sampleID_refGenome_TR_methylation_summary.tsv  # Summary of VCF file fields, but reformatted 
                      ├─ sampleID.log                                   # log of analysis
````

- Additional fields in the methylation VCF file: 

| Field          | Description                                                         |
| ------------------- | ------------------------------------------------------------------- |
| `TR_LEN`             | Lengths of the tandem repeat alleles in bp, one per haplotype (HP1,HP2)                              |
| `TR_N_CPG`             | Number of CpG sites in the TR sequence, one per haplotype (HP1,HP2)                        |
| `TR_PATTERN`               | TR allele motif pattern, decomposed using [uTR](https://github.com/morisUtokyo/uTR) (motif HP1,motif HP2)                |
| `TR_AM` | Average methylation percentage at the TR region, one per haplotype (HP1,HP2) |
| `TR_N_METH_VALID`            | Number of valid CpG sites used for average methylation calculation across the TR region, one per haplotype, as reported by modkit (HP1,HP2)                                            |
| `UPSTREAM_TR_AM` | Average methylation percentage in upstream region of TR (length of flanking region specified in config), one per haplotype (HP1,HP2) |
| `UPSTREAM_TR_N_METH_VALID`               |  Number of valid CpG sites used for average methylation calculation across the upstream TR region, one per haplotype, as reported by modkit (HP1,HP2)                |
| `DOWNSTREAM_TR_AM` | Average methylation percentage in downstream region of TR (length of flanking region specified in config), one per haplotype (HP1,HP2) |
| `DOWNSTREAM_TR_N_METH_VALID`            |  Number of valid CpG sites used for average methylation calculation across the downstream TR region, one per haplotype, as reported by modkit (HP1,HP2)             |
| `TR_CPG_METH_HP1`  | Methylation percentages of individual CpGs in TR region of HP1, ordered by CpG position within the allele            |
| `TR_CPG_DEPTH_HP1` | Read depths for individual CpG positions in TR region of HP1, ordered by CpG position within the allele |
| `TR_CPG_METH_HP2`            | Methylation percentages of individual CpGs in TR region of HP2, ordered by CpG position within the allele                                            |
| `TR_CPG_DEPTH_HP2`  | Read depths for individual CpG positions in TR region of HP2, ordered by CpG position within the allele            |

  
  
  
**Note:** for detailed analysis of selected TRs, you can use a follow-up script [recreate-selected-TR-alignments.sh](workflow/scripts/recreate-selected-TR-alignments.sh), which will generate a local alignment file, fasta file and bed file of the TR alleles of interest, which you can use for IGV visualisation.

Details on how to use that script:  
```
Usage: bash recreate-selected-TR-alignments.sh -v <methylated.vcf.gz> -r <reference.fa> -i <phased.bam> -o <output_dir> -s <sample_id> [-e <extend>] [-h <haploid_chroms>] (-d <id1,id2,...> | -l <id_list.txt>)
  -v  VCF file (LongTR or methylated), bgzipped and indexed (.tbi)
  -r  Reference FASTA (indexed with .fai)
  -i  Phased BAM (indexed with .bai), contains HP tags
  -o  Output base directory
  -s  Sample ID (used in filenames)
  -e  Flank size to extend around TR allele when building mini-reference (default: 1000)
  -h  Comma-separated haploid chromosomes (e.g. "chrX,chrY,chrY") (default: none)
  -d  Comma-separated TR IDs to process (e.g. "chr1_15235-15322,chr7_
  -l  File with TR IDs (one per line)
```
  

**8. TE methylation analysis**
```
projectID/
├── 08_TE_methylation_calling/
    ├── sampleID/
        ├── RefGenome/
            ├── non_ref_TE/                                                  # output folder for non-reference TE methylation analysis 
                ├── sampleID_RefGenome.table.pass.summary.meth.phased.txt    # summary file, extended with haplotype-aware methylation information, as calculated by modkit (see below)
                ├── sampleID_RefGenome.table.pass.summary.meth.unphased.txt  # summary file, extended with haplotype-unaware methylation information, as calculated by modkit (see below)
                ├── modkit/                                                  # detailed modkit output information

            ├── ref_TE/                                                      # output folder for methylation analysis of annotated TEs in reference genome 
                ├── TE_type/                                                 # class of TE to be analysed, as specified in the config
                      ├── sampleID_RefGenome_methylation_summary.tsv         # summary table containing methylation information per annotated TE (see below)
                      ├── cpg_resolution/                                    # directory containing TE-specific modkit pileups
                      ├── variants/                                          # directory containing sequence variants (based on phased VCFs) that intersect with the annotated TEs
 ```
  
- Non-reference TE insertion methylation summary table:  

| Column          | Description                                                         |
| ------------------- | ------------------------------------------------------------------- |
| `UUID`             | Unique identifier of TE insertion                              |
| `Chrom`             | Chromosome                        |
| `Start`               | Start position of TE insertion                |
| `End` | End position of TE insertion |
| `Strand`            | Strand information                                            |
| `Family:Subfamily`             | Classification of TE                        |
| `LengthIns`               | Length of insertion (in bp)                |
| `Consensus` | Consensus sequence of TE insertion, with lower case corresponding to the TE sequence and upper case to flanking regions |
| `Phasing`            | Phasing information in format sampleID_phased_alignment|phaseblock:HP:supporting_reads                                         |
| `TEAverageMeth`             | Average methylation across TE region, as calculated by modkit stats                              |
| `TE_Nvalid`             | Number of valid CpGs used for average TE methylation calculation, as reported by modkit stats                        |
| `Upstream_AverageMeth`               | Average methylation across upstream TE region, as calculated by modkit stats  (flanking region length specified in config) |
| `UpstreamN_valid` | Number of valid CpGs used for average upstream TE methylation calculation, as reported by modkit stats |
| `Downstream_AverageMeth`               | Average methylation across downstream TE region, as calculated by modkit stat (flanking region length specified in config) |
| `DownstreamN_valid` | Number of valid CpGs used for average downstream TE methylation calculation, as reported by modkit stats |


- reference TEs summary table description of columns:

| Column                    | Description                                                                                       |
|---------------------------|-------------------------------------------------------------------------------------------------|
| `#chr`                    | Chromosome                                                                                       |
| `start`                   | Start position of annotated TE                                                                   |
| `end`                     | End position of annotated TE                                                                     |
| `family`                  | RepeatMasker family classification of TE                                                        |
| `strand`                  | Strand information                                                                               |
| `TE_avgMeth_phased`       | Haplotype-aware average methylation across TE region, as calculated by modkit stats             |
| `TE_Nvalid_phased`        | Number of valid CpGs used for haplotype-aware average TE methylation calculation, reported by modkit stats |
| `upstream_avgMeth_phased` | Haplotype-aware average methylation across upstream TE region (flanking region length in config)|
| `upstream_Nvalid_phased`  | Number of valid CpGs used for haplotype-aware average upstream TE methylation calculation       |
| `downstream_avgMeth_phased`| Haplotype-aware average methylation across downstream TE region (flanking region length in config) |
| `downstream_Nvalid_phased`| Number of valid CpGs used for haplotype-aware average downstream TE methylation calculation     |
| `TE_avgMeth_unphased`     | Haplotype-unaware average methylation across TE region                                          |
| `TE_Nvalid_unphased`      | Number of valid CpGs used for haplotype-unaware average TE methylation calculation               |
| `upstream_avgMeth_unphased`| Haplotype-unaware average methylation across upstream TE region (flanking region length in config) |
| `upstream_Nvalid_unphased`| Number of valid CpGs used for haplotype-unaware average upstream TE methylation calculation     |
| `downstream_avgMeth_unphased`| Haplotype-unaware average methylation across downstream TE region (flanking region length in config) |
| `downstream_Nvalid_unphased`| Number of valid CpGs used for haplotype-unaware average downstream TE methylation calculation   |
| `total_SNP`               | Total number of called SNVs intersecting with annotated TE regions (SNP_count + INDEL_count)     |
| `SNP_count`               | Total number of intersecting SNPs                                                               |
| `INDEL_count`             | Total number of intersecting INDELs                                                             |
| `SV_count`                | Total number of SVs intersecting with annotated TE regions                                      |
| `SV_types`                | Type(s) of intersecting SV                                                                      |
| `SV_IDs`                  | ID(s) of intersecting SV                                                                        |


**QC evaluation**

```
projectID/
    ├── qc/
        ├── basecalling/
              ├── sampleID/
                    ├── sampleID_summary.tsv      # summary file generated by dorado
        ├── fastq/
              ├── pre-filtering/    
                    ├── sampleID/                 # directory containing nanoplot output for raw fastq files
              ├── post-filtering/
                    ├── sampleID/                 # directory containing nanoplot output for fastq files post-filtering
        ├── phased_bam/
                    ├── sampleID/                 # directory containing nanoplot output for fastq files post-filtering
                        ├── RefGenome/            # directory containing cramino and nanoplot output for phased bam file
```

You can find more details on the QC metrics generated by the specific tools on their github pages: [cramino](https://github.com/wdecoster/cramino) and  [nanoplot](https://github.com/wdecoster/NanoPlot).


**logs**

```
projectID/
    ├── logs/
          ├── slurm/                            # directory containing logs associated to SLURM jobs
          ├── snakemake-rules/                  # directory containing logs associated with each rule in the snakemake workflow

```
---

## QUESTIONS?
Please leave any feedback, issue or question on the [Issues section](https://github.com/leenput/repeatome_pipeline/issues).   
