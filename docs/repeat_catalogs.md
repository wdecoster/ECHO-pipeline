# BUNDLED REPEAT CATALOGS - echoDB

The ECHO pipeline includes multiple repeat catalogs specifically designed to capture the repetitive elements of interest, both tandem repeats (TRs) and transposable elements (TEs), 
for the selected reference genome (GRCh38 or T2T CHM13v2). These catalogs, which are hosted within this [sandbox zenodo repository](https://sandbox.zenodo.org/records/430049), are compiled from published resources and adapted to ensure compatibility with the pipeline. 
They define the genomic loci where genotyping and/or methylation profiling is performed, enabling analysis of the human repeatome. 

**Version:** echoDB v1  
**Release date:** January 2026

---

## Supported reference assemblies

echoDB currently supports the following assemblies:

- **GRCh38**
- **T2T-CHM13v2**


Default catalogs used by the ECHO pipeline:

| Assembly        | TR genotyping            | TR methylation        | TE characterization |
|-----------------|--------------------------|------------------------|----------------|
| GRCh38          | adotto_genome_wide       | adotto_cpg_str        | te_all         |
| T2T-CHM13v2     | strchive_disease_loci    | strchive_disease_loci | te_all         |


---

## Tandem repeat (TR) catalogs

For TR analysis, the catalogs specify the 1-based genomic coordinates of repeat loci to be profiled in the following format:
```bash
chr  start  end  TRmotif  TR_ID
```

### Genome-wide (GRCh38)

#### Adotto genome-wide TRs
- **ID:** `adotto_genome_wide`
- **Loci:** 1,784,804
- **Description:**  
  Genome-wide TR catalog derived from the Adotto repeat catalog, consolidating multiple repeat discovery sources.  
  Provides broad coverage of STRs and VNTRs, though locus boundaries may be heterogeneous due to merged annotations.  
  **Recommended default** for sensitive, genome-wide TR genotyping in ECHO. For TR regions with multiple repeat annotations,  
  the longest-spanning TRF annotation was chosen as the canonical motif (TRmotif) and used for any downstream motif-based filtering. 
- **Source:** [project adotto catalog v1.2.1](https://zenodo.org/records/13987414), released as part of the GIAB tandem repeat benchmark variant set ([English et al. 2025](https://www.nature.com/articles/s41587-024-02225-z)).

#### HipSTR reference TRs
- **ID:** `hipstr_reference_genome_wide`
- **Loci:** 1,638,945
- **Description:**  
  Conservative STR catalog based on the HipSTR reference set, originally designed for short-read genotyping.  
  Restricted to canonical STRs with simple structures and well-defined boundaries.
- **Source:** [HipSTR GitHub] (https://github.com/HipSTR-Tool/HipSTR-references/blob/master/human/hg38.hipstr_reference.bed.gz)

---

### CpG-containing TR subsets (GRCh38)

These catalogs contain loci with **CpGs in the canonical repeat motif (TRmotif)**, recommended for methylation profiling.
CpGs may occur within the motif body and/or at motif junctions.

#### CpG STRs (1-6 bp motifs)
- **adotto_cpg_str**: 8,424 loci  
- **hipstr_reference_cpg_str**: 9,318 loci  

#### CpG VNTRs (7–100 bp motifs)
- **adotto_cpg_vntr**
- Derived from the Adotto genome-wide catalog.

---

### Targeted TR panels

#### Pathogenic loci
- **ID:** `strchive_disease_loci`
- **Loci:** 73
- **Description:**  
  Disease-associated STR loci curated by [STRchive](https://strchive.org/loci/).
- **Available for:** GRCh38 and T2T-CHM13v2

#### Forensic loci (GRCh38)
- **ID:** `strbase_loci`
- **Loci:** 77
- **Description:**  
  STRbase loci formatted for LongTR-compatible input.
- **Source:** [NIST STRbase](https://strbase.nist.gov/)

---

---

## Transposable element (TE) catalogs

### Genome-wide TE annotations

echoDB provides genome-wide and class-specific transposable element (TE) annotations
for each supported reference assembly. All TE catalogs are derived from UCSC
RepeatMasker tracks and filtered to retain transposable elements only.

- **ID:** `te_all`
- **Classes included:** LINE, SINE, LTR, DNA, retroposon, helitron

Available for:
- GRCh38 (source [here](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz))
- T2T-CHM13v2 (source [here](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz))

---

### TE class-specific catalogs

Class-specific BED files are provided using the following pattern:
echoDBv1/TEs/refGenome/TE_classes/refGenome_TEs_<class>.bed

- **Classes:** `LINE`, `SINE`, `LTR`, `DNA`, `retroposon`, `helitron`
- **Description:**  
Class-specific BED files filtered by RepeatMasker class annotation.

---
## notes:
  - All BED files use chr-prefixed contig names.
  - TR catalogs are formatted for LongTR input (as provided in this DB; coordinate conventions follow the catalog source/format).
  - CpG STR catalogs are recommended for downstream methylation profiling to reduce computational burden
  - Catalogs are provided as annotation resources and do not include methylation or activity estimates


---

# CUSTOM REPEAT CATALOGS 

In addition to the bundled **echoDB** catalogs, users may supply **custom repeat catalogs** to analyze specific loci or repeat sets of interest. This allows targeted analyses such as:

- candidate loci identified in previous studies  
- custom panels for targeted applications  
- reduced catalogs for faster exploratory analyses  

Custom catalogs can be provided independently for **TR (epi)genotyping** or **TE characterization**.  

## Custom tandem repeat catalogs

Custom TR catalogs must follow the same **tab-separated 1-based BED-like format** used by ECHO:

```bash
chr   start   end   TRmotif   TR_ID
```
**Requirements**:  
- Coordinates must correspond to the reference genome used for read alignment  
- Chromosome naming must match the reference assembly (e.g. `chr1`, `chr2`, etc.) 
- TR_ID values must be unique
- Motifs should represent the canonical repeat unit where possible  
- Overlapping loci are allowed but may increase runtime and complicate interpretation

For methylation profiling of TR loci, users may optionally restrict analysis to **TRs containing CpGs within the repeat motif**.  

This can be achieved by generating a custom subset of TR loci where:  

- the repeat motif itself contains a **CG dinucleotide**, or  
- **CpGs are formed across repeat boundaries** when motifs concatenate.  

Restricting catalogs to CpG-containing TRs can substantially **reduce computational burden** and focus methylation analyses on informative loci.  

These catalogs can be supplied to ECHO via the standard `--tr_catalog` parameter.  

## Custom transposable element catalogs  

Users may also provide **custom TE catalogs** to restrict analysis to specific transposable element insertions or TE families of interest.

This can be useful for:  
- focusing on **specific TE families or subfamilies** (e.g. *L1HS*, *AluYb8*, *SVA_F*)
- analyzing **candidate polymorphic insertions**
- restricting analysis to **experimentally validated loci**
- building **reduced catalogs for faster exploratory analyses**

Custom TE catalogs should be provided in **BED format**, specifying the genomic coordinates of TE loci to analyze.

```bash
chr   start   end   TEfamily   TE_ID
```

Example catalog entries:
```
chr1 3456789 3457234 LINE/L1 L1HS_1
chr3 8765432 8765667 SINE/Alu AluYb8_23
chr7 10234567 10234999 LTR/ERV ERVK_5
```

**Requirements**:  
- Coordinates must correspond to the **reference genome used for read alignment**  
- **TE_ID values must be unique**  
- Chromosome naming must match the reference assembly (e.g. `chr1`, `chr2`, etc.)  
- Overlapping TE loci are allowed but may increase computational runtime

Custom TE catalogs can be supplied to the ECHO pipeline using the `--te_catalog` parameter.

---

