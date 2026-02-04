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
- Optional: TR locus filtering for downstream methylation anaysis (default: YES)   
- Using a custom [script](workflow/scripts/TR-longTR-methylation_v4.sh), methylation information is incorporated for each TR allele in the VCF.
- Motif decomposition using [uTR](https://github.com/morisUtokyo/uTR).

**TE characterization**  
- non-ref TE analysis: identification of TE insertions not present in the reference genome using [TLDR](https://github.com/adamewing/tldr), followed by methylation information extraction using [script](workflow/scripts/TLDR-methylation_v2.sh)  
- ref TE analysis: two scripts ([script1](workflow/scripts/ref_TE_avg_meth.sh) and [script2](workflow/scripts/ref_TE_cpg_res.sh)) to analyse and summarize sequence variants and methylation information across annotated TEs in the reference genome.
