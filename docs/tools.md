### **Overview of tools**
**Raw read manipulation**  
- Basecalling with [dorado](https://github.com/nanoporetech/dorado) (v1.1.1) (*optional, GPU availability needed*)   
- Conversion from ubam to fastq format using [samtools](https://www.htslib.org/) (v1.22.1) (*optional*)     
- Read filtering using [chopper](https://github.com/wdecoster/chopper) (v0.11.0) (*optional*)   

**QC**  
- Pre- and post-filtered read QC and evaluation of mapping and phasing metrics with [cramino](https://github.com/wdecoster/cramino) (v1.1.0) and [nanoplot](https://github.com/wdecoster/NanoPlot) (v1.46.1) 

**Alignment**  
- Read alignment to human reference genome ([GRCh38](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/)/[T2T-CHM13v2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/)) with [minimap2](https://github.com/lh3/minimap2) (v2.30).      

**Variant calling**  
- Small variant calling (SNVs and Indels) using [Clair3](https://github.com/HKU-BAL/Clair3) (v1.2.0)  
- Structural variant (SV) calling using [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) (v2.6.3)   
- Filter variants with [BCFtools](https://samtools.github.io/bcftools/bcftools.html) (v1.22)   

**Phasing**    
- Phase and haplotag reads with [LongPhase](https://github.com/twolinin/longphase) (v1.7.3) (`--snp-file <SNP.vcf> --sv-file <SV.vcf> --mod-file <modcall.vcf>`)

**Methylation**  
- Generate phased/unphased methylation pileups at CpG sites with [modkit](https://github.com/nanoporetech/modkit) (v0.5.0) (`--cpg --ignore h --combine-strands`)   

**TR characterization**  
- TR genotyping using [LongTR](https://github.com/gymrek-lab/LongTR) (v1.2) (`--phased-bam`) 
- Optional: TR locus filtering for downstream methylation anaysis (in default mode with bundled catalogs: YES)   
- Using a custom [script](workflow/scripts/TR-longTR-methylation.sh), methylation information is incorporated for each TR allele in the (filtered) VCF.
- Motif decomposition using [uTR](https://github.com/morisUtokyo/uTR).

**TE characterization**  
- non-ref TE analysis: identification of TE insertions not present in the reference genome using [TLDR](https://github.com/adamewing/tldr) (v1.2.2), followed by methylation information extraction using [script](workflow/scripts/TLDR-methylation_v2.sh)  
- ref TE analysis: two scripts ([script1](workflow/scripts/ref_TE_avg_meth.sh) and [script2](workflow/scripts/ref_TE_cpg_res.sh)) to analyse and summarize sequence variants and methylation information across annotated TEs in the reference genome.
