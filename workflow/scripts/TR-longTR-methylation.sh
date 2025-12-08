#!/usr/bin/env bash
#
# -----------------------------------------------------------------------------
# Script Name:    TR-longTR-methylation.sh
# Description:    This script calculate DNA methylation (DNAm) levels at Tandem Repeats (TR) sites.
#                 The script takes as input a vcf file obtained from LongTR tool. Then it 
#                 calculates haplotype-specific DNAm level at cpg sites and averages 
#                 for TR and their upstream and downstream region. Finally, it recreates 
#                 a vcf file.
# Author:         Leena Putzeys, Brando Poggiali, Nicole Flack, Henry Barton
# Date Created:   2025-03-02
# Last Modified:  2025-11-24
# Version:        5.1.0
# License:        MIT
# Dependencies:   [modkit, bgzip, tabix, samtools, minimap2, bedtools, bcftools, awk]
# Usage: ./TR-longTR-methylation_v4.sh -v <vcf_file> -r <reference_fasta> -i <phased_bam> 
#        -o <output_dir> -s <sample_id> -e <extend> -f <flanking_bases> 
#        -h <haploid_chromosomes>
# -----------------------------------------------------------------------------


set -euo pipefail
# source functions
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/utils_TR-methylation.sh"

# Print usage instructions
usage() {
    echo "Usage: $0 -v <vcf_file> -r <reference_fasta> -i <phased_bam> -o <output_dir> -s <sample_id> -e <extend_consensus> -f <flanking_bases> -h <haploid_chromosomes>"
    echo "  -v <vcf_file>         Input VCF file with variant data and corrected headers"
    echo "  -r <reference_fasta>  Reference genome in FASTA format"
    echo "  -i <phased_bam>       Input phased BAM file"
    echo "  -o <output_dir>       Temporary output directory for generated intermediate files"
    echo "  -s <sample_id>        SampleID to be used in file names"
    echo "  -t <threads>          Number of threads"
    echo "  -e <extend>           Number of flanking bases to include for remapping reads"
    echo "  -f <flanking_bases>   Number of flanking bases to include for up- and downstream methylation analysis"
    echo "  -h <haploid_chromosomes> Comma-separated list of haploid chromosomes (e.g. chrX,chrY)"

    exit 1
}

# Parse command-line arguments
while getopts "v:r:i:o:s:t:e:f:h:" opt; do
    case $opt in
        v) INPUT_VCF="$OPTARG" ;;
        r) REFERENCE_FASTA="$OPTARG" ;;
        i) PHASED_BAM="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        s) SAMPLE_ID="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        e) EXTEND="$OPTARG" ;;
        f) FLANKING_BASES="$OPTARG" ;;
        h) HAPLOID_CHROMOSOMES="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check required arguments
if [[ -z "$INPUT_VCF" || -z "$REFERENCE_FASTA" || -z "$PHASED_BAM" || -z "$OUTPUT_DIR"  || -z "$SAMPLE_ID" ]]; then
    usage
fi

# Verify tools are installed
MISSING_TOOLS=()
for tool in samtools minimap2 awk modkit bgzip tabix bedtools bcftools uTR; do
    if ! command -v "$tool" &> /dev/null; then
        MISSING_TOOLS+=("$tool")
    fi
done

# Define uTR tool path
# UTR_TOOL_DIR="/ifs/software/research/unique/pipeline_tools/uTR/uTR"
UTR_EXE="${UTR_EXE:-uTR}" 

if [ ${#MISSING_TOOLS[@]} -ne 0 ]; then
    echo "Error: The following required tools are not installed: ${MISSING_TOOLS[*]}"
    exit 1
fi

#Check if files are present
[[ ! -f "$INPUT_VCF" ]] && echo "VCF file not found!" && exit 1
[[ ! -f "$REFERENCE_FASTA" ]] && echo "Reference FASTA not found!" && exit 1
[[ ! -f "$PHASED_BAM" ]] && echo "Phased BAM not found!" && exit 1

#Function to process each line of a vcf file
process_line(){
    local LINE="$1"
    MODIFIED_LINE=$(modify_gt_based_on_pdp "$LINE")
    
    # Parse VCF fields
    CHROM=$(echo "$MODIFIED_LINE" | cut -f1)
    POS=$(echo "$MODIFIED_LINE" | cut -f2)
    TR_ID=$(echo "$MODIFIED_LINE" | cut -f3)
    REF=$(echo "$MODIFIED_LINE" | cut -f4)
    ALT=$(echo "$MODIFIED_LINE" | cut -f5)
    TR_START=$(echo "$MODIFIED_LINE" | grep -o "START=[0-9]*" | cut -d'=' -f2)
    TR_END=$(echo "$MODIFIED_LINE" | grep -o "END=[0-9]*" | cut -d"=" -f2)
    GENOTYPE=$(echo "$MODIFIED_LINE" | cut -f10 | cut -d':' -f1)
    
    # Parse alternate genotypes and alleles
    IFS='|' read -ra GENO_ARRAY <<< "$GENOTYPE"
    IFS=',' read -ra ALT_ALLELES <<< "$ALT"
    ALL_ALLELES=("$REF" "${ALT_ALLELES[@]}")
    
    echo "$LINE" >> "$MAIN_LOG"
    echo "" >> "$MAIN_LOG"
    echo "Genotype of the locus: ${GENO_ARRAY[@]}" >> "$MAIN_LOG"
    echo "List of all alleles: ${ALL_ALLELES[@]}"  >> "$MAIN_LOG"
    
    #Declare variables with empty values
    HP1_AVG_TR_VALUE="."
    HP1_TR_NCOV="."
    HP1_AVG_upTR_VALUE="."
    HP1_upTR_NCOV="."
    HP1_AVG_downTR_VALUE="."
    HP1_downTR_NCOV="."
    HP1_TR_LEN="."
    HP1_TR_N_CPG="."
    HP2_AVG_TR_VALUE="."
    HP2_TR_NCOV="."
    HP2_AVG_upTR_VALUE="."
    HP2_upTR_NCOV="."
    HP2_AVG_downTR_VALUE="."
    HP2_downTR_NCOV="."
    
    TR_LEN="."
    TR_N_CPG="."
    HP1_TR_LEN="."
    HP1_TR_N_CPG="."
    HP2_TR_LEN="."
    HP2_TR_N_CPG="."
    
    TR_CPG_METH_HP1="."
    TR_CPG_DEPTH_HP1="."
    TR_CPG_METH_HP2="."
    TR_CPG_DEPTH_HP2="." 
    
    HP1_TR_AVG_VALUE="."
    HP1_TR_NCOV="."
    HP1_upTR_AVG_VALUE="."
    HP1_upTR_NCOV="."
    HP1_downTR_AVG_VALUE="."
    HP1_downTR_NCOV="."
    HP2_TR_N_CPG="."
    HP2_TR_AVG_VALUE="."
    HP2_TR_NCOV="."
    HP2_upTR_AVG_VALUE="."
    HP2_upTR_NCOV="."
    HP2_downTR_AVG_VALUE="."
    HP2_downTR_NCOV="."
    
    TR_PATTERN="."
    HP1_TR_PATTERN="."
    HP2_TR_PATTERN="."
    
    # Determine repeat extension boundaries
    UPSTREAM_START=$((TR_START - EXTEND))
    UPSTREAM_END=$((TR_START - 1))
    DOWNSTREAM_START=$((TR_END + 1))
    DOWNSTREAM_END=$((TR_END + EXTEND))
    
    # Extract flanking sequences, they are the same for every allele
    UPSTREAM_SEQ=$(samtools faidx "$REFERENCE_FASTA" "${CHROM}:${UPSTREAM_START}-${UPSTREAM_END}" | tail -n +2 | tr -d '\n' | tr '[:upper:]' '[:lower:]')
    DOWNSTREAM_SEQ=$(samtools faidx "$REFERENCE_FASTA" "${CHROM}:${DOWNSTREAM_START}-${DOWNSTREAM_END}" | tail -n +2 | tr -d '\n' | tr '[:upper:]' '[:lower:]')
    
    
    for i in "${!GENO_ARRAY[@]}"; do #iterate over indices to work with haplotypes [0 is HP1 and 1 is HP2]
        allele=${GENO_ARRAY[${i}]}
        
        if [[ "$allele" == "." ]]; then
            if [[ "$i" == 0 ]]; then
                echo "Haplotype 1 is not processed because allele is ." >> "$MAIN_LOG"
            elif [[ "$i" == 1 ]]; then
                echo "Haplotype 2 is not processed because allele is ." >> "$MAIN_LOG"
            fi
            continue
        fi

        allele_seq=${ALL_ALLELES[${allele}]}
        
        if [[ "$allele_seq" == "<DEL>" ]]; then
            if [[ "$i" == 0 ]]; then
                echo "Haplotype 1 is not processed because allele is <DEL>" >> "$MAIN_LOG"
            elif [[ "$i" == 1 ]]; then
                echo "Haplotype 2 is not processed because allele is <DEL>" >> "$MAIN_LOG"
            fi
            continue
        fi 
        
        #Count length of allele, number of CpG sites (CG) in the allele sequence (a final C is not included in the count)
        if is_haploid "$CHROM"; then
            HAPLOTYPE="haploid"
            TR_LEN=${#allele_seq}
            TR_N_CPG=$(echo "${allele_seq}" | grep -o "CG" | wc -l)
            echo "Process allele: "${HAPLOTYPE}""  >> "$MAIN_LOG"
            echo -e "Allele ${allele}: ${allele_seq} \n  with length: ${TR_LEN}" >> "$MAIN_LOG"
        elif [[ "$i" -eq 0 ]]; then
            HAPLOTYPE=1
            HP1_TR_LEN=${#allele_seq}
            HP1_TR_N_CPG=$(echo "${allele_seq}" | grep -o "CG" | wc -l)
            echo "Process allele: "${HAPLOTYPE}"" >> "$MAIN_LOG"
            echo -e "Allele ${allele}: ${allele_seq} \n  with length: ${HP1_TR_LEN}" >> "$MAIN_LOG"
        elif [[ "$i" -eq 1 ]]; then
            HAPLOTYPE=2
            HP2_TR_LEN=${#allele_seq}
            HP2_TR_N_CPG=$(echo "${allele_seq}" | grep -o "CG" | wc -l)
            echo "Process allele: "${HAPLOTYPE}"" >> "$MAIN_LOG"
            echo -e "Allele ${allele}: ${allele_seq} \n  with length: ${HP2_TR_LEN}" >> "$MAIN_LOG"
        fi
        
        #Create fasta file
        SEQ="${UPSTREAM_SEQ}$(echo "${allele_seq}" | tr '[:lower:]' '[:upper:]')${DOWNSTREAM_SEQ}"
        HEADER="${CHROM}_${POS}_${TR_ID}_${HAPLOTYPE}"
        REGION_FASTA="${OUTPUT_DIR}/${HEADER}_region.fasta"
        STR_ALLELE_FASTA="${OUTPUT_DIR}/${HEADER}_STR_allele.fasta"
        uTR_out="${OUTPUT_DIR}/${HEADER}_uTR.out"
        TR_START_REL=$((EXTEND + 1))
        TR_END_REL=$((EXTEND + ${#REF}))
        echo "**** Start Analysis for: $HEADER ****" >> "$MAIN_LOG"
        echo ">$HEADER" > "$REGION_FASTA"
        echo ">$HEADER" > "$STR_ALLELE_FASTA"
        echo "$SEQ" >> "$REGION_FASTA"
        echo "$allele_seq" | tr -d '\r' | tr -cd 'ACGTNacgtn' >> "$STR_ALLELE_FASTA"
        #echo "$allele_seq" >> "$STR_ALLELE_FASTA"
        samtools faidx "$REGION_FASTA"
         
        # Extract reads overlapping the region and convert to FASTQ
        REGION_BAM="${OUTPUT_DIR}/${TR_ID}_extracted.bam"
        samtools view -h "$PHASED_BAM" "$CHROM:$TR_START-$TR_END" | samtools addreplacerg -r "ID:${TR_ID}" - | samtools view -b - > "$REGION_BAM"
        
        # Filter reads by HP tag and run uTR to get pattern of TR allele
        if [[ "$HAPLOTYPE" == 1 ]]; then
            echo "Filter HP1 reads" >> "$MAIN_LOG"
            FILTERED_BAM="${OUTPUT_DIR}/${CHROM}_${TR_ID}_extracted_HP1.bam"
            samtools view -h -d HP:1 -b -o "$FILTERED_BAM" "$REGION_BAM"        
        elif [[ "$HAPLOTYPE" == 2 ]]; then
            echo "Filter HP2 reads" >> "$MAIN_LOG"
            FILTERED_BAM="${OUTPUT_DIR}/${CHROM}_${TR_ID}_extracted_HP2.bam"
            samtools view -h -d HP:2 -b -o "$FILTERED_BAM" "$REGION_BAM"        
        else
            echo "No filtering since chromosome is haploid" >> "$MAIN_LOG"
            FILTERED_BAM="$REGION_BAM"
        fi
        
        # Validate the bam file
        if [[ ! -s "$FILTERED_BAM" ]]; then
            echo "Error: no reads extracted for $TR_ID, $HAPLOTYPE" >> "$MAIN_LOG"
            continue
        fi
        
        # Convert to FASTQ
         echo "Converting to fast $HEADER" >> "$MAIN_LOG"
        REGION_FASTQ="${OUTPUT_DIR}/${CHROM}_${TR_ID}_${HAPLOTYPE}.extracted.fastq"
        samtools fastq -t -T MM,ML,HP,PS "$FILTERED_BAM" > "$REGION_FASTQ" 2>/dev/null
        
        # Map with minimap2
        echo "Remapping to TR $HEADER sequence" >> "$MAIN_LOG"
        REGION_SAM="${OUTPUT_DIR}/${CHROM}_${TR_ID}_${HAPLOTYPE}_mapped.sam"
        minimap2 -ax map-ont -y --sam-hit-only "$REGION_FASTA" "$REGION_FASTQ" > "$REGION_SAM" 2>> "$MAIN_LOG"
        
        # Convert, sort, and index BAM
        REGION_REALIGNED_BAM="${ALIGNMENTS}/${CHROM}_${TR_ID}_${HAPLOTYPE}_mapped.sorted.bam"
        samtools view -bS "$REGION_SAM" | samtools sort -o "$REGION_REALIGNED_BAM" 2>/dev/null
        samtools index "$REGION_REALIGNED_BAM"
        
        #If bam file does not have any read stop the analysis of the allele
        if [[ $(samtools view -c "$REGION_REALIGNED_BAM") -eq 0 ]]; then
            echo "Skipping $BAM_FILE: no reads"
            continue
        fi
        
        # Perform modkit pileup on the generated BAM file
        echo "Running modkit pileup for ${HEADER}" >> "$MAIN_LOG"
        PILEUP_OUTPUT="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_region_modkit_pileup.bed"
        modkit pileup "${REGION_REALIGNED_BAM}" "${PILEUP_OUTPUT}" --cpg --ref "${REGION_FASTA}" --ignore h --combine-strands --mod-threshold m:0.8 2>> "$MAIN_LOG"
        # Brief pause to ensure modkit output is fully written before bgzip
        timeout 10 bash -c "while [ ! -s '$PILEUP_OUTPUT' ]; do sleep 0.1; done"
        echo "Modkit pileup completed successfully" >> "$MAIN_LOG"
        # Compress and index the phased pileup output
        bgzip "$PILEUP_OUTPUT"
        tabix "${PILEUP_OUTPUT}.gz"
        
        # Perform modkit stats for the TR region for each phased output
        TR_ALL_START=$((EXTEND + 1))
        TR_ALL_END=$((EXTEND + ${#ALL_ALLELES[${allele}]}))
        TR_REGION_END=$((EXTEND + ${#ALL_ALLELES[${allele}]} + ${EXTEND}))
         
        REGION_BED="${OUTPUT_DIR}/${HEADER}.bed"
        REGION_UPSTREAM_BED="${OUTPUT_DIR}/${HEADER}_upstream.bed"
        REGION_DOWNSTREAM_BED="${OUTPUT_DIR}/${HEADER}_downstream.bed"
        GENOME_COORDINATES="${OUTPUT_DIR}/${HEADER}_coordinate.txt"
        echo -e "$HEADER\t$TR_ALL_START\t$TR_ALL_END\t$TR_ID" > "$REGION_BED"
        
        #Create genome fasta file to use bedtools flank
        echo -e "${HEADER}\t${TR_REGION_END}" > "${GENOME_COORDINATES}" 
        
        bedtools flank -i "$REGION_BED" -l "$FLANKING_BASES" -r 0 -g "${GENOME_COORDINATES}" > "$REGION_UPSTREAM_BED"
        bedtools flank -i "$REGION_BED" -l 0 -r "$FLANKING_BASES" -g "${GENOME_COORDINATES}" > "$REGION_DOWNSTREAM_BED"
        
        STAT_OUTPUT_TR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_modkit_stats.tsv"
        STAT_OUTPUT_upTR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_upstream_modkit_stats.tsv"
        STAT_OUTPUT_downTR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_downstream_modkit_stats.tsv"
        
        echo "Running modkit stats for TR ${HEADER}" >> "$MAIN_LOG"
        modkit stats --regions "$REGION_BED" --min-coverage 3 -o "${STAT_OUTPUT_TR}" "${PILEUP_OUTPUT}.gz" 2>/dev/null
        modkit stats --regions "$REGION_UPSTREAM_BED" --min-coverage 3 -o "${STAT_OUTPUT_upTR}" "${PILEUP_OUTPUT}.gz" 2>/dev/null
        modkit stats --regions "$REGION_DOWNSTREAM_BED" --min-coverage 3 -o "${STAT_OUTPUT_downTR}" "${PILEUP_OUTPUT}.gz" 2>/dev/null
        
        #Extraction of methylation values and run uTR to get allele pattern
        echo "Extraction of average methylation values" >> "$MAIN_LOG"
        
        TR_AVG_METHYLATION=$(awk 'NR==2 {if ($8 == "") print "."; else print $8}' "$STAT_OUTPUT_TR")
        TR_COV_METH=$(awk 'NR==2 {if ($7 == "") print "."; else print $7}' "$STAT_OUTPUT_TR")
        upTR_AVG_METHYLATION=$(awk 'NR==2 {if ($8 == "") print "."; else print $8}' "$STAT_OUTPUT_upTR")
        upTR_COV_METH=$(awk 'NR==2 {if ($7 == "") print "."; else print $7}' "$STAT_OUTPUT_upTR")
        downTR_AVG_METHYLATION=$(awk 'NR==2 {if ($8 == "") print "."; else print $8}' "$STAT_OUTPUT_downTR")
        downTR_COV_METH=$(awk 'NR==2 {if ($7 == "") print "."; else print $7}' "$STAT_OUTPUT_downTR")
        
        echo "Extraction of CpG methylation values" >> "$MAIN_LOG"
        PILEUP_TR_OUTPUT="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_TR_modkit_pileup.bed"
        
        
        if [[ -s "${PILEUP_OUTPUT}.gz" ]] &&  bedtools intersect -a "${PILEUP_OUTPUT}.gz" -b "$REGION_BED" > "$PILEUP_TR_OUTPUT" 2>/dev/null; then
            if [[ -s "$PILEUP_TR_OUTPUT" ]]; then
                echo "Intersect succeeded, continuing with next step..." >> "$MAIN_LOG"
                TR_CPG_METH=$(cut -f11 "$PILEUP_TR_OUTPUT" | paste -sd, -)
                TR_CPG_DEPTH=$(cut -f10 "$PILEUP_TR_OUTPUT" | paste -sd, -)
            else
                TR_CPG_METH="."
                TR_CPG_DEPTH="."
                echo "No CpG sites in the allele" >> "$MAIN_LOG"
            fi
        else
            TR_CPG_METH="."
            TR_CPG_DEPTH="."
            echo "No CpG sites in the region" >> "$MAIN_LOG"
        fi
        
        if "$UTR_EXE" -f "$STR_ALLELE_FASTA" -y -o "$uTR_out" 2>/dev/null; then
            TR_PATTERN=$(extract_uTR_features "$uTR_out")
            echo -e " Pattern: ${TR_PATTERN}" >> "$MAIN_LOG" 
        else
            echo "uTR failed to decompose the DNA sequence. Continuing with the rest of the pipeline..."
            echo "uTR failed to decompose the DNA sequence" >> "$MAIN_LOG"
        fi
        
        
        if [[ "$HAPLOTYPE" == "haploid" ]]; then
            TR_AVG_VALUE="$TR_AVG_METHYLATION"
            TR_NCOV="$TR_COV_METH"
            upTR_AVG_VALUE="$upTR_AVG_METHYLATION"
            upTR_NCOV="$upTR_COV_METH"
            downTR_AVG_VALUE="$downTR_AVG_METHYLATION"
            downTR_NCOV="$downTR_COV_METH"
            TR_CPG_METH_HP1="$TR_CPG_METH"
            TR_CPG_DEPTH_HP1="$TR_CPG_DEPTH"
            TR_PATTERN="$TR_PATTERN"
        elif [[ "$HAPLOTYPE" == 1 ]]; then
            HP1_TR_AVG_VALUE="$TR_AVG_METHYLATION"
            HP1_TR_NCOV="$TR_COV_METH"
            HP1_upTR_AVG_VALUE="$upTR_AVG_METHYLATION"
            HP1_upTR_NCOV="$upTR_COV_METH"
            HP1_downTR_AVG_VALUE="$downTR_AVG_METHYLATION"
            HP1_downTR_NCOV="$downTR_COV_METH"
            TR_CPG_METH_HP1="$TR_CPG_METH"
            TR_CPG_DEPTH_HP1="$TR_CPG_DEPTH"
            HP1_TR_PATTERN="$TR_PATTERN"
        elif [[ "$HAPLOTYPE" == 2 ]]; then
            HP2_TR_AVG_VALUE="$TR_AVG_METHYLATION"
            HP2_TR_NCOV="$TR_COV_METH"
            HP2_upTR_AVG_VALUE="$upTR_AVG_METHYLATION"
            HP2_upTR_NCOV="$upTR_COV_METH"
            HP2_downTR_AVG_VALUE="$downTR_AVG_METHYLATION"
            HP2_downTR_NCOV="$downTR_COV_METH"
            TR_CPG_METH_HP2="$TR_CPG_METH"
            TR_CPG_DEPTH_HP2="$TR_CPG_DEPTH"
            HP2_TR_PATTERN="$TR_PATTERN"
        fi
        
        
        # Cleanup temporary files
        rm -f "$REGION_FASTA"
        rm -f "${REGION_FASTA}.fai"
        rm -f "$STR_ALLELE_FASTA"
        rm -f "$REGION_FASTQ"
        rm -f "$REGION_SAM"
        rm -f "$REGION_BAM"
        rm -f "$REGION_REALIGNED_BAM"
        rm -f "${REGION_REALIGNED_BAM}.bai"
        rm -f "$FILTERED_BAM"
        rm -f "$REGION_BED"
        rm -f "$REGION_UPSTREAM_BED"
        rm -f "$REGION_DOWNSTREAM_BED"
        rm -f "$GENOME_COORDINATES"
        
        rm -f "$STAT_OUTPUT_TR"
        rm -f "$STAT_OUTPUT_upTR"
        rm -f "$STAT_OUTPUT_downTR"
        rm -f "$uTR_out"
        
        rm -f "$PILEUP_TR_OUTPUT"
        rm -f "${PILEUP_OUTPUT}.gz" "${PILEUP_OUTPUT}.gz.tbi"
    
    done   
    
    # Create vcf line
    echo "Recreate vcf line with methylation information" >> "$MAIN_LOG"    
    echo "" >> "$MAIN_LOG"
    # Extract the information from  the modified VCF line
    BASE_VCF_FIELDS=$(echo "$MODIFIED_LINE" | cut -f1-8)
    OLD_FORMAT=$(echo "$MODIFIED_LINE" | cut -f9)
    OLD_SAMPLE=$(echo "$MODIFIED_LINE" | cut -f10)
    
    #Add information to the new VCF line 
    NEW_FORMAT="${OLD_FORMAT}:TR_LEN:TR_N_CPG:TR_PATTERN:TR_AM:TR_N_METH_VALID:UPSTREAM_TR_AM:UPSTREAM_TR_N_METH_VALID:DOWNSTREAM_TR_AM:DOWNSTREAM_TR_N_METH_VALID:TR_CPG_METH_HP1:TR_CPG_DEPTH_HP1:TR_CPG_METH_HP2:TR_CPG_DEPTH_HP2"
    if is_haploid "$CHROM"; then  
    NEW_SAMPLE="${OLD_SAMPLE}:${TR_LEN}:${TR_N_CPG}:${TR_PATTERN}:${TR_AVG_VALUE}:${TR_NCOV}:${upTR_AVG_VALUE}:${upTR_NCOV}:${downTR_AVG_VALUE}:${downTR_NCOV}:${TR_CPG_METH_HP1}:${TR_CPG_DEPTH_HP1}:${TR_CPG_METH_HP2}:${TR_CPG_DEPTH_HP2}"
    else
    NEW_SAMPLE="${OLD_SAMPLE}:${HP1_TR_LEN},${HP2_TR_LEN}:${HP1_TR_N_CPG},${HP2_TR_N_CPG}:${HP1_TR_PATTERN},${HP2_TR_PATTERN}:${HP1_TR_AVG_VALUE},${HP2_TR_AVG_VALUE}:${HP1_TR_NCOV},${HP2_TR_NCOV}:${HP1_upTR_AVG_VALUE},${HP2_upTR_AVG_VALUE}:${HP1_upTR_NCOV},${HP2_upTR_NCOV}:${HP1_downTR_AVG_VALUE},${HP2_downTR_AVG_VALUE}:${HP1_downTR_NCOV},${HP2_downTR_NCOV}:${TR_CPG_METH_HP1}:${TR_CPG_DEPTH_HP1}:${TR_CPG_METH_HP2}:${TR_CPG_DEPTH_HP2}"
    fi   

    # Rebuild the VCF line
    NEW_LINE="${BASE_VCF_FIELDS}\t${NEW_FORMAT}\t${NEW_SAMPLE}"
    # append to output vcf
    echo -e "$NEW_LINE" >> "$OUTPUT_VCF"
    #"${OUTPUT_DIR}/${CHROM}_${POS}_${TR_ID}.line.tsv"
}

# ******* START ANALYSIS *******
echo ""
echo "########## START ANALYSIS ##########"
echo ""

#Activate cleanup function when exit script
# trap cleanup EXIT

# Define paths inside the output directory
OUTPUT_VCF="${INPUT_VCF%.vcf.gz}_methylated.vcf"
ALIGNMENTS="${OUTPUT_DIR}/alignments"
METH="${OUTPUT_DIR}/methylation"
MAIN_LOG="${OUTPUT_VCF%.vcf}.log"

#Create output directories
mkdir -p "$OUTPUT_DIR" "$ALIGNMENTS" "$METH" 

# Copy existing headers from the input VCF
zgrep '^##' "$INPUT_VCF" | grep -v '^##bcftools'> "$OUTPUT_VCF"

# Add new FORMAT fields for allele and haplotype-specific methylation info in vcf file
cat <<EOF >> "$OUTPUT_VCF"
##FORMAT=<ID=TR_LEN,Number=2,Type=Integer,Description="Lengths of the tandem repeat alleles in bp, one per haplotype">
##FORMAT=<ID=TR_N_CPG,Number=2,Type=Integer,Description="Number of CpG sites in the TR sequence, one per haplotype">
##FORMAT=<ID=TR_PATTERN,Number=2,Type=String,Description="Decomposed TR allele pattern using uTR tool, one per haplotype">
##FORMAT=<ID=TR_AM,Number=2,Type=Float,Description="Average methylation percentage at the TR region, one per haplotype">
##FORMAT=<ID=TR_N_METH_VALID,Number=2,Type=Integer,Description="Number of valid CpG sites used for TR methylation, one per haplotype">
##FORMAT=<ID=UPSTREAM_TR_AM,Number=2,Type=Float,Description="Average methylation percentage upstream of TR, one per haplotype">
##FORMAT=<ID=UPSTREAM_TR_N_METH_VALID,Number=2,Type=Integer,Description="Number of valid CpGs upstream of TR, one per haplotype">
##FORMAT=<ID=DOWNSTREAM_TR_AM,Number=2,Type=Float,Description="Average methylation percentage downstream of TR, one per haplotype">
##FORMAT=<ID=DOWNSTREAM_TR_N_METH_VALID,Number=2,Type=Integer,Description="Number of valid CpGs downstream of TR, one per haplotype">
##FORMAT=<ID=TR_CPG_METH_HP1,Number=.,Type=Float,Description="CpG methylation percentages for haplotype 1, ordered by CpG position within the allele">
##FORMAT=<ID=TR_CPG_DEPTH_HP1,Number=.,Type=Integer,Description="CpG coverage depths for haplotype 1, ordered by CpG position within the allele">
##FORMAT=<ID=TR_CPG_METH_HP2,Number=.,Type=Float,Description="CpG methylation percentages for haplotype 2, ordered by CpG position within the allele">
##FORMAT=<ID=TR_CPG_DEPTH_HP2,Number=.,Type=Integer,Description="CpG coverage depths for haplotype 2, ordered by CpG position within the allele">
EOF

# Append column names and the body to the header
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$SAMPLE_ID" >> "$OUTPUT_VCF"

# Export functions and variables in order to be seen by the process_line function run as subshell
export -f process_line
export -f modify_gt_based_on_pdp
export -f is_haploid
export -f extract_uTR_features

export REFERENCE_FASTA
export PHASED_BAM
export OUTPUT_DIR
export OUTPUT_VCF
export SAMPLE_ID
export EXTEND
export FLANKING_BASES
export HAPLOID_CHROMOSOMES
export ALIGNMENTS
export METH
export UTR_EXE
export MAIN_LOG

# Run process_line on headless input VCF
echo "Processing lines in ${INPUT_VCF}"

# set -x  # trace commands
# exec 3>&1 4>&2   # save stdout/stderr

# xargs could be dropped entirely here (could avoid writting so many tmp files) and speed could be adjusted in snakemake with chunking
xargs -P "$THREADS" -n 1 -d '\n' bash -c 'set -e; process_line "$1"' _ < <(bcftools view -H "$INPUT_VCF")
wait
# Append results and clean up
rm -r "${OUTPUT_DIR}"

echo "Processing complete" 
