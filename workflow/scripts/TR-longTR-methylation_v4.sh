#!/bin/bash
# -----------------------------------------------------------------------------
# Script Name:    ref_TE_cpg_res.sh
# Description:    A script to analyse TRs genotyped using LongTR, adding
#                 allele-specific methylation information.
# Author:         Leena Putzeys, Brando Poggiali
# Date Created:   2025-03-02
# Last Modified:  2025-07-08
# Version:        2.0.0
# License:        MIT
# Dependencies:   [modkit, bgzip, tabix, samtools, minimap2, bedtools, bcftools, awk]
# Usage: ./script.sh -v <vcf_file> -r <reference_fasta> -i <phased_bam> 
#        -o <output_dir> -s <sample_id> -e <extend> -f <flanking_bases> 
#        -h <haploid_chromosomes>
# -----------------------------------------------------------------------------


set -e

# Print usage instructions
usage() {
    echo "Usage: $0 -v <vcf_file> -r <reference_fasta> -i <phased_bam> -o <output_dir> -s <sample_id> -e <extend_consensus> -f <flanking_bases> -h <haploid_chromosomes>"
    echo "  -v <vcf_file>         Input VCF file with variant data"
    echo "  -r <reference_fasta>  Reference genome in FASTA format"
    echo "  -i <phased_bam>       Input phased BAM file"
    echo "  -o <output_dir>       Output directory for generated files"
    echo "  -s <sample_id>        SampleID to be used in file names"
    echo "  -e <extend>           Number of flanking bases to include for remapping reads (default: 1000)"
    echo "  -f <flanking_bases>   Number of flanking bases to include for up- and downstream methylation analysis (default: 250)"
    echo "  -h <haploid_chromosomes> Comma-separated list of haploid chromosomes (e.g. chrX,chrY)"
    exit 1
}

# Defaults
FLANKING_BASES=200
HAPLOID_CHROMOSOMES=""
EXTEND=1000

# Parse command-line arguments
while getopts "v:r:i:o:s:e:f:h:" opt; do
    case $opt in
        v) VCF_FILE="$OPTARG" ;;
        r) REFERENCE_FASTA="$OPTARG" ;;
        i) PHASED_BAM="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        s) SAMPLE_ID="$OPTARG" ;;
        e) EXTEND="$OPTARG" ;;
        f) FLANKING_BASES="$OPTARG" ;;
        h) HAPLOID_CHROMOSOMES="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check required arguments
if [[ -z "$VCF_FILE" || -z "$REFERENCE_FASTA" || -z "$PHASED_BAM" || -z "$OUTPUT_DIR"  || -z "$SAMPLE_ID" ]]; then
    usage
fi

# Verify tools are installed
MISSING_TOOLS=()
for tool in samtools minimap2 awk modkit bgzip tabix bedtools bcftools; do
    if ! command -v "$tool" &> /dev/null; then
        MISSING_TOOLS+=("$tool")
    fi
done

# Define uTR tool path
UTR_TOOL_DIR="/ifs/software/research/unique/pipeline_tools/uTR/uTR"


if [ ${#MISSING_TOOLS[@]} -ne 0 ]; then
    echo "Error: The following required tools are not installed: ${MISSING_TOOLS[*]}"
    exit 1
fi

#Check if files are present
[[ ! -f "$VCF_FILE" ]] && echo "VCF file not found!" && exit 1
[[ ! -f "$REFERENCE_FASTA" ]] && echo "Reference FASTA not found!" && exit 1
[[ ! -f "$PHASED_BAM" ]] && echo "Phased BAM not found!" && exit 1

# Check if a chromosome is haploid
#iIFS=',' read -r -a HAPLOID_ARRAY <<< "$HAPLOID_CHROMOSOMES"


#++ Define functions ++

#Function to check if a chromosome is haploid
is_haploid() {
    local chrom="$1"
    IFS=',' read -r -a HAPLOID_ARRAY <<< "$HAPLOID_CHROMOSOMES"

    for haploid_chrom in "${HAPLOID_ARRAY[@]}"; do
        if [[ "$chrom" == "$haploid_chrom" ]]; then
            return 0
        fi
    done
    return 1
}

# Function to check if a VCF entry lacks phasing information
# and modify the GT field based on PDP values and handle haploid chromosomes
modify_gt_based_on_pdp() {
    local LINE FORMAT SAMPLE_INFO GT_INDEX PDP_INDEX GT_VALUE PDP_VALUE MODIFIED_GT CHROM

    LINE="$1"
    FORMAT=$(echo "$LINE" | cut -f9)
    SAMPLE_INFO=$(echo "$LINE" | cut -f10)
    CHROM=$(echo "$LINE" | cut -f1)  # Extract chromosome name

    # Find the indexes of GT and PDP in the FORMAT field
    GT_INDEX=$(echo "$FORMAT" | awk -F':' '{for (i=1; i<=NF; i++) if ($i=="GT") print i}')
    PDP_INDEX=$(echo "$FORMAT" | awk -F':' '{for (i=1; i<=NF; i++) if ($i=="PDP") print i}')

    # Extract GT value
    GT_VALUE=$(echo "$SAMPLE_INFO" | cut -d':' -f$GT_INDEX)

    # Check PDP value
    if is_haploid "$CHROM"; then
        MODIFIED_GT="${GT_VALUE}"
    else
        # Extract PDP value
        PDP_VALUE=$(echo "$SAMPLE_INFO" | cut -d':' -f$PDP_INDEX)

        # Split GT and PDP into alleles
        IFS='|' read -r GT1 GT2 <<< "$GT_VALUE"
        IFS='|' read -r PDP1 PDP2 <<< "$PDP_VALUE"

        # Modify GT based on PDP values
        if [[ "$PDP1" != "0" ]]; then
            GT1="$GT1"
        else
            GT1="."
        fi

        if [[ "$PDP2" != "0" ]]; then
            GT2="$GT2"
        else
            GT2="."
        fi

        # Construct modified GT
        MODIFIED_GT="$GT1|$GT2"
    fi

    # Replace the old GT value with the new one
    MODIFIED_INFO=$(echo "$SAMPLE_INFO" | awk -v gt_index="$GT_INDEX" -v new_gt="$MODIFIED_GT" 'BEGIN{FS=OFS=":"}{$gt_index=new_gt; print}')

    # Construct the modified VCF line
    echo -e "$(echo "$LINE" | cut -f1-9)\t$MODIFIED_INFO"
}

# Function to extract uTR features
extract_uTR_features() {
    local utr_output="$1"

    # Default values (if no valid data is found)
    local pat_value="."

    # Extract values using awk
    if [[ -s "$utr_output" ]]; then
        pat_value=$(awk '/#Pat / {sub(/#Decomp.*/, "", $0); print substr($0, index($0, "<"))}' "$utr_output" | sort -u | paste -sd "," - | sed 's/[[:space:]]*$//')
    fi

    # Ensure extracted values are non-empty
    echo "${pat_value}"
}


# Define paths inside the output directory
INPUT_VCF="${OUTPUT_DIR}/${SAMPLE_ID}_input_sorted.vcf"
OUTPUT_VCF="${OUTPUT_DIR}/${SAMPLE_ID}_methylated.vcf"
OUTPUT_SORTED_VCF="${OUTPUT_DIR}/${SAMPLE_ID}_methylated_sorted.vcf.gz"
TMP_DIR="${OUTPUT_DIR}/temp"
ALIGNMENTS="${OUTPUT_DIR}/alignments"
METH="${OUTPUT_DIR}/methylation"
LOGS="${OUTPUT_DIR}/logs"
OUTPUT_SUMMARY="${OUTPUT_DIR}/${SAMPLE_ID}_TR_summary.tsv"
CORRECT_HEADER_VCF_FILE="${OUTPUT_DIR}/${SAMPLE_ID}_corrected_header.vcf"

#Create output directory
mkdir -p "$OUTPUT_DIR" "$TMP_DIR" "$ALIGNMENTS" "$METH" "$LOGS"

#Checkif body of the vcf file is not empty.
if ! zgrep -v '^#' "$VCF_FILE" | grep -q .; then
    echo "Compressed VCF has no body — exiting."
    exit 1
fi



#VCF file outputed by longTR has an issue in the header formatting so it is necessary to modify the header
bcftools annotate   --header-lines <(echo '##FORMAT=<ID=DFLANKINDEL,Number=1,Type=Integer,Description="Total number of reads with an indel in the regions flanking the STR">') -Oz -o "$CORRECT_HEADER_VCF_FILE" "$VCF_FILE" 2>/dev/null

#Sort vcf file
bcftools sort "$CORRECT_HEADER_VCF_FILE" -Oz -o "$INPUT_VCF" 2>/dev/null

# Copy existing headers from the input VCF
zgrep '^##' "$CORRECT_HEADER_VCF_FILE" | grep -v '^##bcftools'> "$OUTPUT_VCF"

#bcftools view -h "$OUTPUT_VCF".tmp | grep -Ev '^##INFO=<ID=(NSKIP|NFILT|INEXACT_ALLELE|BPDIFFS|DP|DSNP|DFLANKINDEL|REFAC|AC)' > "$OUTPUT_VCF"


# Add FORMAT fields for individual-specific methylation info
# Add new FORMAT lines
cat <<EOF >> "$OUTPUT_VCF"
##FORMAT=<ID=TR_LEN,Number=G,Type=Integer,Description="Lengths of the tandem repeat alleles in bp, one per haplotype">
##FORMAT=<ID=TR_N_CPG,Number=G,Type=Integer,Description="Number of CpG sites in the TR sequence, one per haplotype">
##FORMAT=<ID=TR_PATTERN,Number=G,Type=String,Description="Decomposed TR allele pattern using uTR tool, one per haplotype">
##FORMAT=<ID=TR_AM,Number=G,Type=Float,Description="Average methylation percentage at the TR region, one per haplotype">
##FORMAT=<ID=TR_N_METH_VALID,Number=G,Type=Integer,Description="Number of valid CpG sites used for TR methylation, one per haplotype">
##FORMAT=<ID=UPSTREAM_TR_AM,Number=G,Type=Float,Description="Average methylation percentage upstream of TR, one per haplotype">
##FORMAT=<ID=UPSTREAM_TR_N_METH_VALID,Number=G,Type=Integer,Description="Number of valid CpGs upstream of TR, one per haplotype">
##FORMAT=<ID=DOWNSTREAM_TR_AM,Number=G,Type=Float,Description="Average methylation percentage downstream of TR, one per haplotype">
##FORMAT=<ID=DOWNSTREAM_TR_N_METH_VALID,Number=G,Type=Integer,Description="Number of valid CpGs downstream of TR, one per haplotype">
##FORMAT=<ID=TR_CPG_METH_HP1,Number=.,Type=Float,Description="CpG methylation percentages for haplotype 1, ordered by CpG position within the allele">
##FORMAT=<ID=TR_CPG_DEPTH_HP1,Number=.,Type=Integer,Description="CpG coverage depths for haplotype 1, ordered by CpG position within the allele">
##FORMAT=<ID=TR_CPG_METH_HP2,Number=.,Type=Float,Description="CpG methylation percentages for haplotype 2, ordered by CpG position within the allele">
##FORMAT=<ID=TR_CPG_DEPTH_HP2,Number=.,Type=Integer,Description="CpG coverage depths for haplotype 2, ordered by CpG position within the allele">
EOF


# Append column names and the body to the header
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$SAMPLE_ID" >> "$OUTPUT_VCF"

# First processing loop: read VCF file and generate TR-specific FASTA and BED entries
echo ""
echo "START ANALYSIS"
echo ""


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

#Define log output
LOG_file="${LOGS}/${CHROM}_${POS}_${TR_ID}.log"

# Parse alternate genotypes and alleles
IFS='|' read -ra GENO_ARRAY <<< "$GENOTYPE"
IFS=',' read -ra ALT_ALLELES <<< "$ALT"
ALL_ALLELES=("$REF" "${ALT_ALLELES[@]}")

echo "$LINE" > "$LOG_file"
echo "" >> "$LOG_file"
echo "Genotype of the locus: ${GENO_ARRAY[@]}" >> "$LOG_file"
echo "List of all alleles: ${ALL_ALLELES[@]}"  >> "$LOG_file"
echo "" >> "$LOG_file"

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
genotype=${GENO_ARRAY[${i}]}

if [[ "$genotype" == "." && "$i" == 0 ]]; then
    echo "Haplotype 1 is not processed because allele is ." >> "$LOG_file"
    continue
elif [[ "$genotype" == "." &&  "$i" == 1 ]]; then
    echo "Haplotype 2 is not processed because allele is ." >> "$LOG_file"
    continue
fi

allele_seq=${ALL_ALLELES[${genotype}]}

#Count length of allele, number of CpG sites (CG) in the allele sequence (a final C is not included in the count)
if is_haploid "$CHROM"; then
    HAPLOTYPE="haploid"
    TR_LEN=${#allele_seq}
    TR_N_CPG=$(echo "${allele_seq}" | grep -o "CG" | wc -l)
    echo "Process allele: "${HAPLOTYPE}""  >> "$LOG_file"
    echo -e "Allele ${genotype}: ${allele_seq} \n  with length: ${TR_LEN}" >> "$LOG_file"
elif [[ "$i" -eq 0 ]]; then
    HAPLOTYPE=1
    HP1_TR_LEN=${#allele_seq}
    HP1_TR_N_CPG=$(echo "${allele_seq}" | grep -o "CG" | wc -l)
    echo "Process allele: "${HAPLOTYPE}"" >> "$LOG_file"
    echo -e "Allele ${genotype}: ${allele_seq} \n  with length: ${HP1_TR_LEN}" >> "$LOG_file"
elif [[ "$i" -eq 1 ]]; then
    HAPLOTYPE=2
    HP2_TR_LEN=${#allele_seq}
    HP2_TR_N_CPG=$(echo "${allele_seq}" | grep -o "CG" | wc -l)
    echo "Process allele: "${HAPLOTYPE}"" >> "$LOG_file"
    echo -e "Allele ${genotype}: ${allele_seq} \n  with length: ${HP2_TR_LEN}" >> "$LOG_file"
fi

#Create fasta file
SEQ="${UPSTREAM_SEQ}$(echo "${allele_seq}" | tr '[:lower:]' '[:upper:]')${DOWNSTREAM_SEQ}"
HEADER="${CHROM}_${POS}_${TR_ID}_${HAPLOTYPE}"
REGION_FASTA="${TMP_DIR}/${HEADER}_region.fasta"
STR_ALLELE_FASTA="${TMP_DIR}/${HEADER}_STR_allele.fasta"
uTR_out="${TMP_DIR}/${HEADER}_uTR.out"
TR_START_REL=$((EXTEND + 1))
TR_END_REL=$((EXTEND + ${#REF}))
echo "**** Start Analysis for: $HEADER ****" >> "$LOG_file"
echo ">$HEADER" > "$REGION_FASTA"
echo ">$HEADER" > "$STR_ALLELE_FASTA"
echo "$SEQ" >> "$REGION_FASTA"
echo "$allele_seq" >> "$STR_ALLELE_FASTA"
samtools faidx "$REGION_FASTA"
 
# Extract reads overlapping the region and convert to FASTQ
REGION_BAM="${TMP_DIR}/${TR_ID}_extracted.bam"
samtools view -h "$PHASED_BAM" "$CHROM:$TR_START-$TR_END" | samtools addreplacerg -r "ID:${TR_ID}" - | samtools view -b - > "$REGION_BAM"

# Filter reads by HP tag and run uTR to get pattern of TR allele
if [[ "$HAPLOTYPE" == 1 ]]; then
    echo "Filter HP1 reads" >> "$LOG_file"
    FILTERED_BAM="${TMP_DIR}/${CHROM}_${TR_ID}_extracted_HP1.bam"
    samtools view -h "$REGION_BAM" | awk '/^@/ || /HP:i:1/' | samtools view -b - > "$FILTERED_BAM"
elif [[ "$HAPLOTYPE" == 2 ]]; then
    echo "Filter HP2 reads" >> "$LOG_file"
    FILTERED_BAM="${TMP_DIR}/${CHROM}_${TR_ID}_extracted_HP2.bam"
    samtools view -h "$REGION_BAM" | awk '/^@/ || /HP:i:2/' | samtools view -b - > "$FILTERED_BAM"
else
    echo "No filtering since chromosome is haploid" >> "$LOG_file"
    FILTERED_BAM="$REGION_BAM"
fi

# Validate the bam file
if [[ ! -s "$FILTERED_BAM" ]]; then
    echo "Error: no reads extracted for $TR_ID, $HAPLOTYPE" >> "$LOG_file"
    continue
fi

# Convert to FASTQ
 echo "Converting to fast $HEADER" >> "$LOG_file"
REGION_FASTQ="${TMP_DIR}/${CHROM}_${TR_ID}_${HAPLOTYPE}.extracted.fastq"
samtools fastq -t -T MM,ML,HP,PS "$FILTERED_BAM" > "$REGION_FASTQ" 2>/dev/null

# Map with minimap2
echo "Remapping to TR $HEADER sequence" >> "$LOG_file"
REGION_SAM="${TMP_DIR}/${CHROM}_${TR_ID}_${HAPLOTYPE}_mapped.sam"
minimap2 -ax map-ont -y --sam-hit-only "$REGION_FASTA" "$REGION_FASTQ" > "$REGION_SAM" 2>/dev/null

# Convert, sort, and index BAM
REGION_REALIGNED_BAM="${ALIGNMENTS}/${CHROM}_${TR_ID}_${HAPLOTYPE}_mapped.sorted.bam"
samtools view -bS "$REGION_SAM" | samtools sort -o "$REGION_REALIGNED_BAM" 2>/dev/null
samtools index "$REGION_REALIGNED_BAM"

# Perform modkit pileup on the generated BAM file
echo "Running modkit pileup for ${HEADER}" >> "$LOG_file"
PILEUP_OUTPUT="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_region_modkit_pileup.bed"
modkit pileup "${REGION_REALIGNED_BAM}" "${PILEUP_OUTPUT}" --cpg --ref "${REGION_FASTA}" --ignore h --combine-strands --mod-threshold m:0.8 2>/dev/null

echo "Modkit pileup completed successfully" >> "$LOG_file"
# Compress and index the phased pileup output
bgzip "$PILEUP_OUTPUT"
tabix "${PILEUP_OUTPUT}.gz"

# Perform modkit stats for the TR region for each phased output
TR_ALL_START=$((EXTEND + 1))
TR_ALL_END=$((EXTEND + ${#ALL_ALLELES[${genotype}]}))
TR_REGION_END=$((EXTEND + ${#ALL_ALLELES[${genotype}]} + ${EXTEND}))
 
REGION_BED="${TMP_DIR}/${HEADER}.bed"
REGION_UPSTREAM_BED="${TMP_DIR}/${HEADER}_upstream.bed"
REGION_DOWNSTREAM_BED="${TMP_DIR}/${HEADER}_downstream.bed"
GENOME_COORDINATES="${TMP_DIR}/${HEADER}_coordinate.txt"
echo -e "$HEADER\t$TR_ALL_START\t$TR_ALL_END\t$TR_ID" > "$REGION_BED"

#Create genome fasta file to use bedtools flank
echo -e "${HEADER}\t${TR_REGION_END}" > "${GENOME_COORDINATES}" 

bedtools flank -i "$REGION_BED" -l "$FLANKING_BASES" -r 0 -g "${GENOME_COORDINATES}" > "$REGION_UPSTREAM_BED"
bedtools flank -i "$REGION_BED" -l 0 -r "$FLANKING_BASES" -g "${GENOME_COORDINATES}" > "$REGION_DOWNSTREAM_BED"

STAT_OUTPUT_TR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_modkit_stats.tsv"
STAT_OUTPUT_upTR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_upstream_modkit_stats.tsv"
STAT_OUTPUT_downTR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_downstream_modkit_stats.tsv"

echo "Running modkit stats for TR ${HEADER}" >> "$LOG_file"
modkit stats --regions "$REGION_BED" --min-coverage 3 -o "${STAT_OUTPUT_TR}" "${PILEUP_OUTPUT}.gz" 2>/dev/null
modkit stats --regions "$REGION_UPSTREAM_BED" --min-coverage 3 -o "${STAT_OUTPUT_upTR}" "${PILEUP_OUTPUT}.gz" 2>/dev/null
modkit stats --regions "$REGION_DOWNSTREAM_BED" --min-coverage 3 -o "${STAT_OUTPUT_downTR}" "${PILEUP_OUTPUT}.gz" 2>/dev/null

#Extraction of methylation values and run uTR to get allele pattern
echo "Extraction of average methylation values" >> "$LOG_file"

TR_AVG_METHYLATION=$(awk 'NR==2 {if ($8 == "") print "."; else print $8}' "$STAT_OUTPUT_TR")
TR_COV_METH=$(awk 'NR==2 {if ($7 == "") print "."; else print $7}' "$STAT_OUTPUT_TR")
upTR_AVG_METHYLATION=$(awk 'NR==2 {if ($8 == "") print "."; else print $8}' "$STAT_OUTPUT_upTR")
upTR_COV_METH=$(awk 'NR==2 {if ($7 == "") print "."; else print $7}' "$STAT_OUTPUT_upTR")
downTR_AVG_METHYLATION=$(awk 'NR==2 {if ($8 == "") print "."; else print $8}' "$STAT_OUTPUT_downTR")
downTR_COV_METH=$(awk 'NR==2 {if ($7 == "") print "."; else print $7}' "$STAT_OUTPUT_downTR")


echo "Extraction of CpG methylation values" >> "$LOG_file"
PILEUP_TR_OUTPUT="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_TR_modkit_pileup.bed"

if bedtools intersect -a "${PILEUP_OUTPUT}.gz" -b "$REGION_BED" > "$PILEUP_TR_OUTPUT" 2>/dev/null; then
    if [[ -s "$PILEUP_TR_OUTPUT" ]]; then
        echo "Intersect succeeded, continuing with next step..." >> "$LOG_file"
        TR_CPG_METH=$(cut -f11 "$PILEUP_TR_OUTPUT" | paste -sd, -)
        TR_CPG_DEPTH=$(cut -f10 "$PILEUP_TR_OUTPUT" | paste -sd, -)
    else
        TR_CPG_METH="."
        TR_CPG_DEPTH="."
        echo "No CpG sites in the allele" >> "$LOG_file"
    fi
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
    "$UTR_TOOL_DIR" -f "$STR_ALLELE_FASTA" -y -o "$uTR_out"
    TR_PATTERN=$(extract_uTR_features "$uTR_out")
    echo -e "${TR_PATTERN}" >> "$LOG_file"
elif [[ "$HAPLOTYPE" == 1 ]]; then
    HP1_TR_AVG_VALUE="$TR_AVG_METHYLATION"
    HP1_TR_NCOV="$TR_COV_METH"
    HP1_upTR_AVG_VALUE="$upTR_AVG_METHYLATION"
    HP1_upTR_NCOV="$upTR_COV_METH"
    HP1_downTR_AVG_VALUE="$downTR_AVG_METHYLATION"
    HP1_downTR_NCOV="$downTR_COV_METH"
    TR_CPG_METH_HP1="$TR_CPG_METH"
    TR_CPG_DEPTH_HP1="$TR_CPG_DEPTH"
    "$UTR_TOOL_DIR" -f "$STR_ALLELE_FASTA" -y -o "$uTR_out"
    HP1_TR_PATTERN=$(extract_uTR_features "$uTR_out")
    echo -e "${HP1_TR_PATTERN}" >> "$LOG_file"
elif [[ "$HAPLOTYPE" == 2 ]]; then
    HP2_TR_AVG_VALUE="$TR_AVG_METHYLATION"
    HP2_TR_NCOV="$TR_COV_METH"
    HP2_upTR_AVG_VALUE="$upTR_AVG_METHYLATION"
    HP2_upTR_NCOV="$upTR_COV_METH"
    HP2_downTR_AVG_VALUE="$downTR_AVG_METHYLATION"
    HP2_downTR_NCOV="$downTR_COV_METH"
    TR_CPG_METH_HP2="$TR_CPG_METH"
    TR_CPG_DEPTH_HP2="$TR_CPG_DEPTH"
    "$UTR_TOOL_DIR" -f "$STR_ALLELE_FASTA" -y -o "$uTR_out"
    HP2_TR_PATTERN=$(extract_uTR_features "$uTR_out")
    echo -e "${HP2_TR_PATTERN}" >> "$LOG_file"
fi


# Cleanup temporary files
#rm -f "$REGION_FASTQ"
#rm -f "$REGION_SAM"
#rm -f "$REGION_BAM"
#rm -f "$REGION_BED"
#rm -f "$REGION_BAM"
#rm -f "${PILEUP_OUTPUT}.gz" "${PILEUP_OUTPUT}.gz.tbi"


done   

# Create vcf line
echo "Recreate vcf line with methylation information" >> "$LOG_file"    
echo "" >> "$LOG_file"
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

# Append to output VCF
echo -e "$NEW_LINE" > "${TMP_DIR}/${HEADER}.line.tsv"
	 
}

# Export functions and variables
export -f process_line
export -f modify_gt_based_on_pdp
export -f is_haploid
export -f extract_uTR_features

export REFERENCE_FASTA
export PHASED_BAM
export OUTPUT_DIR
export SAMPLE_ID
export EXTEND
export FLANKING_BASES
export HAPLOID_CHROMOSOMES
export TMP_DIR
export ALIGNMENTS
export METH
export UTR_TOOL_DIR
export LOGS

#NUM_JOBS=$(nproc)
#echo "Starting processing with ${NUM_JOBS} parallel jobs..."


bcftools view -H "$INPUT_VCF" | xargs -P 8 -n 1 -d '\n' bash -c 'process_line "$0"' 

cat "${TMP_DIR}"/*.line.tsv >> "$OUTPUT_VCF"

bcftools annotate   --remove 'INFO/NSKIP,INFO/NFILT,INFO/INEXACT_ALLELE,INFO/BPDIFFS,INFO/DP,INFO/DSNP,INFO/DFLANKINDEL,INFO/REFAC,INFO/AC' -Oz -o "$OUTPUT_VCF".gz "$OUTPUT_VCF" 2>/dev/null

bcftools sort "$OUTPUT_VCF".gz -Oz -o "$OUTPUT_SORTED_VCF" 2>/dev/null

rm "$OUTPUT_VCF" "$OUTPUT_VCF".gz
echo "Pipeline complete"


#echo ""
#echo "STEP5: GENERATE SUMMARY FILE"
#echo ""
#
## Print header for summary file
#echo -e "CHROM\tPOS\tID\tREF_MOTIF\tTR_REF_LENGTH\t${SAMPLE_ID}_GT\t${SAMPLE_ID}_TR_LEN\tTR_PATTERN\tTR_AM\tTR_N_METH_VALID\tUPSTREAM_TR_AM\tUPSTREAM_TR_N_METH_VALID\tDOWNSTREAM_TR_AM\tDOWNSTREAM_TR_N_METH_VALID" > "$OUTPUT_SUMMARY"
#
## Process each non-header line in the output VCF
#awk -F'\t' '
#BEGIN { OFS="\t" }
#!/^#/ {
#    chrom=$1;
#    pos=$2;
#    id=$3;
#    ref=$4;
#    info=$8;
#    tr_am=$12;
#    tr_n_meth_valid=$13;
#    upstream_tr_am=$14;
#    upstream_tr_n_meth_valid=$15;
#    downstream_tr_am=$16;
#    downstream_tr_n_meth_valid=$17;
#    tr_len=$18;
#    tr_pat=$20;
#
#    # Compute REF allele length
#    ref_length = length(ref);
#
#    # Extract values from INFO field
#    motif = "."; 
#    if (match(info, /MOTIF=([^;]+)/, arr)) motif = arr[1];
#
#    # Extract GT (first field of sample column)
#    gt = $10;
#    split(gt, gt_fields, ":");
#    gt = gt_fields[1];
#
#    # Print extracted values
#    print chrom, pos, id, motif, ref_length, gt, tr_len, tr_pat, tr_am, tr_n_meth_valid, upstream_tr_am, upstream_tr_n_meth_valid, downstream_tr_am, downstream_tr_n_meth_valid;
#}' "$OUTPUT_VCF" >> "$OUTPUT_SUMMARY"
#
#echo ""
#echo "Output files are in $OUTPUT_DIR. Remove temporary files stored in: $TMP_DIR"
## Uncomment the following line to automatically remove temporary files after the run
#rm -rf "$TMP_DIR"
#echo ""
#echo "--------------------------------"
#echo "Pipeline completed successfully!"
#echo "--------------------------------"
#echo ""
