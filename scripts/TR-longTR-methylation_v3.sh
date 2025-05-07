#!/bin/bash

# A script to analyse TRs genotyped using LongTR, adding allele-specific methylation information.

# Usage:
# ./script.sh -v <vcf_file> -r <reference_fasta> -i <phased_bam> -o <output_dir> -s <sample_id> -e <extend> -f <flanking_bases> -h <haploid_chromosomes>

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
for tool in samtools minimap2 awk modkit; do
    if ! command -v "$tool" &> /dev/null; then
        echo "Error: $tool is not installed. Please install it and try again."
        exit 1
    fi
done

# Check if a chromosome is haploid
IFS=',' read -r -a HAPLOID_ARRAY <<< "$HAPLOID_CHROMOSOMES"

is_haploid() {
    local chrom="$1"
    for haploid_chrom in "${HAPLOID_ARRAY[@]}"; do
        if [[ "$chrom" == "$haploid_chrom" ]]; then
            return 0
        fi
    done
    return 1
}

# Create output directory if it doesnt exist
mkdir -p "$OUTPUT_DIR"

# Define paths inside the output directory
INPUT_VCF="${OUTPUT_DIR}/${SAMPLE_ID}_input_sorted.vcf"
MULTIFASTA="${OUTPUT_DIR}/${SAMPLE_ID}_alleles.fasta"
OUTPUT_BED="${OUTPUT_DIR}/${SAMPLE_ID}_alleles.bed"
OUTPUT_UPSTREAM_BED="${OUTPUT_DIR}/${SAMPLE_ID}_upstream_alleles.bed"
OUTPUT_DOWNSTREAM_BED="${OUTPUT_DIR}/${SAMPLE_ID}_downstream_alleles.bed"
OUTPUT_VCF="${OUTPUT_DIR}/${SAMPLE_ID}_methylated.vcf"
TMP_DIR="${OUTPUT_DIR}/temp"
ALIGNMENTS="${OUTPUT_DIR}/alignments"
METH="${OUTPUT_DIR}/methylation"
OUTPUT_SUMMARY="${OUTPUT_DIR}/${SAMPLE_ID}_TR_summary.tsv"

mkdir -p "$TMP_DIR" "$ALIGNMENTS" "$METH"

if [[ "$VCF_FILE" == *.gz ]]; then
    gunzip -c "$VCF_FILE" > "${OUTPUT_DIR}/${SAMPLE_ID}_input.vcf"
    VCF_FILE="${OUTPUT_DIR}/${SAMPLE_ID}_input.vcf"
fi

grep "^#" "$VCF_FILE" > "$INPUT_VCF" && grep -v "^#" "$VCF_FILE" | sort -V -k1,1 -k2,2 >> "$INPUT_VCF"

# Prepare output files
> "$MULTIFASTA"
> "$OUTPUT_BED"
> "$OUTPUT_VCF"

# Copy existing headers from the input VCF
grep "^##" "$VCF_FILE" > "$OUTPUT_VCF"

# Add header information to VCF output file
echo "##METHYLATION: TR_AM - Description=\"Average methylation percentage for TR alleles (ref,alt)\">" >> "$OUTPUT_VCF"
echo "##METHYLATION: TR_N_METH_VALID - Description=\"Count of valid sites used for AM calculation for TR alleles (ref,alt)\">" >> "$OUTPUT_VCF"
echo "##METHYLATION: UPSTREAM_TR_AM - Description=\"Average methylation percentage for upstream sequence (region specified by flanking bp) of TR alleles (ref,alt)\">" >> "$OUTPUT_VCF"
echo "##METHYLATION: UPSTREAM_TR_N_METH_VALID - Description=\"Count of valid sites used for UPSTREAM AM calculation for TR alleles (ref,alt)\">" >> "$OUTPUT_VCF"
echo "##METHYLATION: DOWNSTREAM_TR_AM - Description=\"Average methylation percentage for downstream sequence (region specified by flanking bp) of TR alleles (ref,alt)\">" >> "$OUTPUT_VCF"
echo "##METHYLATION: DOWNSTREAM_TR_N_METH_VALID - Description=\"Count of valid sites used for downstream AM calculation for TR alleles (ref,alt)\">" >> "$OUTPUT_VCF"
echo "##TR SEQUENCE COMPOSITION: LENGTH OF TR ALLELE (HP1|HP2)\">" >> "$OUTPUT_VCF"
echo "##TR SEQUENCE COMPOSITION: uTR info output per allele (HP1;HP2). Format: (A,B,C), where A shows the length of the DNA string, B means the number of types detected, and C is the dissimilarity ratio between the tandem repeat pattern and the string\">" >> "$OUTPUT_VCF"
echo "##TR SEQUENCE COMPOSITION: TR pattern as defined by uTR (HP1,HP2)\">" >> "$OUTPUT_VCF"
echo "##TR SEQUENCE COMPOSITION: TR motif decomposition as defined by uTR (HP1,HP2)\">" >> "$OUTPUT_VCF"
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$SAMPLE_ID\t$TR_ID\tTR_AM\tTR_N_METH_VALID\tUPSTREAM_TR_AM\tUPSTREAM_TR_N_METH_VALID\tDOWNSTREAM_TR_AM\tDOWNSTREAM_TR_N_METH_VALID\tTR_LEN\tTR_INFO_VALUE\tTR_PAT\tTR_DECOMP" >> "$OUTPUT_VCF"

# Function to check if a VCF entry lacks phasing information
# Function to modify the GT field based on PDP values and handle haploid chromosomes
modify_gt_based_on_pdp() {
    local LINE FORMAT INFO GT_INDEX PDP_INDEX GT_VALUE PDP_VALUE MODIFIED_GT CHROM

    LINE="$1"
    FORMAT=$(echo "$LINE" | cut -f9)
    INFO=$(echo "$LINE" | cut -f10)
    CHROM=$(echo "$LINE" | cut -f1)  # Extract chromosome name

    # Find the indexes of GT and PDP in the FORMAT field
    GT_INDEX=$(echo "$FORMAT" | awk -F':' '{for (i=1; i<=NF; i++) if ($i=="GT") print i}')
    PDP_INDEX=$(echo "$FORMAT" | awk -F':' '{for (i=1; i<=NF; i++) if ($i=="PDP") print i}')

    # Extract GT value
    GT_VALUE=$(echo "$INFO" | cut -d':' -f$GT_INDEX)

    # Check PDP value
    if is_haploid "$CHROM"; then
        # Haploid chromosome case (no PDP field)
        MODIFIED_GT="${GT_VALUE}|."
    else
        # Extract PDP value
        PDP_VALUE=$(echo "$INFO" | cut -d':' -f$PDP_INDEX)

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
    MODIFIED_INFO=$(echo "$INFO" | awk -v gt_index="$GT_INDEX" -v new_gt="$MODIFIED_GT" 'BEGIN{FS=OFS=":"}{$gt_index=new_gt; print}')

    # Construct the modified VCF line
    echo -e "$(echo "$LINE" | cut -f1-9)\t$MODIFIED_INFO"
}


# First processing loop: read VCF file and generate TR-specific FASTA and BED entries
echo ""
echo "STEP1: PREPARING INPUT FILES"
echo ""

while read -r LINE; do
    # Skip VCF header lines
    if [[ $LINE == \#* ]]; then
        continue
    fi

    # Modify GT field based on PDP values
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

    # Parse alternate alleles
    IFS=',' read -ra ALT_ALLELES <<< "$ALT"

    # Determine repeat extension boundaries
    UPSTREAM_START=$((TR_START - EXTEND))
    UPSTREAM_END=$((TR_START - 1))
    DOWNSTREAM_START=$((TR_END + 1))
    DOWNSTREAM_END=$((TR_END + EXTEND))

    # Extract flanking sequences
    UPSTREAM_SEQ=$(samtools faidx "$REFERENCE_FASTA" "${CHROM}:${UPSTREAM_START}-${UPSTREAM_END}" | tail -n +2 | tr -d '\n' | tr '[:upper:]' '[:lower:]')
    DOWNSTREAM_SEQ=$(samtools faidx "$REFERENCE_FASTA" "${CHROM}:${DOWNSTREAM_START}-${DOWNSTREAM_END}" | tail -n +2 | tr -d '\n' | tr '[:upper:]' '[:lower:]')


    # Generate FASTA and BED entries based on genotype
    if is_haploid "$CHROM"; then
        # Haploid chromosome logic
        case "$GENOTYPE" in
            "0|.")
                REF_SEQ="${UPSTREAM_SEQ}$(echo "${REF}" | tr '[:lower:]' '[:upper:]')${DOWNSTREAM_SEQ}"
                HEADER="${CHROM}_${POS}_${TR_ID}_HP1_REF"
                TR_START_REL=$((EXTEND + 1))
                TR_END_REL=$((EXTEND + ${#REF}))
                echo ">$HEADER" >> "$MULTIFASTA"
                echo "$REF_SEQ" >> "$MULTIFASTA"
                echo -e ">$HEADER\n$REF" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP1.fasta"
                echo -e "$HEADER\t$TR_START_REL\t$TR_END_REL\t$TR_ID" >> "$OUTPUT_BED"
                ;;
            "1|.")
               ALT_SEQ="${UPSTREAM_SEQ}$(echo "${ALT_ALLELES[0]}" | tr '[:lower:]' '[:upper:]')${DOWNSTREAM_SEQ}"
               HEADER="${CHROM}_${POS}_${TR_ID}_HP1_ALT"
               TR_START_ALT=$((EXTEND + 1))
               TR_END_ALT=$((EXTEND + ${#ALT_ALLELES[0]}))
               echo ">$HEADER" >> "$MULTIFASTA"
               echo "$ALT_SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER\n${ALT_ALLELES[0]}" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP1.fasta"
               echo -e "$HEADER\t$TR_START_ALT\t$TR_END_ALT\t$TR_ID" >> "$OUTPUT_BED"
               ;;
        esac
    else
        # Diploid chromosome logic
        case "$GENOTYPE" in
            "1|2")
               # Heterozygous two alternate alleles
               ALT_SEQ1="${UPSTREAM_SEQ}$(echo "${ALT_ALLELES[0]}" | tr '[:lower:]' '[:upper:]')${DOWNSTREAM_SEQ}"
               ALT_SEQ2="${UPSTREAM_SEQ}$(echo "${ALT_ALLELES[1]}" | tr '[:lower:]' '[:upper:]')${DOWNSTREAM_SEQ}"
               HEADER_ALT1="${CHROM}_${POS}_${TR_ID}_HP1_ALT1"
               HEADER_ALT2="${CHROM}_${POS}_${TR_ID}_HP2_ALT2"
               TR_START_ALT1=$((EXTEND + 1))
               TR_END_ALT1=$((EXTEND + ${#ALT_ALLELES[0]}))
               TR_START_ALT2=$((EXTEND + 1))
               TR_END_ALT2=$((EXTEND + ${#ALT_ALLELES[1]}))
               echo ">$HEADER_ALT1" >> "$MULTIFASTA"
               echo "$ALT_SEQ1" >> "$MULTIFASTA"
               echo ">$HEADER_ALT2" >> "$MULTIFASTA"
               echo -e ">$HEADER_ALT1\n${ALT_ALLELES[0]}" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP1.fasta"
               echo "$ALT_SEQ2" >> "$MULTIFASTA"
               echo -e ">$HEADER_ALT2\n${ALT_ALLELES[1]}" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP2.fasta"
               echo -e "$HEADER_ALT1\t$TR_START_ALT1\t$TR_END_ALT1\t$TR_ID" >> "$OUTPUT_BED"
               echo -e "$HEADER_ALT2\t$TR_START_ALT2\t$TR_END_ALT2\t$TR_ID" >> "$OUTPUT_BED"
               ;;
           "2|1")
               # Heterozygous two alternate alleles
               ALT_SEQ1="${UPSTREAM_SEQ}$(echo "${ALT_ALLELES[1]}" | tr '[:lower:]' '[:upper:]')${DOWNSTREAM_SEQ}"
               ALT_SEQ2="${UPSTREAM_SEQ}$(echo "${ALT_ALLELES[0]}" | tr '[:lower:]' '[:upper:]')${DOWNSTREAM_SEQ}"
               HEADER_ALT1="${CHROM}_${POS}_${TR_ID}_HP1_ALT2"
               HEADER_ALT2="${CHROM}_${POS}_${TR_ID}_HP2_ALT1"
               TR_START_ALT1=$((EXTEND + 1))
               TR_END_ALT1=$((EXTEND + ${#ALT_ALLELES[1]}))
               TR_START_ALT2=$((EXTEND + 1))
               TR_END_ALT2=$((EXTEND + ${#ALT_ALLELES[0]}))
               echo ">$HEADER_ALT1" >> "$MULTIFASTA"
               echo "$ALT_SEQ1" >> "$MULTIFASTA"
               echo -e ">$HEADER_ALT1\n${ALT_ALLELES[1]}" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP1.fasta"
               echo ">$HEADER_ALT2" >> "$MULTIFASTA"
               echo "$ALT_SEQ2" >> "$MULTIFASTA"
               echo -e ">$HEADER_ALT2\n${ALT_ALLELES[0]}" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP2.fasta"
               echo -e "$HEADER_ALT1\t$TR_START_ALT1\t$TR_END_ALT1\t$TR_ID" >> "$OUTPUT_BED"
               echo -e "$HEADER_ALT2\t$TR_START_ALT2\t$TR_END_ALT2\t$TR_ID" >> "$OUTPUT_BED"
               ;;
           "0|0")
               # Homozygous reference
               REF_UPPER=$(echo "$REF" | tr '[:lower:]' '[:upper:]')
               SEQ="${UPSTREAM_SEQ}${REF_UPPER}${DOWNSTREAM_SEQ}"
               TR_START_REL=$((EXTEND + 1))
               TR_END_REL=$((EXTEND + ${#REF}))
               HEADER_HP1="${CHROM}_${POS}_${TR_ID}_HP1_REF"
               HEADER_HP2="${CHROM}_${POS}_${TR_ID}_HP2_REF"
               echo ">$HEADER_HP1" >> "$MULTIFASTA"
               echo "$SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP1\n$REF" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP1.fasta"
               echo ">$HEADER_HP2" >> "$MULTIFASTA"
               echo "$SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP2\n$REF" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP2.fasta"
               echo -e "$HEADER_HP1\t$TR_START_REL\t$TR_END_REL\t$TR_ID" >> "$OUTPUT_BED"
               echo -e "$HEADER_HP2\t$TR_START_REL\t$TR_END_REL\t$TR_ID" >> "$OUTPUT_BED"
               ;;
           "0|1")
               # Heterozygous reference and first alternate
               REF_UPPER=$(echo "$REF" | tr '[:lower:]' '[:upper:]')
               ALT_UPPER=$(echo "${ALT_ALLELES[0]}" | tr '[:lower:]' '[:upper:]')
               REF_SEQ="${UPSTREAM_SEQ}${REF_UPPER}${DOWNSTREAM_SEQ}"
               ALT_SEQ="${UPSTREAM_SEQ}${ALT_UPPER}${DOWNSTREAM_SEQ}"
               TR_START_REF=$((EXTEND + 1))
               TR_END_REF=$((EXTEND + ${#REF}))
               TR_START_ALT=$((EXTEND + 1))
               TR_END_ALT=$((EXTEND + ${#ALT_ALLELES[0]}))
               HEADER_REF="${CHROM}_${POS}_${TR_ID}_HP1_REF"
               echo ">$HEADER_REF" >> "$MULTIFASTA"
               echo "$REF_SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP1\n$REF" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP1.fasta"
               HEADER_ALT="${CHROM}_${POS}_${TR_ID}_HP2_ALT1"
               echo ">$HEADER_ALT" >> "$MULTIFASTA"
               echo "$ALT_SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP2\n${ALT_ALLELES[0]}" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP2.fasta"
               echo -e "$HEADER_REF\t$TR_START_REF\t$TR_END_REF\t$TR_ID" >> "$OUTPUT_BED"
               echo -e "$HEADER_ALT\t$TR_START_ALT\t$TR_END_ALT\t$TR_ID" >> "$OUTPUT_BED"
               ;;
           "1|0")
               # Heterozygous reference and first alternate
               REF_UPPER=$(echo "$REF" | tr '[:lower:]' '[:upper:]')
               ALT_UPPER=$(echo "${ALT_ALLELES[0]}" | tr '[:lower:]' '[:upper:]')
               REF_SEQ="${UPSTREAM_SEQ}${REF_UPPER}${DOWNSTREAM_SEQ}"
               ALT_SEQ="${UPSTREAM_SEQ}${ALT_UPPER}${DOWNSTREAM_SEQ}"
               TR_START_REF=$((EXTEND + 1))
               TR_END_REF=$((EXTEND + ${#REF}))
               TR_START_ALT=$((EXTEND + 1))
               TR_END_ALT=$((EXTEND + ${#ALT_ALLELES[0]}))
               HEADER_ALT="${CHROM}_${POS}_${TR_ID}_HP1_ALT1"
               HEADER_REF="${CHROM}_${POS}_${TR_ID}_HP2_REF"
               echo ">$HEADER_ALT" >> "$MULTIFASTA"
               echo "$ALT_SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP1\n${ALT_ALLELES[0]}" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP1.fasta"
               echo ">$HEADER_REF" >> "$MULTIFASTA"
               echo "$REF_SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP2\n$REF" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP2.fasta"
               echo -e "$HEADER_ALT\t$TR_START_ALT\t$TR_END_ALT\t$TR_ID" >> "$OUTPUT_BED"
               echo -e "$HEADER_REF\t$TR_START_REF\t$TR_END_REF\t$TR_ID" >> "$OUTPUT_BED"
               ;;
           "1|1")
               # Homozygous first alternate
               ALT_SEQ="${UPSTREAM_SEQ}$(echo "${ALT_ALLELES[0]}" | tr '[:lower:]' '[:upper:]')${DOWNSTREAM_SEQ}"
               HEADER_HP1="${CHROM}_${POS}_${TR_ID}_HP1_ALT1"
               HEADER_HP2="${CHROM}_${POS}_${TR_ID}_HP2_ALT1"
               TR_START_REL=$((EXTEND + 1))
               TR_END_REL=$((EXTEND + ${#ALT_ALLELES[0]}))
               echo ">$HEADER_HP1" >> "$MULTIFASTA"
               echo "$ALT_SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP1\n${ALT_ALLELES[0]}" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP1.fasta"
               echo ">$HEADER_HP2" >> "$MULTIFASTA"
               echo "$ALT_SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP2\n${ALT_ALLELES[0]}" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP2.fasta"
               echo -e "$HEADER_HP1\t$TR_START_REL\t$TR_END_REL\t$TR_ID" >> "$OUTPUT_BED"
               echo -e "$HEADER_HP2\t$TR_START_REL\t$TR_END_REL\t$TR_ID" >> "$OUTPUT_BED"
               ;;
           ".|.")
               # Unphased or missing data for both alleles
               HEADER_UNPHASED_HP1="${CHROM}_${POS}_${TR_ID}_unphased1"
               HEADER_UNPHASED_HP2="${CHROM}_${POS}_${TR_ID}_unphased2"
               UNPHASED_SEQ="${UPSTREAM_SEQ}$(echo "$REF" | tr '[:lower:]' '[:upper:]')${DOWNSTREAM_SEQ}"
               TR_START_UNPHASED=$((EXTEND + 1))
               TR_END_UNPHASED=$((EXTEND + ${#REF}))
               echo ">$HEADER_UNPHASED_HP1" >> "$MULTIFASTA"
               echo "$UNPHASED_SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_UNPHASED_HP1\n$REF" > "${TMP_DIR}/${CHROM}_${TR_ID}_unphased1.fasta"
               echo ">$HEADER_UNPHASED_HP2" >> "$MULTIFASTA"
               echo "$UNPHASED_SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_UNPHASED_HP2\n$REF" > "${TMP_DIR}/${CHROM}_${TR_ID}_unphased2.fasta"
               echo -e "$HEADER_UNPHASED_HP1\t$TR_START_UNPHASED\t$TR_END_UNPHASED\t$TR_ID" >> "$OUTPUT_BED"
               echo -e "$HEADER_UNPHASED_HP2\t$TR_START_UNPHASED\t$TR_END_UNPHASED\t$TR_ID" >> "$OUTPUT_BED"
               ;;
           ".|0")
               # Unphased and reference
               REF_UPPER=$(echo "$REF" | tr '[:lower:]' '[:upper:]')
               SEQ="${UPSTREAM_SEQ}${REF_UPPER}${DOWNSTREAM_SEQ}"
               TR_START_REL=$((EXTEND + 1))
               TR_END_REL=$((EXTEND + ${#REF}))
               HEADER_HP1="${CHROM}_${POS}_${TR_ID}_unphased1"
               HEADER_HP2="${CHROM}_${POS}_${TR_ID}_HP2_REF"
               echo ">$HEADER_HP1" >> "$MULTIFASTA"
               echo "$SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP1\n$REF" > "${TMP_DIR}/${CHROM}_${TR_ID}_unphased1.fasta"
               echo ">$HEADER_HP2" >> "$MULTIFASTA"
               echo "$SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP2\n$REF" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP2.fasta"
               echo -e "$HEADER_HP1\t$TR_START_REL\t$TR_END_REL\t$TR_ID" >> "$OUTPUT_BED"
               echo -e "$HEADER_HP2\t$TR_START_REL\t$TR_END_REL\t$TR_ID" >> "$OUTPUT_BED"
               ;;
            "0|.")
               # Reference allele and unphased
               REF_UPPER=$(echo "$REF" | tr '[:lower:]' '[:upper:]')
               SEQ="${UPSTREAM_SEQ}${REF_UPPER}${DOWNSTREAM_SEQ}"
               TR_START_REL=$((EXTEND + 1))
               TR_END_REL=$((EXTEND + ${#REF}))
               HEADER_HP1="${CHROM}_${POS}_${TR_ID}_HP1_REF"
               HEADER_HP2="${CHROM}_${POS}_${TR_ID}_unphased2"
               echo ">$HEADER_HP1" >> "$MULTIFASTA"
               echo "$SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP1\n$REF" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP1.fasta"
               echo ">$HEADER_HP2" >> "$MULTIFASTA"
               echo "$SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP2\n$REF" > "${TMP_DIR}/${CHROM}_${TR_ID}_unphased2.fasta"
               echo -e "$HEADER_HP1\t$TR_START_REL\t$TR_END_REL\t$TR_ID" >> "$OUTPUT_BED"
               echo -e "$HEADER_HP2\t$TR_START_REL\t$TR_END_REL\t$TR_ID" >> "$OUTPUT_BED"
               ;;
            ".|1")
               # Heterozygous reference and first alternate
               REF_UPPER=$(echo "$REF" | tr '[:lower:]' '[:upper:]')
               ALT_UPPER=$(echo "${ALT_ALLELES[0]}" | tr '[:lower:]' '[:upper:]')
               REF_SEQ="${UPSTREAM_SEQ}${REF_UPPER}${DOWNSTREAM_SEQ}"
               ALT_SEQ="${UPSTREAM_SEQ}${ALT_UPPER}${DOWNSTREAM_SEQ}"
               TR_START_REF=$((EXTEND + 1))
               TR_END_REF=$((EXTEND + ${#REF}))
               TR_START_ALT=$((EXTEND + 1))
               TR_END_ALT=$((EXTEND + ${#ALT_ALLELES[0]}))
               HEADER_HP1="${CHROM}_${POS}_${TR_ID}_unphased1"
               echo ">$HEADER_HP1" >> "$MULTIFASTA"
               echo "$REF_SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP1\n$REF" > "${TMP_DIR}/${CHROM}_${TR_ID}_unphased1.fasta"
               HEADER_HP2="${CHROM}_${POS}_${TR_ID}_HP2_ALT1"
               echo ">$HEADER_HP2" >> "$MULTIFASTA"
               echo "$ALT_SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP2\n${ALT_ALLELES[0]}" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP2.fasta"
               echo -e "$HEADER_HP1\t$TR_START_REF\t$TR_END_REF\t$TR_ID" >> "$OUTPUT_BED"
               echo -e "$HEADER_HP2\t$TR_START_ALT\t$TR_END_ALT\t$TR_ID" >> "$OUTPUT_BED"
               ;;
           "1|.")
               # Heterozygous reference and first alternate
               REF_UPPER=$(echo "$REF" | tr '[:lower:]' '[:upper:]')
               ALT_UPPER=$(echo "${ALT_ALLELES[0]}" | tr '[:lower:]' '[:upper:]')
               REF_SEQ="${UPSTREAM_SEQ}${REF_UPPER}${DOWNSTREAM_SEQ}"
               ALT_SEQ="${UPSTREAM_SEQ}${ALT_UPPER}${DOWNSTREAM_SEQ}"
               TR_START_REF=$((EXTEND + 1))
               TR_END_REF=$((EXTEND + ${#REF}))
               TR_START_ALT=$((EXTEND + 1))
               TR_END_ALT=$((EXTEND + ${#ALT_ALLELES[0]}))
               HEADER_HP1="${CHROM}_${POS}_${TR_ID}_HP1_ALT1"
               echo ">$HEADER_HP1" >> "$MULTIFASTA"
               echo "$ALT_SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP1\n${ALT_ALLELES[0]}" > "${TMP_DIR}/${CHROM}_${TR_ID}_HP1.fasta"
               HEADER_HP2="${CHROM}_${POS}_${TR_ID}_unphased2"
               echo ">$HEADER_HP2" >> "$MULTIFASTA"
               echo "$REF_SEQ" >> "$MULTIFASTA"
               echo -e ">$HEADER_HP2\n$REF" > "${TMP_DIR}/${CHROM}_${TR_ID}_unphased2.fasta"
               echo -e "$HEADER_HP1\t$TR_START_ALT\t$TR_END_ALT\t$TR_ID" >> "$OUTPUT_BED"
               echo -e "$HEADER_HP2\t$TR_START_REF\t$TR_END_REF\t$TR_ID" >> "$OUTPUT_BED"
               ;;
        esac
    fi
done < "$INPUT_VCF"

# Create an upstream-shifted BED file (flanking bp from start)
awk -v flank="$FLANKING_BASES" '{OFS="\t"} {start=$2-flank; if (start < 0) start=0; print $1, start, $2}' "$OUTPUT_BED" > "$OUTPUT_UPSTREAM_BED" 

# Downstream: from end to (end + flank)
awk -v flank="$FLANKING_BASES" '{OFS="\t"} {print $1, $3, $3+flank}' "$OUTPUT_BED" > "$OUTPUT_DOWNSTREAM_BED"

# Second processing loop where we extract the reads that span each TR region
# and remap them against newly created reference TR sequences.

echo ""
echo "STEP2: REMAPPING AND METHYLATION ANALYSIS"
echo ""

while read -r LINE; do
    # Skip VCF header lines
    if [[ $LINE == \#* ]]; then
        continue
    fi

    MODIFIED_LINE=$(modify_gt_based_on_pdp "$LINE")

    # Parse and extract same parameters than in first loop
    CHROM=$(echo "$MODIFIED_LINE" | cut -f1)
    POS=$(echo "$MODIFIED_LINE" | cut -f2)
    TR_ID=$(echo "$MODIFIED_LINE" | cut -f3)
    TR_START=$(echo "$MODIFIED_LINE" | grep -o "START=[0-9]*" | cut -d'=' -f2)
    TR_END=$(echo "$MODIFIED_LINE" | grep -o "END=[0-9]*" | cut -d"=" -f2)

    # Log progress
    echo -e ".............Processing ${CHROM}:${TR_ID}............."
    echo ""


    echo "Extracting overlapping reads"

    # Check again if haploid
    if is_haploid "$CHROM"; then
        # Haploid chromosome logic: process only one sequence (ALT or REF)
        for HAPLOTYPE in ALT REF; do
            # Extract reads overlapping the region
            REGION_BAM="${TMP_DIR}/${TR_ID}_extracted.bam"
            samtools view -h "$PHASED_BAM" "$CHROM:$TR_START-$TR_END" | samtools addreplacerg -r "ID:${TR_ID}" - | samtools view -b - > "$REGION_BAM"

            # Check if haplotype exists in MULTIFASTA
            if ! grep -q -E ">${CHROM}_${POS}_${TR_ID}_HP1_${HAPLOTYPE}$" "$MULTIFASTA"; then
               continue
            fi

            # Convert to FASTQ
            REGION_FASTQ="${TMP_DIR}/${TR_ID}.extracted.fastq"
            samtools fastq -t -T MM,ML,HP,PS "$REGION_BAM" > "$REGION_FASTQ" 2>/dev/null

            # Proceed with processing the haploid allele sequence
            REGION_FASTA="${TMP_DIR}/${TR_ID}_${HAPLOTYPE}_reference.fasta"
            grep -A 1 -E ">${CHROM}_${POS}_${TR_ID}_HP1_${HAPLOTYPE}" "$MULTIFASTA" > "$REGION_FASTA"
            if [[ ! -s "$REGION_FASTA" ]]; then
                echo "Error: no FASTA sequence found for ${CHROM}:${TR_ID}:${HAPLOTYPE}" >&2
                continue
            fi
            samtools faidx "$REGION_FASTA"

            # Map with minimap2
            echo "Remapping to TR ${HAPLOTYPE} sequence"
            REGION_SAM="${TMP_DIR}/${TR_ID}_mapped.sam"
            minimap2 -ax map-ont -y --sam-hit-only "$REGION_FASTA" "$REGION_FASTQ" > "$REGION_SAM" 2>/dev/null

            # Convert, sort, and index BAM
            REGION_SORTED_BAM="${ALIGNMENTS}/${CHROM}_${TR_ID}_mapped.sorted.bam"
            samtools view -bS "$REGION_SAM" | samtools sort -o "$REGION_SORTED_BAM"
            samtools index "$REGION_SORTED_BAM"

            # Perform modkit pileup on the generated BAM file
            echo "Running modkit pileup for ${HAPLOTYPE}"
            PILEUP_OUTPUT="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_modkit_pileup.bed"
            modkit pileup "${REGION_SORTED_BAM}" "${PILEUP_OUTPUT}" --cpg --ref "${REGION_FASTA}" --ignore h --combine-strands --mod-threshold m:0.8 2>/dev/null
            bgzip "$PILEUP_OUTPUT"
            tabix "${PILEUP_OUTPUT}.gz"

            # Perform modkit stats for the TR region for each phased output
            STAT_OUTPUT_TR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_modkit_stats.tsv"
            STAT_OUTPUT_upTR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_upstream_modkit_stats.tsv"
            STAT_OUTPUT_downTR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_downstream_modkit_stats.tsv"
            echo "Running modkit stats for ${HAPLOTYPE}"
            modkit stats --regions "$OUTPUT_BED" --min-coverage 3 -o "${STAT_OUTPUT_TR}" "${PILEUP_OUTPUT}.gz" 2>/dev/null
            modkit stats --regions "$OUTPUT_UPSTREAM_BED" --min-coverage 3 -o "${STAT_OUTPUT_upTR}" "${PILEUP_OUTPUT}.gz" 2>/dev/null
            modkit stats --regions "$OUTPUT_DOWNSTREAM_BED" --min-coverage 3 -o "${STAT_OUTPUT_downTR}" "${PILEUP_OUTPUT}.gz" 2>/dev/null
            echo "Sucessfully processed ${CHROM}:${TR_ID}:${HAPLOTYPE}!"
            echo ""
        done
    else
        for HAPLOTYPE in HP1_REF HP2_REF HP1_ALT1 HP1_ALT2 HP2_ALT1 HP2_ALT2 unphased1 unphased2; do
            # Extract reads overlapping the region
            REGION_BAM="${TMP_DIR}/${TR_ID}_extracted.bam"
            samtools view -h "$PHASED_BAM" "$CHROM:$TR_START-$TR_END" | samtools addreplacerg -r "ID:${TR_ID}" - | samtools view -b - > "$REGION_BAM"

            # Check if haplotype exists in MULTIFASTA
            if ! grep -q -E ">${CHROM}_${POS}_${TR_ID}_${HAPLOTYPE}" "$MULTIFASTA"; then
               continue
            fi

            # Filter reads by HP tag
            if [[ "$HAPLOTYPE" == HP1* ]]; then
                echo "Filter HP1 reads"
                FILTERED_BAM="${TMP_DIR}/${CHROM}_${TR_ID}_extracted_HP1.bam"
                samtools view -h "$REGION_BAM" | awk '/^@/ || /HP:i:1/' | samtools view -b - > "$FILTERED_BAM"
            elif [[ "$HAPLOTYPE" == HP2* ]]; then
                echo "Filter HP2 reads"
                FILTERED_BAM="${TMP_DIR}/${CHROM}_${TR_ID}_extracted_HP2.bam"
                samtools view -h "$REGION_BAM" | awk '/^@/ || /HP:i:2/' | samtools view -b - > "$FILTERED_BAM"
            else
                echo "Filter unphased reads"
                FILTERED_BAM="${TMP_DIR}/${CHROM}_${TR_ID}_extracted_UNPHASED.bam"
                samtools view -h "$REGION_BAM" | awk '/^@/ || !/HP:i:/' | samtools view -b - > "$FILTERED_BAM"
            fi

            # Validate the bam file
            if [[ ! -s "$FILTERED_BAM" ]]; then
                echo "Error: no reads extracted for $TR_ID, $HAPLOTYPE"
                continue
            fi

            # Convert to FASTQ
            REGION_FASTQ="${TMP_DIR}/${CHROM}_${TR_ID}_${HAPLOTYPE}.extracted.fastq"
            samtools fastq -t -T MM,ML,HP,PS "$FILTERED_BAM" > "$REGION_FASTQ" 2>/dev/null

            # Proceed with processing the haplotype
            REGION_FASTA="${TMP_DIR}/${CHROM}_${TR_ID}_${HAPLOTYPE}_reference.fasta"
            grep -A 1 -E ">${CHROM}_${POS}_${TR_ID}_${HAPLOTYPE}" "$MULTIFASTA" > "$REGION_FASTA"
            samtools faidx "$REGION_FASTA"

            # Map with minimap2
            echo "Remapping to TR $HAPLOTYPE sequence"
            REGION_SAM="${TMP_DIR}/${CHROM}_${TR_ID}_${HAPLOTYPE}_mapped.sam"
            minimap2 -ax map-ont -y --sam-hit-only "$REGION_FASTA" "$REGION_FASTQ" > "$REGION_SAM" 2>/dev/null

            # Convert, sort, and index BAM
            REGION_SORTED_BAM="${ALIGNMENTS}/${CHROM}_${TR_ID}_${HAPLOTYPE}_mapped.sorted.bam"
            samtools view -bS "$REGION_SAM" | samtools sort -o "$REGION_SORTED_BAM" 2>/dev/null
            samtools index "$REGION_SORTED_BAM"

            # Perform modkit pileup on the generated BAM file
            echo "Running modkit pileup for ${HAPLOTYPE}"
            PILEUP_OUTPUT="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_modkit_pileup.bed"
            modkit pileup "${REGION_SORTED_BAM}" "${PILEUP_OUTPUT}" --cpg --ref "${REGION_FASTA}" --ignore h --combine-strands --mod-threshold m:0.8 2>/dev/null

            # Check for phased output files
            if  [[ -e "$PILEUP_OUTPUT" && -s "$PILEUP_OUTPUT" ]]; then
                 echo "Modkit pileup completed successfully"
                 # Compress and index the phased pileup output
                 bgzip "$PILEUP_OUTPUT"
                 tabix "${PILEUP_OUTPUT}.gz"

                 # Perform modkit stats for the TR region for each phased output
                 STAT_OUTPUT_TR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_modkit_stats.tsv"
                 STAT_OUTPUT_upTR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_upstream_modkit_stats.tsv"
                 STAT_OUTPUT_downTR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_downstream_modkit_stats.tsv"
                 echo "Running modkit stats for TR ${HAPLOTYPE}"
                 modkit stats --regions "$OUTPUT_BED" --min-coverage 3 -o "${STAT_OUTPUT_TR}" "${PILEUP_OUTPUT}.gz" 2>/dev/null
                 modkit stats --regions "$OUTPUT_UPSTREAM_BED" --min-coverage 3 -o "${STAT_OUTPUT_upTR}" "${PILEUP_OUTPUT}.gz" 2>/dev/null
                 modkit stats --regions "$OUTPUT_DOWNSTREAM_BED" --min-coverage 3 -o "${STAT_OUTPUT_downTR}" "${PILEUP_OUTPUT}.gz" 2>/dev/null
            else
                 # Log error and continue to next TR entry
                 echo "Error: modkit pileup output is empty or missing for ${TR_ID}:${HAPLOTYPE}" >&2
                 continue
            fi

            echo "Sucessfully processed ${CHROM}:${TR_ID}:${HAPLOTYPE}!"
            echo ""
        done
    fi
done < "$INPUT_VCF"

# Third processing loop: assign average allele-specific TR methylation values to VCF

echo "STEP3: UPDATE VCF OUTPUT WITH METHYLATION DATA"
echo ""

declare -A METHYLATION_VALUES

while read -r LINE; do
    if [[ $LINE == \#* ]]; then
        continue
    fi

    MODIFIED_LINE=$(modify_gt_based_on_pdp "$LINE")

    # Parse relevant VCF fields (same als in loop 1 and 2)
    CHROM=$(echo "$MODIFIED_LINE" | cut -f1)
    TR_ID=$(echo "$MODIFIED_LINE" | cut -f3)

    # Initialize methylation values with defaults for unphased
    TR_METHYLATION=".,."
    TR_NVALID=".,."
    upTR_METHYLATION=".,."
    upTR_NVALID=".,."
    downTR_METHYLATION=".,."
    downTR_NVALID=".,."

    echo "processing ${CHROM}:${TR_ID}"

    if is_haploid "$CHROM"; then
        # Haploid chromosome: assign methylation for REF or ALT
        if [[ "$GT" == "1|." ]]; then
            HAPLOTYPE="ALT"
        elif [[ "$GT" == "0|." ]]; then
            HAPLOTYPE="REF"
        else
            echo "Skipping haploid TR $CHROM:$TR_ID (unsupported GT=$GT)"
            continue
        fi

        STAT_FILE_TR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_modkit_stats.tsv"
        STAT_FILE_upTR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_upstream_modkit_stats.tsv"
        STAT_FILE_downTR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_downstream_modkit_stats.tsv"

        # Check if haplotype's stats file exists
        if [[ -e "$STAT_FILE_TR" && -s "$STAT_FILE_TR" && -e "$STAT_FILE_upTR" && -s "$STAT_FILE_upTR" && -e "$STAT_FILE_downTR" && -s "$STAT_FILE_downTR" ]]; then
            TR_AVG_METHYLATION=$(awk 'NR==2 {if ($8 == "") print "."; else print $8}' "$STAT_FILE_TR")
            TR_COV_METH=$(awk 'NR==2 {if ($7 == "") print "."; else print $7}' "$STAT_FILE_TR")
            upTR_AVG_METHYLATION=$(awk 'NR==2 {if ($8 == "") print "."; else print $8}' "$STAT_FILE_upTR")
            upTR_COV_METH=$(awk 'NR==2 {if ($7 == "") print "."; else print $7}' "$STAT_FILE_upTR")
            downTR_AVG_METHYLATION=$(awk 'NR==2 {if ($8 == "") print "."; else print $8}' "$STAT_FILE_downTR")
            downTR_COV_METH=$(awk 'NR==2 {if ($7 == "") print "."; else print $7}' "$STAT_FILE_downTR")
            TR_METHYLATION="${TR_AVG_METHYLATION},."
            TR_NVALID="${TR_COV_METH},."
            upTR_METHYLATION="${upTR_AVG_METHYLATION},."
            upTR_NVALID="${upTR_COV_METH},."
            downTR_METHYLATION="${downTR_AVG_METHYLATION},."
            downTR_NVALID="${downTR_COV_METH},."
        else
            echo "Warning: missing or empty stats file for haploid TR $TR_ID"
        fi

    else
        # Diploid chromosomes
        HP1_TR_VALUE="."
        HP2_TR_VALUE="."
        HP1_TR_NCOV="."
        HP2_TR_NCOV="."
        HP1_upTR_VALUE="."
        HP2_upTR_VALUE="."
        HP1_upTR_NCOV="."
        HP2_upTR_NCOV="."
        HP1_downTR_VALUE="."
        HP2_downTR_VALUE="."
        HP1_downTR_NCOV="."
        HP2_downTR_NCOV="."


        for HAPLOTYPE in HP1_REF HP1_ALT1 HP1_ALT2 HP2_REF HP2_ALT1 HP2_ALT2; do
            STAT_FILE_TR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_modkit_stats.tsv"
            STAT_FILE_upTR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_upstream_modkit_stats.tsv"
            STAT_FILE_downTR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_downstream_modkit_stats.tsv"
            if [[ -e "$STAT_FILE_TR" && -s "$STAT_FILE_TR" && -e "$STAT_FILE_upTR" && -s "$STAT_FILE_upTR" && -e "$STAT_FILE_downTR" && -s "$STAT_FILE_downTR" ]]; then
                TR_AVG_METHYLATION=$(awk 'NR==2 {if ($8 == "") print "."; else print $8}' "$STAT_FILE_TR")
                TR_COV_METH=$(awk 'NR==2 {if ($7 == "") print "."; else print $7}' "$STAT_FILE_TR")
                upTR_AVG_METHYLATION=$(awk 'NR==2 {if ($8 == "") print "."; else print $8}' "$STAT_FILE_upTR")
                upTR_COV_METH=$(awk 'NR==2 {if ($7 == "") print "."; else print $7}' "$STAT_FILE_upTR")
                downTR_AVG_METHYLATION=$(awk 'NR==2 {if ($8 == "") print "."; else print $8}' "$STAT_FILE_downTR")
                downTR_COV_METH=$(awk 'NR==2 {if ($7 == "") print "."; else print $7}' "$STAT_FILE_downTR")
                if [[ "$HAPLOTYPE" == HP1* ]]; then
                    HP1_TR_VALUE="$TR_AVG_METHYLATION"
                    HP1_TR_NCOV="$TR_COV_METH"
                    HP1_upTR_VALUE="$upTR_AVG_METHYLATION"
                    HP1_upTR_NCOV="$upTR_COV_METH"
                    HP1_downTR_VALUE="$downTR_AVG_METHYLATION"
                    HP1_downTR_NCOV="$downTR_COV_METH"
                elif [[ "$HAPLOTYPE" == HP2* ]]; then
                    HP2_TR_VALUE="$TR_AVG_METHYLATION"
                    HP2_TR_NCOV="$TR_COV_METH"
                    HP2_upTR_VALUE="$upTR_AVG_METHYLATION"
                    HP2_upTR_NCOV="$upTR_COV_METH"
                    HP2_downTR_VALUE="$downTR_AVG_METHYLATION"
                    HP2_downTR_NCOV="$downTR_COV_METH"
                fi
            else
                continue
            fi
        done

        TR_METHYLATION="${HP1_TR_VALUE},${HP2_TR_VALUE}"
        TR_NVALID="${HP1_TR_NCOV},${HP2_TR_NCOV}"
        upTR_METHYLATION="${HP1_upTR_VALUE},${HP2_upTR_VALUE}"
        upTR_NVALID="${HP1_upTR_NCOV},${HP2_upTR_NCOV}"
        downTR_METHYLATION="${HP1_downTR_VALUE},${HP2_downTR_VALUE}"
        downTR_NVALID="${HP1_downTR_NCOV},${HP2_downTR_NCOV}"
    fi

    # Append methylation values to the VCF entry
    #echo -e "$MODIFIED_LINE\t$TR_ID\t$TR_METHYLATION\t$TR_NVALID\t$upTR_METHYLATION\t$upTR_NVALID\t$downTR_METHYLATION\t$downTR_NVALID" >> "$OUTPUT_VCF"

    #Store values in an associative array**
    METHYLATION_VALUES["$CHROM:$TR_ID"]="$TR_METHYLATION $TR_NVALID $upTR_METHYLATION $upTR_NVALID $downTR_METHYLATION $downTR_NVALID"

done < "$INPUT_VCF"

echo ""
echo "STEP4: RUNNING uTR AND SUMMARY OF TR SEQUENCE FEATURES"
echo ""

# Define uTR tool path
UTR_TOOL_DIR="/ifs/software/research/unique/pipeline_tools/uTR/uTR"

# Function to extract uTR features
extract_uTR_features() {
    local utr_output="$1"
    
    # Default values (if no valid data is found)
    local info_value="."
    local TR_len="."
    local pat_value="."
    local decomp_value="."

    # Extract values using awk
    if [[ -s "$utr_output" ]]; then
        info_value=$(awk -F'[()]' '/#Info / {print $2; exit}' "$utr_output" | tr -d '[:space:]')
        TR_len=$(awk -F'[()]' '/#Info / {print $2; exit}' "$utr_output" | cut -d  "," -f1) # Extract length
        pat_value=$(awk '/#Pat / {sub(/#Decomp.*/, "", $0); print substr($0, index($0, "<"))}' "$utr_output" | sort -u | paste -sd "," - | sed 's/[[:space:]]*$//')
        decomp_value=$(awk '/#Decomp / {print substr($0, index($0, "["))}' "$utr_output" | sort -u | paste -sd "," - | sed 's/[[:space:]]*$//')
    fi

    # Ensure extracted values are non-empty
    echo -e "${info_value}\t${TR_len}\t${pat_value}\t${decomp_value}"
}

# Process input file to extract motif characteristics
while read -r LINE; do
    # Skip VCF header lines
    if [[ $LINE == \#* ]]; then
        continue
    fi

    MODIFIED_LINE=$(modify_gt_based_on_pdp "$LINE")

    # Parse relevant VCF fields (same als in loop 1 and 2)
    CHROM=$(echo "$MODIFIED_LINE" | cut -f1)
    TR_ID=$(echo "$MODIFIED_LINE" | cut -f3)

    # Retrieve stored methylation values 🔥
    if [[ -n "${METHYLATION_VALUES["$CHROM:$TR_ID"]}" ]]; then
        read TR_METHYLATION TR_NVALID upTR_METHYLATION upTR_NVALID downTR_METHYLATION downTR_NVALID <<< "${METHYLATION_VALUES["$CHROM:$TR_ID"]}"
    else
        TR_METHYLATION=".,."
        TR_NVALID=".,."
        upTR_METHYLATION=".,."
        upTR_NVALID=".,."
        downTR_METHYLATION=".,."
        downTR_NVALID=".,."
        echo "Warning: No methylation values found for $CHROM:$TR_ID"
    fi

    # Log progress
    echo ""
    echo -e ".............Processing ${CHROM}:${TR_ID} - motif analysis .............."
    echo ""

    # Initiate values
    HP1_info_value="."
    HP1_TR_len="."
    HP1_pat_value="."
    HP1_decomp_value="."
    HP2_info_value="."
    HP2_TR_len="."
    HP2_pat_value="."
    HP2_decomp_value="."

    # Loop over haplotypes
    for HAPLOTYPE in HP1 HP2 unphased1 unphased2; do
    # Motif analysis with uTR
        uTR_in="${TMP_DIR}/${CHROM}_${TR_ID}_${HAPLOTYPE}.fasta"
        uTR_out="${TMP_DIR}/${CHROM}_${TR_ID}_${HAPLOTYPE}_uTR.out"

        # Check if the input FASTA file exists
        if [[ ! -f "$uTR_in" ]]; then
            continue
        fi
        
        # Run uTR tool
        echo "Running uTR for ${CHROM}_${TR_ID}_${HAPLOTYPE}..."
        "$UTR_TOOL_DIR" -f "$uTR_in" -y -o "$uTR_out" 

        # Check if uTR output file was created
        if [[ ! -s "$uTR_out" ]]; then
            continue
        fi

        # Extract and log features
        read info_value TR_len pat_value decomp_value <<< $(extract_uTR_features "$uTR_out")
 
        if [[ "$HAPLOTYPE" == HP1* ]]; then
            HP1_info_value="$info_value"
            HP1_TR_len="$TR_len"
            HP1_pat_value="$pat_value"
            HP1_decomp_value="$decomp_value"
        elif [[ "$HAPLOTYPE" == HP2* ]]; then
            HP2_info_value="$info_value"
            HP2_TR_len="$TR_len"
            HP2_pat_value="$pat_value"
            HP2_decomp_value="$decomp_value"
        else
            continue
        fi
    done

    INFO_VALUE="${HP1_info_value};${HP2_info_value}"
    TR_LEN="${HP1_TR_len}|${HP2_TR_len}"
    PAT="${HP1_pat_value},${HP2_pat_value}"
    DECOMP="${HP1_decomp_value},${HP2_decomp_value}"

    # Append motif composition values to the VCF entry
    echo -e "$MODIFIED_LINE\t$TR_ID\t$TR_METHYLATION\t$TR_NVALID\t$upTR_METHYLATION\t$upTR_NVALID\t$downTR_METHYLATION\t$downTR_NVALID\t$TR_LEN\t$INFO_VALUE\t$PAT\t$DECOMP" >> "$OUTPUT_VCF"

done < "$INPUT_VCF"

echo ""
echo "STEP5: GENERATE SUMMARY FILE"
echo ""

# Print header for summary file
echo -e "CHROM\tPOS\tID\tREF_MOTIF\tTR_REF_LENGTH\t${SAMPLE_ID}_GT\t${SAMPLE_ID}_TR_LEN\tTR_PATTERN\tTR_AM\tTR_N_METH_VALID\tUPSTREAM_TR_AM\tUPSTREAM_TR_N_METH_VALID\tDOWNSTREAM_TR_AM\tDOWNSTREAM_TR_N_METH_VALID" > "$OUTPUT_SUMMARY"

# Process each non-header line in the output VCF
awk -F'\t' '
BEGIN { OFS="\t" }
!/^#/ {
    chrom=$1;
    pos=$2;
    id=$3;
    ref=$4;
    info=$8;
    tr_am=$12;
    tr_n_meth_valid=$13;
    upstream_tr_am=$14;
    upstream_tr_n_meth_valid=$15;
    downstream_tr_am=$16;
    downstream_tr_n_meth_valid=$17;
    tr_len=$18;
    tr_pat=$20;

    # Compute REF allele length
    ref_length = length(ref);

    # Extract values from INFO field
    motif = "."; 
    if (match(info, /MOTIF=([^;]+)/, arr)) motif = arr[1];

    # Extract GT (first field of sample column)
    gt = $10;
    split(gt, gt_fields, ":");
    gt = gt_fields[1];

    # Print extracted values
    print chrom, pos, id, motif, ref_length, gt, tr_len, tr_pat, tr_am, tr_n_meth_valid, upstream_tr_am, upstream_tr_n_meth_valid, downstream_tr_am, downstream_tr_n_meth_valid;
}' "$OUTPUT_VCF" >> "$OUTPUT_SUMMARY"

echo ""
echo "Output files are in $OUTPUT_DIR. Remove temporary files stored in: $TMP_DIR"
# Uncomment the following line to automatically remove temporary files after the run
#rm -rf "$TMP_DIR"
echo ""
echo "--------------------------------"
echo "Pipeline completed successfully!"
echo "--------------------------------"
echo ""
