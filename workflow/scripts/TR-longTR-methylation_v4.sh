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

mkdir -p "$OUTPUT_DIR" "$TMP_DIR" "$ALIGNMENTS" "$METH"

#Sort vcf file
bcftools sort "$VCF_FILE" -Oz -o "$INPUT_VCF"


# Copy existing headers from the input VCF
bcftools view -h "$VCF_FILE" > "$OUTPUT_VCF"

# Add FORMAT fields for individual-specific methylation info
echo '##FORMAT=<ID=TR_AM,Number=1,Type=Float,Description="Average methylation percentage for this individual at the TR region">' >> "$OUTPUT_VCF"
echo '##FORMAT=<ID=TR_N_METH_VALID,Number=1,Type=Integer,Description="Number of valid CpG sites used for TR methylation calculation in this individual">' >> "$OUTPUT_VCF"
echo '##FORMAT=<ID=UPSTREAM_TR_AM,Number=1,Type=Float,Description="Average methylation percentage upstream of TR for this individual">' >> "$OUTPUT_VCF"
echo '##FORMAT=<ID=UPSTREAM_TR_N_METH_VALID,Number=1,Type=Integer,Description="Valid CpGs upstream of TR in this individual">' >> "$OUTPUT_VCF"
echo '##FORMAT=<ID=DOWNSTREAM_TR_AM,Number=1,Type=Float,Description="Average methylation percentage downstream of TR for this individual">' >> "$OUTPUT_VCF"
echo '##FORMAT=<ID=DOWNSTREAM_TR_N_METH_VALID,Number=1,Type=Integer,Description="Valid CpGs downstream of TR in this individual">' >> "$OUTPUT_VCF"
echo '##FORMAT=<ID=TR_LEN,Number=1,Type=String,Description="Length of TR alleles for this individual, format HP1|HP2">' >> "$OUTPUT_VCF"
echo '##FORMAT=<ID=TR_INFO_VALUE,Number=1,Type=String,Description="uTR info output per allele (HP1;HP2), format (A,B,C)">' >> "$OUTPUT_VCF"
echo '##FORMAT=<ID=TR_PAT,Number=1,Type=String,Description="TR pattern per allele, format HP1,HP2">' >> "$OUTPUT_VCF"
echo '##FORMAT=<ID=TR_DECOMP,Number=1,Type=String,Description="TR motif decomposition per allele, format HP1,HP2">' >> "$OUTPUT_VCF"

# Add VCF column headers
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$SAMPLE_ID" >> "$OUTPUT_VCF"

# Function to check if a VCF entry lacks phasing information
# Function to modify the GT field based on PDP values and handle haploid chromosomes
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


# First processing loop: read VCF file and generate TR-specific FASTA and BED entries
echo ""
echo "START ANALYSIS"
echo ""


bcftools view "$INPUT_VCF" | while read -r LINE; do
    
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

    # Parse alternate genotypes and alleles
    IFS='|' read -ra GENO_ARRAY <<< "$GENOTYPE"
    IFS=',' read -ra ALT_ALLELES <<< "$ALT"
    ALL_ALLELES=("$REF" "${ALT_ARRAY[@]}")

    # Determine repeat extension boundaries
    UPSTREAM_START=$((TR_START - EXTEND))
    UPSTREAM_END=$((TR_START - 1))
    DOWNSTREAM_START=$((TR_END + 1))
    DOWNSTREAM_END=$((TR_END + EXTEND))

    # Extract flanking sequences
    UPSTREAM_SEQ=$(samtools faidx "$REFERENCE_FASTA" "${CHROM}:${UPSTREAM_START}-${UPSTREAM_END}" | tail -n +2 | tr -d '\n' | tr '[:upper:]' '[:lower:]')
    DOWNSTREAM_SEQ=$(samtools faidx "$REFERENCE_FASTA" "${CHROM}:${DOWNSTREAM_START}-${DOWNSTREAM_END}" | tail -n +2 | tr -d '\n' | tr '[:upper:]' '[:lower:]')


    for i in "${!GENO_ARRAY[@]}"; do #iterate over indices to work with haplotypes [0 is HP1 and 1 is HP2]
	allele=${GENO_ARRAY[$i]}

        if [[ "$allele" == "." ]]; then
            TR_VALUE="."
            TR_NCOV="."
            upTR_VALUE="."
            upTR_NCOV="."
            downTR_VALUE="."
            downTR_NCOV="." 
        fi
       
        if is_haploid "$CHROM"; then
            HAPLOTYPE="haploid"  
        elif [[ "$i" -eq 0 ]]; then
            HAPLOTYPE=1
        elif [[ "$i" -eq 1 ]]; then
            HAPLOTYPE=2
        fi
        
        #Create fasta file
        SEQ="${UPSTREAM_SEQ}$(echo "${ALL_ALLELES[${i}]}" | tr '[:lower:]' '[:upper:]')${DOWNSTREAM_SEQ}"
        HEADER="${CHROM}_${POS}_${TR_ID}_${HAPLOTYPE}"
        REGION_FASTA="${HEADER}.fasta"
        TR_START_REL=$((EXTEND + 1))
        TR_END_REL=$((EXTEND + ${#REF}))
        echo ">$HEADER" >> "$REGION_FASTA"
        echo "$SEQ" >> "$REGION_FASTA"
        samtools faidx "$REGION_FASTA"
         
        # Extract reads overlapping the region and convert to FASTQ
        REGION_BAM="${TMP_DIR}/${TR_ID}_extracted.bam"
        samtools view -h "$PHASED_BAM" "$CHROM:$TR_START-$TR_END" | samtools addreplacerg -r "ID:${TR_ID}" - | samtools view -b - > "$REGION_BAM"
        
        # Filter reads by HP tag
        if [[ "$HAPLOTYPE" == 1 ]]; then
            echo "Filter HP1 reads"
            FILTERED_BAM="${TMP_DIR}/${CHROM}_${TR_ID}_extracted_HP1.bam"
            samtools view -h "$REGION_BAM" | awk '/^@/ || /HP:i:1/' | samtools view -b - > "$FILTERED_BAM"
        elif [[ "$HAPLOTYPE" == 2 ]]; then
            echo "Filter HP2 reads"
            FILTERED_BAM="${TMP_DIR}/${CHROM}_${TR_ID}_extracted_HP2.bam"
            samtools view -h "$REGION_BAM" | awk '/^@/ || /HP:i:2/' | samtools view -b - > "$FILTERED_BAM"
        else
            echo "No filtering since chromosome is haploid"
            FILTERED_BAM="$REGION_BAM"
        fi

        # Validate the bam file
        if [[ ! -s "$FILTERED_BAM" ]]; then
            echo "Error: no reads extracted for $TR_ID, $HAPLOTYPE"
            continue
        fi


        # Convert to FASTQ
        REGION_FASTQ="${TMP_DIR}/${CHROM}_${TR_ID}_${HAPLOTYPE}.extracted.fastq"
        samtools fastq -t -T MM,ML,HP,PS "$FILTERED_BAM" > "$REGION_FASTQ" 2>/dev/null

        # Map with minimap2
        echo "Remapping to TR $HAPLOTYPE sequence"
        REGION_SAM="${TMP_DIR}/${CHROM}_${TR_ID}_${HAPLOTYPE}_mapped.sam"
        minimap2 -ax map-ont -y --sam-hit-only "$REGION_FASTA" "$REGION_FASTQ" > "$REGION_SAM" 2>/dev/null

        # Convert, sort, and index BAM
        REGION_REALIGNED_BAM="${ALIGNMENTS}/${CHROM}_${TR_ID}_${HAPLOTYPE}_mapped.sorted.bam"
        samtools view -bS "$REGION_SAM" | samtools sort -o "$REGION_REALIGNED_BAM" 2>/dev/null
        samtools index "$REGION_REALIGNED_BAM"

        # Perform modkit pileup on the generated BAM file
        echo "Running modkit pileup for ${allele}"
        PILEUP_OUTPUT="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_modkit_pileup.bed"
        modkit pileup "${REGION_REALIGNED_BAM}" "${PILEUP_OUTPUT}" --cpg --ref "${REGION_FASTA}" --ignore h --combine-strands --mod-threshold m:0.8 2>/dev/null

        # Check for phased output files
        if  [[ -e "$PILEUP_OUTPUT" && -s "$PILEUP_OUTPUT" ]]; then
            echo "Modkit pileup completed successfully"
            # Compress and index the phased pileup output
            bgzip "$PILEUP_OUTPUT"
            tabix "${PILEUP_OUTPUT}.gz"

            # Perform modkit stats for the TR region for each phased output
            TR_ALL_START=$((EXTEND + 1))
            TR_ALL_END=$((EXTEND + ${#ALL_ALLELES[${i}]}))
            REGION_BED="${HEADER}.bed"
            REGION_UPSTREAM_BED="${HEADER}_upstream.bed"
            REGION_DOWNSTREAM_BED="${HEADER}_downstream.bed"
        
            echo -e "$HEADER\t$TR_ALL_START\t$TR_ALL_END\t$TR_ID" > "$REGION_BED"
            bedtools flank -i "$OUTPUT_BED" -l "$FLANKING_BASES" -r 0 > "$OUTPUT_UPSTREAM_BED"
            bedtools flank -i "$OUTPUT_BED" -l 0 -r "$FLANKING_BASES" > "$OUTPUT_DOWNSTREAM_BED"

            STAT_OUTPUT_TR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_modkit_stats.tsv"
            STAT_OUTPUT_upTR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_upstream_modkit_stats.tsv"
            STAT_OUTPUT_downTR="${METH}/${CHROM}_${TR_ID}_${HAPLOTYPE}_downstream_modkit_stats.tsv"
            echo "Running modkit stats for TR ${HAPLOTYPE}"
            modkit stats --regions "$REGION_BED" --min-coverage 3 -o "${STAT_OUTPUT_TR}" "${PILEUP_OUTPUT}.gz" 2>/dev/null
            modkit stats --regions "$OUTPUT_UPSTREAM_BED" --min-coverage 3 -o "${STAT_OUTPUT_upTR}" "${PILEUP_OUTPUT}.gz" 2>/dev/null
            modkit stats --regions "$OUTPUT_DOWNSTREAM_BED" --min-coverage 3 -o "${STAT_OUTPUT_downTR}" "${PILEUP_OUTPUT}.gz" 2>/dev/null
            else
                 # Log error and continue to next TR entry
                 echo "Error: modkit pileup output is empty or missing for ${TR_ID}:${HAPLOTYPE}" >&2
                 continue
            fi


            # Cleanup temporary files
            rm -f "$REGION_FASTQ"
            rm -f "$REGION_SAM"
            rm -f "$REGION_BAM"
            rm -f "$REGION_BED"
            rm -f "$REGION_BAM"
            rm -f "${PILEUP_OUTPUT}.gz" "${PILEUP_OUTPUT}.gz.tbi"


    done
done

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
    GT=$(echo "$MODIFIED_LINE" | cut -f10 | cut -d':' -f1)

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
rm -rf "$TMP_DIR"
echo ""
echo "--------------------------------"
echo "Pipeline completed successfully!"
echo "--------------------------------"
echo ""
