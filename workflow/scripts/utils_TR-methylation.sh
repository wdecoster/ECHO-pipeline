#!/bin/bash

#++++ Define functions ++++
#Clean up when exit script
cleanup() {
    echo "Performing cleanup..."
    rm -rf "${TMP_DIR}" "${METH}" "${ALIGNMENTS}" "${LOGS}"
    rm -f "$INPUT_VCF" "$CORRECT_HEADER_VCF_FILE" "$OUTPUT_VCF" "$OUTPUT_VCF.gz"
}


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