#!/bin/bash
# -----------------------------------------------------------------------------
# Script Name:    summarise_ref_TE_avg_meth.sh
# Description:    Aggregate modification statistics. 
# Author:         Leena Putzeys, Brando Poggiali, Nicole Flack
# Date Created:   2025-11-12
# Last Modified:  2025-11-19
# Version:        1.0.0
# License:        MIT
# Dependencies:   [modkit, bgzip, tabix, bedtools, bcftools, awk]
# -----------------------------------------------------------------------------

# for debugging in bash pipeline, x for logging executed commands
set -euo pipefail
#set -x
set -m # for jobmanagement and listing within the bash shell

# Usage and help function
usage() {
    echo "Usage: $0 -c <TE_catalog> -s <sample_id> -o <working_dir> -x <threads> "
    exit 1
}

# Check if no arguments were provided
if [ $# -eq 0 ]; then
    usage
fi
# check if help argument is provided
if [[ "${1:-}" == "--help" ]]; then
    usage
fi
# ------------ Required tools ------------
REQUIRED_TOOLS=(bedtools awk)
for tool in "${REQUIRED_TOOLS[@]}"; do
    if ! command -v "$tool" &>/dev/null; then
        echo "Error: Required tool '$tool' not found in PATH."
        exit 1
    fi
done

# ------------ Default Values ------------
TE_CATALOG=""
SAMPLE_ID=""
OUTDIR=""

# ------------ Parse Args ------------ 
# assings input arguments of each --option to the parameter.
while getopts "p:u:v:t:c:s:o:x:f:" opt; do
    case $opt in
        c) TE_CATALOG="$OPTARG" ;;
        s) SAMPLE_ID="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        *) usage ;;
    esac
done

# ------------ Check Inputs ------------
## logical operator string matching: -z True if the length of string = zero
if [[ -z "$TE_CATALOG" || -z "$SAMPLE_ID" || -z "$OUTDIR" ]]; then
    usage
fi

# -------- Define inputs ------------
#Set up input files
# preserving original scrip definitions
SNP_INTERSECT="${OUTDIR}/variants/${SAMPLE_ID}_SNPs_intersect.bed"
SV_INTERSECT="${OUTDIR}/variants/${SAMPLE_ID}_SV_intersect.bed"
SNP_INTERSECT_COUNT="${OUTDIR}/variants/${SAMPLE_ID}_SNPs_intersect_count.bed"
SV_INTERSECT_COUNT="${OUTDIR}/variants/${SAMPLE_ID}_SV_intersect_count.bed"



# ------------ Previous step (5) Create summary files -----------------------------------------------------------
echo "-- (5) -- : Create summary file..."
#++ Time measuring ++#
start_t=$(date +%s)

#- (4.1) Intersect TE catalog with intersected SNV, SV count files and SV codes -
(
  paste "${TE_CATALOG}" "${SNP_INTERSECT_COUNT}" | awk 'BEGIN { OFS="\t" } $2 == $9 && $3 == $10 {print $1, $2, $3, $4, $5, $6, $7, $15}' > "${OUTDIR}/TE_cat_SNP_count.bed"
  
  paste "${OUTDIR}/TE_cat_SNP_count.bed" "${SV_INTERSECT_COUNT}" | awk 'BEGIN { OFS="\t" } $2 == $10 && $3 == $11 { print $1, $2, $3, $4, $5, $6, $7, $8, $16 }' > "${OUTDIR}/TE_cat_SNP_SV_count.bed"
  
  awk '{
    key = $1 FS $2 FS $3;
    chr[$1,$2,$3] = $1;
    start[$1,$2,$3] = $2;
    end[$1,$2,$3] = $3;

    id = $10;
    if (id != "" && id != ".") {
        sv[key] = (key in sv ? sv[key] "," id : id);
    }
  } END {
    for (k in sv) {
      split(k, fields, FS);
      print fields[1], fields[2], fields[3], sv[k];
    }
  }' OFS="\t" "${SV_INTERSECT}" > "${OUTDIR}/interesected_SV_collapsed.bed"
  
  sort -k1,1V -k2,2n -k3,3n "${OUTDIR}/interesected_SV_collapsed.bed" > "${OUTDIR}/interesected_SV_collapsed_sorted.bed"
  
  bedtools intersect -a "${OUTDIR}/TE_cat_SNP_SV_count.bed" -b "${OUTDIR}/interesected_SV_collapsed_sorted.bed"  -wao -f 1.0 -r | cut -f1-9,13 > "${OUTDIR}/TE_cat_SNP_SV_count_SV_code.bed"
  
  rm "${OUTDIR}/TE_cat_SNP_count.bed" "${OUTDIR}/TE_cat_SNP_SV_count.bed" "${OUTDIR}/interesected_SV_collapsed.bed" "${OUTDIR}/interesected_SV_collapsed_sorted.bed"
) &
VARIANT_JOB=$!

#- (4.2) Merge unphased methylation stats from TE, and upstream/downstrem regions -
(

paste "${OUTDIR}/mod_unphased/${SAMPLE_ID}_full_stats_unphased.tsv" "${OUTDIR}/mod_unphased/${SAMPLE_ID}_upstream_stats_unphased.tsv" "${OUTDIR}/mod_unphased/${SAMPLE_ID}_downstream_stats_unphased.tsv" | awk 'BEGIN { OFS="\t" } $2 == $11 && $3 == $18 { print $1, $2, $3, $4, $5, $7, $8, $15, $16, $23, $24 }' > "${OUTDIR}/unphased_stats_TE_and_flanking.bed"

) &
UNPHASED_JOB=$!


#- (4.3) Merge phased methylation stats from TE, and upstream/downstrem regions -
(
echo "Empty flanking" > "${OUTDIR}/phased_stats_upstream.bed"

paste "${OUTDIR}/mod_phased/${SAMPLE_ID}_full_stats_1.tsv" "${OUTDIR}/mod_phased/${SAMPLE_ID}_full_stats_2.tsv" | awk 'BEGIN { OFS="\t" } $2 == $10 && $3 == $11 { print $1, $2, $3, $4, $5, $7 "," $15, $8 "," $16 }' > "${OUTDIR}/phased_stats.bed"

paste "${OUTDIR}/phased_stats.bed" "${OUTDIR}/mod_phased/${SAMPLE_ID}_upstream_stats_1.tsv" "${OUTDIR}/mod_phased/${SAMPLE_ID}_upstream_stats_2.tsv" | awk 'BEGIN { OFS="\t" } $2 == $10 && $2 == $18 { print $1, $2, $3, $4, $5, $6, $7, $14 "," $22, $15 "," $23 }' >> "${OUTDIR}/phased_stats_upstream.bed" 

paste "${OUTDIR}/phased_stats_upstream.bed" "${OUTDIR}/mod_phased/${SAMPLE_ID}_downstream_stats_1.tsv" "${OUTDIR}/mod_phased/${SAMPLE_ID}_downstream_stats_2.tsv" | awk 'BEGIN { OFS="\t" } $3 == $11 && $3 == $19 { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $16 "," $24, $17 "," $25 }' > "${OUTDIR}/phased_stats_TE_and_flanking.bed"

rm "${OUTDIR}/phased_stats.bed" "${OUTDIR}/phased_stats_upstream.bed"
) &
PHASED_JOB=$!

wait ${VARIANT_JOB} ${UNPHASED_JOB} ${PHASED_JOB}

#Create summary file merging variant TE file with unphased methylation stat file and then the phased one.
bedtools intersect -a "${OUTDIR}/TE_cat_SNP_SV_count_SV_code.bed" -b "${OUTDIR}/unphased_stats_TE_and_flanking.bed" -wao -f 1.0 -r | awk 'BEGIN { OFS="\t" } { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $16, $17,$18, $19,$20, $21 }' > "${OUTDIR}/summary_variants_unphased.tsv"

#Print header summary file
echo -e "#chr\tstart\tend\tfamily\t.\tstrand\tID\tTE_length\tTE_avgMeth_phased\tTE_Nvalid_phased\tupstream_avgMeth_phased\tupstream_Nvalid_phased\tdownstream_avgMeth_phased\tdownstream_Nvalid_phased\tTE_avgMeth_unphased\tTE_Nvalid_unphased\tupstream_avgMeth_unphased\tupstream_Nvalid_unphased\tdownstream_avgMeth_unphased\tdownstream_Nvalid_unphased\ttotal_SNV\tSV_count\tSV_IDs" > "${OUTDIR}/${SAMPLE_ID}_methylation_summary.bed"

#
bedtools intersect -a "${OUTDIR}/summary_variants_unphased.tsv" -b "${OUTDIR}/phased_stats_TE_and_flanking.bed"  -wao -f 1.0 -r | awk 'BEGIN { OFS="\t" } { print $1, $2, $3, $4, $5, $6, $7, $3 - $2, $23, $22, $25, $24, $27, $26, $12, $11, $14, $13, $16, $15, $8, $9, $10 }' >> "${OUTDIR}/${SAMPLE_ID}_methylation_summary.bed"

#Remove intermediary files
rm "${OUTDIR}/phased_stats_TE_and_flanking.bed" "${OUTDIR}/summary_variants_unphased.tsv" "${OUTDIR}/unphased_stats_TE_and_flanking.bed" "${OUTDIR}/TE_cat_SNP_SV_count_SV_code.bed"

end_t=$(date +%s)
elapsed=$(( (end_t - start_t) / 60 ))

echo "time for (5) Summary creation : $elapsed minutes"
echo "-- (5) -- completed."
