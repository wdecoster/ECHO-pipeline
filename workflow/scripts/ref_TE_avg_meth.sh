#!/bin/bash
# -----------------------------------------------------------------------------
# Script Name:    ref_TE_cpg_res.sh
# Description:    Intersecting bed files containing methylation information of
#                 all cpg sites with TE catalogs.
# Author:         Leena Putzeys, Brando Poggiali
# Date Created:   2025-03-02
# Last Modified:  2025-07-08
# Version:        2.0.0
# License:        MIT
# Dependencies:   [modkit, bgzip, tabix, bedtools, bcftools, awk]
# -----------------------------------------------------------------------------


# for debugging in bash pipeline, x for logging executed commands
set -euo pipefail
#set -x
set -m # for jobmanagement and listing within the bash shell

# Usage and help function
usage() {
    echo "Usage: $0 -p <phased_dir> -u <unphased_dir> -v <phased_variation_dir> \\
                -t <TE_type> -c <TE_catalog> -s <sample_id> -o <output_dir> -x <threads> [-f <flank_bp>] "
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
REQUIRED_TOOLS=(modkit bgzip tabix bcftools bedtools awk)
for tool in "${REQUIRED_TOOLS[@]}"; do
    if ! command -v "$tool" &>/dev/null; then
        echo "Error: Required tool '$tool' not found in PATH."
        exit 1
    fi
done


# =========================================== 
# First half of Script summarizes methylation and variation 
# =========================================== 

# ------------ Default Values ------------
PHASED_PILEUP=""
UNPHASED_PILEUP=""
VARIATION=""
TE=""
TE_CATALOG=""
SAMPLE_ID=""
OUTDIR=""
TOTAL_THREADS=2
FLANK=250

# ------------ Parse Args ------------ 
# assings input arguments of each --option to the parameter.
while getopts "p:u:v:t:c:s:o:x:f:" opt; do
    case $opt in
        p) PHASED_PILEUP="$OPTARG" ;;
        u) UNPHASED_PILEUP="$OPTARG" ;;
        v) VARIATION="$OPTARG" ;;
        t) TE="$OPTARG" ;; #Options: all, all_cpg, DNA, DNA_cpg, LINE, LINE_cpg, helitron, helitron_cpg, SINE, SINE_cpg, LTR, LTR_cpg, retroposon, retroposon_cpg
        c) TE_CATALOG="$OPTARG" ;;
        s) SAMPLE_ID="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        x) TOTAL_THREADS="$OPTARG" ;; # must be even ! (n/2)
        f) FLANK="$OPTARG" ;; # optional
        *) usage ;;
    esac
done

# ------------ Check Inputs ------------
## logical operator string matching: -z True if the length of string = zero
if [[ -z "$PHASED_PILEUP" || -z "$UNPHASED_PILEUP" || -z "$TE" || -z "$TE_CATALOG" || -z "$SAMPLE_ID" || -z "$OUTDIR" ]]; then
    usage
fi

# create directory structure 
mkdir -p "$OUTDIR/mod_phased" "$OUTDIR/mod_unphased" "$OUTDIR/variants"

# -------- Define inputs ------------
#Set up input files
phased_pileup_1="${PHASED_PILEUP}/${SAMPLE_ID}_haplotype_1.bed.gz"
phased_pileup_2="${PHASED_PILEUP}/${SAMPLE_ID}_haplotype_2.bed.gz"
unphased_pileup="${UNPHASED_PILEUP}/${SAMPLE_ID}_unphased.bed.gz"
SNP_data="${VARIATION}/${SAMPLE_ID}_phased.vcf.gz"
SV_data="${VARIATION}/${SAMPLE_ID}_phased_SV.vcf.gz"
SNP_filt="${VARIATION}/${SAMPLE_ID}_SNP_filt.vcf"
SV_filt="${VARIATION}/${SAMPLE_ID}_SV_filt.vcf"

# ------------ (1) Generate shifted bed files --------------------------------------------------------------
# Generate shifted bed file for upstream TE regions (default -250 bp, or number of bases specified by flank variable)

UPSTREAM_TE_CATALOG="${OUTDIR}/${TE}_upstream.bed"
DOWNSTREAM_TE_CATALOG="${OUTDIR}/${TE}_downstream.bed"
# coarse-grained" parallelization for two separate output files. with the background parameter '&'
#++ Time measuring ++#
start_t=$(date +%s.%N)
awk -v flank="$FLANK" '{OFS="\t"} {start=$2-flank; if (start < 0) start=0; print $1, start, $2, $7}' "$TE_CATALOG" > "$UPSTREAM_TE_CATALOG" & 
UPstream_PID=$!
awk -v flank="$FLANK" '{OFS="\t"} {print $1, $3, $3+flank, $7}' "$TE_CATALOG" > "$DOWNSTREAM_TE_CATALOG" &
DOWNstream_PID=$!
# wait $UPstream_PID $DOWNstream_PID

#++ time ++#
end_t=$(date +%s.%N)
elapsed=$(echo "$end_t - $start_t" | bc)
echo "Executed time for (1) generating shifted bed files: $elapsed seconds"

# ------------ (2) Filter VCF ------------------------------------------------------------------------------------
# Now filtering for at least 5 covered reads used for variant calling

#++ Time measuring ++#
start_t=$(date +%s.%N)
# is there actually a benefit of --threads ??
echo "Thread count: $((TOTAL_THREADS))"
bcftools view "$SNP_data" -i 'FILTER="PASS" & FORMAT/DP>5' --threads $((TOTAL_THREADS/4)) --output-file "$SNP_filt" & 
VCF_SNP_PID=$!
bcftools view "$SV_data" -i 'FILTER="PASS" & INFO/SUPPORT>5' --threads $((TOTAL_THREADS/4)) --output-file "$SV_filt" &
VCF_SV_PID=$!
wait $UPstream_PID $DOWNstream_PID $VCF_SNP_PID $VCF_SV_PID


# Check and process SNP_filt, skip if already exists
if [[ ! -f "$SNP_filt.gz" ]]; then
    bgzip -@ 8 "$SNP_filt" && tabix "$SNP_filt.gz" &
    bg1_PID=$!
else
    echo "SNP_filt already compressed, skipping bgzip"
    if [[ ! -f "$SNP_filt.gz.tbi" ]]; then
        tabix "$SNP_filt.gz" &
        bg1_PID=$!
    else
        bg1_PID=""
    fi
fi

# Check and process SV_filt, skip if already exists
if [[ ! -f "$SV_filt.gz" ]]; then
    bgzip -@ 8 "$SV_filt" && tabix "$SV_filt.gz" &
    bg2_PID=$!
else
    echo "SV_filt already compressed, skipping bgzip"
    if [[ ! -f "$SV_filt.gz.tbi" ]]; then
        tabix "$SV_filt.gz" &
        bg2_PID=$!
    else
        bg2_PID=""
    fi
fi

# Only wait if background processes were started
[[ -n "$bg1_PID" ]] || [[ -n "$bg2_PID" ]] && wait $bg1_PID $bg2_PID

end_t=$(date +%s.%N)
elapsed=$(echo "$end_t - $start_t" | bc)
echo "Executed time for (2) filtering vcf files : $elapsed seconds"


# ------------ (3) Modkit stats analysis -----------------------------------------------------------
echo "-- (3) -- : Running modkit stats operations..."
#++ Time measuring ++#
start_t=$(date +%s)

# reassigning max threads manually for parallelization
max_parallel_jobs=2 #9 for retroposon
# goal is to keep always 2 jobs in parallel with max possible threads for each
threads_per_job=$((TOTAL_THREADS / 2)) #8 for retroposon
# making sure, at least one job is running in parallel
[[ $max_parallel_jobs -lt 1 ]] && max_parallel_jobs=1

echo "-- Starting modkit stats calculation with $threads_per_job threads per jobs and $max_parallel_jobs maximal paralel executed jobs --"


rm -f "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_stats_2.tsv" "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_stats_1.tsv" "${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_stats_unphased.tsv" "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_upstream_stats_1.tsv" "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_upstream_stats_2.tsv" "${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_upstream_stats_unphased.tsv" "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_downstream_stats_1.tsv" "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_downstream_stats_2.tsv" "${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_downstream_stats_unphased.tsv"


# including job array struccture
cmds=(
  "modkit stats -t $threads_per_job --regions \"$TE_CATALOG\"           -c m -o \"${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_stats_1.tsv\"        \"$phased_pileup_1\""
  "modkit stats -t $threads_per_job --regions \"$TE_CATALOG\"           -c m -o \"${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_stats_2.tsv\"        \"$phased_pileup_2\""
  "modkit stats -t $threads_per_job --regions \"$TE_CATALOG\"           -c m -o \"${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_stats_unphased.tsv\" \"$unphased_pileup\""
  "modkit stats -t $threads_per_job --regions \"$UPSTREAM_TE_CATALOG\"  -c m -o \"${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_upstream_stats_1.tsv\"  \"$phased_pileup_1\""
  "modkit stats -t $threads_per_job --regions \"$UPSTREAM_TE_CATALOG\"  -c m -o \"${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_upstream_stats_2.tsv\"  \"$phased_pileup_2\""
  "modkit stats -t $threads_per_job --regions \"$UPSTREAM_TE_CATALOG\"  -c m -o \"${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_upstream_stats_unphased.tsv\" \"$unphased_pileup\""
  "modkit stats -t $threads_per_job --regions \"$DOWNSTREAM_TE_CATALOG\" -c m -o \"${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_downstream_stats_1.tsv\" \"$phased_pileup_1\""
  "modkit stats -t $threads_per_job --regions \"$DOWNSTREAM_TE_CATALOG\" -c m -o \"${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_downstream_stats_2.tsv\" \"$phased_pileup_2\""
  "modkit stats -t $threads_per_job --regions \"$DOWNSTREAM_TE_CATALOG\" -c m -o \"${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_downstream_stats_unphased.tsv\" \"$unphased_pileup\""
)

# starting to loop thirough commands
for cmd in "${cmds[@]}"; do
    # counting all jobs in current shell
    while (( $(jobs | wc -l) >= max_parallel_jobs )); do
            # sleep for 1 second, if max number of jobs is reached (continue if a jobs slot is free)
            sleep 1
    done
    # start command in background
    {
        # echo "START: $cmd"
        eval "$cmd"
    } &
done

# wait for all instances
# wait
end_t=$(date +%s)
elapsed=$(( (end_t - start_t) / 60 ))

echo "time for (3) modkit stats: $elapsed minutes"
echo "-- (3) -- completed."


# ------------ (4) Variant Intersections ----------------------------------------------------------------
echo "-- (4) -- : Running variant intersections..."
#++ Time measuring ++#

# Report SNPs and SVs that intersect with TE elements of interest
# define paths and max parallel jobs
SNP_INTERSECT="${OUTDIR}/variants/${SAMPLE_ID}_${TE}_SNPs_intersect.bed"
SV_INTERSECT="${OUTDIR}/variants/${SAMPLE_ID}_${TE}_SV_intersect.bed"
SNP_INTERSECT_COUNT="${OUTDIR}/variants/${SAMPLE_ID}_${TE}_SNPs_intersect_count.bed"
SV_INTERSECT_COUNT="${OUTDIR}/variants/${SAMPLE_ID}_${TE}_SV_intersect_count.bed"

# execute them all parrallel in background
max_parallel_jobs=3 #9 for retroposon
# Command array for variant intersections
variant_intersect_cmds=(
    "bedtools intersect -a \"${TE_CATALOG}\" -b \"${SNP_filt}.gz\" -wa -wb > \"${SNP_INTERSECT}\""
    "bedtools intersect -a \"${TE_CATALOG}\" -b \"${SV_filt}.gz\" -wa -wb > \"${SV_INTERSECT}\""
    "bedtools intersect -a \"${TE_CATALOG}\" -b \"${SNP_filt}.gz\" -wa -wb -C > \"${SNP_INTERSECT_COUNT}\""
    "bedtools intersect -a \"${TE_CATALOG}\" -b \"${SV_filt}.gz\" -wa -wb -C > \"${SV_INTERSECT_COUNT}\""
)

# starting to loop thourgh commands
for v_cmd in "${variant_intersect_cmds[@]}"; do
    # counting all jobs in current shell
    while (( $(jobs | wc -l) >= max_parallel_jobs )); do
            sleep 1
    done
    # start command in background
    {
        eval "$v_cmd"
    } &
done
# wait on intersect jobs to finish
wait

end_t=$(date +%s)
elapsed=$(( (end_t - start_t) / 60 ))

echo "time for (3) + (4) intersection variant : $elapsed minutes"
echo "-- (4) -- completed."


# ------------ (5) Create summary files -----------------------------------------------------------
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
    sv[key] = (key in sv ? sv[key] "," $10 : $10);
  } END {
    for (k in sv) {
      split(k, fields, FS);
      print fields[1], fields[2], fields[3], sv[k];
    }
  }' OFS="\t" "${SV_INTERSECT}" > "${OUTDIR}/interesected_SV_collapsed.bed"
  
  sort -k1,1 -k2,2n -k3,3n "${OUTDIR}/interesected_SV_collapsed.bed" > "${OUTDIR}/interesected_SV_collapsed_sorted.bed"
  
  bedtools intersect -a "${OUTDIR}/TE_cat_SNP_SV_count.bed" -b "${OUTDIR}/interesected_SV_collapsed_sorted.bed"  -wao -f 1.0 -r | cut -f1-9,13 > "${OUTDIR}/TE_cat_SNP_SV_count_SV_code.bed"
  
  rm "${OUTDIR}/TE_cat_SNP_count.bed" "${OUTDIR}/TE_cat_SNP_SV_count.bed" "${OUTDIR}/interesected_SV_collapsed.bed" "${OUTDIR}/interesected_SV_collapsed_sorted.bed"
) &
VARIANT_JOB=$!

#- (4.2) Merge unphased methylation stats from TE, and upstream/downstrem regions -
(

paste "${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_stats_unphased.tsv" "${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_upstream_stats_unphased.tsv" "${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_downstream_stats_unphased.tsv" | awk 'BEGIN { OFS="\t" } $2 == $11 && $3 == $18 { print $1, $2, $3, $4, $5, $7, $8, $15, $16, $23, $24 }' > "${OUTDIR}/unphased_stats_TE_and_flanking.bed"

) &
UNPHASED_JOB=$!


#- (4.3) Merge phased methylation stats from TE, and upstream/downstrem regions -
(
echo "Empty flanking" > "${OUTDIR}/phased_stats_upstream.bed"

paste "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_stats_1.tsv" "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_stats_2.tsv" | awk 'BEGIN { OFS="\t" } $2 == $10 && $3 == $11 { print $1, $2, $3, $4, $5, $7 "," $15, $8 "," $16 }' > "${OUTDIR}/phased_stats.bed"

paste "${OUTDIR}/phased_stats.bed" "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_upstream_stats_1.tsv" "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_upstream_stats_2.tsv" | awk 'BEGIN { OFS="\t" } $2 == $10 && $2 == $18 { print $1, $2, $3, $4, $5, $6, $7, $14 "," $22, $15 "," $23 }' >> "${OUTDIR}/phased_stats_upstream.bed" 

paste "${OUTDIR}/phased_stats_upstream.bed" "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_downstream_stats_1.tsv" "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_downstream_stats_2.tsv" | awk 'BEGIN { OFS="\t" } $3 == $11 && $3 == $19 { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $16 "," $24, $17 "," $25 }' > "${OUTDIR}/phased_stats_TE_and_flanking.bed"

rm "${OUTDIR}/phased_stats.bed" "${OUTDIR}/phased_stats_upstream.bed"
) &
PHASED_JOB=$!

wait ${VARIANT_JOB} ${UNPHASED_JOB} ${PHASED_JOB}

#Create summary file merging variant TE file with unphased methylation stat file and then the phased one.
bedtools intersect -a "${OUTDIR}/TE_cat_SNP_SV_count_SV_code.bed" -b "${OUTDIR}/unphased_stats_TE_and_flanking.bed" -wao -f 1.0 -r | awk 'BEGIN { OFS="\t" } { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $16, $17,$18, $19,$20, $21 }' > "${OUTDIR}/summary_variants_unphased.tsv"

#Print header summary file
echo -e "#chr\tstart\tend\tfamily\t.\tstrand\tID\tTE_length\tTE_avgMeth_phased\tTE_Nvalid_phased\tupstream_avgMeth_phased\tupstream_Nvalid_phased\tdownstream_avgMeth_phased\tdownstream_Nvalid_phased\tTE_avgMeth_unphased\tTE_Nvalid_unphased\tupstream_avgMeth_unphased\tupstream_Nvalid_unphased\tdownstream_avgMeth_unphased\tdownstream_Nvalid_unphased\ttotal_SNP\tSNP_count\tINDEL_count\tSV_count\tSV_types\tSV_IDs" > "${OUTDIR}/${SAMPLE_ID}_${TE}_methylation_summary.bed"

#
bedtools intersect -a "${OUTDIR}/summary_variants_unphased.tsv" -b "${OUTDIR}/phased_stats_TE_and_flanking.bed"  -wao -f 1.0 -r | awk 'BEGIN { OFS="\t" } { print $1, $2, $3, $4, $5, $6, $7, $3 - $2, $23, $22, $25, $24, $27, $26, $12, $11, $14, $13, $16, $15, $8, $9, $10 }' >> "${OUTDIR}/${SAMPLE_ID}_${TE}_methylation_summary.bed"

#Remove intermediary files
rm "${OUTDIR}/phased_stats_TE_and_flanking.bed" "${OUTDIR}/summary_variants_unphased.tsv" "${OUTDIR}/unphased_stats_TE_and_flanking.bed"

end_t=$(date +%s)
elapsed=$(( (end_t - start_t) / 60 ))

echo "time for (5) Summary creation : $elapsed minutes"
echo "-- (5) -- completed."


