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
                -t <TE_type> -c <TE_catalog> -s <sample_id> -o <output_dir> -x <threads ==n/2> [-f <flank_bp>] "
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

# ------------ (1) Generate shifted bed files ------------
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

# ------------ (2) Filter VCF ------------
# Now filtering for at least 5 covered reads used for variant calling

#++ Time measuring ++#
start_t=$(date +%s.%N)
# is there actually a benefit of --threads ??
echo "Thread count: $((TOTAL_THREADS/2))"
bcftools view "$SNP_data" -i 'FILTER="PASS" & FORMAT/DP>5' --threads "$((TOTAL_THREADS/2))" --output-file "$SNP_filt" & 
VCF_SNP_PID=$!
bcftools view "$SV_data" -i 'FILTER="PASS" & INFO/SUPPORT>5' --threads "$((TOTAL_THREADS/2))" --output-file "$SV_filt" &
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


# ------------ (3) Modkit stats Operationen ------------
echo "-- (3) -- : Running modkit stats operations..."
#++ Time measuring ++#
start_t=$(date +%M)

# reassigning max threads manually for parallelization
max_parallel_jobs=4 #9 for retroposon
# goal is to keep always 2 jobs in parallel with max possible threads for each
threads_per_job=$((TOTAL_THREADS / 4)) #8 for retroposon
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
##################################
elapsed=$(echo "$end_t - $start_t" | bc)

echo "time for (3) modkit stats: $elapsed seconds"
echo "-- (3) -- completed."


# ------------ (4) Variant Intersections ------------
echo "-- (4) -- : Running variant intersections..."
#++ Time measuring ++#
start_t=$(date +%M)

# Report SNPs and SVs that intersect with TE elements of interest
# define paths and max parallel jobs
SNP_INTERSECT="${OUTDIR}/variants/${SAMPLE_ID}_${TE}_SNPs_intersect.bed"
SV_INTERSECT="${OUTDIR}/variants/${SAMPLE_ID}_${TE}_SV_intersect.bed"
SNP_INTERSECT_COUNT="${OUTDIR}/variants/${SAMPLE_ID}_${TE}_SNPs_intersect_count.bed"
SV_INTERSECT_COUNT="${OUTDIR}/variants/${SAMPLE_ID}_${TE}_SV_intersect_count.bed"

# execute them all parrallel in background
max_parallel_jobs=4 #9 for retroposon
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

end_t=$(date +%M)
elapsed=$(echo "$end_t - $start_t" | bc)
echo "time for (4) intersection varaint : $elapsed minutes"
echo "-- (4) -- completed."

echo "All operations completed successfully!"

# report about the generated files
echo ""
echo "Generated files:"
find "${OUTDIR}" -name "${SAMPLE_ID}_${TE}_*" -type f | wc -l
echo "Generated files list: "
find "${OUTDIR}" -name "${SAMPLE_ID}_${TE}_*" -type f | sort

# ------------END OF FIRST HALF OF THE SCRIPT------------ #
echo "END of first half of the script, continuing with python while loop script..."
