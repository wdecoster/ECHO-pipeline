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
# Dependencies:   [modkit, bgzip, tabix, bedtools, awk]
# -----------------------------------------------------------------------------

# for debugging in bash pipeline, x for logging executed commands
set -euo pipefail

# Usage and help function
usage() {
    echo "Usage: $0 -p <phased_dir> -u <unphased_dir> -t <TE_type> -c <TE_catalog>  \\
                -s <sample_id> -o <output_dir> [-f <flank_bp>] "
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

# ------------ Default Values ------------
PHASED_PILEUP=""
UNPHASED_PILEUP=""
TE=""
TE_CATALOG=""
SAMPLE_ID=""
OUTDIR=""
FLANK=250

# ------------ Parse Args ------------ 
# assings input arguments of each --option to the parameter.
while getopts "p:u:t:c:s:o:f:" opt; do
    case $opt in
        p) PHASED_PILEUP="$OPTARG" ;;
        u) UNPHASED_PILEUP="$OPTARG" ;;
        t) TE="$OPTARG" ;; #Options: all, all_cpg, DNA, DNA_cpg, LINE, LINE_cpg, helitron, helitron_cpg, SINE, SINE_cpg, LTR, LTR_cpg, retroposon, retroposon_cpg
        c) TE_CATALOG="$OPTARG" ;;
        s) SAMPLE_ID="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
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
mkdir -p "$OUTDIR/mod_phased" "$OUTDIR/mod_unphased"

# -------- Define inputs ------------
#Set up input files
phased_pileup_1="${PHASED_PILEUP}/${SAMPLE_ID}_haplotype_1.bed.gz"
phased_pileup_2="${PHASED_PILEUP}/${SAMPLE_ID}_haplotype_2.bed.gz"
unphased_pileup="${UNPHASED_PILEUP}/${SAMPLE_ID}_unphased.bed.gz"

# ------------ (1) Generate shifted bed files ------------
# Generate shifted bed file for upstream TE regions (default -250 bp, or number of bases specified by flank variable)
mkdir -p "${OUTDIR}/TMP"
UPSTREAM_TE_CATALOG="${OUTDIR}/TMP/${TE}_upstream.bed"
DOWNSTREAM_TE_CATALOG="${OUTDIR}/TMP/${TE}_downstream.bed"
# coarse-grained" parallelization for two separate output files. with the background parameter '&'
awk -v flank="$FLANK" '{OFS="\t"} {start=$2-flank; if (start < 0) start=0; print $1, start, $2, $7}' "$TE_CATALOG" > "$UPSTREAM_TE_CATALOG" & 
UPstream_PID=$!
awk -v flank="$FLANK" '{OFS="\t"} {print $1, $3, $3+flank, $7}' "$TE_CATALOG" > "$DOWNSTREAM_TE_CATALOG" &
DOWNstream_PID=$!
wait $UPstream_PID $DOWNstream_PID

echo "Upstream and Downstream bed files created"
echo
# ------------ (2) Bedtools intersect Operationen ------------ #
# Intersect modkit pileup data on (upstream/downstream) TE catalog to obtain methylation of individual CpGs that are embedded in TE elements or their flanking regions.

# bash-command arrays for launching bedtools intersect 
bedtools_intersect_cmds=(
    "zcat \"$phased_pileup_1\" | bedtools intersect -a - -b \"$TE_CATALOG\" -wa -wb > \"${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_pileup_1.bed\""
    "zcat \"$phased_pileup_2\" | bedtools intersect -a - -b \"$TE_CATALOG\" -wa -wb > \"${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_pileup_2.bed\""
    "zcat \"$unphased_pileup\" | bedtools intersect -a - -b \"$TE_CATALOG\" -wa -wb > \"${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_pileup_unphased.bed\""
    "zcat \"$phased_pileup_1\" | bedtools intersect -a - -b \"$UPSTREAM_TE_CATALOG\" -wa -wb > \"${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_upstream_pileup_1.bed\""
    "zcat \"$phased_pileup_2\" | bedtools intersect -a - -b \"$UPSTREAM_TE_CATALOG\" -wa -wb > \"${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_upstream_pileup_2.bed\""
    "zcat \"$unphased_pileup\" | bedtools intersect -a - -b \"$UPSTREAM_TE_CATALOG\" -wa -wb > \"${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_upstream_pileup_unphased.bed\""
    "zcat \"$phased_pileup_1\" | bedtools intersect -a - -b \"$DOWNSTREAM_TE_CATALOG\" -wa -wb > \"${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_downstream_pileup_1.bed\""
    "zcat \"$phased_pileup_2\" | bedtools intersect -a - -b \"$DOWNSTREAM_TE_CATALOG\" -wa -wb > \"${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_downstream_pileup_2.bed\""
    "zcat \"$unphased_pileup\" | bedtools intersect -a - -b \"$DOWNSTREAM_TE_CATALOG\" -wa -wb > \"${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_downstream_pileup_unphased.bed\""
)


echo "Running bedtools intersect operations..."

# 9 == Number of all commands, process them all in parallel
max_parallel_jobs=9

# starting to loop thourgh commands
for b_cmd in "${bedtools_intersect_cmds[@]}"; do
    # counting all background jobs in current shell
    while (( $(jobs | wc -l) >= max_parallel_jobs )); do
            sleep 1
    done
    # start command in background
    {
        # echo "START: $b_cmd"
        eval "$b_cmd"
    } &
done
# wait on all background jobs from bedtools intersect
wait

rm -rd "${OUTDIR}/TMP"
echo " ---- completed ----"

