#!/bin/bash
# Prepare shifted bed files and filtered VCFs for processing with modkit stats
# Steps 1 and 2 from ref_TE_avg_meth.sh

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

# ------------ Initialise Flag Values ------------
## hard-coded values set here will override rule-defined values
PHASED_PILEUP=""
UNPHASED_PILEUP=""
VARIATION=""
TE=""
TE_CATALOG=""
SAMPLE_ID=""
OUTDIR=""
TOTAL_THREADS=""
FLANK=""

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
# _filt are outputs
# original script writes them to 03_phasing, write to current outdir instead to avoid overwriting comparison run files
SNP_filt="${OUTDIR}/variants/${SAMPLE_ID}_SNP_filt.vcf"
SV_filt="${OUTDIR}/variants/${SAMPLE_ID}_SV_filt.vcf"


# ------------ (1) Generate shifted bed files --------------------------------------------------------------
# Generate shifted bed file for upstream TE regions

UPSTREAM_TE_CATALOG="${OUTDIR}/${TE}_upstream.bed"
DOWNSTREAM_TE_CATALOG="${OUTDIR}/${TE}_downstream.bed"
# coarse-grained" parallelization for two separate output files. with the background parameter '&'
#++ Time measuring ++#
start_t=$(date +%s.%N)
awk -v flank="$FLANK" '{OFS="\t"} {start=$2-flank; if (start < 0) start=0; print $1, start, $2, $7}' "$TE_CATALOG" > "$UPSTREAM_TE_CATALOG" & 
UPstream_PID=$!
awk -v flank="$FLANK" '{OFS="\t"} {print $1, $3, $3+flank, $7}' "$TE_CATALOG" > "$DOWNSTREAM_TE_CATALOG" &
DOWNstream_PID=$!

#++ time ++#
end_t=$(date +%s.%N)
elapsed=$(echo "$end_t - $start_t" | bc)
echo "Executed time for (1) generating shifted bed files: $elapsed seconds"

# ------------ (2) Filter VCF ------------------------------------------------------------------------------------
# Now filtering for at least 5 covered reads used for variant calling

#++ Time measuring ++#
start_t=$(date +%s.%N)
# --threads is only used for compression/decompression when using --output-type/-O
bcftools view "${SNP_data}" \
    -i 'FILTER="PASS" & FORMAT/DP>5' \
    --threads "${TOTAL_THREADS}" \
    --output-file "${SNP_filt}.gz" \
    -O z && tabix "${SNP_filt}.gz" &
VCF_SNP_PID=$!

bcftools view "${SV_data}" \
     -i 'FILTER="PASS" & INFO/SUPPORT>5' \
     --threads "${TOTAL_THREADS}" \
     --output-file "${SV_filt}.gz" \
     -O z && tabix "${SV_filt}.gz" &
VCF_SV_PID=$!
wait $UPstream_PID $DOWNstream_PID $VCF_SNP_PID $VCF_SV_PID

end_t=$(date +%s.%N)
elapsed=$(echo "$end_t - $start_t" | bc)
echo "Executed time for (2) filtering vcf files : $elapsed seconds"

# --------- Variant Intersections (previously step 4) --------------------------

# Report SNPs and SVs that intersect with TE elements of interest
# define paths and max parallel jobs
SNP_INTERSECT="${OUTDIR}/variants/${SAMPLE_ID}_SNPs_intersect.bed"
SV_INTERSECT="${OUTDIR}/variants/${SAMPLE_ID}_SV_intersect.bed"
SNP_INTERSECT_COUNT="${OUTDIR}/variants/${SAMPLE_ID}_SNPs_intersect_count.bed"
SV_INTERSECT_COUNT="${OUTDIR}/variants/${SAMPLE_ID}_SV_intersect_count.bed"

echo "Starting serial bedtools intersects..."
echo "call 1"
bedtools intersect -a "${TE_CATALOG}" -b ${SNP_filt}.gz -wa -wb > "${SNP_INTERSECT}"
echo "call 2"
bedtools intersect -a "${TE_CATALOG}" -b "${SV_filt}.gz" -wa -wb > "${SV_INTERSECT}"
echo "call 3"
bedtools intersect -a "${TE_CATALOG}" -b "${SNP_filt}.gz" -wa -wb -C > "${SNP_INTERSECT_COUNT}"
echo "call 4"
bedtools intersect -a "${TE_CATALOG}" -b "${SV_filt}.gz" -wa -wb -C > "${SV_INTERSECT_COUNT}"

echo "Variant intersections & preparation for modkit stats complete"
