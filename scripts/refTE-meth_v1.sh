#!/bin/bash

#Set up directories
PHASED_PILEUP="/home/leena/UNIQUE/GIAB/modkit/HG002/phased"
UNPHASED_PILEUP="/home/leena/UNIQUE/GIAB/modkit/HG002/unphased"
VARIATION="/home/leena/UNIQUE/GIAB/phased-variation"
TE="retroposon" #Options: all, all_cpg, DNA, DNA_cpg, LINE, LINE_cpg, helitron, helitron_cpg, SINE, SINE_cpg, LTR, LTR_cpg, retroposon, retroposon_cpg
TE_CATALOG="/home/leena/UNIQUE/repeat-catalogues/TEs/GRCh38/TE_classes/GRCh38_TEs_${TE}.bed"
OUTDIR="/home/leena/UNIQUE/GIAB/TEs/HG002/GRCh38/refTE"
SAMPLE_ID="HG002"

#Set up input files
phased_pileup_1="${PHASED_PILEUP}/GRCh38_HG002_sup_modkit_m0.8_1.bed.gz"
phased_pileup_2="${PHASED_PILEUP}/GRCh38_HG002_sup_modkit_m0.8_2.bed.gz"
unphased_pileup="${UNPHASED_PILEUP}/GRCh38_HG002_sup_modkit_m0.8.bed.gz"
SNP_data="${VARIATION}/${SAMPLE_ID}/phased_sup_SNPs.vcf.gz"
SV_data="${VARIATION}/${SAMPLE_ID}/phased_sup_SV.vcf.gz"
SNP_filt="${VARIATION}/${SAMPLE_ID}/phased_sup_SNPs_filt.vcf"
SV_filt="${VARIATION}/${SAMPLE_ID}/phased_sup_SV_filt.vcf"

mkdir -p "${OUTDIR}/mod_phased" "${OUTDIR}/mod_unphased" "${OUTDIR}/variants" 

# Generate shifted bed file for upstream TE regions (-250 bp)
UPSTREAM_TE_CATALOG="${OUTDIR}/${TE}_upstream.bed"
awk '{OFS="\t"} {start=$2-250; if (start < 0) start=0; print $1, start, $2, $4}' "$TE_CATALOG" > "$UPSTREAM_TE_CATALOG"

# Filter variant data
bcftools view "$SNP_data" -i 'FILTER="PASS" & FORMAT/DP>5' --threads 8 --output "$SNP_filt"
bcftools view "$SV_data" -i 'FILTER="PASS" & INFO/SUPPORT>5' --threads 8 --output "$SV_filt"
bgzip "$SNP_filt"
bgzip "$SV_filt"
tabix "$SNP_filt.gz"
tabix "$SV_filt.gz"

# Intersect modkit pileup data on (upstream) TE catalog to obtain methylation of individual CpGs that are embedded in TE elements or their upstream regions.
zcat "$phased_pileup_1" | bedtools intersect -a - -b "${TE_CATALOG}" -wa -wb > "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_pileup_1.bed"
zcat "$phased_pileup_2" | bedtools intersect -a - -b "${TE_CATALOG}" -wa -wb > "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_pileup_2.bed"
zcat "$unphased_pileup" | bedtools intersect -a - -b "${TE_CATALOG}" -wa -wb > "${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_pileup_unphased.bed"
zcat "$phased_pileup_1" | bedtools intersect -a - -b "${UPSTREAM_TE_CATALOG}" -wa -wb > "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_upstream_pileup_1.bed"
zcat "$phased_pileup_2" | bedtools intersect -a - -b "${UPSTREAM_TE_CATALOG}" -wa -wb > "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_upstream_pileup_2.bed"
zcat "$unphased_pileup" | bedtools intersect -a - -b "${UPSTREAM_TE_CATALOG}" -wa -wb > "${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_upstream_pileup_unphased.bed"

# Modkit stats across specified TE regions
modkit stats -t 16 --regions "$TE_CATALOG" -c m -o "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_stats_1.tsv" "$phased_pileup_1" 2>/dev/null
modkit stats -t 16 --regions "$TE_CATALOG" -c m -o "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_stats_2.tsv" "$phased_pileup_2" 2>/dev/null
modkit stats -t 16 --regions "$TE_CATALOG" -c m -o "${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_stats_unphased.tsv" "$unphased_pileup" 2>/dev/null
modkit stats -t 16 --regions "$UPSTREAM_TE_CATALOG" -c m -o "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_upstream_stats_1.tsv" "$phased_pileup_1" 2>/dev/null
modkit stats -t 16 --regions "$UPSTREAM_TE_CATALOG" -c m -o "${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_upstream_stats_2.tsv" "$phased_pileup_2" 2>/dev/null
modkit stats -t 16 --regions "$UPSTREAM_TE_CATALOG" -c m -o "${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_upstream_stats_unphased.tsv" "$unphased_pileup" 2>/dev/null

# Report and SNPs and SVs that intersect with TE elements of interest
SNP_INTERSECT="${OUTDIR}/variants/${SAMPLE_ID}_${TE}_SNPs_intersect.bed"
SV_INTERSECT="${OUTDIR}/variants/${SAMPLE_ID}_${TE}_SV_intersect.bed"
SNP_INTERSECT_COUNT="${OUTDIR}/variants/${SAMPLE_ID}_${TE}_SNPs_intersect_count.bed"
SV_INTERSECT_COUNT="${OUTDIR}/variants/${SAMPLE_ID}_${TE}_SV_intersect_count.bed"

bedtools intersect -a "${TE_CATALOG}" -b "${SNP_filt}.gz" -wa -wb > "$SNP_INTERSECT"
bedtools intersect -a "${TE_CATALOG}" -b "${SV_filt}.gz" -wa -wb > "$SV_INTERSECT"
bedtools intersect -a "${TE_CATALOG}" -b "${SNP_filt}.gz" -wa -wb -C > "$SNP_INTERSECT_COUNT"
bedtools intersect -a "${TE_CATALOG}" -b "${SV_filt}.gz" -wa -wb -C > "$SV_INTERSECT_COUNT"

# Make summarizing table of output
OUTPUT="${OUTDIR}/${SAMPLE_ID}_summary_per_ref_${TE}.txt"

echo "making the final summary file.."
echo -e "chr\tstart\tend\tfamily\t.\tstrand\tID\tTE_length\tTE_avgMeth_phased\tTE_Nvalid_phased\tupstream_avgMeth_phased\tupstream_Nvalid_phased\tTE_avgMeth_unphased\tTE_Nvalid_unphased\tupstream_avgMeth_unphased\tupstream_Nvalid_unphased\ttotal_SNP\tSNP_count\tINDEL_count\tSV_count\tSV_types\tSV_IDs" > "$OUTPUT"

# Define methylation stats files
stats_1="${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_stats_1.tsv"
stats_2="${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_stats_2.tsv"
stats_unph="${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_stats_unphased.tsv"
up_1="${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_upstream_stats_1.tsv"
up_2="${OUTDIR}/mod_phased/${SAMPLE_ID}_${TE}_upstream_stats_2.tsv"
up_unph="${OUTDIR}/mod_unphased/${SAMPLE_ID}_${TE}_upstream_stats_unphased.tsv"

# Function to get methylation percent,Nvalid by coordinates
read_meth_stats() {
    local stats_file=$1
    local chr=$2
    local start=$3
    local end=$4

    awk -v chr="$chr" -v start="$start" -v end="$end" -F'\t' \
        '$1 == chr && $2 == start && $3 == end {print $8 "," $7; exit}' "$stats_file"
}


# Main loop: Read BED line by line
while IFS=$'\t' read -r chr start end family dot strand id; do
    region="${chr}:${start}-${end}"

    # Extract length
    length=$((end - start))

    # Extract SNPs in region (PASS only) and count total
    total_snp=$(bcftools view -r "$region" "${SNP_filt}.gz" | grep -v "^#" | wc -l)
    snp_count=$(bcftools view -r "$region" -v snps "${SNP_filt}.gz" | grep -v "^#" | wc -l)
    indel_count=$(bcftools view -r "$region" -v indels "${SNP_filt}.gz" | grep -v "^#" | wc -l)

    # Methylation values (TE and upstream)
    TE_hp1=$(read_meth_stats "$stats_1" "$chr" "$start" "$end")
    TE_hp2=$(read_meth_stats "$stats_2" "$chr" "$start" "$end")
    TE_unph=$(read_meth_stats "$stats_unph" "$chr" "$start" "$end")

    up_start=$(( start - 250 )); [[ $up_start -lt 0 ]] && up_start=0
    up_end=$start

    UP_hp1=$(read_meth_stats "$up_1" "$chr" "$up_start" "$up_end")
    UP_hp2=$(read_meth_stats "$up_2" "$chr" "$up_start" "$up_end")
    UP_unph=$(read_meth_stats "$up_unph" "$chr" "$up_start" "$up_end")

    # Fallbacks if missing 
    [ -z "$TE_hp1" ] && TE_hp1="NA,NA"
    [ -z "$TE_hp2" ] && TE_hp2="NA,NA"
    [ -z "$TE_unph" ] && TE_unph="NA,NA"
    [ -z "$UP_hp1" ] && UP_hp1="NA,NA"
    [ -z "$UP_hp2" ] && UP_hp2="NA,NA"
    [ -z "$UP_unph" ] && UP_unph="NA,NA"

    # Format outputs
    TE_avgMeth="$(echo "$TE_hp1" | cut -d',' -f1),$(echo "$TE_hp2" | cut -d',' -f1)"
    TE_Nvalid="$(echo "$TE_hp1" | cut -d',' -f2),$(echo "$TE_hp2" | cut -d',' -f2)"
    UP_avgMeth="$(echo "$UP_hp1" | cut -d',' -f1),$(echo "$UP_hp2" | cut -d',' -f1)"
    UP_Nvalid="$(echo "$UP_hp1" | cut -d',' -f2),$(echo "$UP_hp2" | cut -d',' -f2)"
    TE_avgMeth_unphased="$(echo "$TE_unph" | cut -d',' -f1)"
    TE_Nvalid_unphased="$(echo "$TE_unph" | cut -d',' -f2)"
    UP_avgMeth_unphased="$(echo "$UP_unph" | cut -d',' -f1)"
    UP_Nvalid_unphased="$(echo "$UP_unph" | cut -d',' -f2)"

    # Get overlapping SVs for this TE region from SV_INTERSECT
    sv_matches=$(awk -v chr="$chr" -v start="$start" -v end="$end" \
        '$1 == chr && $2 == start && $3 == end' "$SV_INTERSECT")

    sv_count=0
    svtypes="NA"
    sv_ids="NA"

    if [[ -n "$sv_matches" ]]; then
        sv_count=$(echo "$sv_matches" | wc -l)
        svtypes=$(echo "$sv_matches" | awk -F'\t' '{
            n = split($15, info, ";");
            for (i = 1; i <= n; i++) {
                split(info[i], kv, "=");
                if (kv[1] == "SVTYPE") {
                    print kv[2];
                }
            }
        }' | sort | uniq -c | awk '{printf "%s:%s,", $2, $1}' | sed 's/,$//')
        sv_ids=$(echo "$sv_matches" | awk '{print $10}' | paste -sd "," -)
    fi

    # Output line
    echo -e "${chr}\t${start}\t${end}\t${family}\t${dot}\t${strand}\t${id}\t${length}\t${TE_avgMeth}\t${TE_Nvalid}\t${UP_avgMeth}\t${UP_Nvalid}\t${TE_avgMeth_unphased}\t${TE_Nvalid_unphased}\t${UP_avgMeth_unphased}\t${UP_Nvalid_unphased}\t${total_snp}\t${snp_count}\t${indel_count}\t${sv_count}\t${svtypes}\t${sv_ids}" >> "$OUTPUT"

done < "$TE_CATALOG"

echo "Summary written to: $OUTPUT"





