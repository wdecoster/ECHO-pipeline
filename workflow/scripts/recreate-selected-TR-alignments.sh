#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Script Name:    recreate-selected-TR-alignments.sh
# Description:    Recreate per-haplotype alignment BAMs for selected TR IDs
#                 from a LongTR (or methylated) VCF + a phased BAM, suitable
#                 for IGV visualization. Outputs sorted/indexed BAMs plus the
#                 per-TR mini reference FASTA used for remapping.
#
# Dependencies:   [samtools, minimap2, bedtools, bcftools, awk, bgzip, tabix]
#
# Example:
#   bash recreate-selected-TR-alignments.sh \
#     -v sample_TR_methylation.vcf.gz \
#     -r GRCh38.fa \
#     -i sample.phased.bam \
#     -o out/recreated_alignments \
#     -s SAMPLE1 \
#     -d "chr1_123456,chr3_789012" \
#     -e 1000 \
#     -h "chrX,chrY,chrM"
#
# Or with a file:
#   bash recreate-selected-TR-alignments.sh ... -l ids.txt
#   (ids.txt: one TR ID per line; IDs are the VCF ID column)
# -----------------------------------------------------------------------------

set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 -v <vcf.gz> -r <reference.fa> -i <phased.bam> -o <output_dir> -s <sample_id> [-e <extend>] [-h <haploid_chroms>] (-d <id1,id2,...> | -l <id_list.txt>)
  -v  VCF file (LongTR or methylated), bgzipped and indexed (.tbi)
  -r  Reference FASTA (indexed with .fai)
  -i  Phased BAM (indexed with .bai), contains HP tags
  -o  Output base directory
  -s  Sample ID (used in filenames)
  -e  Flank size to extend around TR allele when building mini-reference (default: 1000)
  -h  Comma-separated haploid chromosomes (e.g. "chrX,chrY,chrY") (default: none)
  -d  Comma-separated TR IDs to process
  -l  File with TR IDs (one per line)
EOF
  exit 1
}

# Defaults
EXTEND=1000
HAPLOID_CHROMOSOMES=""

VCF=""
REF=""
PHASED_BAM=""
OUTDIR=""
SAMPLE=""
ID_CSV=""
ID_FILE=""

while getopts "v:r:i:o:s:e:h:d:l:" opt; do
  case $opt in
    v) VCF="$OPTARG" ;;
    r) REF="$OPTARG" ;;
    i) PHASED_BAM="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    s) SAMPLE="$OPTARG" ;;
    e) EXTEND="$OPTARG" ;;
    h) HAPLOID_CHROMOSOMES="$OPTARG" ;;
    d) ID_CSV="$OPTARG" ;;
    l) ID_FILE="$OPTARG" ;;
    *) usage ;;
  esac
done

# Basic checks
[[ -z "$VCF" || -z "$REF" || -z "$PHASED_BAM" || -z "$OUTDIR" || -z "$SAMPLE" ]] && usage
if [[ -z "${ID_CSV}" && -z "${ID_FILE}" ]]; then
  echo "Error: provide TR IDs via -d or -l" >&2; exit 1
fi

for tool in samtools minimap2 bedtools bcftools awk bgzip tabix; do
  command -v "$tool" >/dev/null 2>&1 || { echo "Missing tool: $tool" >&2; exit 1; }
done

[[ -f "$VCF" ]] || { echo "VCF not found: $VCF" >&2; exit 1; }
[[ -f "${VCF}.tbi" ]] || { echo "VCF index (.tbi) missing for: $VCF" >&2; exit 1; }
[[ -f "$REF" ]] || { echo "Reference FASTA not found: $REF" >&2; exit 1; }
[[ -f "${REF}.fai" ]] || { echo "Reference FASTA index (.fai) missing: $REF.fai" >&2; exit 1; }
[[ -f "$PHASED_BAM" ]] || { echo "Phased BAM not found: $PHASED_BAM" >&2; exit 1; }
[[ -f "${PHASED_BAM}.bai" ]] || { echo "BAM index (.bai) missing: ${PHASED_BAM}.bai" >&2; exit 1; }

mkdir -p "$OUTDIR"
LOG="${OUTDIR}/recreate_alignments_${SAMPLE}.log"
: > "$LOG"

# Helpers
is_haploid() {
  local chrom="$1"
  IFS=',' read -r -a arr <<< "$HAPLOID_CHROMOSOMES"
  for c in "${arr[@]}"; do
    [[ "$chrom" == "$c" ]] && return 0
  done
  return 1
}

# If entries have PDP in FORMAT and you want to zero-out haplotypes with PDP=0 (like your main),
# we can derive a masked GT on the fly. Otherwise we just trust GT as-is.
mask_gt_with_pdp_if_present() {
  # Input: full VCF line (no header)
  # Output: GT string with '.' where PDP allele support is 0 (diploid only)
  local line="$1"
  local fmt sample chrom
  fmt=$(echo "$line" | cut -f9)
  sample=$(echo "$line" | cut -f10)
  chrom=$(echo "$line" | cut -f1)

  # indices
  local gt_idx pdp_idx
  gt_idx=$(echo "$fmt" | awk -F':' '{for (i=1;i<=NF;i++) if($i=="GT") print i}')
  [[ -z "$gt_idx" ]] && { echo "$sample" | cut -d':' -f1; return 0; }

  local gt
  gt=$(echo "$sample" | cut -d':' -f"$gt_idx")
  # Haploid chromosomes: keep GT as is
  if is_haploid "$chrom"; then
    echo "$gt"; return 0
  fi

  pdp_idx=$(echo "$fmt" | awk -F':' '{for (i=1;i<=NF;i++) if($i=="PDP") print i}')
  if [[ -z "$pdp_idx" ]]; then
    echo "$gt"; return 0
  fi

  local pdp gt1 gt2 p1 p2
  pdp=$(echo "$sample" | cut -d':' -f"$pdp_idx")
  IFS='|' read -r gt1 gt2 <<< "$gt"
  IFS='|' read -r p1 p2 <<< "$pdp"

  [[ "$p1" == "0" ]] && gt1="."
  [[ "$p2" == "0" ]] && gt2="."
  echo "${gt1}|${gt2}"
}

# Collect IDs
TMP_IDS="${OUTDIR}/.ids.tmp"
: > "$TMP_IDS"
if [[ -n "$ID_CSV" ]]; then
  echo "$ID_CSV" | tr ',' '\n' | sed '/^[[:space:]]*$/d' >> "$TMP_IDS"
fi
if [[ -n "$ID_FILE" ]]; then
  sed 's/\r$//' "$ID_FILE" | sed '/^[[:space:]]*$/d' >> "$TMP_IDS"
fi
sort -u "$TMP_IDS" -o "$TMP_IDS"

echo "Processing $(wc -l < "$TMP_IDS") TR ID(s)..." | tee -a "$LOG"

# Per-ID loop
while read -r TRID; do
  [[ -z "$TRID" ]] && continue

  echo "== ID: ${TRID} ==" | tee -a "$LOG"

  # Pull the VCF line for this ID
  # Note: we assume unique IDs. If not unique, we take all matches (loop).
  bcftools view -H -i "ID==\"${TRID}\"" "$VCF" | while IFS= read -r LINE; do
    [[ -z "$LINE" ]] && continue

    CHROM=$(echo "$LINE" | cut -f1)
    POS=$(echo "$LINE" | cut -f2)
    ID_FIELD=$(echo "$LINE" | cut -f3)
    REF_ALLELE=$(echo "$LINE" | cut -f4)
    ALT_ALLELES=$(echo "$LINE" | cut -f5)
    INFO=$(echo "$LINE" | cut -f8)

    # START/END from INFO
    TR_START=$(echo "$INFO" | tr ';' '\n' | awk -F'=' '$1=="START"{print $2}')
    TR_END=$(echo "$INFO" | tr ';' '\n' | awk -F'=' '$1=="END"{print $2}')

    [[ -z "$TR_START" || -z "$TR_END" ]] && {
      echo "  ! Missing START/END for ${TRID}; skipping" | tee -a "$LOG"; continue; }

    # Build allele list
    IFS=',' read -r -a ALT_ARR <<< "$ALT_ALLELES"
    ALL_ALLELES=("$REF_ALLELE" "${ALT_ARR[@]}")

    # Determine GT (mask with PDP where available for diploid)
    GT=$(mask_gt_with_pdp_if_present "$LINE")
    # For haploids: we'll use a single "haploid" pass; for diploid, iterate over two HPs
    HAP_ARRAY=()
    if is_haploid "$CHROM"; then
      HAP_ARRAY=("haploid")
    else
      HAP_ARRAY=("1" "2")
    fi

    # Flanks from reference
    UPSTREAM_START=$((TR_START - EXTEND))
    [[ $UPSTREAM_START -lt 1 ]] && UPSTREAM_START=1
    UPSTREAM_END=$((TR_START - 1))
    DOWNSTREAM_START=$((TR_END + 1))
    DOWNSTREAM_END=$((TR_END + EXTEND))

    UPSTREAM_SEQ=$(samtools faidx "$REF" "${CHROM}:${UPSTREAM_START}-${UPSTREAM_END}" | tail -n +2 | tr -d '\n' | tr '[:upper:]' '[:lower:]' || true)
    DOWNSTREAM_SEQ=$(samtools faidx "$REF" "${CHROM}:${DOWNSTREAM_START}-${DOWNSTREAM_END}" | tail -n +2 | tr -d '\n' | tr '[:upper:]' '[:lower:]' || true)
    [[ -z "$UPSTREAM_SEQ" ]] && UPSTREAM_SEQ=""
    [[ -z "$DOWNSTREAM_SEQ" ]] && DOWNSTREAM_SEQ=""

    # Output folder for this TR
    TRDIR="${OUTDIR}/${CHROM}_${POS}_${ID_FIELD}"
    mkdir -p "$TRDIR"

    # Loop haplotypes
    for idx in "${!HAP_ARRAY[@]}"; do
      HP_LABEL="${HAP_ARRAY[$idx]}"

      # Get allele index for diploid; for haploid we use REF (0) or the called allele?
      if [[ "$HP_LABEL" == "haploid" ]]; then
        # Use GT if it is a single value like "0" or "1", else fall back to 0
        ALLELE_IDX=$(echo "$GT" | sed 's/|.*//' )
        [[ -z "$ALLELE_IDX" || "$ALLELE_IDX" == "." ]] && ALLELE_IDX="0"
      else
        # HP1 is left allele, HP2 is right allele in "a|b"
        if [[ "$HP_LABEL" == "1" ]]; then
          ALLELE_IDX=$(echo "$GT" | cut -d'|' -f1)
        else
          ALLELE_IDX=$(echo "$GT" | cut -d'|' -f2)
        fi
      fi

      if [[ -z "$ALLELE_IDX" || "$ALLELE_IDX" == "." ]]; then
        echo "  - HP${HP_LABEL}: no supported allele (masked/not present); skipping" | tee -a "$LOG"
        continue
      fi

      ALLELE_SEQ="${ALL_ALLELES[$ALLELE_IDX]}"
      if [[ "$ALLELE_SEQ" == "<DEL>" ]]; then
        echo "  - HP${HP_LABEL}: allele is <DEL>; skipping" | tee -a "$LOG"
        continue
      fi

      # Build mini-reference (flanks + UPPERCASE allele + flanks)
      HEADER="${CHROM}_${POS}_${ID_FIELD}_${HP_LABEL}"
      REGION_FASTA="${TRDIR}/${HEADER}_region.fasta"
      echo "  - Building mini-reference: ${HEADER}" | tee -a "$LOG"

      SEQ="${UPSTREAM_SEQ}$(echo "${ALLELE_SEQ}" | tr '[:lower:]' '[:upper:]')${DOWNSTREAM_SEQ}"
      echo ">$HEADER" > "$REGION_FASTA"
      echo "$SEQ" >> "$REGION_FASTA"
      samtools faidx "$REGION_FASTA"


      # Write BEDs within the mini-reference (0-based, half-open) ---
      UPLEN=${#UPSTREAM_SEQ}
      ALLEN=${#ALLELE_SEQ}
      DOWNLEN=${#DOWNSTREAM_SEQ}
      CONTIG="${HEADER}"
      TR_START_0=${UPLEN}
      TR_END_0=$((UPLEN + ALLEN))
      CONTIG_LEN=$((UPLEN + ALLEN + DOWNLEN))
      TR_BED="${TRDIR}/${HEADER}_TR.bed"
      echo -e "${CONTIG}\t${TR_START_0}\t${TR_END_0}\tTR" > "$TR_BED"

      echo "    * Wrote BEDs: $(basename "$TR_BED")" | tee -a "$LOG"

      # Extract reads overlapping the original TR locus from phased BAM
      REGION_BAM="${TRDIR}/${HEADER}_locus_reads.bam"
      samtools view -h "$PHASED_BAM" "${CHROM}:${TR_START}-${TR_END}" | samtools view -b -o "$REGION_BAM"

      # Filter by HP tag if diploid
      if [[ "$HP_LABEL" == "1" ]]; then
        FILTERED_BAM="${TRDIR}/${HEADER}_HP1.bam"
        samtools view -h -d HP:1 -b -o "$FILTERED_BAM" "$REGION_BAM"
      elif [[ "$HP_LABEL" == "2" ]]; then
        FILTERED_BAM="${TRDIR}/${HEADER}_HP2.bam"
        samtools view -h -d HP:2 -b -o "$FILTERED_BAM" "$REGION_BAM"
      else
        FILTERED_BAM="$REGION_BAM"
      fi

      if [[ ! -s "$FILTERED_BAM" ]]; then
        echo "    * No reads for HP${HP_LABEL}; skipping remap" | tee -a "$LOG"
        rm -f "$REGION_BAM"
        continue
      fi

      # Convert to FASTQ and remap to the mini-reference
      REGION_FASTQ="${TRDIR}/${HEADER}.fastq"
      samtools fastq -t -T MM,ML,HP,PS "$FILTERED_BAM" > "$REGION_FASTQ" 2>/dev/null || true
      if [[ ! -s "$REGION_FASTQ" ]]; then
        echo "    * No FASTQ reads after extraction; skipping" | tee -a "$LOG"
        rm -f "$REGION_BAM" "$FILTERED_BAM"
        continue
      fi

      REGION_SAM="${TRDIR}/${HEADER}_mapped.sam"
      minimap2 -ax map-ont -y --sam-hit-only "$REGION_FASTA" "$REGION_FASTQ" > "$REGION_SAM" 2>>"$LOG"

      # Sort/index BAM for IGV
      REALIGNED_BAM="${TRDIR}/${HEADER}_mapped.sorted.bam"
      samtools view -bS "$REGION_SAM" | samtools sort -o "$REALIGNED_BAM"
      samtools index "$REALIGNED_BAM"

      echo "    * Wrote: ${REALIGNED_BAM} (+ .bai) and ${REGION_FASTA}" | tee -a "$LOG"

      # Keep intermediates that help reproducibility; remove bulky ones
      rm -f "$REGION_SAM" "$REGION_FASTQ" "$REGION_BAM"
      [[ "$FILTERED_BAM" != "$REGION_BAM" ]] && rm -f "$FILTERED_BAM"

    done # haplotypes
  done

done < "$TMP_IDS"

echo "Done. BAMs ready under: $OUTDIR"

