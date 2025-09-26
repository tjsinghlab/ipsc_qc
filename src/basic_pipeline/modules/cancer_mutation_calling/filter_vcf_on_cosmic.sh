#!/bin/bash
set -euo pipefail

# ----------------------------
# Defaults
# ----------------------------
SAMPLE=""
COSMIC_DIR=""
OUTPUT_DIR=""
REF_DIR=""

# ----------------------------
# Parse arguments
# ----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2 ;;
    --cosmic_dir) COSMIC_DIR="$2"; shift 2 ;;
    --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
    --ref_dir) REF_DIR="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 --sample SAMPLE --cosmic_dir COSMIC_DIR --output_dir OUTDIR --ref_dir REF_DIR"
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

# ----------------------------
# Check required arguments
# ----------------------------
if [[ -z "$SAMPLE" || -z "$COSMIC_DIR" || -z "$OUTPUT_DIR" || -z "$REF_DIR" ]]; then
  echo "[ERROR] Missing required arguments." >&2
  echo "Run with --help for usage."
  exit 1
fi

echo "[INFO][${SAMPLE}] Looking for COSMIC data in $COSMIC_DIR"

# ----------------------------
# Locate and extract COSMIC
# ----------------------------
COSMIC_TAR="$COSMIC_DIR/Cosmic_CancerGeneCensus_Tsv_v101_GRCh38.tar"
COSMIC_EXTRACT_DIR="$COSMIC_DIR/unpacked"
COSMIC_FILE="$COSMIC_EXTRACT_DIR/Cosmic_CancerGeneCensus_v101_GRCh38.tsv.gz"

if [[ ! -f "$COSMIC_TAR" ]]; then
  echo "[ERROR][${SAMPLE}] COSMIC tarball not found: $COSMIC_TAR" >&2
  exit 1
fi

mkdir -p "$COSMIC_EXTRACT_DIR"

# Only extract if the file is missing
if [[ ! -f "$COSMIC_FILE" ]]; then
  echo "[INFO][${SAMPLE}] Extracting COSMIC tarball..."
  tar -xf "$COSMIC_TAR" -C "$COSMIC_EXTRACT_DIR"
fi

if [[ ! -f "$COSMIC_FILE" ]]; then
  echo "[ERROR][${SAMPLE}] Could not find Cosmic_CancerGeneCensus_v101_GRCh38.tsv.gz after extraction" >&2
  exit 1
fi

echo "[INFO][${SAMPLE}] Found COSMIC file: $COSMIC_FILE"

# ----------------------------
# Validate other inputs
# ----------------------------
VCF_GZ_FILE="${OUTPUT_DIR}/variant_calling/${SAMPLE}${SAMPLE}.variant_filtered.vcf.gz"
GENE_LIST="$REF_DIR/genes.txt"

for f in "$VCF_GZ_FILE" "$GENE_LIST"; do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR][${SAMPLE}] File not found: $f" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"
OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE}_cosmic.tsv"

# ----------------------------
# Preview VCF
# ----------------------------
echo "[INFO][${SAMPLE}] VCF preview (first 5 variant lines):"
zcat "$VCF_GZ_FILE" | grep -v '^#' | head -n 5 | \
  awk '{print "CHR="$1, "POS="$2, "REF="$4, "ALT="$5, "INFO="$8}'

# ----------------------------
# Preview COSMIC
# ----------------------------
echo "[INFO][${SAMPLE}] COSMIC preview (first 5 lines):"
if file "$COSMIC_FILE" | grep -q "gzip"; then
  zcat "$COSMIC_FILE" | head -n 5 | \
    awk -F'\t' '{print "GENE="$1, "CHR="$4, "START="$5, "END="$6}'
else
  head -n 5 "$COSMIC_FILE" | \
    awk -F'\t' '{print "GENE="$1, "CHR="$4, "START="$5, "END="$6}'
fi
# echo "[INFO][${SAMPLE}] COSMIC preview (first 5 lines):"
# zcat "$COSMIC_FILE" | head -n 5 | \
#   awk -F'\t' '{print "GENE="$1, "CHR="$4, "START="$5, "END="$6}'

# ----------------------------
# Filter COSMIC for genes
# ----------------------------
FILTERED_COSMIC="${COSMIC_EXTRACT_DIR}/filtered_${SAMPLE}.tsv"
echo "[INFO][${SAMPLE}] Filtering COSMIC for genes in $GENE_LIST..."
zcat "$COSMIC_FILE" | awk -F'\t' 'NR==1 || FNR==NR {genes[$1]; next} $1 in genes' "$GENE_LIST" - > "$FILTERED_COSMIC"
echo "[INFO][${SAMPLE}] Filtered COSMIC contains $(($(wc -l < "$FILTERED_COSMIC") - 1)) entries."

# ----------------------------
# Extract positions
# ----------------------------
COSMIC_POS="${COSMIC_EXTRACT_DIR}/cosmic_positions_${SAMPLE}.tmp"
VCF_POS="${COSMIC_EXTRACT_DIR}/vcf_positions_${SAMPLE}.tmp"

# COSMIC columns: GENE_SYMBOL(1), CHROMOSOME(4), GENOME_START(5), GENOME_STOP(6)
awk -F'\t' 'NR>1 {chr=$4; sub(/^chr/,"",chr); print chr"\t"$5"\t"$6"\t"$1}' "$FILTERED_COSMIC" > "$COSMIC_POS"

# VCF: chrom(1), pos(2), ref(4), alt(5), info(8)
zcat "$VCF_GZ_FILE" | grep -v '^#' | \
  awk -F'\t' '{chr=$1; sub(/^chr/,"",chr); print chr"\t"$2"\t"$2"\t"$4"\t"$5"\t"$8}' > "$VCF_POS"

# ----------------------------
# Check chromosome formats
# ----------------------------
echo "[INFO][${SAMPLE}] Chromosome format check:"
echo "COSMIC unique chromosomes:"; cut -f1 "$COSMIC_POS" | sort -u | head
echo "VCF unique chromosomes:"; cut -f1 "$VCF_POS" | sort -u | head

# ----------------------------
# Cross-reference
# ----------------------------
echo "[INFO][${SAMPLE}] Cross-referencing VCF with COSMIC..."
awk 'NR==FNR {cosmic[$1":"$2]=$0; next} 
     {key=$1":"$2; if (key in cosmic) print cosmic[key]"\t"$0}' \
  "$COSMIC_POS" "$VCF_POS" > "$OUTPUT_FILE" || true

NUM_MATCHES=$(wc -l < "$OUTPUT_FILE")
echo "[INFO][${SAMPLE}] Final output contains $NUM_MATCHES matched variants."
echo "[INFO][${SAMPLE}] Filtered results saved to: $OUTPUT_FILE"
