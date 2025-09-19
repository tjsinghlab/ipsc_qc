#!/bin/bash
set -euo pipefail

# Defaults
SAMPLE=""
COSMIC_TAR=""
OUTPUT_DIR=""
REF_DIR=""


# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2 ;;
    --cosmic) COSMIC_TAR="$2"; shift 2 ;;
    --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
    --ref_dir) REF_DIR="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 --sample Sample name --vcf input.vcf.gz --cosmic user-downloaded COSMIC reference file --output_dir directory to sample-specific output directory"
      exit 0
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

# Check required arguments
if [[ -z "$SAMPLE" || -z "$COSMIC_TAR" || -z "$OUTPUT_DIR" || -z "$REF_DIR"]]; then
  echo "Error: Missing required arguments."
  echo "Run with --help for usage."
  exit 1
fi

# Validate input files
VCF_GZ_FILE="${OUTPUT_DIR}/${SAMPLE}_variants.vcf.gz"

for f in "$VCF_GZ_FILE" "$COSMIC_TAR" "$GENE_LIST"; do
  if [[ ! -f "$f" ]]; then
    echo "Error: File not found - $f"
    exit 1
  fi
done

# Make sure output directory exists
mkdir -p "$(dirname "$OUTPUT_DIR")"

# Create sample-specific temp directory
TEMP_DIR=$(mktemp -d -t cosmic_${SAMPLE}_XXXX)

# Gene list is fixed inside Docker image
GENE_LIST="$REF_DIR/gene_list.txt"

# Extract COSMIC tar file
echo "[INFO][$SAMPLE] Extracting COSMIC data..."
tar -xf "$COSMIC_TAR" -C "$TEMP_DIR"

COSMIC_GZ_FILE=$(find "$TEMP_DIR" -name "*.tsv.gz" | head -n 1)
if [[ -z "$COSMIC_GZ_FILE" ]]; then
  echo "[ERROR][$SAMPLE] No TSV.GZ file found in COSMIC tarball!"
  exit 1
fi

echo "[INFO][$SAMPLE] Using COSMIC data from: $COSMIC_GZ_FILE"
COSMIC_FILE="${COSMIC_GZ_FILE%.gz}"
gunzip -c "$COSMIC_GZ_FILE" > "$COSMIC_FILE"

# Decompress VCF
TEMP_VCF=$(mktemp -t vcf_${SAMPLE}_XXXX)
echo "[INFO][$SAMPLE] Decompressing VCF..."
zcat "$VCF_GZ_FILE" | grep -v '^#' > "$TEMP_VCF"

# Filter COSMIC for selected genes (from /ref/gene_list.txt inside image)
echo "[INFO][$SAMPLE] Filtering COSMIC for genes in $GENE_LIST..."
grep -F -f "$GENE_LIST" "$COSMIC_FILE" > "${TEMP_DIR}/filtered_cosmic.tmp"

# Extract COSMIC positions
awk -F'\t' '{print $15"\t"$16"\t"$17"\t"$1}' "${TEMP_DIR}/filtered_cosmic.tmp" > "${TEMP_DIR}/cosmic_positions.tmp"

# Extract VCF positions
awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8}' "$TEMP_VCF" > "${TEMP_DIR}/vcf_positions.tmp"

OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE}_cosmic.tsv"
# Cross-reference
echo "[INFO][$SAMPLE] Cross-referencing VCF with COSMIC..."
awk 'NR==FNR {cosmic[$1,$2]=$0; next} ($1,$2) in cosmic {print cosmic[$1,$2]"\t"$0}' \
  "${TEMP_DIR}/cosmic_positions.tmp" \
  "${TEMP_DIR}/vcf_positions.tmp" \
  > "$OUTPUT_FILE"

echo "[INFO][$SAMPLE] Final output: $(wc -l < "$OUTPUT_FILE") matched variants"
echo "[INFO][$SAMPLE] Results saved to: $OUTPUT_FILE"

# Cleanup
rm -rf "$TEMP_DIR" "$TEMP_VCF"
