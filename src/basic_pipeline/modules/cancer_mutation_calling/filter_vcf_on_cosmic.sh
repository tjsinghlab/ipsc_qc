#!/bin/bash

#Usage: ./filter_vcf_cosmic.sh input.vcf.gz cosmic_data.tar gene_list.txt output.tsv

VCF_GZ_FILE="$1"
COSMIC_TAR="$2" #/Cosmic_CancerGeneCensus_Tsv_v101_GRCh37.tar
GENE_LIST="$3"
OUTPUT_FILE="$4"

# Check if input files exist
if [[ ! -f "$VCF_GZ_FILE" || ! -f "$COSMIC_TAR" || ! -f "$GENE_LIST" ]]; then
    echo "Error: One or more input files not found!"
    exit 1
fi

# Extract COSMIC tar file
echo "Extracting COSMIC data..."
TEMP_DIR=$(mktemp -d)
tar -xf "$COSMIC_TAR" -C "$TEMP_DIR"

# Find the COSMIC TSV.GZ file (assuming there is only one)
COSMIC_GZ_FILE=$(find "$TEMP_DIR" -name "*.tsv.gz" | head -n 1)

if [[ -z "$COSMIC_GZ_FILE" ]]; then
    echo "Error: No TSV.GZ file found in the extracted COSMIC data!"
    exit 1
fi

echo "Using COSMIC data from: $COSMIC_GZ_FILE"

# Decompress the COSMIC TSV.GZ file
COSMIC_FILE="${COSMIC_GZ_FILE%.gz}"  # Remove .gz extension
gunzip -c "$COSMIC_GZ_FILE" > "$COSMIC_FILE"

# Check COSMIC file size
echo "COSMIC file contains $(wc -l < "$COSMIC_FILE") lines."

# Decompress the VCF file
TEMP_VCF="filtered_variants.tmp"
echo "Decompressing VCF file..."
zcat "$VCF_GZ_FILE" | grep -v '^#' > "$TEMP_VCF"

# Check VCF file size
echo "VCF file contains $(wc -l < "$TEMP_VCF") variants."

# Read the gene list into an array
echo "Reading gene list from $GENE_LIST..."
mapfile -t GENES < "$GENE_LIST"

# Print each gene name to verify
echo "Genes being filtered:"
for gene in "${GENES[@]}"; do
    echo " - $gene"
done

# Extract COSMIC entries for the selected genes
echo "Filtering COSMIC data for selected genes..."
grep -F -f "$GENE_LIST" "$COSMIC_FILE" > "filtered_cosmic.tmp"

# Check filtered COSMIC file size
echo "Filtered COSMIC contains $(wc -l < "filtered_cosmic.tmp") lines."

# Extract chromosome, start, and stop positions from the filtered COSMIC file
awk -F'\t' '{print $15"\t"$16"\t"$17"\t"$1}' "filtered_cosmic.tmp" > "cosmic_positions.tmp"

# Extract chromosome and position from the VCF file
awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8}' "$TEMP_VCF" > "vcf_positions.tmp"

# Cross-reference VCF positions with COSMIC positions
echo "Cross-referencing VCF positions with COSMIC..."
awk 'NR==FNR {cosmic[$1,$2]=$0; next} ($1,$2) in cosmic {print cosmic[$1,$2]"\t"$0}' "cosmic_positions.tmp" "vcf_positions.tmp" > "$OUTPUT_FILE"

# Check final output size
echo "Final output contains $(wc -l < "$OUTPUT_FILE") matched variants."

# Cleanup
rm -rf "$TEMP_DIR" "$TEMP_VCF" "filtered_cosmic.tmp" "cosmic_positions.tmp" "vcf_positions.tmp"

echo "Filtered results saved to: $OUTPUT_FILE"