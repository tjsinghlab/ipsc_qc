#!/usr/bin/env bash
set -euo pipefail

# ----------------------
# Directories
# ----------------------
REF_DIR="${1:-/refs}"  # default to /refs if not provided
mkdir -p "$REF_DIR"

CHR_DIR="$REF_DIR/chr"

# ----------------------
# Tool checks
# ----------------------
TOOLS=("bowtie2" "samtools" "bcftools" "vcftools")

for tool in "${TOOLS[@]}"; do
    if command -v "$tool" >/dev/null 2>&1; then
        echo "[INFO] $tool found in PATH"
    elif command -v module >/dev/null 2>&1; then
        echo "[INFO] Loading $tool module..."
        module load "$tool"
        if ! command -v "$tool" >/dev/null 2>&1; then
            echo "[ERROR] Failed to load $tool via module"
            exit 1
        fi
    else
        echo "[ERROR] $tool not found. Install it, use Docker, or load HPC module."
        exit 1
    fi
done

# ----------------------
# Reference files
# ----------------------
# Files that need to be supplied by user
USER_FILES=("Cosmic_CancerGeneCensus_Tsv_v101_GRCh37.tar")

for f in "${USER_FILES[@]}"; do
    if [ ! -f "$REF_DIR/$f" ]; then
        echo "[ERROR] Required user-supplied reference missing: $REF_DIR/$f"
        echo "Please download this file from COSMIC and place it in $REF_DIR"
        exit 1
    else
        echo "[INFO] Found user-supplied reference: $REF_DIR/$f"
    fi
done

# Files that can be downloaded automatically
declare -A DOWNLOAD_FILES
DOWNLOAD_FILES["GCF_000027325.1_ASM2732v1_genomic.fna.gz"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz"
DOWNLOAD_FILES["training_file1"]="https://raw.githubusercontent.com/CahanLab/CellNet/master/data/hsTFs.rda"
DOWNLOAD_FILES["training_file2"]="https://raw.githubusercontent.com/CahanLab/CellNet/master/data/mmTFs.rda"

download_if_missing() {
    local file_path="$1"
    local url="$2"
    if [ -f "$file_path" ]; then
        echo "[INFO] Found existing reference: $file_path"
    else
        echo "[INFO] Downloading $file_path from $url..."
        wget -O "$file_path" "$url"
        echo "[INFO] Download complete: $file_path"
    fi
}

for file in "${!DOWNLOAD_FILES[@]}"; do
    download_if_missing "$REF_DIR/$file" "${DOWNLOAD_FILES[$file]}"
done

# ----------------------
# Check dbSNP chr folder
# ----------------------
if [ ! -d "$CHR_DIR" ]; then
    echo "[ERROR] Missing 'chr' folder in reference directory: $CHR_DIR"
    echo "Please create this folder and add dbSNP build 142 common SNP GTF files for each chromosome."
    exit 1
else
    echo "[INFO] Found 'chr' folder: $CHR_DIR"
fi

# Expected chromosomes 1-22, X=23, Y=24
CHRS=( {1..22} 23 24 )

for chr in "${CHRS[@]}"; do
    file="$CHR_DIR/snp142common_$chr.gtf"
    if [ ! -f "$file" ]; then
        echo "[ERROR] Missing SNP file for chromosome $chr: $file"
        echo "Please download it from UCSC Table Browser (snp142Common) and place it in $CHR_DIR"
        exit 1
    else
        echo "[INFO] Found SNP file: $file"
    fi
done

echo "[INFO] All reference files and tools are present."

