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
        echo "[ERROR] $tool not found. Install it or load HPC module."
        exit 1
    fi
done

# ----------------------
# User-supplied references (required)
# ----------------------
COSMIC_FILE="Cosmic_CancerGeneCensus_Tsc_v101_GRCh38.tar"
if [ ! -f "$REF_DIR/$COSMIC_FILE" ]; then
    echo "[ERROR] Required COSMIC reference missing: $REF_DIR/$COSMIC_FILE"
    echo "Please download this file from COSMIC and place it in $REF_DIR"
    exit 1
else
    echo "[INFO] Found COSMIC reference: $REF_DIR/$COSMIC_FILE"
fi

# Optional user-supplied file
# ----------------------
# Check or create genes.txt
# ----------------------
GENES_FILE="$REF_DIR/genes.txt"

if [ -f "$GENES_FILE" ]; then
    echo "[INFO] Found genes.txt: $GENES_FILE"
else
    echo "[WARN] genes.txt not found. Creating default genes.txt..."
    cat <<EOL > "$GENES_FILE"
BRCA1
BRCA2
TP53
BCOR
EGFR
EOL
    echo "[INFO] Created default genes.txt with BRCA1, BRCA2, TP53, BCOR, EGFR"
fi


# ----------------------
# PACNet training files
# ----------------------
PACNET_FILES=("Hs_expTrain_Jun-20-2017.rda" "Hs_stTrain-Jun-20-2017.rda")

for f in "${PACNET_FILES[@]}"; do
    if [ -f "$REF_DIR/$f" ]; then
        echo "[INFO] Found PACNet training file: $REF_DIR/$f"
    else
        echo "[WARN] PACNet training file missing: $REF_DIR/$f"
        echo "[INFO] Attempting to download from S3..."
        if command -v aws >/dev/null 2>&1; then
            aws s3 cp "s3://cellnet-rnaseq/ref/cnproc/HS/$f" "$REF_DIR/" --no-sign-request
            if [ -f "$REF_DIR/$f" ]; then
                echo "[INFO] Successfully downloaded $f"
            else
                echo "[ERROR] Failed to download PACNet file: $f"
                exit 1
            fi
        else
            echo "[ERROR] aws CLI not found. Please install AWS CLI v2 to download PACNet training files."
            exit 1
        fi
    fi
done

# ----------------------
# Auto-downloadable references
# ----------------------
FNA_FILE="GCF_000027235.1_ASM273v1_genomic.fna.gz"
FNA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/235/GCF_000027235.1_ASM273v1/GCF_000027235.1_ASM273v1_genomic.fna.gz"

if [ -f "$REF_DIR/$FNA_FILE" ]; then
    echo "[INFO] Found reference genome: $REF_DIR/$FNA_FILE"
else
    echo "[INFO] Reference genome not found, downloading from NCBI..."
    wget -O "$REF_DIR/$FNA_FILE" "$FNA_URL"
    echo "[INFO] Download complete: $REF_DIR/$FNA_FILE"
fi

# ----------------------
# Check dbSNP chr folder
# ----------------------
if [ ! -d "$CHR_DIR" ]; then
    echo "[ERROR] Missing 'chr' folder in reference directory: $CHR_DIR"
    echo "Please create this folder and add chromosome files named chr1..chr24"
    exit 1
else
    echo "[INFO] Found 'chr' folder: $CHR_DIR"
fi

# Expect chromosome files: chr1 .. chr24
CHRS=( {1..24} )

for chr in "${CHRS[@]}"; do
    file="$CHR_DIR/chr$chr"
    if [ ! -f "$file" ]; then
        echo "[ERROR] Missing chromosome file: $file"
        exit 1
    else
        echo "[INFO] Found chromosome file: $file"
    fi
done

echo "[INFO] All reference files and tools are present."
