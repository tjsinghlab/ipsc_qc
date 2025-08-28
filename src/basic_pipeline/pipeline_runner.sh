#!/usr/bin/env bash
set -euo pipefail

# ----------------------
# Defaults
# ----------------------
PROJECT="default_project"
OUTPUT_DIR="$(pwd)/ipsc_qc_outputs"
REF_DIR="$(pwd)/ref"
VCF_DIR=""
RSEM_DIR=""
BAM_DIR=""
FASTQ_DIR=""

# ----------------------
# Usage
# ----------------------
usage() {
    echo "Usage: pipeline_runner.sh [options]"
    echo ""
    echo "Required arguments:"
    echo "  --vcf_dir PATH       Directory containing VCF files"
    echo "  --rsem_dir PATH      Directory containing RSEM outputs"
    echo "  --bam_dir PATH       Directory containing BAM files"
    echo "  --fastq_dir PATH     Directory containing FASTQ files"
    echo ""
    echo "Optional arguments:"
    echo "  --ref_dir PATH       Reference directory (default: ./ref)"
    echo "  --project NAME       Project name (default: ${PROJECT})"
    echo "  --output_dir PATH    Output directory (default: ./ipsc_qc_outputs)"
    echo ""
    exit 1
}

# ----------------------
# Parse arguments
# ----------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --project)    PROJECT="$2"; shift ;;
        --ref_dir)    REF_DIR="$2"; shift ;;
        --vcf_dir)    VCF_DIR="$2"; shift ;;
        --rsem_dir)   RSEM_DIR="$2"; shift ;;
        --bam_dir)    BAM_DIR="$2"; shift ;;
        --fastq_dir)  FASTQ_DIR="$2"; shift ;;
        --input_dir)  INPUT_DIR="$2"; shift ;;
        --output_dir) OUTPUT_DIR="$2"; shift ;;
        -h|--help)    usage ;;
        *) echo "[ERROR] Unknown option: $1"; usage ;;
    esac
    shift
done

# ----------------------
# Validate required args
# ----------------------
if [[ -z "$REF_DIR" || -z "$VCF_DIR" || -z "$RSEM_DIR" || -z "$BAM_DIR" || -z "$FASTQ_DIR" ]]; then
    echo "[ERROR] Missing one or more required arguments."
    usage
fi

# ----------------------
# Create output directory if not exists
# ----------------------
mkdir -p "$OUTPUT_DIR"

# ----------------------
# Check tools and references
# ----------------------
bash scripts/bash/reference_check.sh "$REF_DIR"

# ----------------------
# RSEM and BAM file names must match
# ----------------------
RSEM_FILES=($(ls "$RSEM_DIR" | sort))
BAM_FILES=($(ls "$BAM_DIR" | sort))

if [ ${#RSEM_FILES[@]} -eq 0 ] || [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "[ERROR] One or both directories ($RSEM_DIR, $BAM_DIR) are empty."
    exit 1
fi

if [ ${#RSEM_FILES[@]} -ne ${#BAM_FILES[@]} ]; then
    echo "[ERROR] Mismatch in number of files between RSEM and BAM directories."
    echo "  RSEM count: ${#RSEM_FILES[@]}"
    echo "  BAM count:  ${#BAM_FILES[@]}"
    exit 1
fi

for i in "${!RSEM_FILES[@]}"; do
    base_rsem=$(basename "${RSEM_FILES[$i]}" .txt | sed 's/\..*//')
    base_bam=$(basename "${BAM_FILES[$i]}" .bam | sed 's/\..*//')

    if [ "$base_rsem" != "$base_bam" ]; then
        echo "[ERROR] File mismatch at index $i:"
        echo "  RSEM: ${RSEM_FILES[$i]}"
        echo "  BAM:  ${BAM_FILES[$i]}"
        exit 1
    fi
done

echo "[INFO] RSEM and BAM files validated. Names and order match."

# ----------------------
# Loop over VCF inputs
# ----------------------
VCF_FILES=("$VCF_DIR"/*.vcf)
if [ ${#VCF_FILES[@]} -eq 0 ]; then
    echo "[ERROR] No VCF files found in $VCF_DIR"
    exit 1
fi

# ----------------------
# Run modules
# ----------------------
for vcf in "${VCF_FILES[@]}"; do
    echo "[INFO] Processing: $vcf (Project: $PROJECT)"

    echo "[INFO] Running cancer mutation calling..."
    bash modules/cancer_mutation_calling/step1_call_mutations.sh "$vcf" "$REF_DIR" "$OUTPUT_DIR"
    Rscript modules/cancer_mutation_calling/step2_process_mutations.R "$vcf" "$REF_DIR" "$OUTPUT_DIR"

    echo "[INFO] Running mycoplasma detection..."
    bash modules/mycoplasma_detection/detect_mycoplasma.sh "$vcf" "$REF_DIR" "$OUTPUT_DIR"

    echo "[INFO] Running PacNet..."
    Rscript modules/pacnet/run_pacnet.R "$vcf" "$REF_DIR" "$RSEM_DIR" "$OUTPUT_DIR"

    echo "[INFO] Running eKaryo..."
    Rscript modules/eKaryo/run_eKaryo.R "$vcf" "$REF_DIR" "$BAM_DIR" "$OUTPUT_DIR"

    echo "[INFO] Running outlier detection..."
    Rscript modules/outliers/detect_outliers.R "$vcf" "$REF_DIR" "$FASTQ_DIR" "$OUTPUT_DIR"
done

echo "[INFO] Pipeline completed successfully."

