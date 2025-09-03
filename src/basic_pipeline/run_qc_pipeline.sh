#!/usr/bin/env bash
set -euo pipefail

# -------------------------------
# Usage
# -------------------------------
usage() {
    echo "Usage: $0 -f FASTQ_DIR -o OUTPUT_DIR -p PROJECT_NAME -c COSMIC_TAR"
    echo
    echo "Arguments:"
    echo "  -f   Path to FASTQ directory (required)"
    echo "  -o   Path to output directory (default: ./outputs)"
    echo "  -p   Project name (default: 'RNAseq_Project')"
    echo "  -c   Path to Cosmic_CancerGeneCensus_Tsv_v101_GRCh37.tar (required)"
    exit 1
}

# -------------------------------
# Parse arguments
# -------------------------------
FASTQ_DIR=""
OUTPUT_DIR="./outputs"
PROJECT="RNAseq_Project"
COSMIC_TAR=""

while getopts "f:o:p:c:" opt; do
    case ${opt} in
        f) FASTQ_DIR=$OPTARG ;;
        o) OUTPUT_DIR=$OPTARG ;;
        p) PROJECT=$OPTARG ;;
        c) COSMIC_TAR=$OPTARG ;;
        *) usage ;;
    esac
done

if [[ -z "$FASTQ_DIR" ]] || [[ -z "$COSMIC_TAR" ]]; then
    echo "[ERROR] FASTQ directory and COSMIC tarball must be provided."
    usage
fi

# Normalize paths
FASTQ_DIR=$(realpath "$FASTQ_DIR")
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")
COSMIC_TAR=$(realpath "$COSMIC_TAR")

# -------------------------------
# Run Docker
# -------------------------------
IMAGE="rna_pipeline:latest"

mkdir -p "$OUTPUT_DIR"

docker run --rm -it \
    -v "$FASTQ_DIR":/workspace/fastq_dir \
    -v "$OUTPUT_DIR":/workspace/outputs \
    -v "$COSMIC_TAR":/opt/ref/Cosmic_CancerGeneCensus_Tsv_v101_GRCh37.tar \
    -e PROJECT="$PROJECT" \
    -e REF_DIR=/opt/ref \
    "$IMAGE" \
    bash pipeline_runner.sh \
        --fastq_dir /workspace/fastq_dir \
        --output_dir /workspace/outputs \
        --project "$PROJECT"
