#!/usr/bin/env bash
set -euo pipefail

# ----------------------
# Defaults
# ----------------------
PROJECT="Project [$(date '+%a %b %d %Y %H:%M')]"
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
    echo "  --project NAME       Project name (default: ${PROJECT} and time/date of pipeline run)"
    echo "  --output_dir PATH    Output directory (default: ./ipsc_qc_outputs)"
    echo ""
    echo "See requirements.txt and README.md for more details on arguments and inputs."
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
        --output_dir) OUTPUT_DIR="$2"; shift ;;
        -h|--help)    usage ;;
        *) echo "[ERROR] Unknown option: $1"; usage ;;
    esac
    shift
done

# ----------------------
# Validate required args
# ----------------------
if [[ -z "$VCF_DIR" || -z "$RSEM_DIR" || -z "$BAM_DIR" || -z "$FASTQ_DIR" ]]; then
    echo "[ERROR] Missing one or more required arguments."
    usage
fi

mkdir -p "$OUTPUT_DIR"
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$LOG_DIR"

# ----------------------
# Extract unique sample names from BAM files
# ----------------------
SAMPLES=($(ls "$BAM_DIR"/*.bam | xargs -n1 basename | sed 's/\..*//' | sort -u))
if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "[ERROR] No BAM files found in $BAM_DIR"
    exit 1
fi
echo "[INFO] Found ${#SAMPLES[@]} samples: ${SAMPLES[*]}"

# ----------------------
# Process each sample
# ----------------------
for sample in "${SAMPLES[@]}"; do
    echo "[INFO] Processing sample: $sample"

    bam_file=$(ls "$BAM_DIR"/"$sample"*.bam 2>/dev/null || true)
    vcf_file=$(ls "$VCF_DIR"/"$sample"*.vcf* 2>/dev/null || true)
    rsem_file=$(ls "$RSEM_DIR"/"$sample"* 2>/dev/null || true)

    # Determine base name before the last underscore
    base_name=$(basename "$sample")
    base_name="${base_name%_*}"  # removes everything after the last underscore

    # Collect FASTQ files for this sample
    fastq_files=($(ls "$FASTQ_DIR"/"$base_name"_*.fastq* 2>/dev/null || true))

    if [[ -z "$bam_file" || -z "$vcf_file" || -z "$rsem_file" || ${#fastq_files[@]} -eq 0 ]]; then
        echo "[WARN] Missing files for sample $sample. Skipping..."
        continue
    fi

    sample_outdir="$OUTPUT_DIR/$sample"
    mkdir -p "$sample_outdir"

    # ----------------------
    # Run PacNet
    # ----------------------
    echo "[INFO] Assessing pluripotency for $sample..."
    Rscript modules/PACNet/run_pacnet.R \
        --vcf "$vcf_file" --rsem "$rsem_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/pacnet_${sample}.log" 2>&1

    # ----------------------
    # Run cancer mutation calling
    # ----------------------
    echo "[INFO] Checking for COSMIC mutations in $sample..."
    bash modules/cancer_mutation_calling/filter_vcf_on_cosmic.sh \
        --vcf "$vcf_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/filtered_vcf_for_cosmic_${sample}.log" 2>&1

    Rscript modules/cancer_mutation_calling/cancer_mutation_mapping.R \
        --vcf "$vcf_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/cancer_mutation_mapping_${sample}.log" 2>&1

    # ----------------------
    # Run eKaryo
    # ----------------------
    echo "[INFO] Assessing chromosomal integrity for $sample..."
    Rscript modules/eSNPKaryotyping/run_eSNPKaryotyping.R \
        --bam "$bam_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/ekaryo_${sample}.log" 2>&1
    
    # ----------------------
    # Run Mycoplasma detection
    # ----------------------
    echo "[INFO] Detecting mycoplasma in $sample..."
    bash modules/mycoplasma_detection/detect_mycoplasma.sh \
        --bam "$bam_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/mycoplasma_${sample}.log" 2>&1

    # ----------------------
    # Run outlier detection (PCA / pacnet scores)
    # ----------------------
    echo "[INFO] Detecting outliers in $sample..."
    Rscript modules/outlier_detection/outlier_detection.R \
        --fastq "${fastq_files[@]}" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/outliers_${sample}.log" 2>&1

    # ----------------------
    # Compile HTML summary
    # ----------------------
    echo "[INFO] Generating summary for $sample..."
    Rscript modules/report_builder/generate_html_summary.R \
        --output_dir "$OUTPUT_DIR" --project "$PROJECT" \
        > "$LOG_DIR/html_summary.log" 2>&1

    echo "[INFO] Finished processing $sample. Results in $sample_outdir"
done

echo "[INFO] Pipeline completed. All results in $OUTPUT_DIR"
