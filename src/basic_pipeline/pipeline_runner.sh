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
    fastq_files=($(ls "$FASTQ_DIR"/"$sample"*.fastq* 2>/dev/null || true))

    if [[ -z "$bam_file" || -z "$vcf_file" || -z "$rsem_file" || ${#fastq_files[@]} -eq 0 ]]; then
        echo "[WARN] Missing files for sample $sample. Skipping..."
        continue
    fi

    sample_outdir="$OUTPUT_DIR/$sample"
    mkdir -p "$sample_outdir"

    # ----------------------
    # Run Mycoplasma detection
    # ----------------------
    echo "[INFO] Running Mycoplasma detection for $sample..."
    bash modules/mycoplasma_detection/detect_mycoplasma.sh \
        --bam "$bam_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/mycoplasma_${sample}.log" 2>&1

    # ----------------------
    # Run cancer mutation calling
    # ----------------------
    echo "[INFO] Running cancer mutation calling for $sample..."
    bash modules/cancer_mutation_calling/step1_call_mutations.sh \
        --vcf "$vcf_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/mutation_calling_step1_${sample}.log" 2>&1

    Rscript modules/cancer_mutation_calling/step2_process_mutations.R \
        --vcf "$vcf_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/mutation_calling_step2_${sample}.log" 2>&1

    # ----------------------
    # Run PacNet
    # ----------------------
    echo "[INFO] Running PacNet for $sample..."
    Rscript modules/pacnet/run_pacnet.R \
        --vcf "$vcf_file" --rsem "$rsem_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/pacnet_${sample}.log" 2>&1

    # ----------------------
    # Run eKaryo
    # ----------------------
    echo "[INFO] Running eKaryo for $sample..."
    Rscript modules/eKaryo/run_eKaryo.R \
        --bam "$bam_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/ekaryo_${sample}.log" 2>&1

    # ----------------------
    # Run outlier detection (PCA / pacnet scores)
    # ----------------------
    echo "[INFO] Running outlier detection for $sample..."
    Rscript modules/outliers/detect_outliers.R \
        --fastq "${fastq_files[@]}" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/outliers_${sample}.log" 2>&1

    echo "[INFO] Finished processing $sample. Results in $sample_outdir"
done

# ----------------------
# Compile HTML summary
# ----------------------
echo "[INFO] Generating HTML summary..."
Rscript modules/report_builder/generate_html_summary.R \
    --output_dir "$OUTPUT_DIR" --project "$PROJECT" \
    > "$LOG_DIR/html_summary.log" 2>&1

echo "[INFO] Pipeline completed successfully. Summary report: $OUTPUT_DIR/${PROJECT}_summary.html"
