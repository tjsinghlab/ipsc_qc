#!/usr/bin/env bash
set -euo pipefail

# ===========================================
# CONFIG / DEFAULTS
# ===========================================
PROJECT="Project [$(date '+%a %b %d %Y %H:%M')]"
OUTPUT_DIR="$(pwd)/ipsc_qc_outputs"
REF_DIR="$(pwd)/ref"
FASTQ_DIR=""

# Paths to WDL runners
PY_RUNNER1="runner_wdl_stage1.py"   # generates BAMs
PY_RUNNER2="runner_wdl_stage2.py"   # generates VCF + RSEM

# ===========================================
# USAGE
# ===========================================
usage() {
    echo "Usage: pipeline_runner.sh [options]"
    echo ""
    echo "Required arguments:"
    echo "  --fastq_dir PATH     Directory containing FASTQ files"
    echo ""
    echo "Optional arguments:"
    echo "  --ref_dir PATH       Reference directory (default: ./ref)"
    echo "  --project NAME       Project name (default: ${PROJECT})"
    echo "  --output_dir PATH    Output directory (default: ./ipsc_qc_outputs)"
    echo ""
    exit 1
}

# ===========================================
# ARG PARSING
# ===========================================
while [[ $# -gt 0 ]]; do
    case $1 in
        --project)    PROJECT="$2"; shift ;;
        --ref_dir)    REF_DIR="$2"; shift ;;
        --fastq_dir)  FASTQ_DIR="$2"; shift ;;
        --output_dir) OUTPUT_DIR="$2"; shift ;;
        -h|--help)    usage ;;
        *) echo "[ERROR] Unknown option: $1"; usage ;;
    esac
    shift
done

# ===========================================
# VALIDATION
# ===========================================
if [[ -z "$FASTQ_DIR" ]]; then
    echo "[ERROR] Missing required argument: --fastq_dir"
    usage
fi

mkdir -p "$OUTPUT_DIR"
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$LOG_DIR"

# ===========================================
# COLLECT SAMPLES (based on FASTQ _R1/_R2 pairs)
# ===========================================
SAMPLES=()
for fq in "${FASTQ_DIR}"/*_R1.fastq.gz; do
    base=$(basename "$fq")
    sample="${base%_R1.fastq.gz}"   # strip suffix
    SAMPLES+=("$sample")
done
echo "[INFO] Found ${#SAMPLES[@]} samples: ${SAMPLES[*]}"

# ===========================================
# MAIN LOOP
# ===========================================
for sample in "${SAMPLES[@]}"; do
    echo "==========================================="
    echo "[INFO] Processing sample: $sample"
    echo "==========================================="

    sample_outdir="$OUTPUT_DIR/$sample"
    mkdir -p "$sample_outdir"

    fq1="${FASTQ_DIR}/${sample}_R1.fastq.gz"
    fq2="${FASTQ_DIR}/${sample}_R2.fastq.gz"

    # ------------------------------------------------
    # STEP 1: Run WDL stage 1 (BAM generation)
    # ------------------------------------------------
    echo "[STEP] Running WDL stage 1 for $sample..."
    ( python3 "$PY_RUNNER1" --fastq1 "$fq1" --fastq2 "$fq2" --outdir "$sample_outdir" > "$LOG_DIR/${sample}_wdl1.log" 2>&1 )
    bam_file="$sample_outdir/${sample}.bam"
    bai_file="${bam_file}.bai"

    if [[ ! -f "$bam_file" ]]; then
        echo "[ERROR] BAM not found for $sample. Skipping..."
        continue
    fi

    # ------------------------------------------------
    # STEP 2: Run WDL stage 2 (VCF + RSEM generation)
    # ------------------------------------------------
    echo "[STEP] Running WDL stage 2 for $sample..."
    ( python3 "$PY_RUNNER2" --bam "$bam_file" --bai "$bai_file" --outdir "$sample_outdir" > "$LOG_DIR/${sample}_wdl2.log" 2>&1 )
    vcf_file="$sample_outdir/${sample}.vcf"
    rsem_file="$sample_outdir/${sample}.rsem"

    if [[ ! -f "$vcf_file" || ! -f "$rsem_file" ]]; then
        echo "[ERROR] VCF or RSEM genes file not found for $sample. Skipping..."
        continue
    fi

    # ------------------------------------------------
    # STEP 3: Run downstream modules
    # ------------------------------------------------
    echo "[STEP] Running PACNet for $sample..."
    Rscript modules/PACNet/run_pacnet.R \
        --vcf "$vcf_file" --rsem "$rsem_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/pacnet_${sample}.log" 2>&1

    echo "[STEP] Running COSMIC mutation calling for $sample..."
    bash modules/cancer_mutation_calling/filter_vcf_on_cosmic.sh \
        --vcf "$vcf_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/filtered_vcf_for_cosmic_${sample}.log" 2>&1

    Rscript modules/cancer_mutation_calling/cancer_mutation_mapping.R \
        --vcf "$vcf_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/cancer_mutation_mapping_${sample}.log" 2>&1

    echo "[STEP] Running eSNPKaryotyping for $sample..."
    Rscript modules/eSNPKaryotyping/run_eSNPKaryotyping.R \
        --bam "$bam_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/ekaryo_${sample}.log" 2>&1

    echo "[STEP] Running Mycoplasma detection for $sample..."
    bash modules/mycoplasma_detection/detect_mycoplasma.sh \
        --bam "$bam_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/mycoplasma_${sample}.log" 2>&1

    echo "[STEP] Running Outlier detection for $sample..."
    Rscript modules/outlier_detection/outlier_detection.R \
        --fastq "$fq1" "$fq2" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" \
        > "$LOG_DIR/outliers_${sample}.log" 2>&1

    # ------------------------------------------------
    # STEP 4: Compile HTML summary
    # ------------------------------------------------
    echo "[STEP] Generating HTML summary for $sample..."
    Rscript modules/report_builder/generate_html_summary.R \
        --output_dir "$OUTPUT_DIR" --project "$PROJECT" \
        > "$LOG_DIR/html_summary_${sample}.log" 2>&1

    echo "[DONE] Finished processing $sample. Results in $sample_outdir"
done

echo "[INFO] Pipeline completed. All results in $OUTPUT_DIR"
