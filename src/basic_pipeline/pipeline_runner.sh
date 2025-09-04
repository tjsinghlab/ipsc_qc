#!/usr/bin/env bash
set -euo pipefail

# ===========================================
# Configuration
# ===========================================
PROJECT="Project [$(date '+%a %b %d %Y %H:%M')]"
REF_DIR="/ref"   # fixed inside Docker image
FASTQ_DIR=""
OUTPUT_DIR=""
COSMIC_DIR=""

# Paths to WDL runners
PY_RUNNER1="/pipeline/modules/preprocessing/bulk_RNAseq_preprocess/run_wdl.py" 
PY_RUNNER2="/pipeline/modules/preprocessing/GATK_variant_calling/gatk4-rna-germline-calling_run.py" 

# ===========================================
# Usage
# ===========================================
usage() {
    echo "Usage: pipeline_runner.sh --fastq_dir PATH --output_dir PATH [--cosmic_dir PATH]"
    echo ""
    echo "Required arguments:"
    echo "  --fastq_dir PATH     Directory containing FASTQ files"
    echo "  --output_dir PATH    Output directory for pipeline results"
    echo ""
    echo "Optional arguments:"
    echo "  --cosmic_dir PATH    Directory containing COSMIC references (if provided, runs mutation calling)"
    echo "  --project NAME       Project name (default: ${PROJECT})"
    echo ""
    exit 1
}

# ===========================================
# Arguments
# ===========================================
while [[ $# -gt 0 ]]; do
    case $1 in
        --project)    PROJECT="$2"; shift ;;
        --fastq_dir)  FASTQ_DIR="$2"; shift ;;
        --output_dir) OUTPUT_DIR="$2"; shift ;;
        --cosmic_dir) COSMIC_DIR="$2"; shift ;;
        -h|--help)    usage ;;
        *) echo "[ERROR] Unknown option: $1"; usage ;;
    esac
    shift
done

# ===========================================
# Input check
# ===========================================
if [[ -z "$FASTQ_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "[ERROR] --fastq_dir and --output_dir must be provided"
    usage
fi

mkdir -p "$OUTPUT_DIR"
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$LOG_DIR"

# ===========================================
# Identify samples
# ===========================================
SAMPLES=()
for fq in "${FASTQ_DIR}"/*_R1.fastq.gz; do
    base=$(basename "$fq")
    sample="${base%_R1.fastq.gz}"
    SAMPLES+=("$sample")
done
echo "[INFO] Found ${#SAMPLES[@]} samples: ${SAMPLES[*]}"

# ===========================================
# Loop over samples
# ===========================================
for sample in "${SAMPLES[@]}"; do
    echo "[INFO] Processing sample: $sample"

    sample_outdir="$OUTPUT_DIR/$sample"
    mkdir -p "$sample_outdir"

    fq1="${FASTQ_DIR}/${sample}_R1.fastq.gz"
    fq2="${FASTQ_DIR}/${sample}_R2.fastq.gz"

    # ------------------------------------------------
    # STEP 1: Preprocessing
    # ------------------------------------------------
    echo "[STEP] Running bulk RNAseq pre-processing for $sample..."
    ( python3 "$PY_RUNNER1" --fastq1 "$fq1" --fastq2 "$fq2" --outdir "$sample_outdir" > "$LOG_DIR/${sample}_wdl1.log" 2>&1 )
    bam_file="$sample_outdir/${sample}.bam"
    bai_file="${bam_file}.bai"
    rsem_file="$sample_outdir/RSEM_outputs/${sample}.rsem.genes.results.gz"

    # ------------------------------------------------
    # STEP 2: Variant Calling
    # ------------------------------------------------
    echo "[STEP] Calling germline variants for $sample..."
    ( python3 "$PY_RUNNER2" --bam "$bam_file" --bai "$bai_file" --outdir "$sample_outdir" > "$LOG_DIR/${sample}_wdl2.log" 2>&1 )
    vcf_file="$sample_outdir/${sample}.vcf"
    
    if [[ ! -f "$vcf_file" ]]; then
        echo "[ERROR] VCF file not found for $sample. Skipping..."
        continue
    fi

    # ------------------------------------------------
    # STEP 3: QC Modules
    # ------------------------------------------------
    echo "[STEP] Running PACNet for $sample..."
    Rscript /pipeline/modules/PACNet/run_pacnet.R \
        --rsem "$rsem_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" --sample "$sample" --project "$PROJECT" \
        > "$LOG_DIR/pacnet_${sample}.log" 2>&1

    if [[ -n "$COSMIC_DIR" ]]; then
        echo "[STEP] Running COSMIC mutation calling for $sample..."
        bash /pipeline/modules/cancer_mutation_calling/filter_vcf_on_cosmic.sh \
            --vcf "$vcf_file" --ref_dir "$REF_DIR" --cosmic_dir "$COSMIC_DIR" --output_dir "$sample_outdir" --sample "$sample" \
            > "$LOG_DIR/filtered_vcf_for_cosmic_${sample}.log" 2>&1

        Rscript /pipeline/modules/cancer_mutation_calling/cancer_mutation_mapping.R \
            --vcf "$vcf_file" --ref_dir "$REF_DIR" --cosmic_dir "$COSMIC_DIR" --output_dir "$sample_outdir" --sample "$sample" \
            > "$LOG_DIR/cancer_mutation_mapping_${sample}.log" 2>&1
    else
        echo "[INFO] Skipping COSMIC mutation calling (no --cosmic_dir provided)."
    fi

    echo "[STEP] Running eSNPKaryotyping for $sample..."
    Rscript /pipeline/modules/eSNPKaryotyping/run_eSNPKaryotyping.R \
        --vcf "$vcf_file" --bam "$bam_file" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" --sample "$sample" \
        > "$LOG_DIR/ekaryo_${sample}.log" 2>&1

    echo "[STEP] Running Mycoplasma detection for $sample..."
    bash /pipeline/modules/mycoplasma_detection/detect_mycoplasma.sh \
        --fastq "$fq1" "$fq2" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" --sample "$sample" \
        > "$LOG_DIR/mycoplasma_${sample}.log" 2>&1

    echo "[STEP] Running Outlier detection for $sample..."
    Rscript /pipeline/modules/outlier_detection/outlier_detection.R \
        --output_dir "$sample_outdir" --sample "$sample" --project "$PROJECT" \
        > "$LOG_DIR/outliers_${sample}.log" 2>&1

    # ------------------------------------------------
    # STEP 4: Compile HTML summary
    # ------------------------------------------------
    echo "[STEP] Generating HTML summary for $sample..."
    Rscript /pipeline/modules/report_builder/generate_html_summary.R \
        --output_dir "$OUTPUT_DIR" --project "$PROJECT" --sample "$sample" \
        > "$LOG_DIR/html_summary_${sample}.log" 2>&1

    echo "[DONE] Finished processing $sample. Results in $sample_outdir"
done

echo "[INFO] Pipeline completed. All results written to $OUTPUT_DIR"
