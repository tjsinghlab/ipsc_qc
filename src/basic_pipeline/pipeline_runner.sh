#!/usr/bin/env bash
set -euo pipefail

# ===========================================
# Configuration
# ===========================================
PROJECT="Project [$(date '+%a %b %d %Y %H:%M')]"
REF_DIR="/ref" #mounted to image by user
FASTQ_DIR="/data" #mounted to image by user
OUTPUT_DIR="/output" #mounted to image by user
COSMIC_DIR="/cosmic" #mounted to image by user

# Paths to WDL runners
PY_RUNNER1="/pipeline/modules/preprocessing/wdlplay/warp-pipelines/bulk_RNAseq_preprocess/run_wdl.py" 
PY_RUNNER2="/pipeline/modules/preprocessing/wdlplay/warp-pipelines/GATK_variant_calling/gatk4-rna-germline-calling_run.py" 

# ===========================================
# Usage
# ===========================================
usage() {
    echo "Usage: pipeline_runner.sh --fastq_dir PATH --output_dir PATH --ref_dir PATH [--cosmic_dir PATH] [--project NAME]"
    echo ""
    echo "Required arguments:"
    echo "  --fastq_dir PATH     Directory containing FASTQ files"
    echo "  --output_dir PATH    Output directory for pipeline results"
    echo "  --ref_dir PATH       Directory containing reference files (genome, annotations, etc.)"
    echo ""
    echo "Optional arguments:"
    echo "  --cosmic_dir PATH    Directory containing COSMIC references"
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
        --ref_dir)    REF_DIR="$2"; shift ;;
        --cosmic_dir) COSMIC_DIR="$2"; shift ;;
        -h|--help)    usage ;;
        *) echo "[ERROR] Unknown option: $1"; usage ;;
    esac
    shift
done

# ===========================================
# Input check
# ===========================================
# if [[ -z "$FASTQ_DIR" || -z "$OUTPUT_DIR" || -z "$REF_DIR" ]]; then
#     echo "[ERROR] --fastq_dir, --output_dir, and --ref_dir must be provided"
#     usage
# fi

mkdir -p "$OUTPUT_DIR"
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$LOG_DIR"

# ------------------------------------------------
# Reference resolution
# ------------------------------------------------
echo "[STEP] Resolving reference files..."

# Helper: check for a file/dir in REF_DIR or /ref
resolve_ref() {
    local name="$1"
    local path_rel="$2"
    local fallback="$3"  # optional: if not found in either, use fallback (download/build/etc.)

    if [[ -e "$REF_DIR/$path_rel" ]]; then
        echo "[INFO] Found $name in user-provided REF_DIR."
        echo "$REF_DIR/$path_rel"
    elif [[ -e "/ref/$path_rel" ]]; then
        echo "[INFO] Using $name from Docker image."
        echo "/ref/$path_rel"
    else
        if [[ -n "$fallback" ]]; then
            echo "[WARN] $name not found. Running fallback step..."
            eval "$fallback"
            echo "$REF_DIR/$path_rel"
        else
            echo "[ERROR] Missing required reference: $name ($path_rel)"
            exit 1
        fi
    fi
}

# STAR index (special case: build if missing)
STAR_INDEX=$(resolve_ref "STAR index" "star_index_oh75" \
"mkdir -p \"$REF_DIR/star_index_oh75\" && \
singularity exec /ref/images/gtex_rnaseq_v10.sif \
    STAR --runMode genomeGenerate \
         --genomeDir \"$REF_DIR/star_index_oh75\" \
         --genomeFastaFiles \"$REF_DIR/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta\" \
         --sjdbGTFfile \"$REF_DIR/gencode.v39.GRCh38.annotation.gtf\" \
         --sjdbOverhang 75 \
         --runThreadN 4")

# RSEM reference (special case: build if missing)
RSEM_REF=$(resolve_ref "RSEM reference" "rsem_reference" \
"mkdir -p \"$REF_DIR/rsem_reference\" && \
singularity exec /ref/images/gtex_rnaseq_v10.sif \
    rsem-prepare-reference \
        \"$REF_DIR/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta\" \
        \"$REF_DIR/rsem_reference/rsem_reference\" \
        --gtf \"$REF_DIR/gencode.v39.GRCh38.annotation.gtf\" \
        --num-threads 4")

# GTF annotation
GENES_GTF=$(resolve_ref "GTF annotation" "gencode.v39.GRCh38.genes.collapsed_only.gtf")

# Reference FASTA + index + dict
REF_FASTA=$(resolve_ref "Reference FASTA" "Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta")
REF_FASTA_FAI=$(resolve_ref "FASTA index" "Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta.fai")
REF_DICT=$(resolve_ref "FASTA dict" "Homo_sapiens_assembly38_noALT_noHLA_noDecoy.dict")

# dbSNP
DBSNP_VCF=$(resolve_ref "dbSNP VCF" "Homo_sapiens_assembly38.dbsnp138.vcf")
DBSNP_IDX=$(resolve_ref "dbSNP index" "Homo_sapiens_assembly38.dbsnp138.vcf.idx")

# Known indels
KNOWN_VCF1=$(resolve_ref "Known indels 1" "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
KNOWN_VCF2=$(resolve_ref "Known indels 2" "Homo_sapiens_assembly38.known_indels.vcf.gz")
KNOWN_VCF1_IDX=$(resolve_ref "Known indels index 1" "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi")
KNOWN_VCF2_IDX=$(resolve_ref "Known indels index 2" "Homo_sapiens_assembly38.known_indels.vcf.gz.tbi")

# GATK image
GATK_SIF=$(resolve_ref "GATK 4.6.1.0 SIF" "images/gatk_4.6.1.0.sif")

# Gene list
GENE_LIST=$(resolve_ref "Gene list" "genes.txt")

# Mycoplasma genome (fallback: download if missing)
MYCO_REF=$(resolve_ref "Mycoplasma genome" "GCF_000027325.1_ASM2732v1_genomic.fna.gz" \
"wget -O \"$REF_DIR/GCF_000027325.1_ASM2732v1_genomic.fna.gz\" \
 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz")

# PACNet training data (fallback: AWS S3 download if missing)
PACNET_EXP=$(resolve_ref "PACNet expTrain" "Hs_expTrain_Jun-20-2017.rda" \
"aws s3 cp s3://cellnet-rnaseq/ref/cnproc/HS/Hs_expTrain_Jun-20-2017.rda \"$REF_DIR/\" --no-sign-request")
PACNET_ST=$(resolve_ref "PACNet stTrain" "Hs_stTrain_Jun-20-2017.rda" \
"aws s3 cp s3://cellnet-rnaseq/ref/cnproc/HS/Hs_stTrain_Jun-20-2017.rda \"$REF_DIR/\" --no-sign-request")

# UCSC SNP142 common SNPs directory
SNP142_DIR=$(resolve_ref "dbSNP142 directory" "chr")

# NCBI reference genome FASTA
NCBI_FASTA=$(resolve_ref "NCBI GRCh38 FASTA" "GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna" \
"wget -O \"$REF_DIR/GCF_000001405.26_GRCh38_genomic.fna\" \
 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38.p13")

# ===========================================
# Identify samples
# ===========================================
SAMPLES=()
shopt -s nullglob

for fq in "${FASTQ_DIR}"/*_R1*.fastq.gz; do
    base=$(basename "$fq")
    # Remove everything from _R1 onward
    sample="${base%%_R1*}"
    SAMPLES+=("$sample")
done

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "[ERROR] No FASTQ files found in $FASTQ_DIR"
    exit 1
fi

echo "[INFO] Found ${#SAMPLES[@]} samples: ${SAMPLES[*]}"


# ===========================================
# Loop over samples
# ===========================================
for sample in "${SAMPLES[@]}"; do
    echo "[INFO] Processing sample: $sample"

    sample_outdir="$OUTPUT_DIR/$sample"
    mkdir -p "$sample_outdir"

    fq1="${FASTQ_DIR}/${sample}_R1_001.fastq.gz"
    fq2="${FASTQ_DIR}/${sample}_R2_001.fastq.gz"

    # ------------------------------------------------
    # STEP 1: Preprocessing
    # ------------------------------------------------
    echo "[STEP] Running bulk RNAseq pre-processing for $sample..."
    ( python3 "$PY_RUNNER1" \--fastq1 "$fq1" --fastq2 "$fq2" --output_dir "$sample_outdir" --sample "$sample" > "$LOG_DIR/${sample}_bulk_preprocess.log" 2>&1 )
    bam_file="$sample_outdir/${sample}.bam"
    bai_file="${bam_file}.bai"
    rsem_file="$sample_outdir/RSEM_outputs/${sample}.rsem.genes.results.gz"

    # ------------------------------------------------
    # STEP 2: Variant Calling
    # ------------------------------------------------
    echo "[STEP] Calling germline variants for $sample..."
    ( python3 "$PY_RUNNER2" --bam "$bam_file" --bai "$bai_file" --outdir "$sample_outdir" --sample "$sample" > "$LOG_DIR/${sample}_wdl2.log" 2>&1 )
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
