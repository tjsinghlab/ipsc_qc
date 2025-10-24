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

mkdir -p "$OUTPUT_DIR"
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$LOG_DIR"

# ------------------------------------------------
# Reference resolution
# ------------------------------------------------
echo "[STEP] Resolving reference files..."

resolve_ref() {
    local name=$1 path=$2 cmd=$3
    local f="$REF_DIR/$path"
    if [ -s "$f" ]; then
        echo "[INFO] Found $name at $f"
    else
        echo "[WARN] Missing $name, fetching..."
        eval "$cmd"
        [ -s "$f" ] || { echo "[ERROR] Failed to fetch $name"; exit 1; }
    fi
    echo "$f"
}

# ----------------------
# Define resources
# ----------------------
resources=(
  

  "ANNOTATION_GTF|gencode.v39.GRCh38.annotation.gtf|wget -O - http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz | gunzip -c > $REF_DIR/gencode.v39.GRCh38.annotation.gtf"

# Requester pays enabled for these resources and cannot be downloaded without credentials
#   "REF_FASTA|Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta|wget -O $REF_DIR/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta https://storage.googleapis.com/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta"
#   "REF_FASTA_FAI|Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta.fai|wget -O $REF_DIR/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta.fai https://storage.googleapis.com/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta.fai"
#   "REF_DICT|Homo_sapiens_assembly38_noALT_noHLA_noDecoy.dict|wget -O $REF_DIR/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.dict https://storage.googleapis.com/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.dict"
#   "GENES_GTF|gencode.v39.GRCh38.genes.collapsed_only.gtf|wget -O - http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.GRCh38.genes.collapsed_only.gtf.gz | gunzip -c > $REF_DIR/gencode.v39.GRCh38.genes.collapsed_only.gtf"

  "DBSNP_VCF|Homo_sapiens_assembly38.dbsnp138.vcf|wget -O $REF_DIR/Homo_sapiens_assembly38.dbsnp138.vcf https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
  "DBSNP_IDX|Homo_sapiens_assembly38.dbsnp138.vcf.idx|wget -O $REF_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.idx https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"

  "KNOWN_VCF1|Mills_and_1000G_gold_standard.indels.hg38.vcf.gz|wget -O $REF_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  "KNOWN_VCF2|Homo_sapiens_assembly38.known_indels.vcf.gz|wget -O $REF_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
  "KNOWN_VCF1_IDX|Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi|wget -O $REF_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
  "KNOWN_VCF2_IDX|Homo_sapiens_assembly38.known_indels.vcf.gz.tbi|wget -O $REF_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"

  "GATK_SIF|gatk_4.6.1.0.sif|singularity pull $REF_DIR/gatk_4.6.1.0.sif library://biocontainers/gatk4:4.6.1.0--py310hdfd78af_0"
  "GTEX_SIF|gtex_rnaseq_v10.sif|singularity pull $REF_DIR/gtex_rnaseq_v10.sif docker://broadinstitute/gtex_rnaseq:V10"

  "GENE_LIST|genes.txt|wget -O $REF_DIR/genes.txt https://raw.githubusercontent.com/tjsinghlab/ipsc_qc/main/src/basic_pipeline/ref_files/genes.txt"

  "MYCO_ORALE_REF|GCF_000420105.1_ASM42010v1_genomic.fna.gz|wget -O $REF_DIR/GCF_000420105.1_ASM42010v1_genomic.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/420/105/GCF_000420105.1_ASM42010v1/GCF_000420105.1_ASM42010v1_genomic.fna.gz"
  "MYCO_FERMENTANS_REF|GCF_003704055.1_ASM370405v1_genomic.fna.gz|wget -O $REF_DIR/GCF_003704055.1_ASM370405v1_genomic.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/704/055/GCF_003704055.1_ASM370405v1/GCF_003704055.1_ASM370405v1_genomic.fna.gz"

  "PACNET_EXP|Hs_expTrain_Jun-20-2017.rda|aws s3 cp s3://cellnet-rnaseq/ref/cnproc/HS/Hs_expTrain_Jun-20-2017.rda $REF_DIR/ --no-sign-request"
  "PACNET_ST|Hs_stTrain_Jun-20-2017.rda|aws s3 cp s3://cellnet-rnaseq/ref/cnproc/HS/Hs_stTrain_Jun-20-2017.rda $REF_DIR/ --no-sign-request"

  "NCBI_FASTA|GCF_000001405.26_GRCh38_genomic.fna|wget -O - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz | gunzip -c > $REF_DIR/GCF_000001405.26_GRCh38_genomic.fna"

)

# ----------------------
# Resolve in one loop
# ----------------------
for res in "${resources[@]}"; do
  IFS="|" read -r var file cmd <<< "$res"
  declare "$var=$(resolve_ref "$var" "$file" "$cmd")"
done

# Special cases: STAR + RSEM need building
STAR_INDEX=$(resolve_ref "STAR index" "star_index_oh75" \
"mkdir -p $REF_DIR/star_index_oh75 && \
 singularity exec $REF_DIR/gtex_rnaseq_v10.sif STAR --runMode genomeGenerate \
   --genomeDir $REF_DIR/star_index_oh75 \
   --genomeFastaFiles $REF_DIR/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
   --sjdbGTFfile $REF_DIR/gencode.v39.GRCh38.annotation.gtf \
   --sjdbOverhang 75 --runThreadN 4")

RSEM_REF=$(resolve_ref "RSEM reference" "rsem_reference" \
"mkdir -p $REF_DIR/rsem_reference && \
 singularity exec $REF_DIR/gtex_rnaseq_v10.sif rsem-prepare-reference \
   $REF_DIR/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
   $REF_DIR/rsem_reference/rsem_reference \
   --gtf $REF_DIR/gencode.v39.GRCh38.annotation.gtf \
   --num-threads 4")

# UCSC SNP142 common SNPs directory
#SNP142_DIR=$(resolve_ref "dbSNP142 directory" "chr") #how to get from UCSC without using their interactive web table?? may need to be downloaded by user beforehand and included in provided /ref directory
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
# ===========================================
# Phase 1: Run WDL pipelines for all samples
# ===========================================
# for sample in "${SAMPLES[@]}"; do
#     echo "[INFO] Running WDLs for sample: $sample"
#     sample_outdir="$OUTPUT_DIR/$sample"
#     mkdir -p "$sample_outdir"

#     fq1="${FASTQ_DIR}/${sample}_R1_001.fastq.gz"
#     fq2="${FASTQ_DIR}/${sample}_R2_001.fastq.gz"

#     bam_file="$sample_outdir/Mark_duplicates_outputs/${sample}.Aligned.sortedByCoord.out.md.bam"

#     bai_file="${bam_file}.bai"
#     rsem_file="$sample_outdir/RSEM_outputs/${sample}.rsem.genes.results.gz"

#     echo "Looking for $bam_file and $rsem_file"
#     # Step 1: Preprocessing
#     if [ -f "$bam_file" ] && [ -f "$rsem_file" ]; then
#         echo "[SKIP] $bam_file and $rsem_file already present. Skipping preprocessing for $sample"
#     else
#         echo "[RUN] Running preprocessing for $sample..."
#         python3 "$PY_RUNNER1" \
#             --fastq1 "$fq1" \
#             --fastq2 "$fq2" \
#             --output_dir "$sample_outdir" \
#             --sample "$sample" \
#             > "$LOG_DIR/${sample}_bulk_preprocess.log" 2>&1
#     fi

#     # Step 2: Variant Calling
#     if [ -d "$sample_outdir/variant_calling" ] && \
#     compgen -G "$sample_outdir/variant_calling/*.vcf.gz" > /dev/null; then
#         echo "[SKIP] Variant calling outputs already present at $sample_outdir/variant_calling/. Skipping variant calling for $sample"
#     else
#         echo "[RUN] Running variant calling for $sample..."
#         python3 "$PY_RUNNER2" \
#             --output_dir "$sample_outdir" \
#             --sample "$sample" \
#             > "$LOG_DIR/${sample}_wdl2.log" 2>&1
#     fi
# done

#!/usr/bin/env bash
set -euo pipefail

# ===========================================
# Function: run_sample
# ===========================================
run_sample() {
    local sample="$1"

    echo "[INFO] Beginning processing for $sample"
    local sample_outdir="$OUTPUT_DIR/$sample"
    mkdir -p "$sample_outdir"

    local fq1="${FASTQ_DIR}/${sample}_R1_001.fastq.gz"
    local fq2="${FASTQ_DIR}/${sample}_R2_001.fastq.gz"
    local bam_file="$sample_outdir/Mark_duplicates_outputs/${sample}.Aligned.sortedByCoord.out.md.bam"
    local rsem_file="$sample_outdir/RSEM_outputs/${sample}.rsem.genes.results.gz"

    # Optional: guard against memory overrun (soft cap)
    if [[ -n "${USE_MEM_MB:-}" ]]; then
        ulimit -v $(( USE_MEM_MB * 1024 )) 2>/dev/null || true
    fi

    # ----------------------
    # Step 1: Preprocessing
    # ----------------------
    if [ -s "$bam_file" ] && [ -s "$rsem_file" ]; then
        echo "[SKIP] Preprocessing already complete for $sample"
    else
        echo "[RUN] Running preprocessing for $sample..."
        python3 "$PY_RUNNER1" \
            --fastq1 "$fq1" \
            --fastq2 "$fq2" \
            --output_dir "$sample_outdir" \
            --sample "$sample" \
            > "$LOG_DIR/${sample}_bulk_preprocess.log" 2>&1

        # Wait/check loop — ensure both files exist & nonzero
        echo "[WAIT] Checking preprocessing completion for $sample..."
        for i in {1..30}; do
            if [ -s "$bam_file" ] && [ -s "$rsem_file" ]; then
                echo "[OK] Preprocessing complete for $sample"
                break
            fi
            sleep 10
            if [ $i -eq 30 ]; then
                echo "[ERROR] Preprocessing failed or incomplete for $sample after waiting 5 minutes."
                return 1
            fi
        done
    fi

    # Delete broken 'data' symlink silently
    find . -xtype l -name data -delete 2>/dev/null || true

    # ----------------------
    # Step 2: Variant Calling
    # ----------------------
    if [ -d "$sample_outdir/variant_calling" ] && \
       compgen -G "$sample_outdir/variant_calling/*.vcf.gz" > /dev/null; then
        echo "[SKIP] Variant calling already complete for $sample"
    else
        echo "[RUN] Running variant calling for $sample..."
        python3 "$PY_RUNNER2" \
            --output_dir "$sample_outdir" \
            --sample "$sample" \
            > "$LOG_DIR/${sample}_wdl2.log" 2>&1
    fi

    echo "[DONE] Preprocessing and variant calling for sample $sample completed."
}

export -f run_sample
export OUTPUT_DIR FASTQ_DIR LOG_DIR PY_RUNNER1 PY_RUNNER2 REF_DIR COSMIC_DIR

# ===========================================
# Dynamic Resource Allocation
# ===========================================

# Detect total cores and memory
TOTAL_CORES=$(nproc)
TOTAL_MEM_MB=$(awk '/MemTotal/ {print int($2/1024)}' /proc/meminfo)

# Allocate 85% of both
USE_CORES=$(awk -v total="$TOTAL_CORES" 'BEGIN {printf "%d", total * 0.85}')
USE_MEM_MB=$(awk -v total="$TOTAL_MEM_MB" 'BEGIN {printf "%d", total * 0.85}')

# Estimate memory required per job (in MB)
MEM_PER_JOB_MB=50000  # 50GB per sample (adjust as needed)

# Compute max jobs limited by memory and cores
MAX_JOBS_BY_MEM=$(( USE_MEM_MB / MEM_PER_JOB_MB ))
MAX_JOBS=$(( USE_CORES < MAX_JOBS_BY_MEM ? USE_CORES : MAX_JOBS_BY_MEM ))
MAX_JOBS=$(( MAX_JOBS > 0 ? MAX_JOBS : 1 ))

echo "[INFO] Detected $TOTAL_CORES cores and ${TOTAL_MEM_MB}MB total memory."
echo "[INFO] Using ${USE_CORES} cores and ${USE_MEM_MB}MB (≈85%) for processing."
echo "[INFO] Estimated ${MEM_PER_JOB_MB}MB per sample → running up to $MAX_JOBS samples in parallel."

export USE_MEM_MB MAX_JOBS

# ===========================================
# Phase 1: Parallel Preprocessing + Variant Calling
# ===========================================
printf '%s\n' "${SAMPLES[@]}" | xargs -n1 -P "$MAX_JOBS" bash -c 'run_sample "$@"' _

# Wait for all background jobs to complete
wait

# ===========================================
# Verify All Preprocessing Complete
# ===========================================
echo "[STEP] Verifying all samples completed preprocessing..."

incomplete_samples=()
for sample in "${SAMPLES[@]}"; do
    bam="$OUTPUT_DIR/$sample/Mark_duplicates_outputs/${sample}.Aligned.sortedByCoord.out.md.bam"
    rsem="$OUTPUT_DIR/$sample/RSEM_outputs/${sample}.rsem.genes.results.gz"
    if [ ! -s "$bam" ] || [ ! -s "$rsem" ]; then
        incomplete_samples+=("$sample")
    fi
done

if [ ${#incomplete_samples[@]} -gt 0 ]; then
    echo "[WARN] The following samples did not complete preprocessing:"
    printf ' - %s\n' "${incomplete_samples[@]}"
    echo "[WARN] Skipping PACNet and downstream steps."
    exit 1
fi

# ===========================================
# Phase 2: Run PACNet + Downstream Scripts
# ===========================================
echo "[STEP] Running PACNet..."
Rscript /pipeline/modules/PACNet/run_pacnet.R \
    --ref_dir "$REF_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --project "$PROJECT" \
    > "$LOG_DIR/pacnet_all.log" 2>&1

echo "[STEP] Outlier detection..."
if [[ ${#SAMPLES[@]} -ge 4 ]]; then
    Rscript /pipeline/modules/outlier_detection/outlier_detection.R \
        --output_dir "$OUTPUT_DIR" --project "$PROJECT" \
        > "$LOG_DIR/outliers_all_samples.log" 2>&1
else
    echo "[SKIP] Outlier detection skipped: only ${#SAMPLES[@]} sample(s) found (minimum 4 required)."
fi


# # Now loop again for per-sample downstream modules
# echo "[STEP] Running per-sample downstream analyses..."
# for sample in "${SAMPLES[@]}"; do
#     sample_outdir="$OUTPUT_DIR/$sample"
#     echo "[INFO] Finding fastq files..."
#     fq1="${FASTQ_DIR}/${sample}_R1_001.fastq.gz"
#     fq2="${FASTQ_DIR}/${sample}_R2_001.fastq.gz"

#     echo "[STEP] Mycoplasma detection for $sample..."
#     bash /pipeline/modules/mycoplasma_detection/detect_mycoplasma.sh \
#         --fastq1 "$fq1" --fastq2 "$fq2" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" --sample "$sample" \
#         > "$LOG_DIR/mycoplasma_${sample}.log" 2>&1

#     if [[ -n "$COSMIC_DIR" ]]; then
#         echo "[STEP] COSMIC mutation calling for $sample..."
#         Rscript /pipeline/modules/cancer_mutation_calling/COSMIC_cancer_mutation_calling.r \
#             --ref_dir "$REF_DIR" --cosmic_dir "$COSMIC_DIR" --output_dir "$sample_outdir" --sample "$sample" \
#             > "$LOG_DIR/cancer_mutation_mapping_${sample}.log" 2>&1
#     fi

#     echo "[STEP] eSNPKaryotyping for $sample..."
#     Rscript /pipeline/modules/eSNPKaryotyping/run_eSNPKaryotyping.R \
#         --ref_dir "$REF_DIR" --output_dir "$sample_outdir" --sample "$sample" \
#         > "$LOG_DIR/ekaryo_${sample}.log" 2>&1

#     echo "[DONE] Finished processing $sample"
# done

# echo "[STEP] Generating summary PDF..."
# Rscript /pipeline/report_builder.R \
#     --output_dir "$OUTPUT_DIR" --project "$PROJECT" \
#     > "$LOG_DIR/project_summary.log" 2>&1

# ===========================================
# Parallelized execution
# ===========================================

run_downstream() {
    local sample="$1"
    local sample_outdir="$OUTPUT_DIR/$sample"
    mkdir -p "$sample_outdir"

    echo "[INFO] Starting downstream analyses for sample: $sample"

    local fq1="${FASTQ_DIR}/${sample}_R1_001.fastq.gz"
    local fq2="${FASTQ_DIR}/${sample}_R2_001.fastq.gz"

    # -------------------------------------------------
    # Step 1: Mycoplasma detection
    # -------------------------------------------------
    echo "[STEP] Mycoplasma detection for $sample..."
    bash /pipeline/modules/mycoplasma_detection/detect_mycoplasma.sh \
        --fastq1 "$fq1" \
        --fastq2 "$fq2" \
        --ref_dir "$REF_DIR" \
        --output_dir "$sample_outdir" \
        --sample "$sample" \
        > "$LOG_DIR/mycoplasma_${sample}.log" 2>&1

    # -------------------------------------------------
    # Step 2: COSMIC mutation calling (if provided)
    # -------------------------------------------------
    if [[ -n "${COSMIC_DIR:-}" ]]; then
        echo "[STEP] COSMIC mutation calling for $sample..."
        Rscript /pipeline/modules/cancer_mutation_calling/COSMIC_cancer_mutation_calling.r \
            --ref_dir "$REF_DIR" \
            --cosmic_dir "$COSMIC_DIR" \
            --output_dir "$sample_outdir" \
            --sample "$sample" \
            > "$LOG_DIR/cancer_mutation_mapping_${sample}.log" 2>&1
    else
        echo "[SKIP] COSMIC_DIR not provided; skipping COSMIC mutation calling for $sample"
    fi

    # -------------------------------------------------
    # Step 3: eSNPKaryotyping
    # -------------------------------------------------
    echo "[STEP] eSNPKaryotyping for $sample..."
    Rscript /pipeline/modules/eSNPKaryotyping/run_eSNPKaryotyping.R \
        --ref_dir "$REF_DIR" \
        --output_dir "$sample_outdir" \
        --sample "$sample" \
        > "$LOG_DIR/ekaryo_${sample}.log" 2>&1

    echo "[DONE] Finished downstream analyses for $sample"
}

export -f run_downstream
export OUTPUT_DIR FASTQ_DIR LOG_DIR REF_DIR COSMIC_DIR

MAX_JOBS=${MAX_JOBS:-$(nproc)}
printf '%s\n' "${SAMPLES[@]}" | xargs -n1 -P "$MAX_JOBS" bash -c 'run_downstream "$@"' _


echo "[STEP] Organizing output directories..."

for sample in "${SAMPLES[@]}"; do
    sample_outdir="$OUTPUT_DIR/$sample"

    # Organize logs
    mkdir -p "$sample_outdir/logs/log_caper"
    mv -f "$sample_outdir"/caper* "$sample_outdir/logs/log_caper/" 2>/dev/null || true
for d in scripts current rnaseq_pipeline_fastq_workflow log; do
    src="$sample_outdir/$d"
    dest="$sample_outdir/logs/$d"
    if [ -d "$src" ]; then
        rm -rf "$dest"
        mv "$src" "$sample_outdir/logs/"
    fi
done

    
    mkdir -p "$OUTPUT_DIR/plots/cancer_mutations/$sample"
    [ -f "$sample_outdir/cosmic_calling/CancerMutationPlot.png" ] && cp "$sample_outdir/cosmic_calling/CancerMutationPlot.png" "$OUTPUT_DIR/plots/cancer_mutations/$sample/"

    mkdir -p "$OUTPUT_DIR/plots/eSNPKaryotyping/$sample"
    [ -f "$sample_outdir/eSNPKaryotyping/${sample}_PlotGenome.png" ] && cp "$sample_outdir/eSNPKaryotyping/${sample}_PlotGenome.png" "$OUTPUT_DIR/plots/eSNPKaryotyping/$sample/"

    mkdir -p "$OUTPUT_DIR/plots/mycoplasma_detection/$sample"
    [ -f "$sample_outdir/mycoplasma/mycoplasma_alignment_summary.pdf" ] && cp "$sample_outdir/mycoplasma/mycoplasma_alignment_summary.pdf" "$OUTPUT_DIR/plots/mycoplasma_detection/$sample/"

    # # Organize preprocessing
    # mkdir -p "$sample_outdir/preprocessing"
    # for d in \
    #     variant_calling/Mark_duplicates_outputs \
    #     variant_calling/RSEM_outputs \
    #     variant_calling/QC_outputs \
    #     star_out \
    #     fastqc_out \
    #     variant_calling; do
    #     [ -d "$sample_outdir/$d" ] && mv "$sample_outdir/$d" "$sample_outdir/preprocessing/"
    # done
done

mkdir -p "$OUTPUT_DIR/plots/PACNet"
[ -f "$OUTPUT_DIR/pacnet/PACNet_heatmap.png" ] && cp "$OUTPUT_DIR/pacnet/PACNet_heatmap.png" "$OUTPUT_DIR/plots/PACNet/"

mkdir -p "$OUTPUT_DIR/plots/outlier_analysis"
[ -f "$OUTPUT_DIR/outlier_analysis/PCA_pacnet_scores.pdf" ] && cp "$OUTPUT_DIR/outlier_analysis/PCA_pacnet_scores.pdf" "$OUTPUT_DIR/plots/outlier_analysis/"

echo "[INFO] Pipeline completed. Results in $OUTPUT_DIR"