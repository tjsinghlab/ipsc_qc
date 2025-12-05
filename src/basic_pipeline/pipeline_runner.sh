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
KEEP_FILES=true

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
    echo "  --keep_files BOOL    Whether to keep intermediate outputs (default: true)"
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
        --keep_files)
            # Normalize to lowercase
            val=$(echo "$2" | tr '[:upper:]' '[:lower:]')
            if [[ "$val" == "true" ]]; then
                KEEP_FILES=true
            elif [[ "$val" == "false" ]]; then
                KEEP_FILES=false
            else
                echo "[ERROR] --keep_files must be 'true' or 'false'"
                exit 1
            fi
            shift
            ;;
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

  #"GATK_SIF|gatk_4.6.1.0.sif|singularity pull $REF_DIR/gatk_4.6.1.0.sif library://biocontainers/gatk4:4.6.1.0--py310hdfd78af_0"
  #"GTEX_SIF|gtex_rnaseq_V10.sif|singularity pull $REF_DIR/gtex_rnaseq_V10.sif docker://broadinstitute/gtex_rnaseq:V10"

  "GENE_LIST|genes.txt|wget -O $REF_DIR/genes.txt https://raw.githubusercontent.com/tjsinghlab/ipsc_qc/main/src/basic_pipeline/ref_files/genes.txt"

  "MYCO_ORALE_REF|GCF_000420105.1_ASM42010v1_genomic.fna.gz|wget -O $REF_DIR/GCF_000420105.1_ASM42010v1_genomic.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/420/105/GCF_000420105.1_ASM42010v1/GCF_000420105.1_ASM42010v1_genomic.fna.gz"
  "MYCO_FERMENTANS_REF|GCF_003704055.1_ASM370405v1_genomic.fna.gz|wget -O $REF_DIR/GCF_003704055.1_ASM370405v1_genomic.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/704/055/GCF_003704055.1_ASM370405v1/GCF_003704055.1_ASM370405v1_genomic.fna.gz"
  "MYCO_HYORHINIS_REF|GCF_900476065.1_50465_F02_genomic.fna.gz|wget -O $REF_DIR/GCF_900476065.1_50465_F02_genomic.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/476/065/GCF_900476065.1_50465_F02/GCF_900476065.1_50465_F02_genomic.fna.gz"


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
 singularity exec /pipeline/modules/gtex_rnaseq_V10.sif STAR --runMode genomeGenerate \
   --genomeDir $REF_DIR/star_index_oh75 \
   --genomeFastaFiles $REF_DIR/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
   --sjdbGTFfile $REF_DIR/gencode.v39.GRCh38.annotation.gtf \
   --sjdbOverhang 75 --runThreadN 4")

RSEM_REF=$(resolve_ref "RSEM reference" "rsem_reference" \
"mkdir -p $REF_DIR/rsem_reference && \
 singularity exec /pipeline/modules/gtex_rnaseq_V10.sif rsem-prepare-reference \
   $REF_DIR/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
   $REF_DIR/rsem_reference/rsem_reference \
   --gtf $REF_DIR/gencode.v39.GRCh38.annotation.gtf \
   --num-threads 4")

# Initialize maps
declare -A fq1_map fq2_map
SAMPLES=()

for fq in "$FASTQ_DIR"/*.fastq.gz "$FASTQ_DIR"/*.fq.gz; do
    [ -e "$fq" ] || continue
    base=$(basename "$fq")

    # Match:
    #   sample_1.fastq.gz
    #   sample_R1.fastq.gz
    #   sample_R1_001.fastq.gz
    #   sample-S7_L002_R2_001.fastq.gz
    if [[ "$base" =~ ^(.+?)(?:_[0-9]{3})*[._-]([Rr]?[12])(?:_[0-9]{3})?\.f(ast)?q\.gz$ ]]; then
        sample="${BASH_REMATCH[1]}"
        readnum="${BASH_REMATCH[2]}"
        readnum="${readnum//[Rr]/}"   # normalize R1 → 1
    else
        echo "[WARN] Could not parse FASTQ filename: $base" >&2
        continue
    fi

    case "$readnum" in
        1) fq1_map["$sample"]="$fq" ;;
        2) fq2_map["$sample"]="$fq" ;;
        *) echo "[WARN] Unrecognized read number in $base" >&2 ;;
    esac
done

# Build sample list
for s in "${!fq1_map[@]}"; do
    if [[ -n "${fq2_map[$s]:-}" ]]; then
        SAMPLES+=("$s")
    else
        echo "[WARN] Sample '$s' has R1 but no R2" >&2
    fi
done

# Check R2 matches
for s in "${!fq2_map[@]}"; do
    [[ -z "${fq1_map[$s]:-}" ]] && echo "[WARN] Sample '$s' has R2 but no R1" >&2
done

echo "[INFO] Found ${#SAMPLES[@]} paired samples"

# ===========================================
# Parallelized execution
# ===========================================

# Run caper init but only if default.conf file does not already exist
# Ensure .caper directory exists
if [ ! -d "$HOME/.caper" ]; then
  mkdir -p "$HOME/.caper"
fi

# Only initialize if default.conf is missing
if [ ! -f "$HOME/.caper/default.conf" ]; then
  echo "[INFO] No Caper config found — running caper init slurm..."
  caper init slurm --out "$HOME/.caper/default.conf"
else
  echo "[INFO] Existing Caper config found at $HOME/.caper/default.conf. Skipping init."
fi

run_sample() {
    local sample="$1"

    echo "[INFO] Beginning processing for $sample"

    local sample_outdir="$OUTPUT_DIR/$sample"
    mkdir -p "$sample_outdir"

    local fq1="${fq1_map[$sample]}"
    local fq2="${fq2_map[$sample]}"

    if [[ ! -f "$fq1" || ! -f "$fq2" ]]; then
        echo "[ERROR] FASTQs for '$sample' not found" >&2
        return 1
    fi

    echo "[INFO] FASTQs: $(basename "$fq1"), $(basename "$fq2")"

    local bam_file="$sample_outdir/Mark_duplicates_outputs/${sample}.Aligned.sortedByCoord.out.md.bam"
    local rsem_file="$sample_outdir/RSEM_outputs/${sample}.rsem.genes.results.gz"

    # Step 1 — Preprocessing
    if [[ -f "$bam_file" && -f "$rsem_file" ]]; then
        echo "[SKIP] Preprocessing already complete for $sample"
    else
        echo "[RUN] Preprocessing for $sample..."
        python3 "$PY_RUNNER1" \
            --fastq1 "$fq1" \
            --fastq2 "$fq2" \
            --output_dir "$sample_outdir" \
            --sample "$sample" \
            > "$LOG_DIR/${sample}_bulk_preprocess.log" 2>&1
    fi

    # Clean bad symlinks
    find "$sample_outdir" -xtype l -name data -delete 2>/dev/null

    # Step 2 — Variant Calling
    if compgen -G "$sample_outdir/variant_calling/*.vcf.gz" > /dev/null; then
        echo "[SKIP] Variant calling already complete for $sample"
    else
        echo "[RUN] Variant calling for $sample..."
        python3 "$PY_RUNNER2" \
            --output_dir "$sample_outdir" \
            --sample "$sample" \
            > "$LOG_DIR/${sample}_wdl2.log" 2>&1
    fi

    if [[ "$KEEP_FILES" == false ]]; then
        echo "[CLEANUP] Removing intermediate artifacts for $sample..."

        rm -f "$sample_outdir"/RNAseq/*/call-SplitNCigarReads/execution/*.bam
        rm -f "$sample_outdir"/RNAseq/*/call-SplitNCigarReads/execution/*.bai
        rm -rf "$sample_outdir"/RNAseq/*/call-ScatterIntervalList/execution/out/
        rm -rf "$sample_outdir"/RNAseq/*/call-HaplotypeCaller/shard-*
        rm -f "$sample_outdir"/RNAseq/*/call-AddReadGroups/execution/*.bam
        rm -f "$sample_outdir"/RNAseq/*/call-AddReadGroups/execution/*.bai
        rm -f "$sample_outdir"/star_out/*.bai
        rm -f "$sample_outdir"/star_out/*.bam
    else
        echo "[CLEANUP] Skipped (keep_files=true)"
    fi

    echo "[DONE] Completed sample $sample"
}

# Detect memory and cores (with or without slurm)
get_total_memory_gb() {
    local total_gb

    # If running under SLURM, try to detect memory allocation
    if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
        total_gb=$(( SLURM_MEM_PER_NODE / 1024 ))  # MB → GB
    elif [[ -n "${SLURM_MEM_PER_CPU:-}" && -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
        # Sometimes SLURM defines memory per CPU instead
        total_gb=$(( (SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK) / 1024 ))
    else
        # Fallback: detect actual available system memory (in GB)
        # Use /proc/meminfo to get MemAvailable if possible, else MemTotal
        local available_kb
        available_kb=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
        if [[ -z "$available_kb" ]]; then
            available_kb=$(awk '/MemTotal/ {print $2}' /proc/meminfo)
        fi
        # Use 60% of available memory for safety when running outside SLURM
        total_gb=$(awk -v kb="$available_kb" 'BEGIN {printf "%.0f", (kb / 1024 / 1024) * 0.8}')
    fi

    echo "$total_gb"
}

get_total_cores() {
    local cores

    # Prefer SLURM, if present
    if [[ -n "${SLURM_CPUS_ON_NODE:-}" ]]; then
        cores="$SLURM_CPUS_ON_NODE"
    elif [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
        cores="$SLURM_CPUS_PER_TASK"
    else
        # Fall back to the system’s available cores
        # Use nproc if available, otherwise parse /proc/cpuinfo
        if command -v nproc &>/dev/null; then
            cores=$(nproc)
        else
            cores=$(grep -c ^processor /proc/cpuinfo)
        fi
    fi

    echo "$cores"
}

# Each sample needs ~64 GB to be safe
MEM_PER_SAMPLE=64

TOTAL_MEM_GB=$(get_total_memory_gb)
TOTAL_CORES=$(get_total_cores)

# Compute max concurrent jobs (respect memory and cores)
MAX_JOBS=$(( TOTAL_MEM_GB / MEM_PER_SAMPLE ))
(( MAX_JOBS < 1 )) && MAX_JOBS=1
(( MAX_JOBS > TOTAL_CORES )) && MAX_JOBS=$TOTAL_CORES

echo "[INFO] SLURM-aware resource detection:"
echo "       Memory available: ${TOTAL_MEM_GB} GB"
echo "       Cores available:  ${TOTAL_CORES}"
echo "       Max parallel jobs: ${MAX_JOBS}"

# Export the job function and environment
export -f run_sample
export OUTPUT_DIR FASTQ_DIR LOG_DIR PY_RUNNER1 PY_RUNNER2 REF_DIR COSMIC_DIR

# Run samples in parallel with smart throttling
printf '%s\n' "${SAMPLES[@]}" | xargs -n1 -P "$MAX_JOBS" bash -c 'run_sample "$@"' _

# ===========================================
# Step: Run PACNet + downstream scripts
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

# ===========================================
# Parallelized execution
# ===========================================

run_downstream() {
    local sample="$1"
    local sample_outdir="$OUTPUT_DIR/$sample"
    mkdir -p "$sample_outdir"

    echo "[INFO] Starting downstream analyses for sample: $sample"

    local fq1=$(ls "${FASTQ_DIR}/${sample}"*_[Rr]1*.f*q.gz 2>/dev/null || echo "")
    local fq2=$(ls "${FASTQ_DIR}/${sample}"*_[Rr]2*.f*q.gz 2>/dev/null || echo "")

    if [[ -z "$fq1" || -z "$fq2" ]]; then
        echo "[ERROR] Could not find both R1 and R2 files for sample '$sample' in $FASTQ_DIR" >&2
        exit 1
    fi

    # -------------------------------------------------
    # Step: Mycoplasma detection
    # -------------------------------------------------
    local myco_stats="$sample_outdir/mycoplasma/mycoplasma_alignment_stats.tsv"
    if [[ ! -f "$myco_stats" ]]; then
        echo "[STEP] Mycoplasma detection for $sample..."
        bash /pipeline/modules/mycoplasma_detection/detect_mycoplasma.sh \
            --fastq1 "$fq1" \
            --fastq2 "$fq2" \
            --ref_dir "$REF_DIR" \
            --output_dir "$sample_outdir" \
            --sample "$sample" \
            > "$LOG_DIR/mycoplasma_${sample}.log" 2>&1
    else
        echo "[SKIP] Mycoplasma stats exist; skipping for $sample"
    fi

    if [[ "$KEEP_FILES" == false ]]; then
        echo "[CLEANUP] Removing intermediate mycoplasma detection artifacts for $sample..."

        rm -f "$sample_outdir"/mycoplasma/*.bam
        rm -f "$sample_outdir"/mycoplasma/*.bam.bai
        rm -f "$sample_outdir"/mycoplasma/*.bt2
        rm -f "$sample_outdir"/mycoplasma/*.fna
    else
        echo "[CLEANUP] Skipped (keep_files=true)"
    fi

    echo "[DONE] Completed mycoplasma detection for $sample"

    # -------------------------------------------------
    # Step: COSMIC mutation calling (if provided)
    # -------------------------------------------------
    local cosmic_plot="$sample_outdir/cosmic_calling/CancerMutationPlot.png"
    if [[ -n "${COSMIC_DIR:-}" ]]; then
        if [[ ! -f "$cosmic_plot" ]]; then
            echo "[STEP] COSMIC mutation calling for $sample..."
            Rscript /pipeline/modules/cancer_mutation_calling/COSMIC_cancer_mutation_calling.r \
                --ref_dir "$REF_DIR" \
                --cosmic_dir "$COSMIC_DIR" \
                --output_dir "$sample_outdir" \
                --sample "$sample" \
                > "$LOG_DIR/cancer_mutation_mapping_${sample}.log" 2>&1
        else
            echo "[SKIP] COSMIC plot exists; skipping for $sample"
        fi
    else
        echo "[SKIP] COSMIC_DIR not provided; skipping COSMIC mutation calling"
    fi

    echo "[DONE] Completed cancer mutation calling for $sample"

    # -------------------------------------------------
    # Step: eSNPKaryotyping
    # -------------------------------------------------
    local ekaryo_variants=( "$sample_outdir"/eSNPKaryotyping/*_variantTable.csv )
    if [[ ! -f "${ekaryo_variants[0]}" ]]; then
        echo "[STEP] eSNPKaryotyping for $sample..."
        Rscript /pipeline/modules/eSNPKaryotyping/run_eSNPKaryotyping.R \
            --ref_dir "$REF_DIR" \
            --output_dir "$sample_outdir" \
            --sample "$sample" \
            > "$LOG_DIR/ekaryo_${sample}.log" 2>&1
    else
        echo "[SKIP] eSNPKaryotyping variant table exists; skipping for $sample"
    fi

    echo "[DONE] Finished eSNPKaryotyping for $sample"

    if [[ "$KEEP_FILES" == false ]]; then
        echo "[CLEANUP] Removing VCF and BAM files for $sample..."
        rm -f -r "$sample_outdir"/star_out
        rm -f -r "$sample_outdir"/variant_calling
        rm -f -r "$sample_outdir"/Mark_duplicates_outputs
    else
        echo "[CLEANUP] Skipped (keep_files=true)"
    fi
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

done

mkdir -p "$OUTPUT_DIR/plots/PACNet"
[ -f "$OUTPUT_DIR/pacnet/PACNet_heatmap.png" ] && cp "$OUTPUT_DIR/pacnet/PACNet_heatmap.png" "$OUTPUT_DIR/plots/PACNet/"

mkdir -p "$OUTPUT_DIR/plots/outlier_analysis"
[ -f "$OUTPUT_DIR/outlier_analysis/PCA_pacnet_scores.pdf" ] && cp "$OUTPUT_DIR/outlier_analysis/PCA_pacnet_scores.pdf" "$OUTPUT_DIR/plots/outlier_analysis/"
[ -f "$OUTPUT_DIR/outlier_analysis/PCA_counts.pdf" ] && cp "$OUTPUT_DIR/outlier_analysis/PCA_counts.pdf" "$OUTPUT_DIR/plots/outlier_analysis/"

##Construct summary table
echo "[STEP] Writing summary table..."
Rscript /pipeline/modules/summary_writer.R \
    --output_dir "$OUTPUT_DIR" \
    --pacnet_file "$OUTPUT_DIR/pacnet/classification_scores.csv" \
    --samples "$(IFS=,; echo "${SAMPLES[*]}")" \
    --project "$PROJECT"

SAMPLES_STR=$(IFS=, ; echo "${SAMPLES[*]}")

Rscript /pipeline/report_builder.R \
  --output_dir "$OUTPUT_DIR" \
  --project "$PROJECT" \
  --samples "$SAMPLES_STR" \
  --html

echo "[INFO] Summary files written."

echo "[INFO] Pipeline completed. Results in $OUTPUT_DIR"