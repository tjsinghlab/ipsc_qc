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

# ===========================================
# Reference resolution
# ===========================================
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

# ===========================================
# Define reference files
# ===========================================
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

# ===========================================
# Resolve STAR and RSEM refs in one loop
# ===========================================
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

# Parse sample names
# Initialize maps
declare -A fq1_map fq2_map # declare 2 associative arrays (keys are sample names, values are fastq paths)
SAMPLES=()

for fq in "$FASTQ_DIR"/*.fastq.gz "$FASTQ_DIR"/*.fq.gz; do # loop over all files ending in fastq.gz or fq.gz
    [ -e "$fq" ] || continue # check if files exist
    base=$(basename "$fq") # extract base name from fasttq file
    sample="" # initialize variable name for sample
    readnum="" # initialize variable name for read number (1 or 2)

    # --------- Option 1: SAMPLENAME_R1/R2.fastq.gz ---------
    if [[ "$base" =~ ^(.+?)_R([12])\.f(ast)?q\.gz$ ]]; then
        sample="${BASH_REMATCH[1]}"
        readnum="${BASH_REMATCH[2]}"
    # --------- Option 2: SAMPLENAME_1/2.fastq.gz ---------
    elif [[ "$base" =~ ^(.+?)_([12])\.f(ast)?q\.gz$ ]]; then
        sample="${BASH_REMATCH[1]}"
        readnum="${BASH_REMATCH[2]}"
    # --------- Option 3: Complex Illumina-style filenames ---------
    #elif [[ "$base" =~ ^(.+?)(?:_S[0-9]+)?(?:_L[0-9]{3})?[_-]?([Rr]?[12])(?:_[0-9]{3})?\.f(ast)?q\.gz$ ]]; then
    elif [[ "$base" =~ ^(.+?)_S[0-9]+_L[0-9]{3}_[Rr]?([12])_[0-9]{3}\.f(ast)?q\.gz$ ]]; then
        sample="${BASH_REMATCH[1]}"
        readnum="${BASH_REMATCH[2]}"
        readnum="${readnum//[Rr]/}"  # normalize R1/R2 → 1/2
    else
        echo "[WARN] Could not parse FASTQ filename: $base" >&2
        continue
    fi

    # --------- Assign to map ---------
    case "$readnum" in # depending on the read number (1 or 2), store path in fq_map:
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

echo "[CHECK] FASTQ mappings for each sample:"
for s in "${SAMPLES[@]}"; do
    echo "  Sample: $s"
    echo "    R1: ${fq1_map[$s]}"
    echo "    R2: ${fq2_map[$s]}"
done

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
  echo "[INFO] No Caper config found — running caper init..."
  caper init local --out "$HOME/.caper/default.conf"
else
  echo "[INFO] Existing Caper config found at $HOME/.caper/default.conf. Skipping init."
fi

# ===========================================
# Resource detection
# ===========================================
get_total_memory_gb() {
    local total_gb
    if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
        total_gb=$(( SLURM_MEM_PER_NODE / 1024 ))
    elif [[ -n "${SLURM_MEM_PER_CPU:-}" && -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
        total_gb=$(( (SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK) / 1024 ))
    else
        local available_kb
        available_kb=$(awk '/MemAvailable/ {print $2}' /proc/meminfo || true)
        [[ -z "$available_kb" ]] && available_kb=$(awk '/MemTotal/ {print $2}' /proc/meminfo)
        total_gb=$(awk -v kb="$available_kb" 'BEGIN {printf "%.0f", (kb/1024/1024)*0.8}')
    fi
    echo "$total_gb"
}

get_total_cores() {
    if   [[ -n "${SLURM_CPUS_ON_NODE:-}" ]]; then echo "$SLURM_CPUS_ON_NODE"
    elif [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then echo "$SLURM_CPUS_PER_TASK"
    elif command -v nproc >/dev/null;     then nproc
    else grep -c ^processor /proc/cpuinfo
    fi
}

MEM_PER_SAMPLE=96
TOTAL_MEM_GB=$(get_total_memory_gb)
TOTAL_CORES=$(get_total_cores)

MAX_JOBS=$(( TOTAL_MEM_GB / MEM_PER_SAMPLE ))
(( MAX_JOBS < 1 )) && MAX_JOBS=1
(( MAX_JOBS > TOTAL_CORES )) && MAX_JOBS=$TOTAL_CORES

# NEW: hard cap at 5 concurrent jobs
(( MAX_JOBS > 5 )) && MAX_JOBS=5

echo "[INFO] Memory: $TOTAL_MEM_GB GB"
echo "[INFO] Cores : $TOTAL_CORES"
echo "[INFO] Max parallel jobs: $MAX_JOBS"

# ===========================================
# Runner fnctns
# ===========================================
run_preprocess() {
    local sample="$1"
    local fq1="$2"
    local fq2="$3"
    local out="$OUTPUT_DIR/$sample"
    mkdir -p "$out"

    local bam="$out/Mark_duplicates_outputs/${sample}.Aligned.sortedByCoord.out.md.bam"
    local rsem="$out/RSEM_outputs/${sample}.rsem.genes.results.gz"

    if [[ -f "$bam" && -f "$rsem" ]]; then
        echo "[SKIP] Preprocess already complete for $sample"
        return 0
    fi

    echo "[RUN] Preprocessing $sample"

    timeout --preserve-status 8h \
        python3 "$PY_RUNNER1" \
            --fastq1 "$fq1" \
            --fastq2 "$fq2" \
            --output_dir "$out" \
            --sample "$sample" \
        > "$LOG_DIR/${sample}_bulk_preprocess.log" 2>&1

    local rc=$?
    if [[ $rc -eq 124 ]]; then
        echo "[TIMEOUT] Preprocess timed out for $sample"
        return 1
    fi
    return $rc
}

run_variant() {
    local sample="$1"
    local out="$OUTPUT_DIR/$sample"

    if compgen -G "$out/variant_calling/*.vcf.gz" >/dev/null; then
        echo "[SKIP] Variant calling complete for $sample"
        return 0
    fi

    echo "[RUN] Variant calling for $sample"

    timeout --preserve-status 4h \
        python3 "$PY_RUNNER2" \
            --output_dir "$out" \
            --sample "$sample" \
        > "$LOG_DIR/${sample}_variant.log" 2>&1

    local rc=$?
    if [[ $rc -eq 124 ]]; then
        echo "[TIMEOUT] Variant calling timed out for $sample"
        return 1
    fi
    return $rc
}

export -f run_preprocess run_variant
export OUTPUT_DIR LOG_DIR PY_RUNNER1 PY_RUNNER2

# ===========================================
# Parallel execution of WDL pipelines (will retry samples which fail and do not produce populated outputs)
# ===========================================
run_phase_with_retries() {
    local phase="$1"   # "preprocess" or "variant"
    shift
    local samples=("$@")
    local retries=3

    for attempt in $(seq 1 $retries); do
        echo "[PHASE: $phase] ATTEMPT $attempt/$retries"

        # Build argument array for xargs
        RUN_ARGS=()
        for s in "${samples[@]}"; do
            RUN_ARGS+=("$s" "${fq1_map[$s]}" "${fq2_map[$s]}")
        done

        printf '%s\n' "${RUN_ARGS[@]}" \
          | xargs -n 3 -P "$MAX_JOBS" bash -c 'run_'"$phase"' "$@"' _

        # After run, check completeness
        incomplete=()
        for s in "${samples[@]}"; do
            out="$OUTPUT_DIR/$s"
            if [[ "$phase" == "preprocess" ]]; then
                bam="$out/Mark_duplicates_outputs/${s}.Aligned.sortedByCoord.out.md.bam"
                rsem="$out/RSEM_outputs/${s}.rsem.genes.results.gz"
                [[ ! -f "$bam" || ! -f "$rsem" ]] && incomplete+=("$s")
            else
                compgen -G "$out/variant_calling/*.vcf.gz" >/dev/null || incomplete+=("$s")
            fi
        done

        if (( ${#incomplete[@]} == 0 )); then
            echo "[PHASE: $phase] SUCCESS — all samples completed"
            return 0
        fi

        echo "[PHASE: $phase] Incomplete samples: ${incomplete[*]}"
        samples=("${incomplete[@]}")
    done

    echo "[ERROR] Phase '$phase' failed after $retries attempts: ${incomplete[*]}"
    exit 1
}

# ===========================================
# Bulk RNA Preprocessing
# ===========================================
run_phase_with_retries preprocess "${SAMPLES[@]}"

# ===========================================
# Variant Calling
# ===========================================
run_phase_with_retries variant "${SAMPLES[@]}"

# ===========================================
# File Cleanup
# ===========================================
echo "[INFO] File cleanup"

for sample in "${SAMPLES[@]}"; do
    out="$OUTPUT_DIR/$sample"

    if [[ "$KEEP_FILES" == false ]]; then
        echo "[CLEANUP] Removing intermediate artifacts for $sample..."
        rm -f "$out"/RNAseq/*/call-SplitNCigarReads/execution/*.bam
        rm -f "$out"/RNAseq/*/call-SplitNCigarReads/execution/*.bai
        rm -rf "$out"/RNAseq/*/call-ScatterIntervalList/execution/out/
        rm -rf "$out"/RNAseq/*/call-HaplotypeCaller/shard-*
        rm -f "$out"/RNAseq/*/call-AddReadGroups/execution/*.bam
        rm -f "$out"/RNAseq/*/call-AddReadGroups/execution/*.bai
        rm -f "$out"/star_out/*.bai
        rm -f "$out"/star_out/*.bam
    else
        echo "[CLEANUP] keep_files=true — skipping cleanup for $sample"
    fi
done

echo "[INFO] Completed preprocessing and variant calling for all samples."


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
    local fq1="$2"
    local fq2="$3"
    local sample_outdir="$OUTPUT_DIR/$sample"
    mkdir -p "$sample_outdir"

    echo "[INFO] Starting downstream analyses for sample: $sample"
    if [[ -z "$fq1" || -z "$fq2" || ! -f "$fq1" || ! -f "$fq2" ]]; then
        echo "[ERROR] Could not find both R1 and R2 files for sample '$sample' in $FASTQ_DIR" >&2
        return 1   # return, don't exit entire pipeline
    fi

    echo "[INFO] Using FASTQs: $(basename "$fq1"), $(basename "$fq2")"

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
#printf '%s\n' "${SAMPLES[@]}" | xargs -n1 -P "$MAX_JOBS" bash -c 'run_downstream "$@"' _
tmplist=$(mktemp)
for s in "${SAMPLES[@]}"; do
  printf '%s\t%s\t%s\n' "$s" "${fq1_map[$s]}" "${fq2_map[$s]}" >> "$tmplist"
done

# Run downstream analyses in parallel safely
xargs -a "$tmplist" -n1 -P "$MAX_JOBS" -I{} bash -c "IFS=\$'\t'; read -r sample fq1 fq2 <<< \"{}\"; run_downstream \"\$sample\" \"\$fq1\" \"\$fq2\"" _

rm -f "$tmplist"


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

# Construct summary table
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