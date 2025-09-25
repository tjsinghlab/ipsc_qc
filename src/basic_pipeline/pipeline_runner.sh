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
  "GENES_GTF|gencode.v39.GRCh38.genes.collapsed_only.gtf|wget -O $REF_DIR/gencode.v39.GRCh38.genes.collapsed_only.gtf https://storage.googleapis.com/gtex-resources/GENCODE/gencode.v39.GRCh38.genes.collapsed_only.gtf"
  "ANNOTATION_GTF|gencode.v39.GRCh38.annotation.gtf|wget -O $REF_DIR/gencode.v39.GRCh38.annotation.gtf https://storage.googleapis.com/gtex-resources/GENCODE/gencode.v39.GRCh38.annotation.gtf"

  "REF_FASTA|Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta|wget -O $REF_DIR/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta https://storage.googleapis.com/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta"
  "REF_FASTA_FAI|Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta.fai|wget -O $REF_DIR/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta.fai https://storage.googleapis.com/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta.fai"
  "REF_DICT|Homo_sapiens_assembly38_noALT_noHLA_noDecoy.dict|wget -O $REF_DIR/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.dict https://storage.googleapis.com/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.dict"

  "DBSNP_VCF|Homo_sapiens_assembly38.dbsnp138.vcf|wget -O $REF_DIR/Homo_sapiens_assembly38.dbsnp138.vcf https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
  "DBSNP_IDX|Homo_sapiens_assembly38.dbsnp138.vcf.idx|wget -O $REF_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.idx https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"

  "KNOWN_VCF1|Mills_and_1000G_gold_standard.indels.hg38.vcf.gz|wget -O $REF_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  "KNOWN_VCF2|Homo_sapiens_assembly38.known_indels.vcf.gz|wget -O $REF_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
  "KNOWN_VCF1_IDX|Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi|wget -O $REF_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
  "KNOWN_VCF2_IDX|Homo_sapiens_assembly38.known_indels.vcf.gz.tbi|wget -O $REF_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"

  "GATK_SIF|gatk_4.6.1.0.sif|singularity pull $REF_DIR/gatk_4.6.1.0.sif library://biocontainers/gatk4:4.6.1.0--py310hdfd78af_0"
  "GTEX_SIF|gtex_rnaseq_v10.sif|singularity pull $REF_DIR/gtex_rnaseq_v10.sif docker://broadinstitute/gtex_rnaseq:V10"

  "GENE_LIST|genes.txt|wget -O $REF_DIR/genes.txt https://raw.githubusercontent.com/tjsinghlab/ipsc_qc/main/src/basic_pipeline/ref_files/genes.txt"

  "MYCO_REF|GCF_000027325.1_ASM2732v1_genomic.fna.gz|wget -O $REF_DIR/GCF_000027325.1_ASM2732v1_genomic.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz"

  "PACNET_EXP|Hs_expTrain_Jun-20-2017.rda|aws s3 cp s3://cellnet-rnaseq/ref/cnproc/HS/Hs_expTrain_Jun-20-2017.rda $REF_DIR/ --no-sign-request"
  "PACNET_ST|Hs_stTrain_Jun-20-2017.rda|aws s3 cp s3://cellnet-rnaseq/ref/cnproc/HS/Hs_stTrain_Jun-20-2017.rda $REF_DIR/ --no-sign-request"

  "NCBI_FASTA|GCF_000001405.26_GRCh38_genomic.fna|wget -O $REF_DIR/GCF_000001405.26_GRCh38_genomic.fna ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38.p13"
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
for sample in "${SAMPLES[@]}"; do
    echo "[INFO] Running WDLs for sample: $sample"
    sample_outdir="$OUTPUT_DIR/$sample"
    mkdir -p "$sample_outdir"

    fq1="${FASTQ_DIR}/${sample}_R1_001.fastq.gz"
    fq2="${FASTQ_DIR}/${sample}_R2_001.fastq.gz"

    bam_file="$sample_outdir/star_out/${sample}.Aligned.sortedByCoord.out.bam"
    bai_file="${bam_file}.bai"
    rsem_file="$sample_outdir/RSEM_outputs/${sample}.rsem.genes.results.gz"

    # Step 1: Preprocessing
    if [ -f "$bam_file" ] && [ -f "$bai_file" ] && [ -f "$rsem_file" ]; then
        echo "[SKIP] Preprocessing already done for $sample"
    else
        echo "[RUN] Running preprocessing for $sample..."
        python3 "$PY_RUNNER1" \
            --fastq1 "$fq1" \
            --fastq2 "$fq2" \
            --output_dir "$sample_outdir" \
            --sample "$sample" \
            > "$LOG_DIR/${sample}_bulk_preprocess.log" 2>&1
    fi

    # Step 2: Variant Calling
    if ls "$sample_outdir/variant_calling/"*.vcf.gz 1> /dev/null 2>&1; then
        echo "[SKIP] Variant calling already done for $sample"
    else
        echo "[RUN] Running variant calling for $sample..."
        python3 "$PY_RUNNER2" \
            --output_dir "$sample_outdir" \
            --sample "$sample" \
            > "$LOG_DIR/${sample}_wdl2.log" 2>&1
    fi
done

# ===========================================
# Phase 2: Run PACNet + downstream scripts
# ===========================================
echo "[STEP] Running PACNet with all samples..."
Rscript /pipeline/modules/PACNet/run_pacnet.R \
    --ref_dir "$REF_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --project "$PROJECT" \
    > "$LOG_DIR/pacnet_all.log" 2>&1

echo "[STEP] Outlier detection for $sample..."
Rscript /pipeline/modules/outlier_detection/outlier_detection.R \
    --output_dir "$OUTPUT_DIR" --sample "$sample" --project "$PROJECT" \
    > "$LOG_DIR/outliers_${sample}.log" 2>&1

# Now loop again for per-sample downstream modules
for sample in "${SAMPLES[@]}"; do
    sample_outdir="$OUTPUT_DIR/$sample"
    fq1="${FASTQ_DIR}/${sample}_R1_001.fastq.gz"
    fq2="${FASTQ_DIR}/${sample}_R2_001.fastq.gz"

    echo "[STEP] eSNPKaryotyping for $sample..."
    Rscript /pipeline/modules/eSNPKaryotyping/run_eSNPKaryotyping.R \
        --ref_dir "$REF_DIR" --output_dir "$sample_outdir" --sample "$sample" \
        > "$LOG_DIR/ekaryo_${sample}.log" 2>&1

    echo "[STEP] Mycoplasma detection for $sample..."
    bash /pipeline/modules/mycoplasma_detection/detect_mycoplasma.sh \
        --fastq "$fq1" "$fq2" --ref_dir "$REF_DIR" --output_dir "$sample_outdir" --sample "$sample" \
        > "$LOG_DIR/mycoplasma_${sample}.log" 2>&1

    if [[ -n "$COSMIC_DIR" ]]; then
        echo "[STEP] COSMIC mutation calling for $sample..."
        # bash /pipeline/modules/cancer_mutation_calling/filter_vcf_on_cosmic.sh \
        #     --ref_dir "$REF_DIR" --cosmic_dir "$COSMIC_DIR" --output_dir "$sample_outdir" --sample "$sample" \
        #     > "$LOG_DIR/filtered_vcf_for_cosmic_${sample}.log" 2>&1

        # Rscript /pipeline/modules/cancer_mutation_calling/cancer_mutation_mapping.R \
        #     --ref_dir "$REF_DIR" --cosmic_dir "$COSMIC_DIR" --output_dir "$sample_outdir" --sample "$sample" \
        #     > "$LOG_DIR/cancer_mutation_mapping_${sample}.log" 2>&1
        Rscript /pipeline/modules/cancer_mutation_calling/COSMIC_cancer_mutation_calling.r \
            --ref_dir "$REF_DIR" --cosmic_dir "$COSMIC_DIR" --output_dir "$sample_outdir" --sample "$sample" \
            > "$LOG_DIR/cancer_mutation_mapping_${sample}.log" 2>&1
    fi

    echo "[STEP] Generating HTML summary for $sample..."
    Rscript /pipeline/modules/report_builder/generate_html_summary.R \
        --output_dir "$OUTPUT_DIR" --project "$PROJECT" --sample "$sample" \
        > "$LOG_DIR/html_summary_${sample}.log" 2>&1

    echo "[DONE] Finished processing $sample"
done

echo "[INFO] Pipeline completed. Results in $OUTPUT_DIR"