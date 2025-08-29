#!/usr/bin/env bash
#MYCOPLASMA DETECTION
#Use bowtie2 to align reads to mycoplasma genome, and assess alignment
#The alignment will be a proxy for level of contamination

#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz
#--2025-01-24 16:08:05--  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz

set -euo pipefail

# Arguments
REF_DIR="$1"       # reference directory (from pipeline_runner.sh)
FASTQ_DIR="$2"     # fastq directory (from pipeline_runner.sh)
OUTPUT_DIR="$3"    # output directory (from pipeline_runner.sh)

# Check inputs
if [[ ! -d "$REF_DIR" ]]; then
    echo "[ERROR] Reference directory not found: $REF_DIR"
    exit 1
fi
if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "[ERROR] FASTQ directory not found: $FASTQ_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR/mycoplasma"

# Reference genome
GENOME_FILE="GCF_000027325.1_ASM2732v1_genomic.fna.gz"
GENOME_PATH="$REF_DIR/$GENOME_FILE"
GENOME_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/$GENOME_FILE"

if [[ -f "$GENOME_PATH" ]]; then
    echo "[INFO] Using genome file: $GENOME_PATH"
else
    echo "[INFO] Genome file not found, downloading..."
    wget -O "$GENOME_PATH" "$GENOME_URL"
fi

# Unzip if needed
if [[ "$GENOME_PATH" == *.gz ]]; then
    gunzip -c "$GENOME_PATH" > "$OUTPUT_DIR/mycoplasma/mycoplasma_genome.fna"
    GENOME_REF="$OUTPUT_DIR/mycoplasma/mycoplasma_genome.fna"
else
    GENOME_REF="$GENOME_PATH"
fi

# Check bowtie2
if command -v bowtie2 >/dev/null 2>&1; then
    echo "[INFO] bowtie2 found in PATH"
elif command -v module >/dev/null 2>&1; then
    echo "[INFO] Loading bowtie2 via module system..."
    module load bowtie2
else
    echo "[ERROR] bowtie2 not found. Please install or load it."
    exit 1
fi

# Build bowtie2 index
INDEX_PREFIX="$OUTPUT_DIR/mycoplasma/mycoplasma_index"
if [[ ! -f "${INDEX_PREFIX}.1.bt2" ]]; then
    echo "[INFO] Building bowtie2 index..."
    bowtie2-build "$GENOME_REF" "$INDEX_PREFIX"
else
    echo "[INFO] Bowtie2 index already exists: $INDEX_PREFIX"
fi

# Align FASTQ files
for fq in "$FASTQ_DIR"/*fastq.gz "$FASTQ_DIR"/*fq.gz; do
    [[ -e "$fq" ]] || continue  # skip if no files

    sample=$(basename "$fq" | sed 's/_R[12].*//; s/\.fastq.*//; s/\.fq.*//')

    # Paired-end detection
    fq1="$FASTQ_DIR/${sample}_R1.fastq.gz"
    fq2="$FASTQ_DIR/${sample}_R2.fastq.gz"

    out_prefix="$OUTPUT_DIR/mycoplasma/${sample}"

    if [[ -f "$fq1" && -f "$fq2" ]]; then
        echo "[INFO] Running paired-end alignment for $sample"
        bowtie2 -x "$INDEX_PREFIX" -1 "$fq1" -2 "$fq2" \
            -S "${out_prefix}.sam" --no-unal 2> "${out_prefix}_bowtie2.log"
    else
        echo "[INFO] Running single-end alignment for $sample"
        bowtie2 -x "$INDEX_PREFIX" -U "$fq" \
            -S "${out_prefix}.sam" --no-unal 2> "${out_prefix}_bowtie2.log"
    fi

    # Convert to BAM
    samtools view -bS "${out_prefix}.sam" | samtools sort -o "${out_prefix}_sorted.bam"
    samtools index "${out_prefix}_sorted.bam"

    # Clean up SAM to save space
    rm "${out_prefix}.sam"

    echo "[INFO] Finished $sample. Alignment log: ${out_prefix}_bowtie2.log"
done

echo "[INFO] Mycoplasma detection completed. Results in $OUTPUT_DIR/mycoplasma"
