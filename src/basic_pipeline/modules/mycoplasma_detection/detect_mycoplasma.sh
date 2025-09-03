#!/usr/bin/env bash
# Mycoplasma detection per sample
set -euo pipefail

# -------------------------------
# Parse named arguments
# -------------------------------
REF_DIR=""
FASTQ_DIR=""
OUTPUT_DIR=""
SAMPLE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref_dir) REF_DIR="$2"; shift 2 ;;
        --fastq_dir) FASTQ_DIR="$2"; shift 2 ;;
        --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
        --sample) SAMPLE="$2"; shift 2 ;;
        *) echo "[ERROR] Unknown argument: $1"; exit 1 ;;
    esac
done

# Validate inputs
if [[ -z "$SAMPLE" ]]; then
    echo "[ERROR] --sample must be provided"
    exit 1
fi

for dir in "$REF_DIR" "$FASTQ_DIR"; do
    [[ -d "$dir" ]] || { echo "[ERROR] Directory not found: $dir"; exit 1; }
done

SAMPLE_OUT="$OUTPUT_DIR/mycoplasma/$SAMPLE"
mkdir -p "$SAMPLE_OUT"

# -------------------------------
# Reference genome
# -------------------------------
GENOME_FILE="GCF_000027325.1_ASM2732v1_genomic.fna.gz"
GENOME_PATH="$REF_DIR/$GENOME_FILE"
GENOME_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/$GENOME_FILE"

[[ -f "$GENOME_PATH" ]] || wget -O "$GENOME_PATH" "$GENOME_URL"

GENOME_REF="$SAMPLE_OUT/mycoplasma_genome.fna"
gunzip -c "$GENOME_PATH" > "$GENOME_REF"

command -v bowtie2 >/dev/null 2>&1 || { echo "[ERROR] bowtie2 not found"; exit 1; }

INDEX_PREFIX="$OUTPUT_DIR/mycoplasma/mycoplasma_index"
[[ -f "${INDEX_PREFIX}.1.bt2" ]] || bowtie2-build "$GENOME_REF" "$INDEX_PREFIX"

# -------------------------------
# Alignment
# -------------------------------
fq1="$FASTQ_DIR/${SAMPLE}_R1.fastq.gz"
fq2="$FASTQ_DIR/${SAMPLE}_R2.fastq.gz"

[[ -f "$fq1" ]] || { echo "[ERROR] FASTQ file not found: $fq1"; exit 1; }

OUT_PREFIX="$SAMPLE_OUT/${SAMPLE}"

if [[ -f "$fq2" ]]; then
    bowtie2 -x "$INDEX_PREFIX" -1 "$fq1" -2 "$fq2" -S "${OUT_PREFIX}.sam" --no-unal 2> "${OUT_PREFIX}_bowtie2.log"
else
    bowtie2 -x "$INDEX_PREFIX" -U "$fq1" -S "${OUT_PREFIX}.sam" --no-unal 2> "${OUT_PREFIX}_bowtie2.log"
fi

samtools view -bS "${OUT_PREFIX}.sam" | samtools sort -o "${OUT_PREFIX}_sorted.bam"
samtools index "${OUT_PREFIX}_sorted.bam"
rm "${OUT_PREFIX}.sam"

# -------------------------------
# Compute alignment stats
# -------------------------------
TOTAL_READS=$(grep "reads; of these:" "${OUT_PREFIX}_bowtie2.log" | head -n1 | awk '{print $1}')
ALIGNED_READS=$(grep "aligned exactly 1 time" "${OUT_PREFIX}_bowtie2.log" | awk '{sum += $1} END {print sum}')
PERCENT_ALIGNED=$(awk -v a="$ALIGNED_READS" -v t="$TOTAL_READS" 'BEGIN{if(t>0) print (a/t)*100; else print 0}')

ALIGN_STATS="$SAMPLE_OUT/mycoplasma_alignment_stats.tsv"
echo -e "Sample\tTotal_Reads\tPercent_Aligned" > "$ALIGN_STATS"
echo -e "${SAMPLE}\t${TOTAL_READS}\t${PERCENT_ALIGNED}" >> "$ALIGN_STATS"

# -------------------------------
# R plotting
# -------------------------------
export ALIGN_STATS
Rscript --vanilla - <<'EOF'
library(ggplot2)
library(readr)
library(dplyr)

ALIGN_STATS <- Sys.getenv("ALIGN_STATS")
align_stats <- read_tsv(ALIGN_STATS)

p1 <- ggplot(align_stats, aes(x = reorder(Sample, Percent_Aligned), y = Percent_Aligned)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Mycoplasma Contamination per Sample",
       x = "Sample", y = "Percent of Reads Aligned") +
  theme_minimal(base_size = 14)

p2 <- ggplot(align_stats, aes(x = Total_Reads, y = Percent_Aligned)) +
  geom_point(color = "firebrick", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(title = "Correlation: Total Reads vs % Aligned",
       x = "Total Reads", y = "Percent Aligned") +
  theme_minimal(base_size = 14)

pdf(file = file.path(dirname(ALIGN_STATS), "mycoplasma_alignment_summary.pdf"), width = 10, height = 7)
print(p1)
print(p2)
dev.off()
EOF

echo "[INFO] Mycoplasma detection complete for $SAMPLE"
echo "[INFO] Alignment stats: $ALIGN_STATS"
echo "[INFO] PDF plots: $SAMPLE_OUT/mycoplasma_alignment_summary.pdf"
