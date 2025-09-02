#!/usr/bin/env bash
# Mycoplasma detection with alignment % and correlation plot
set -euo pipefail

# -------------------------------
# Parse named arguments
# -------------------------------
REF_DIR=""
FASTQ_DIR=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref_dir) REF_DIR="$2"; shift 2 ;;
        --fastq_dir) FASTQ_DIR="$2"; shift 2 ;;
        --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
        *) echo "[ERROR] Unknown argument: $1"; exit 1 ;;
    esac
done

# Validate inputs
for dir in "$REF_DIR" "$FASTQ_DIR"; do
    [[ -d "$dir" ]] || { echo "[ERROR] Directory not found: $dir"; exit 1; }
done
mkdir -p "$OUTPUT_DIR/mycoplasma"

# -------------------------------
# Reference genome
# -------------------------------
GENOME_FILE="GCF_000027325.1_ASM2732v1_genomic.fna.gz"
GENOME_PATH="$REF_DIR/$GENOME_FILE"
GENOME_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/$GENOME_FILE"

[[ -f "$GENOME_PATH" ]] || wget -O "$GENOME_PATH" "$GENOME_URL"

GENOME_REF="$OUTPUT_DIR/mycoplasma/mycoplasma_genome.fna"
gunzip -c "$GENOME_PATH" > "$GENOME_REF"

command -v bowtie2 >/dev/null 2>&1 || { echo "[ERROR] bowtie2 not found"; exit 1; }

INDEX_PREFIX="$OUTPUT_DIR/mycoplasma/mycoplasma_index"
[[ -f "${INDEX_PREFIX}.1.bt2" ]] || bowtie2-build "$GENOME_REF" "$INDEX_PREFIX"

# -------------------------------
# Alignment stats
# -------------------------------
ALIGN_STATS="$OUTPUT_DIR/mycoplasma/mycoplasma_alignment_stats.tsv"
echo -e "Sample\tTotal_Reads\tPercent_Aligned" > "$ALIGN_STATS"

for fq in "$FASTQ_DIR"/*fastq.gz "$FASTQ_DIR"/*fq.gz; do
    [[ -e "$fq" ]] || continue
    sample=$(basename "$fq" | sed 's/_R[12].*//; s/\.fastq.*//; s/\.fq.*//')
    fq1="$FASTQ_DIR/${sample}_R1.fastq.gz"
    fq2="$FASTQ_DIR/${sample}_R2.fastq.gz"
    out_prefix="$OUTPUT_DIR/mycoplasma/${sample}"

    if [[ -f "$fq1" && -f "$fq2" ]]; then
        bowtie2 -x "$INDEX_PREFIX" -1 "$fq1" -2 "$fq2" -S "${out_prefix}.sam" --no-unal 2> "${out_prefix}_bowtie2.log"
    else
        bowtie2 -x "$INDEX_PREFIX" -U "$fq" -S "${out_prefix}.sam" --no-unal 2> "${out_prefix}_bowtie2.log"
    fi

    samtools view -bS "${out_prefix}.sam" | samtools sort -o "${out_prefix}_sorted.bam"
    samtools index "${out_prefix}_sorted.bam"
    rm "${out_prefix}.sam"

    total_reads=$(grep "reads; of these:" "${out_prefix}_bowtie2.log" | head -n1 | awk '{print $1}')
    aligned_reads=$(grep "aligned exactly 1 time" "${out_prefix}_bowtie2.log" | awk '{sum += $1} END {print sum}')
    percent_aligned=$(awk -v a="$aligned_reads" -v t="$total_reads" 'BEGIN{if(t>0) print (a/t)*100; else print 0}')
    echo -e "${sample}\t${total_reads}\t${percent_aligned}" >> "$ALIGN_STATS"
done

echo "[INFO] Alignment stats saved to $ALIGN_STATS"

# -------------------------------
# Embedded R plotting (bar + correlation)
# -------------------------------
Rscript --vanilla - <<'EOF'
library(ggplot2)
library(readr)
library(dplyr)

ALIGN_STATS <- Sys.getenv("ALIGN_STATS")
align_stats <- read_tsv(ALIGN_STATS)

# Bar plot: % aligned per sample
p1 <- ggplot(align_stats, aes(x = reorder(Sample, Percent_Aligned), y = Percent_Aligned)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Mycoplasma Contamination per Sample",
       x = "Sample", y = "Percent of Reads Aligned") +
  theme_minimal(base_size = 14)

# Correlation plot: total reads vs % aligned
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

echo "[INFO] PDF plots saved: $OUTPUT_DIR/mycoplasma/mycoplasma_alignment_summary.pdf"
