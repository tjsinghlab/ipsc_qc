#!/usr/bin/env bash
# Mycoplasma detection per sample
set -euo pipefail

# -------------------------------
# Parse named arguments
# -------------------------------
FASTQ1_DIR=""
FASTQ2_DIR=""
OUTPUT_DIR=""
SAMPLE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq1) FASTQ1_DIR="$2"; shift 2 ;;
    --fastq2) FASTQ2_DIR="$2"; shift 2 ;;
    --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --ref_dir) REF_DIR="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 --fastq Fastq files for sample --output_dir directory to sample-specific output directory --sample Sample name --ref_dir Reference directory (from Docker image)"
      exit 0
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

# Validate inputs
if [[ -z "$SAMPLE" ]]; then
    echo "[ERROR] --sample must be provided"
    exit 1
fi

# [[ -d "$FASTQ_DIR" ]] || { echo "[ERROR] Directory not found: $FASTQ_DIR"; exit 1; }

SAMPLE_OUT="$OUTPUT_DIR/mycoplasma"
mkdir -p "$SAMPLE_OUT"

# -------------------------------
# Reference genome (from Docker image at /ref/)
# -------------------------------
REF_DIR="/ref"
GENOME_FILE="GCF_000027325.1_ASM2732v1_genomic.fna.gz"
GENOME_PATH="$REF_DIR/$GENOME_FILE"

if [[ ! -f "$GENOME_PATH" ]]; then
    echo "[ERROR] Reference genome not found inside Docker image: $GENOME_PATH"
    exit 1
fi

GENOME_REF="$SAMPLE_OUT/mycoplasma_genome.fna"
gunzip -c "$GENOME_PATH" > "$GENOME_REF"

command -v bowtie2 >/dev/null 2>&1 || { echo "[ERROR] bowtie2 not found"; exit 1; }

INDEX_PREFIX="$OUTPUT_DIR/mycoplasma/mycoplasma_index"
[[ -f "${INDEX_PREFIX}.1.bt2" ]] || bowtie2-build "$GENOME_REF" "$INDEX_PREFIX"

# -------------------------------
# Alignment
# -------------------------------
fq1=$FASTQ1_DIR
fq2=$FASTQ2_DIR

[[ -f "$fq1" ]] || { echo "[ERROR] FASTQ file not found: $fq1"; exit 1; }

OUT_PREFIX="$SAMPLE_OUT/${SAMPLE}"

if [[ -f "$fq2" ]]; then
    bowtie2 -x "$INDEX_PREFIX" -1 "$fq1" -2 "$fq2" -S "${OUT_PREFIX}.sam" 2> "${OUT_PREFIX}_bowtie2.log"
else
    bowtie2 -x "$INDEX_PREFIX" -U "$fq1" -S "${OUT_PREFIX}.sam" 2> "${OUT_PREFIX}_bowtie2.log"
fi

samtools view -bS -@ 4 "${OUT_PREFIX}.sam" -o "${OUT_PREFIX}.bam"
samtools sort -@ 4 -o "${OUT_PREFIX}_sorted.bam" "${OUT_PREFIX}.bam"
samtools index "${OUT_PREFIX}_sorted.bam"
rm "${OUT_PREFIX}.sam" "${OUT_PREFIX}.bam"

# -------------------------------
# Compute alignment stats
# -------------------------------
# TOTAL_READS=$(grep "reads; of these:" "${OUT_PREFIX}_bowtie2.log" | head -n1 | awk '{print $1}')
# ALIGNED_READS=$(grep "aligned exactly 1 time" "${OUT_PREFIX}_bowtie2.log" | awk '{sum += $1} END {print sum}')
# PERCENT_ALIGNED=$(awk -v a="$ALIGNED_READS" -v t="$TOTAL_READS" 'BEGIN{if(t>0) print (a/t)*100; else print 0}')

# ALIGN_STATS="$SAMPLE_OUT/mycoplasma_alignment_stats.tsv"
# echo -e "Sample\tTotal_Reads\tPercent_Aligned" > "$ALIGN_STATS"
# echo -e "${SAMPLE}\t${TOTAL_READS}\t${PERCENT_ALIGNED}" >> "$ALIGN_STATS"

FLAGSTAT_OUT="$SAMPLE_OUT/${SAMPLE}_flagstat.txt"
samtools flagstat "${OUT_PREFIX}_sorted.bam" > "$FLAGSTAT_OUT"

# parse:
# e.g. "53041911 + 0 in total (QC-passed reads + QC-failed reads)"
# and "   123 + 0 mapped (0.00% : N/A)"
TOTAL_READS=$(awk '/in total/ {print $1; exit}' "$FLAGSTAT_OUT")
MAPPED_READS=$(awk '/mapped \(/ {print $1; exit}' "$FLAGSTAT_OUT")

# fallback to zero if parsing failed
TOTAL_READS=${TOTAL_READS:-0}
MAPPED_READS=${MAPPED_READS:-0}

PERCENT_ALIGNED=$(awk -v a="$MAPPED_READS" -v t="$TOTAL_READS" \
  'BEGIN{if(t>0) printf("%.6f", (a/t)*100); else print "0"}')

ALIGN_STATS="$SAMPLE_OUT/mycoplasma_alignment_stats.tsv"
echo -e "Sample\tTotal_Reads\tMapped_Reads\tPercent_Aligned" > "$ALIGN_STATS"
echo -e "${SAMPLE}\t${TOTAL_READS}\t${MAPPED_READS}\t${PERCENT_ALIGNED}" >> "$ALIGN_STATS"

# -------------------------------
# R plotting
# -------------------------------
export ALIGN_STATS
Rscript --vanilla - <<'EOF'
library(ggplot2)
library(readr)

# Read alignment stats
ALIGN_STATS <- Sys.getenv("ALIGN_STATS")
align_stats <- read_tsv(ALIGN_STATS, show_col_types = FALSE)

percent <- as.numeric(align_stats$Percent_Aligned[1])
sample_name <- align_stats$Sample[1]
species_name <- align_stats$Species[1]

# Gradient bar data
df <- data.frame(x = seq(0, 100, length.out = 500), y = 0)

p <- ggplot(df, aes(x, y)) +
  # gradient bar
  geom_tile(aes(fill = x), height = 0.2) +
  scale_fill_gradient(low = "skyblue", high = "red", guide = "none") +
  
  # short vertical line marker
  geom_segment(aes(x = percent, xend = percent, y = 0, yend = 0.35),
               color = "black", size = 1) +
  
  # percent label above the line
  annotate("text", x = percent, y = 0.5,
           label = sprintf("%.4f%%", percent),
           hjust = 0.5, vjust = 0, size = 5) +
  
  labs(title = paste0(sample_name, " — ", species_name),
       x = "Percent Aligned", y = NULL) +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

# Save PDF and PNG
out_pdf <- file.path(dirname(ALIGN_STATS), "mycoplasma_alignment_summary.pdf")
out_png <- file.path(dirname(ALIGN_STATS), "mycoplasma_alignment_summary.png")
ggsave(out_pdf, p, width = 7, height = 2.5)
ggsave(out_png, p, width = 7, height = 2.5)
EOF

echo "[INFO] Mycoplasma detection complete for $SAMPLE"
echo "[INFO] Alignment stats: $ALIGN_STATS"
echo "[INFO] PDF plots: $SAMPLE_OUT/mycoplasma_alignment_summary.pdf"
