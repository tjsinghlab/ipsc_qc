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
# Reference genomes (from Docker image at /ref/)
# -------------------------------
REF_DIR="/ref"
declare -A GENOMES
GENOMES[orale]="GCF_000420105.1_ASM42010v1_genomic.fna.gz"
GENOMES[fermentans]="GCF_003704055.1_ASM370405v1_genomic.fna.gz"

command -v bowtie2 >/dev/null 2>&1 || { echo "[ERROR] bowtie2 not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "[ERROR] samtools not found"; exit 1; }

ALIGN_STATS_FILE="$SAMPLE_OUT/mycoplasma_alignment_stats.tsv"
echo -e "Sample\tSpecies\tTotal_Reads\tMapped_Reads\tPercent_Aligned" > "$ALIGN_STATS_FILE"

for species in "${!GENOMES[@]}"; do
    GENOME_FILE="${GENOMES[$species]}"
    GENOME_PATH="$REF_DIR/$GENOME_FILE"

    if [[ ! -f "$GENOME_PATH" ]]; then
        echo "[ERROR] Reference genome not found: $GENOME_PATH"
        exit 1
    fi

    GENOME_REF="$SAMPLE_OUT/mycoplasma_${species}_genome.fna"
    gunzip -c "$GENOME_PATH" > "$GENOME_REF"

    INDEX_PREFIX="$SAMPLE_OUT/mycoplasma_index_${species}"
    [[ -f "${INDEX_PREFIX}.1.bt2" ]] || bowtie2-build "$GENOME_REF" "$INDEX_PREFIX"

    OUT_PREFIX="$SAMPLE_OUT/${SAMPLE}_${species}"
    if [[ -f "$FASTQ2_DIR" ]]; then
        bowtie2 -x "$INDEX_PREFIX" -1 "$FASTQ1_DIR" -2 "$FASTQ2_DIR" -S "${OUT_PREFIX}.sam" 2> "${OUT_PREFIX}_bowtie2.log"
    else
        bowtie2 -x "$INDEX_PREFIX" -U "$FASTQ1_DIR" -S "${OUT_PREFIX}.sam" 2> "${OUT_PREFIX}_bowtie2.log"
    fi

    samtools view -bS -@ 4 "${OUT_PREFIX}.sam" -o "${OUT_PREFIX}.bam"
    samtools sort -@ 4 -o "${OUT_PREFIX}_sorted.bam" "${OUT_PREFIX}.bam"
    samtools index "${OUT_PREFIX}_sorted.bam"
    rm "${OUT_PREFIX}.sam" "${OUT_PREFIX}.bam"

    FLAGSTAT_OUT="$SAMPLE_OUT/${SAMPLE}_${species}_flagstat.txt"
    samtools flagstat "${OUT_PREFIX}_sorted.bam" > "$FLAGSTAT_OUT"

    TOTAL_READS=$(awk '/in total/ {print $1; exit}' "$FLAGSTAT_OUT")
    MAPPED_READS=$(awk '/mapped \(/ {print $1; exit}' "$FLAGSTAT_OUT")
    TOTAL_READS=${TOTAL_READS:-0}
    MAPPED_READS=${MAPPED_READS:-0}
    PERCENT_ALIGNED=$(awk -v a="$MAPPED_READS" -v t="$TOTAL_READS" 'BEGIN{if(t>0) printf("%.6f", (a/t)*100); else print "0"}')

    echo -e "${SAMPLE}\t${species}\t${TOTAL_READS}\t${MAPPED_READS}\t${PERCENT_ALIGNED}" >> "$ALIGN_STATS_FILE"
done

# -------------------------------
# R plotting
# -------------------------------
export ALIGN_STATS_FILE
Rscript --vanilla - <<'EOF'
library(ggplot2)
library(readr)
library(dplyr)
library(patchwork)

ALIGN_STATS_FILE <- Sys.getenv("ALIGN_STATS_FILE")
align_stats <- read_tsv(ALIGN_STATS_FILE, show_col_types = FALSE)

# Gradient bar data
df <- data.frame(x = seq(0, 100, length.out = 500), y = 0)

plots <- lapply(1:nrow(align_stats), function(i) {
  species_name <- align_stats$Species[i]
  percent <- as.numeric(align_stats$Percent_Aligned[i])
  sample_name <- align_stats$Sample[i]

  p <- ggplot(df, aes(x, y)) +
    # Gradient block
    geom_tile(aes(fill = x), height = 0.05) +
    scale_fill_gradient(low = "skyblue", high = "red", guide = "none") +
    # Marker line (exactly matching the height of the gradient bar)
    geom_segment(aes(x = percent, xend = percent, y = -0.025, yend = 0.025),
                 color = "black", size = 1) +
    # Percentage text (moved down a little, smaller font)
    annotate("text", x = percent, y = 0.08,
             label = sprintf("%.2f%%", percent),
             hjust = 0.5, vjust = 0, size = 3.5) +
    # Shorter, readable title
    labs(title = paste0(sample_name, " ", species_name),
         x = "Percent Aligned", y = NULL) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        vjust = 0,
        size = 14,     # smaller than before (20 → 14)
        face = "bold"
      ),
      plot.margin = margin(10, 10, 10, 10)
    )

  return(p)
})

# Combine plots
final_plot <- wrap_plots(plots, ncol = 1)

# Save outputs
out_pdf <- file.path(dirname(ALIGN_STATS_FILE), "mycoplasma_alignment_summary.pdf")
out_png <- file.path(dirname(ALIGN_STATS_FILE), "mycoplasma_alignment_summary.png")

ggsave(out_pdf, final_plot, width = 7, height = 2.5 * length(plots))
ggsave(out_png, final_plot, width = 7, height = 2.5 * length(plots))

EOF

echo "[INFO] Mycoplasma detection complete for $SAMPLE"
echo "[INFO] Alignment stats: $ALIGN_STATS_FILE"
echo "[INFO] PDF/PNG plots: $SAMPLE_OUT/mycoplasma_alignment_summary.pdf, $SAMPLE_OUT/mycoplasma_alignment_summary.png"