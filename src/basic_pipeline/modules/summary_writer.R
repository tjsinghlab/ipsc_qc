#!/usr/bin/env Rscript

# ===========================================
# Construct final summary table (generalized genes)
# ===========================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

# -------------------------------------------
# Command-line options
# -------------------------------------------
option_list <- list(
  make_option(c("--output_dir"), type="character", help="Output directory containing per-sample results"),
  make_option(c("--pacnet_file"), type="character", help="Path to PACNet classification_scores.csv"),
  make_option(c("--samples"), type="character", help="Comma-separated list of sample names"),
  make_option(c("--project"), type="character", default="Project", help="Project name (optional)")
)

opt <- parse_args(OptionParser(option_list=option_list))
OUTPUT_DIR <- opt$output_dir
PACNET_FILE <- opt$pacnet_file
SAMPLES <- strsplit(opt$samples, ",")[[1]]
#SAMPLES<-c("24246R-06-01_S0_L001", "24246R-06-08_S415_L008", "24246R-06-11_S198_L008")

# -------------------------------------------
# Discover all genes across all samples
# -------------------------------------------
all_genes <- character()
for (sample in SAMPLES) {
  cosmic_file <- file.path(OUTPUT_DIR, sample, "cosmic_calling", paste0(sample, "_CancerMutations_summary.tsv"))
  if (file.exists(cosmic_file)) {
    df <- fread(cosmic_file)
    if (ncol(df) >= 3) {
      setnames(df, c("GENE_SYMBOL", "unique_pos_count", "MUTATION_DESCRIPTION"))
      all_genes <- union(all_genes, unique(df$GENE_SYMBOL))
    }
  }
}

# -------------------------------------------
# Prepare summary table
# -------------------------------------------
summary_list <- list()

for (sample in SAMPLES) {
  sample_dir <- file.path(OUTPUT_DIR, sample)
  
  # --------------------------
  # Mycoplasma %
  # --------------------------
  myco_file <- file.path(sample_dir, "mycoplasma/mycoplasma_alignment_stats.tsv")
  myco_pct <- NA
  if (file.exists(myco_file)) {
    df <- fread(myco_file, skip=1)
    if (nrow(df) > 0) {
      myco_pct <- sum(as.numeric(df[[5]]))
    }
  }
  
  # --------------------------
  # Fraction abnormal karyotype
  # --------------------------
  karyo_file <- file.path(sample_dir, "eSNPKaryotyping", paste0(sample, "_variantTable.csv"))
  frac_karyo <- NA
  if (file.exists(karyo_file)) {
    df <- fread(karyo_file)
    total <- nrow(df)
    abnormal <- sum(df[[6]] != 10)
    if (total > 0) frac_karyo <- round(abnormal/total, 6)
  }
  
  # --------------------------
  # PACNet esc score
  # --------------------------
  pacnet_scores <- fread(PACNET_FILE)
  esc_scores <- subset(pacnet_scores, V1=="esc")

  print(names(esc_scores))
  print(sample)

  esc_score <- NA
  col <- grep(sample, names(esc_scores), value = TRUE)

  if (length(col) == 1) {
      esc_score <- esc_scores[[col]]
  } else if (length(col) > 1) {
      warning(paste("Multiple PACNet columns matched sample:", sample))
      esc_score <- esc_scores[[col[1]]]
  } else {
      warning(paste("No PACNet ESC score column matched sample:", sample))
      esc_score <- NA
  }
  
  # --------------------------
  # COSMIC mutations (dynamic genes)
  # --------------------------
  cosmic_file <- file.path(sample_dir, "cosmic_calling", paste0(sample, "_CancerMutations_summary.tsv"))
  
  frame_counts <- setNames(rep(0, length(all_genes)), all_genes)
  total_counts <- setNames(rep(0, length(all_genes)), all_genes)
  
  if (file.exists(cosmic_file)) {
    df <- fread(cosmic_file)
    if (ncol(df) >= 3) {
      setnames(df, c("GENE_SYMBOL", "MUTATION_DESCRIPTION", "unique_pos_count" ))
      for (gene in all_genes) {
        sub <- df[GENE_SYMBOL == gene]
        total_counts[gene] <- sum(as.numeric(sub$unique_pos_count))
        frame_counts[gene] <- sum(as.numeric(sub$unique_pos_count[grep("frame", sub$MUTATION_DESCRIPTION, ignore.case=TRUE)]))
      }
    }
  }
  
  # --------------------------
  # Construct row
  # --------------------------
  row <- c(
    Sample = sample,
    Total_Mycoplasma_Percent_Alignment = myco_pct,
    Fraction_Abnormal_Karyotype = frac_karyo,
    esc_score = esc_score
  )
  
  # Append dynamic genes
  for (gene in all_genes) {
    row[paste0(gene, "_Frame_Mutations")] <- frame_counts[gene]
    row[paste0(gene, "_Total_Mutations")] <- total_counts[gene]
  }
  
  summary_list[[sample]] <- row
}

# -------------------------------------------
# Write summary
# -------------------------------------------
summary_df <- rbindlist(lapply(summary_list, as.list), fill=TRUE)
output_file <- file.path(OUTPUT_DIR, "final_summary_table.tsv")
fwrite(summary_df, output_file, sep="\t", na="NA", quote=FALSE)

# -------------------------------------------
# Combined heatmap: Mutations (top) + QC Metrics (bottom)
# -------------------------------------------

library(ggplot2)
library(reshape2)
library(patchwork)

df <- summary_df
rnamx <- df$Sample
df$Sample <- NULL
df <- df %>% mutate_if(is.character, as.numeric)

mat <- as.matrix(df)
rownames(mat) <- rnamx

# Split columns by type
mut_cols   <- grep("Mutations$", colnames(mat), value = TRUE)
other_cols <- setdiff(colnames(mat), mut_cols)

mat_mut   <- mat[, mut_cols, drop = FALSE]
mat_other <- mat[, other_cols, drop = FALSE]

# Melt for ggplot2
df_mut   <- melt(mat_mut)
df_other <- melt(mat_other)

colnames(df_mut)   <- c("Sample", "Metric", "Value")
colnames(df_other) <- c("Sample", "Metric", "Value")

# dynamic axis-label scaling
num_samples <- length(SAMPLES)
if (num_samples <= 20) {
  x_cex <- 12
} else if (num_samples <= 50) {
  x_cex <- 10
} else if (num_samples <= 100) {
  x_cex <- 8
} else {
  x_cex <- 6
}

y_cex <- 12

# Color palette
cools <- colorRampPalette(c("black", "limegreen", "yellow"))(100)

# MUTATION HEATMAP (top)
p_mut <- ggplot(df_mut, aes(x = Sample, y = Metric, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = cools, name = "Mutations") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = x_cex),
    axis.text.y = element_text(size = y_cex),
    plot.title  = element_text(size = 16, face = "bold")
  ) +
  labs(title = "Mutation Counts", x = "Samples", y = "Mutation Metrics")

# QC SCORES HEATMAP (bottom, fixed 0–1)
p_other <- ggplot(df_other, aes(x = Sample, y = Metric, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = cools, limits = c(0, 1), name = "QC Score") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = x_cex),
    axis.text.y = element_text(size = y_cex),
    plot.title  = element_text(size = 16, face = "bold")
  ) +
  labs(title = "QC Metrics (0–1)", x = "Samples", y = "QC Metrics")

# STACKED vertically: top = mutations, bottom = QC
combined_plot <- p_mut / p_other + 
  plot_layout(heights = c(1, 1.2))

# Output one PNG
png(
  file = file.path(OUTPUT_DIR, "final_summary_heatmap.png"),
  width = max(12, num_samples * 0.5),
  height = 22,
  units = "in",
  res = 300,
  type = "cairo"
)

print(combined_plot)
dev.off()

message("[INFO] Summary files written to output path.")


