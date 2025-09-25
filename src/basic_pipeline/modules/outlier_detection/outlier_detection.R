#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# -------------------------------
# Command-line options
# -------------------------------
option_list <- list(
  make_option(c("-p", "--project"), type="character", default="default_project",
              help="Project name", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="Sample name", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL,
              help="Path to output directory", metavar="character")
)

opt <- parse_args(OptionParser(option_list = option_list))
project_name <- opt$project
sample <- opt$sample
output_dir <- opt$output_dir
sample_outdir <- file.path(output_dir, "pacnet")

message("[INFO] Processing sample: ", sample)

# -------------------------------
# Load per-sample count matrix
# -------------------------------
counts_file  <- file.path(output_dir, "pacnet/query_matrix.csv")
if (!file.exists(counts_file)) stop("Counts file not found: ", counts_file)
counts <- read.csv(counts_file, row.names = 1)
counts_new<-na.omit(counts)
cpm_counts <- cpm(counts_new)
log_counts <- log2(cpm_counts + 1)
log_counts_filtered <- log_counts[apply(log_counts, 1, var) > 0, ]

# PCA on counts
pca_counts <- prcomp(t(log_counts_filtered), scale. = TRUE)
pc_count<-min(ncol(log_counts_filtered), 5)
pcs_counts <- pca_counts$x[, 1:pc_count]

# Detect outliers using Mahalanobis distance
md_counts <- mahalanobis(pcs_counts, colMeans(pcs_counts), cov(pcs_counts))
cutoff_counts <- qchisq(0.99, df = ncol(pcs_counts))
outliers_counts <- md_counts > cutoff_counts

df_counts <- data.frame(
  PC1 = pcs_counts[,1], PC2 = pcs_counts[,2],
  Outlier = outliers_counts, Sample = rownames(pcs_counts)
)

pca_plot_counts <- ggplot(df_counts, aes(x = PC1, y = PC2, color = Outlier, label = Sample)) +
  geom_point(size = 3) +
  geom_text(hjust = 1.2, vjust = 1.2, size = 3) +
  theme_minimal() +
  labs(title = paste(project_name, "-", sample, "- PCA: Expression Counts"))

# Save counts PCA
counts_pdf <- file.path(output_dir, paste0(project_name, "_", sample, "_PCA_counts.pdf"))
ggsave(filename = counts_pdf, plot = pca_plot_counts, width = 8, height = 6)

# -------------------------------
# PCA on PACNet ESC scores
# -------------------------------
scores_file <- file.path(output_dir, "pacnet/classification_scores.csv")
if (!file.exists(scores_file)) stop("PACNet scores file not found: ", scores_file)
scores <- read.csv(scores_file, row.names = 1)

# Remove random samples and low-variance features
scores <- scores[, !grepl("rand", colnames(scores))]
scores_t <- t(scores)
scores_t <- scores_t[, apply(scores_t, 2, var) > 0, drop = FALSE]

pcs_scores <- prcomp(scores_t, scale. = TRUE)$x[, 1:5]
md_scores <- mahalanobis(pcs_scores, colMeans(pcs_scores), cov(pcs_scores))
cutoff_scores <- qchisq(0.99, df = ncol(pcs_scores))
outliers_scores <- md_scores > cutoff_scores

df_scores <- data.frame(
  PC1 = pcs_scores[,1], PC2 = pcs_scores[,2],
  Outlier = outliers_scores, Sample = rownames(pcs_scores)
)

pca_plot_scores <- ggplot(df_scores, aes(x = PC1, y = PC2, color = Outlier, label = Sample)) +
  geom_point(size = 3) +
  geom_text(hjust = 1.2, vjust = 1.2, size = 3) +
  theme_minimal() +
  labs(title = paste(project_name, "-", sample, "- PCA: PACNet ESC Scores"))

# Save PACNet PCA
scores_pdf <- file.path(output_dir, paste0(project_name, "_", sample, "_PCA_pacnet_scores.pdf"))
ggsave(filename = scores_pdf, plot = pca_plot_scores, width = 8, height = 6)

message("[INFO] PCA plots saved for sample: ", sample)
message("[INFO] Counts PCA: ", counts_pdf)
message("[INFO] PACNet PCA: ", scores_pdf)
