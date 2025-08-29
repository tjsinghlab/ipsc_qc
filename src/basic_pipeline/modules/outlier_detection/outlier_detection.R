library(edgeR)
library(dplyr)
library(ggplot2)
library(patchwork)

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-p", "--project"), type="character", default="default_project",
              help="Project name", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=file.path(cwd, "outputs"),
              help="Path to desired output directory", metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

project_name <- opt$project
output_dir <- opt$output_dir


counts<-read.csv(paste0(output_dir, "bd2_matrix_from_split.csv"),row.names=1)

cpm_counts <- cpm(counts)   # counts per million

# Step 2: Log-transform
log_counts <- log2(cpm_counts + 1)

# Step 3: Remove zero-variance genes
log_counts_filtered <- log_counts[apply(log_counts, 1, var) > 0, ]

# assume counts is a matrix with genes as rows and samples as columns
# transpose so PCA runs on samples
pca <- prcomp(t(log_counts_filtered), scale. = TRUE)

# extract first few PCs for outlier detection
pcs <- pca$x[, 1:5]   # choose how many PCs you want

# compute Mahalanobis distance
md <- mahalanobis(pcs, colMeans(pcs), cov(pcs))

# threshold (e.g. chi-square with df = number of PCs, p < 0.01)
cutoff <- qchisq(0.99, df = ncol(pcs))
outliers_pca <- md > cutoff

# plot PCA with outliers highlighted
library(ggplot2)
df <- data.frame(PC1 = pcs[,1], PC2 = pcs[,2], 
                 Outlier = outliers_pca, Sample = rownames(pcs))
ggplot(df, aes(x = PC1, y = PC2, color = Outlier, label = Sample)) +
  geom_point(size=3) +
  geom_text(hjust=1.2, vjust=1.2, size=3) +
  theme_minimal()


#read in pacnet classification scores
scores<-read.csv(paste0(output_dir, "classification_scores.csv"), row.names=1)

#remove random samples
scores <- scores[, !grepl("rand", colnames(scores))]
scores_t <- t(scores)

scores_t<- scores_t[, apply(scores_t, 2, var) > 0, drop = FALSE]

pcs_scores <- prcomp(scores_t, scale. = TRUE)$x[, 1:5]  # first 5 PCs
md_scores <- mahalanobis(pcs_scores,
                         colMeans(pcs_scores),
                         cov(pcs_scores))

cutoff_scores <- qchisq(0.99, df = ncol(scores_t))
outliers_scores <- md_scores > cutoff_scores

df <- data.frame(PC1 = pcs_scores[,1], PC2 = pcs_scores[,2], 
                 Outlier = outliers_scores, Sample = rownames(pcs_scores))
ggplot(df, aes(x = PC1, y = PC2, color = Outlier, label = Sample)) +
  geom_point(size=3) +
  geom_text(hjust=1.2, vjust=1.2, size=3) +
  theme_minimal()
