#Prepare lognorm counts for PEER analysis
library(Seurat)

setwd("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc")

#counts (RSEM expected_counts) data:
expected_counts<-readRDS("bulk_ipsc_with_timepoint_deconvolution_cuomo.rds")

#Exporting lognorm matrix (log counts+1) for PEER:
lognorm_matrix <- GetAssayData(expected_counts, assay = "RNA", layer = "data")
lognorm_matrix[is.na(lognorm_matrix)] <- 0

#Keep only genes with counts >= 1 in at least 10 samples
keep_genes <- rowSums(counts_matrix >= 1) >= 10
lognorm_matrix <- lognorm_matrix[keep_genes, ]
df_out <- data.frame(GeneID = rownames(lognorm_matrix), lognorm_matrix, check.names = FALSE)

# Write to tab-delimited file
write.table(df_out,
            file = "lognorm_for_peer.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


#use lognorm_for_peer.tsv as input for run_peer.py script

#Analyze results from run_peer.py:
library(stringr)
library(ggplot2)
library(dplyr)

#adjusted counts matrix from peer:
peer_residuals<-read.table("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/peer/iPSC.PEER_residuals.txt", header = T, row.names = 1)

peer_mtx<-as.matrix(peer_residuals)

colnames(peer_mtx)<-str_replace_all(colnames(peer_mtx), "[.]", "-")
peer_assay <- CreateAssayObject(counts = peer_mtx)

expected_counts<-readRDS("bulk_ipsc_with_timepoint_deconvolution_cuomo.rds")
expected_counts@assays$PEER<-peer_assay

DefaultAssay(expected_counts) <- "PEER"

expected_counts@assays$PEER@data<-expected_counts@assays$PEER@counts
expected_counts <- FindVariableFeatures(expected_counts)
expected_counts <- ScaleData(expected_counts)
expected_counts <- RunPCA(expected_counts)

# Replace with your actual Seurat object
seurat_obj <- expected_counts

# Get gene expression matrices
rna_mat <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
peer_mat <- GetAssayData(seurat_obj, assay = "PEER", slot = "counts")

common_genes <- intersect(rownames(rna_mat), rownames(peer_mat))

# Subset both matrices to the shared genes
rna_mat <- rna_mat[common_genes, ]
peer_mat <- peer_mat[common_genes, ]

# Calculate per-gene variance
rna_var <- apply(rna_mat, 1, var)
peer_var <- apply(peer_mat, 1, var)

# Combine into one tidy dataframe
var_df <- data.frame(
  gene = names(rna_var),
  RNA = rna_var,
  PEER = peer_var
) %>%
  tidyr::pivot_longer(cols = c("RNA", "PEER"), names_to = "Assay", values_to = "Variance")

# Plot violin plots of variance distribution
ggplot(var_df, aes(x = Assay, y = Variance, fill = Assay)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  scale_y_log10() +  # Optional: log scale to show variance spread better
  theme_minimal() +
  labs(
    title = "Gene-wise Variance: RNA vs. PEER Residuals",
    x = "Assay",
    y = "Variance (log scale)"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )


#peer factors:
peer_covars<-read.table("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/peer/iPSC.PEER_covariates.txt", header = T, row.names = 1)
peer_covars<-t(peer_covars)
peer_covars<-data.frame(peer_covars)
peer_covars$Sample<-rownames(peer_covars)

metadata<-expected_counts@meta.data
metadata<-left_join(metadata, peer_covars)

rownames(metadata)<-metadata$Sample
expected_counts@meta.data<-metadata

saveRDS(expected_counts, "ipsc_858samples_expected_counts_with_logenorm_peer_factors.rds")

##how cellcycle score is calculated:
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

bulk_ipsc <- CellCycleScoring(bulk_ipsc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(bulk_ipsc, group.by=c("Project","differentiation_status","Phase"), reduction="pca")
saveRDS(bulk_ipsc, "ipsc_858samples_expected_counts_with_logenorm_peer_factors_cellcycle_and_vst_matrix.rds")

#Run regressions with gene expression as response and peer-corrected PCs as predictors
expected_counts<-readRDS("ipsc_858samples_expected_counts_with_logenorm_peer_factors.rds")
meta<-expected_counts@meta.data
peer_mat <- GetAssayData(expected_counts, assay = "PEER", slot = "counts")

pca <- prcomp(t(peer_mat), center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca$x)  # scores (principal components)

pca_df$Project <- meta$Project
#For one gene:
gene_of_interest <- "POU5F1"

expr <- peer_mat[gene_of_interest, ]  # returns a vector of expression per sample
df <- cbind(pca_df, POU5F1_expr = expr[rownames(pca_df)])  # match samples

model <- lm(POU5F1_expr ~ PC1 + PC2 + PC3 + PC4 + PC5, data = df)
summary(model)

ggplot(df, aes(x = PC2, y = POU5F1_expr)) +
  geom_point(color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
  labs(title = "POU5F1 Expression vs. PC2",
       x = "PC2", y = "POU5F1 Expression") +
  theme_minimal()
