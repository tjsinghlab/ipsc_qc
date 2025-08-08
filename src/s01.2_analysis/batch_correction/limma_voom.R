#limma voom
library(Seurat)
library(limma)
library(dplyr)
library(edgeR)
library(ggplot2)

setwd("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc")

expected_counts<-readRDS("ipsc_858samples_expected_counts_with_logenorm_peer_factors.rds")

DefaultAssay(expected_counts)<-"RNA"

meta <- expected_counts@meta.data

#Extract first 10 peer factors
covariates <- paste0("InferredCov", 1:10)

design_formula <- as.formula(paste("~ Project +", paste(covariates, collapse = " + ")))

design <- model.matrix(design_formula, data = meta)
# design <-  model.matrix(~Project + InferredCov1 +InferredCov2 +InferredCov3 + InferredCov4 + InferredCov5 + 
#                           InferredCov6 + InferredCov7 + InferredCov8 + InferredCov9 + InferredCov10)

counts <- as.matrix(GetAssayData(expected_counts, layer = "counts"))
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

v <- voom(dge, design, plot=T)

v$E_corrected <- removeBatchEffect(v$E, batch = meta$Project, design = design)
v$E_corrected<-NULL

fit <- lmFit(v, design)
fit <- eBayes(fit)
topTable(fit, coef = "ProjectNeuroKaire infected lines", adjust = "fdr")

plotSA(fit, main = "Residual SD vs Average Expression")
volcanoplot(fit, coef = "ProjectNeuroKaire infected lines", highlight = 10, names = rownames(fit))

pca <- prcomp(t(v$E), center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca$x)  # scores (principal components)

pca_df$Project <- meta$Project

ggplot(pca_df, aes(x = PC1, y = PC2, color = Project)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA of limma voom Corrected Expression",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")) +
  theme_minimal()

scree_df <- data.frame(
  PC = paste0("PC", 1:10),
  Variance = summary(pca)$importance[2, 1:10]
)

ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_col(fill = "skyblue") +
  labs(title = "limma batch-corrected PCs", y = "Proportion of Variance Explained") +
  theme_minimal()



#Regressions on gene expression by PCs

#For one gene:
gene_of_interest <- "NANOG"

expr <- v$E[gene_of_interest, ]  # returns a vector of expression per sample
df <- cbind(pca_df, POU5F1_expr = expr[rownames(pca_df)])  # match samples

model <- lm(POU5F1_expr ~ PC1 + PC2 + PC3 + PC4 + PC5, data = df)
summary(model)

ggplot(df, aes(x = PC2, y = POU5F1_expr)) +
  geom_point(color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
  labs(title = "POU5F1 Expression vs. PC2",
       x = "PC2", y = "POU5F1 Expression") +
  theme_minimal()
