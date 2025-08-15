#DESeq
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)

setwd("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/")
#abundance (TPM) data:
tpm<-readRDS("ipsc_object_merged_pca_tsne_umap_log2norm_tpm.rds")
tpm_matrix <- GetAssayData(tpm, assay = "RNA", layer = "counts")

#counts (RSEM expected_counts) data:
expected_counts <- readRDS("ipsc_858samples_expected_counts_with_logenorm_peer_factors_cellcycle_and_vst_matrix.rds")

print("Extracting counts assay from object...")
counts_matrix <- GetAssayData(expected_counts, assay = "RNA", layer = "counts")

#deseq wants integers:
counts_matrix <- round(counts_matrix)

print("Extracting PEER-corrected matrix...")
peer_mat <- GetAssayData(expected_counts, assay = "PEER", slot = "counts")

print("Filtering counts matrix on counts > 1 in at least 10 samples...")
common_genes <- intersect(rownames(counts_matrix), rownames(peer_mat))

# Subset both matrices to the shared genes
counts_matrix <- counts_matrix[common_genes, ]

col_data <- expected_counts@meta.data
all(colnames(counts_matrix) == rownames(col_data))  # Should return TRUE

print("Creating SummarizedExperiment object...")
se <- SummarizedExperiment(
  assays = list(counts = counts_matrix),# tpm = tpm_matrix),
  colData = col_data)

print("Creating DESeqDataSet object...")
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = col_data,
  design = ~ Project + InferredCov1 +InferredCov2 +InferredCov3 + InferredCov4 + InferredCov5 + 
    InferredCov6 + InferredCov7 + InferredCov8 + InferredCov9 + InferredCov10 + Phase
)

print("Running DESeq...")
dds <- DESeq(dds)

print("Saving dds object...")
saveRDS(dds, file = "dds_ipsc_858samples_with_10peer_factorsand_cellcycle_phase.rds")

print("Normalizing counts...")
norm_counts <- counts(dds, normalized = TRUE)

print("Calculating variance stabilizing transformation (VST)...")
vsd <- vst(dds, blind = FALSE) # Or use rlog(dds) if dataset is small (<30)

print("Saving VST object...")
saveRDS(vsd, file = "vsd_ipsc_858samples_with_10peer_factors_and_cellcycle_phase.rds")#
print("Done!")

vsd<-readRDS("vsd_ipsc_858samples_with_10peer_factors_and_cellcycle_phase.rds")
dds<-readRDS("dds_ipsc_858samples_with_10peer_factorsand_cellcycle_phase.rds")

mat  <- assay(vsd)

#design to keep for PCA
design_keep <- model.matrix(~ Project, colData(dds))

#covar_names <- c(paste0("InferredCov", 1:10),"Project","Phase")
covar_names <- c("Project","Phase")

#convert covars into numeric matrix for limma
covar_df <- as.data.frame(colData(dds)[, covar_names, drop = FALSE])
covar_mat <- do.call(cbind, lapply(covar_df, function(x) {
  if (is.factor(x) || is.character(x)) {
    model.matrix(~ x - 1)  # one-hot encode
  } else {
    matrix(as.numeric(x), ncol = 1)
  }
}))

#remove batch effects for corrected matrix with limma
library(limma)
mat_adj <- removeBatchEffect(mat,
                             design = design_keep,
                             covariates = covar_mat)

#create a DESeqTransform object with adjusted assay
vs_adj <- vsd
assay(vs_adj) <- mat_adj

#remove mito genes before PCA because they can dominate PCs (but I want to see the underlying genes!!)
mt_genes <- grep("^MT-", rownames(vs_adj), value = TRUE)
mat_noMT <- assay(vs_adj)[!rownames(vs_adj) %in% mt_genes, ]

#run PCA on top 2000 most variable non-mito genes
topVarGenes <- head(order(rowVars(mat_noMT), decreasing = TRUE), 2000)
pca <- prcomp(t(mat_noMT[topVarGenes, ]), center = TRUE)

#extract gene loadings for PC1 and PC2 (rows=genes & cols=PCs)
gene_loadings <- pca$rotation

pc_1_2<-data.frame(gene_loadings)[1:2]
pc_1_2$gene<-rownames(pc_1_2)

bulk_ipsc <- readRDS("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ipsc_858samples_expected_counts_with_logenorm_peer_factors_cellcycle_and_vst_matrix.rds")
#create new assay for corrected matrix
bulk_ipsc[["AdjExpr"]] <- CreateAssayObject(counts = mat_adj)

DefaultAssay(bulk_ipsc) <- "AdjExpr"

#run pca with seurat for object
bulk_ipsc <- NormalizeData(bulk_ipsc)       # optional if not already normalized
bulk_ipsc <- FindVariableFeatures(bulk_ipsc)
bulk_ipsc <- ScaleData(bulk_ipsc, features = rownames(bulk_ipsc))
bulk_ipsc <- RunPCA(bulk_ipsc)
DimPlot(bulk_ipsc, reduction = "pca", group.by=c("Project"))

ElbowPlot(bulk_ipsc)

#run UMAP (did NOT end up liking this)
bulk_ipsc<-FindNeighbors(bulk_ipsc,dims=1:10)
bulk_ipsc<-FindClusters(bulk_ipsc,resolution = 0.05)
bulk_ipsc<-RunUMAP(bulk_ipsc, dims=1:10)
DimPlot(bulk_ipsc, group.by=c("Project","differentiation_status"), reduction="umap")

#Run regressions on differentiation status
#First I'm removing samples that don't have a differentiation status
metadata_status<-subset(bulk_ipsc@meta.data, differentiation_status=="Successful" | differentiation_status=="Failed")

#Metadata with PC1 and PC2 embeddings, subset to just hipsci lines with available diff. status
df<-metadata_status
df$status<-ifelse(df$differentiation_status=="Failed", 0, 1)

df$differentiation_status <- factor(df$differentiation_status, levels = c("Failed", "Successful"))

#Logistic regression
model <- glm(differentiation_status ~ Formative_TPM + Primed_TPM + PC_1 + PC_2,
             data = df,
             family = binomial)

summary(lm(status~Formative_TPM + PC_1 + PC_2, data=df))

summary(model)





###Extra code for PCA and VST plots. Stuff I tried but didn't end up using
#Keeping because may be useful later?



# phase_matrix <- model.matrix(~ Phase - 1, colData(dds))
# covars <- cbind(covars, phase_matrix)

# mat_adj <- limma::removeBatchEffect(mat, design = design_keep, covariates = as.matrix(covars))

# pca <- prcomp(t(mat_adj), center = TRUE, scale. = FALSE)
# plot(pca$x[,1:2])


library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
  #as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst")#, "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

#Calculate Euclidean distances
sampleDists <- dist(t(assay(vsd)))

plotPCA(vsd, intgroup = c("Project"))+ggtitle("PCA on VST normalized counts")
plotPCA(vsd, intgroup = c("Project"), ntop=2000)+ggtitle("PCA on VST normalized counts")


mat<-assay(vsd)

# pick top variable genes
top_var_genes <- head(order(rowVars(mat), decreasing = TRUE), 2000)

# run PCA
pca <- prcomp(t(mat[top_var_genes, ]), scale. = TRUE)

loadings <- pca$rotation  # genes x PCs
loadings_df <- as.data.frame(pca$rotation[, 1:2])
loadings_df$gene <- rownames(loadings_df)

#glmPCA, which I didn't really end up liking
install.packages("glmpca")
library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$Project <- dds$Project
gpca.dat$Sample <- dds$Sample

ggplot(gpca.dat, aes(x = dim1, y = dim2, color = Project)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")



ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Project","Project")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

plotDispEsts(dds)


bulk_ipsc <- readRDS("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ipsc_858samples_expected_counts_with_logenorm_peer_factors.rds")

vst_mat <- assay(vsd)
vst_mat <- vst_mat[, colnames(bulk_ipsc)]
bulk_ipsc[["VST"]] <- CreateAssayObject(counts = vst_mat)

DefaultAssay(bulk_ipsc) <- "VST"

bulk_ipsc<-FindVariableFeatures(bulk_ipsc)
bulk_ipsc <- ScaleData(bulk_ipsc, features = rownames(bulk_ipsc),
                        do.scale = FALSE, do.center = TRUE)
bulk_ipsc<-RunPCA(bulk_ipsc)
DimPlot(bulk_ipsc)
FeaturePlot(bulk_ipsc, features=c("POU5F1","MAP2"), cols=c("grey","red"))


pca_coords <- as.data.frame(pca$x)  # samples x PCs
pca_coords$cell <- rownames(pca_coords)

pca_coords <- pca_coords[match(colnames(bulk_ipsc), pca_coords$cell), ]

bulk_ipsc[["pca_deseq"]] <- CreateDimReducObject(
  embeddings = as.matrix(pca_coords[, grep("^PC", colnames(pca_coords))]),
  key = "PC_",
  assay = DefaultAssay(bulk_ipsc)
)

DimPlot(bulk_ipsc, reduction = "pca_deseq", dims = c(1, 2))

bulk_ipsc <- FindNeighbors(bulk_ipsc, dims = 1:10)
bulk_ipsc <- FindClusters(bulk_ipsc, resolution = 0.05)
bulk_ipsc <- RunUMAP(bulk_ipsc, dims = 1:10)

DimPlot(bulk_ipsc, group.by=c("seurat_clusters","Project"))


#cell cycle accelerators (increased expression -> cancerous)
FeaturePlot(bulk_ipsc, features=c("MYC","CCND1","CCNE1","E2F1","E2F3","AURKA","AURKB","PLK1"), reduction="pca")

#Tumor suppressors (decreased expression -> cancerous)
FeaturePlot(bulk_ipsc, features=c("TP53","CDKN2A","RB1","CDKN1A","CDKN1B"), reduction="pca")

#DNA repair and damage tolerance 
FeaturePlot(bulk_ipsc, features=c("BRCA1","BRCA2","MLH1","MSH2","MSH6","PMS2","RAD51","RAD51C"), reduction="pca")
