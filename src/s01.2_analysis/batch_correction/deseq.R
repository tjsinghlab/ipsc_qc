library(Seurat)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(ggplot2)
setwd("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc")

print("Reading in the expected counts matrix with 10 PEER factors...")
#counts (RSEM expected_counts) data:
#expected_counts<-readRDS("bulk_ipsc_with_timepoint_deconvolution_cuomo.rds")
expected_counts<-readRDS("ipsc_858samples_expected_counts_with_logenorm_peer_factors.rds")

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
  colData = col_data
)

print("Creating DESeqDataSet object...")
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = col_data,
  design = ~ Project + InferredCov1 +InferredCov2 +InferredCov3 + InferredCov4 + InferredCov5 + 
    InferredCov6 + InferredCov7 + InferredCov8 + InferredCov9 + InferredCov10
)

print("Running DESeq...")
dds <- DESeq(dds)

print("Saving dds object...")
saveRDS(dds, file = "dds_ipsc_858samples_with_10peer_factors.rds")

print("Filtering low count genes...")
smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

print("Normalizing counts...")
norm_counts <- counts(dds, normalized = TRUE)

print("Calculating variance stabilizing transformation (VST)...")
vsd <- vst(dds, blind = FALSE) # Or use rlog(dds) if dataset is small (<30)

print("Saving VST object...")
saveRDS(vsd, file = "vsd_ipsc_858samples_with_10peer_factors.rds")

print("Done!")

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
#pcaData <- plotPCA(vsd, intgroup = c("Project"), returnData = TRUE)
#percentVar <- round(100 * attr(pcaData, "percentVar"))

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
 
#Add VST corrected matrix to original Seurat object and re-do clustering/UMAP
bulk_ipsc <- readRDS("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ipsc_858samples_expected_counts_with_logenorm_peer_factors.rds")

vst_mat <- assay(vsd)
vst_mat <- vst_mat[, colnames(bulk_ipsc)]
bulk_ipsc[["VST"]] <- CreateAssayObject(counts = vst_mat)

DefaultAssay(bulk_ipsc) <- "VST"

bulk_ipsc<-FindVariableFeatures(bulk_ipsc)
bulk_ipsc <- ScaleData(bulk_ipsc, features = rownames(bulk_ipsc),
                        do.scale = FALSE, do.center = TRUE) #additional scaling not necessary for VST corrected counts
bulk_ipsc<-RunPCA(bulk_ipsc)
DimPlot(bulk_ipsc)
FeaturePlot(bulk_ipsc, features=c("POU5F1","MAP2"), cols=c("grey","red"))

bulk_ipsc <- FindNeighbors(bulk_ipsc, dims = 1:10)
bulk_ipsc <- FindClusters(bulk_ipsc, resolution = 0.05)
bulk_ipsc <- RunUMAP(bulk_ipsc, dims = 1:10)

DimPlot(bulk_ipsc, group.by=c("seurat_clusters","Project"))

#Pluripotent:
FeaturePlot(bulk_ipsc, features=c("NANOG","PODXL","POU5F1","CNMD","SOX2","B3GALT5","LIN28A","LIN28B","DPPA4","ZFP42"), cols=c("gray","red"), reduction="pca")

#Developing neurons:
FeaturePlot(scrna, features=c("ASCL1","NEUROG2","GFAP","EOMES","PROX1","NCAM1","DCX","CALB2","NEUROD1","SOX10","SOX9","SOX1"), cols=c("grey","red"))

#Mature neurons:
FeaturePlot(scrna, features=c("SLC32A1","RBFOX3","MAP2","NEFL","SYN1","SLC17A7","SLC17A6","PAX6","SNAP25","CHAT","CUX2","ATP1B1"), cols=c("grey","red"))


