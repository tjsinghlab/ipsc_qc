#DESeq
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)

#abundance (TPM) data:
tpm<-readRDS("ipsc_object_merged_pca_tsne_umap_log2norm_tpm.rds")
tpm_matrix <- GetAssayData(tpm, assay = "RNA", layer = "counts")

#counts (RSEM expected_counts) data:
#expected_counts<-readRDS("bulk_ipsc_with_timepoint_deconvolution_cuomo.rds")
expected_counts<-readRDS("ipsc_858samples_expected_counts_with_logenorm_peer_factors.rds")

counts_matrix <- GetAssayData(expected_counts, assay = "RNA", layer = "counts")

#deseq wants integers:
counts_matrix <- round(counts_matrix)

peer_mat <- GetAssayData(seurat_obj, assay = "PEER", slot = "counts")

common_genes <- intersect(rownames(counts_matrix), rownames(peer_mat))

# Subset both matrices to the shared genes
counts_matrix <- counts_matrix[common_genes, ]

col_data <- expected_counts@meta.data
all(colnames(counts_matrix) == rownames(col_data))  # Should return TRUE

se <- SummarizedExperiment(
  assays = list(counts = counts_matrix),# tpm = tpm_matrix),
  colData = col_data
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = col_data,
  design = ~ Project + InferredCov1 +InferredCov2 +InferredCov3 + InferredCov4 + InferredCov5 + 
    InferredCov6 + InferredCov7 + InferredCov8 + InferredCov9 + InferredCov10
)

dds <- DESeq(dds)

smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
nrow(dds)

norm_counts <- counts(dds, normalized = TRUE)

vsd <- vst(dds, blind = FALSE) # Or use rlog(dds) if dataset is small (<30)
head(assay(vsd), 3)

#corrected_mat <- assay(vsd)

# rld <- rlog(dds, blind = FALSE)
# head(assay(rld), 3)


lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
meanSdPlot(cts, ranks = FALSE)

log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)

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


