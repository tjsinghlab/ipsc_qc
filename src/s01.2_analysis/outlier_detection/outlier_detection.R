#Outlier detection
library(e1071)
library(dplyr)
library(tidyverse)
library(stringr)
library(patchwork)
library(ggplot2)
library(Seurat)
base_dir<-"/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/"
setwd(base_dir)

#Read in and format data
#iPSC Seurat objects
merged<-readRDS("ipsc_merged_objects_with_may_data_esc_res_pca_tsne.rds")
#merged<-readRDS("ipsc_and_hart30_merged_tsne_pca.rds")
#merged<-readRDS("merged_ipsc_hipsci_hart_bd2_data_with_esc_scores_and_pca.rds")

# pac_res<-read.csv("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ron_data_may2025/classification_results_2025-05-15.csv")
# colnames(pac_res)<-str_split_i(colnames(pac_res),"X",2)
# colnames(pac_res)<-str_replace_all(colnames(pac_res),"[.]","-")
# colnames(pac_res)[1]<-"X"
# pac_list<-lapply(list.files("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/files_from_pacnet/", full.names = T), function(x){
#   x<-read.csv(x)
#   return(x)
# })
# names(pac_list)<-list.files("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/files_from_pacnet/")
# pac_list$new_ron_data<-pac_res
# 
# pac_list2<-lapply(pac_list, function(df){
#   rownames(df)<-df$X
#   df$X<-NULL
#   df<-t(df)
#   df <- df[1:(nrow(df) - 3), ]
#   df<-data.frame(df)
#   df$sample_name<-rownames(df)
#   return(df)
# })
# 
# obj<-Reduce(rbind, pac_list2)
# obj$Sample<-obj$sample_name
# #merged<-readRDS("ipsc_and_hart30_merged_tsne_pca.rds")
# #merged@meta.data<-left_join(merged@meta.data, obj, "Sample")
# #merged2 <- subset(merged, subset = !is.na(Sample))
# met<-merged@meta.data
# met<-left_join(met, obj, "Sample")
# merged$esc<-met$esc
# merged$neuron<-met$neuron

# merged<-FindNeighbors(merged,dims=1:12,reduction="pca")
# merged<-FindClusters(merged,resolution = 0.2)
# merged<-RunUMAP(merged, dims=1:12, reduction="pca")
# DimPlot(merged, group.by=c("seurat_clusters","SEX","CELL_TYPE","IID"))
#Additional metadata from Puigdevall 2023 publication with failed and successful lines
hipsci_meta<-read.csv("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/hipsci-files.csv", header=T)
hipsci_kilpinen<-read.csv("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/kilpinen_table.csv", header=T)
colnames(hipsci_kilpinen)<-hipsci_kilpinen[1,]
hipsci_kilpinen<-hipsci_kilpinen[-c(1),]
hipsci_meta<-hipsci_meta[-c(2:5,10,11)]
colnames(hipsci_meta)<-c("donor_iPSC","Sample","Sex","Culture","Passage")
hipsci_meta<-full_join(hipsci_meta, hipsci_kilpinen, "donor_iPSC")
hipsci_meta <- as.data.frame(
  apply(hipsci_meta,2, function(x) gsub("\\s+", "", x)))

met<-merged@meta.data
#merged@meta.data<-left_join(merged@meta.data, hipsci_meta, "Sample")
met<-left_join(met, hipsci_meta, "Sample")
merged$dopaminergic_neurons_exp<-met$dopaminergic_neurons_exp
#rownames(merged@meta.data)<-merged@meta.data$Sample
table(merged$dopaminergic_neurons_exp)
DimPlot(merged, group.by="dopaminergic_neurons_exp",reduction="tsne")

merged$macrophage_differentiation<-met$macrophages
merged$sensory_neuron_differentiation<-met$sensory_neurons

#Read in results from pacnet (Ron's 30 samples)
# pac_res<-read.csv("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ron_data_may2025/classification_results_2025-05-15.csv")
# colnames(pac_res)<-str_split_i(colnames(pac_res),"X",2)
# colnames(pac_res)<-str_replace_all(colnames(pac_res),"[.]","-")
# colnames(pac_res)[1]<-"X"
# pac_list<-lapply(list.files("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/files_from_pacnet/", full.names = T), function(x){
#   x<-read.csv(x)
#   return(x)
# })
# #Read in pacnet results from BD2, cancer validation, and hipsci
# #Note: only 500 hipsci lines are accounted for here.
# names(pac_list)<-list.files("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/files_from_pacnet/")
# pac_list$new_ron_data<-pac_res
# pac_list2<-lapply(pac_list, function(df){
#   rownames(df)<-df$X
#   df$X<-NULL
#   df<-t(df)
#   df <- df[1:(nrow(df) - 3), ]
#   df<-data.frame(df)
#   df$sample_name<-rownames(df)
#   return(df)
# })
# 
# obj<-Reduce(rbind, pac_list2)
# obj$Sample<-obj$sample_name
# 
# #Add pacnet info to object
# merged@meta.data<-left_join(merged@meta.data, obj, "Sample")
# #merged2 <- subset(merged, subset = !is.na(Sample))
# met<-merged@meta.data
# rownames(met)<-met$Sample
# merged@meta.data<-met
# 
# #Remove cells that don't have a pacnet score
# valid_cells <- rownames(merged@meta.data[!is.na(merged@meta.data$esc), ])
# merged2 <- merged[, valid_cells]

#Outliers determined visually and confirmed by checking hipsci data
#Includes differentiated neurons from Ron
outliers<-c("ERR1243459","ERR1203456","24246R-11-03","24246R-11-04",
            "ERR2039323","ERR2039333","ERR2278299","ERR2278291")
filtered <- subset(merged2, cells = setdiff(Cells(merged2), outliers))
filtered<-subset(filtered, Sample!="24246R-11-03")
filtered<-subset(filtered, Sample!="24246R-11-04")

#Run one-class SVM on pacnet's esc scores; in progress
seurat_obj<-merged
#train_df<-filtered@meta.data[,c("esc","Group","nFeature_RNA")]
#df <- seurat_obj@meta.data[, c("esc", "Group", "nFeature_RNA")]

#install.packages("solitude")
#library(solitude)
#svm_model <- svm(formula=esc~esc, data=df, type = 'one-classification', kernel = 'radial', nu = 0.05, scale = TRUE, na.omit)
#svm_pred <- predict(svm_model, df)
#seurat_obj$svm_outlier <- !svm_pred

#Run isolation forest analysis on esc scores (this method is batch-agnostic)
library(solitude)
esc_values <- seurat_obj$esc

esc_df <- data.frame(esc = esc_values)
iso_forest <- isolationForest$new()
iso_forest$fit(esc_df)
preds <- iso_forest$predict(esc_df)

seurat_obj@meta.data$iso_outlier_esc <- preds$anomaly_score > 0.62  # threshold adjustable

DimPlot(seurat_obj, group.by="iso_outlier_esc", reduction="pca")
hist(preds$anomaly_score, breaks = 50, col = "plum", main = "Anomaly Scores", xlab = "Score")
abline(v = 0.62, col = "red", lty = 2)


df<-seurat_obj@meta.data
outlier_fraction <- df %>%
  group_by(Project) %>%
  summarise(Fraction_Outliers = mean(iso_outlier_esc))  # since iso_outlier is logical!

ggplot(outlier_fraction, aes(x = reorder(Project, -Fraction_Outliers), y = Fraction_Outliers)) +
  geom_col(fill = "orchid") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Proportion of isolated forest outliers per project",
       y = "Fraction Outliers", x = "Project")


#Run isolation forest on PCs
# seurat_obj2<-CreateSeuratObject(counts=seurat_obj@assays$RNA@layers$counts, meta.data =seurat_obj@meta.data)
# #seurat_obj2<-NormalizeData(seurat_obj2)
# seurat_obj2@assays$RNA$data<-log(seurat_obj2@assays$RNA@layers$counts+1)
# seurat_obj2<-ScaleData(seurat_obj2)
# seurat_obj2<-FindVariableFeatures(seurat_obj2)
# #seurat_obj@reductions$pca=NULL
# seurat_obj2<-RunPCA(seurat_obj2)
# rownames(seurat_obj2@meta.data)<-seurat_obj2@meta.data$Sample
pca_data <- Embeddings(seurat_obj, "pca")[, 1:6]
iso_forest <- isolationForest$new()
iso_forest$fit(as.data.frame(pca_data))
preds<-iso_forest$predict(pca_data)
seurat_obj@meta.data$iso_score_pca <- preds$anomaly_score
seurat_obj@meta.data$iso_outlier_pca <- preds$anomaly_score > 0.625  

DimPlot(seurat_obj, reduction="pca", group.by = "iso_outlier_pca") +
  ggtitle("Isolated Forest Outliers on PCA")
hist(preds$anomaly_score, breaks = 50, col = "plum", main = "Anomaly Scores", xlab = "Score")
abline(v = 0.62, col = "red", lty = 2)
DimPlot(seurat_obj, reduction="pca",group.by="iso_outlier_pca")

#combining differentiation status meta data
ipsc<-seurat_obj
#levels(ipsc$characteristic)<-c("Failed","BD2_thaw","Successful","QC")
ipsc$characteristic[ipsc$characteristic=="bad"]<-"Failed"
ipsc$characteristic[ipsc$characteristic=="neuron"]<-"Successful"

ipsc$differentiation_status<-ipsc$dopaminergic_neurons_exp

meta <- ipsc@meta.data
meta$differentiation_status<-ifelse(is.na(meta$differentiation_status), meta$characteristic, meta$differentiation_status)
#meta$differentiation_status<-as.character(meta$differentiation_status)
#meta$characteristic<-as.character(meta$characteristic)
#meta$differentiation_status[is.na(meta$differentiation_status)] <- 
# meta$characteristic[is.na(meta$differentiation_status)]

ipsc@meta.data <- meta
# Trim the meta.data to match the RNA assay
#ipsc@meta.data <- ipsc@meta.data[colnames(ipsc), ]
# ipsc<-CreateSeuratObject(counts=ipsc@assays$RNA@layers$counts, meta.data = ipsc@meta.data)
# ipsc@assays$RNA$data<-log(ipsc@assays$RNA@layers$counts+1)
# ipsc<-ScaleData(ipsc)
# rownames(ipsc)<-rownames(merged)
# colnames(ipsc)<-colnames(merged)
# ipsc<-FindVariableFeatures(ipsc)
# ipsc<-RunPCA(ipsc)

DimPlot(ipsc, group.by="differentiation_status", reduction="pca")

DimPlot(ipsc, group.by=c("differentiation_status"))

seurat_obj$differentiation_status<-ipsc$differentiation_status


########mahalanobis
#seurat_obj<-seurat_obj2


# Get the ESC metadata
#seurat_obj@reductions$pca=NULL
#seurat_obj<-RunPCA(seurat_obj)
esc_data <- seurat_obj$esc
esc_data<-data.frame(esc_data)
esc_data$Sample<-rownames(esc_data)
# Choose how many PCs to use
n_pcs <- 6

# Extract PCA embeddings
pca_data <- Embeddings(seurat_obj, "pca")[, 1:n_pcs]
pca_data<-data.frame(pca_data)
#pca_data$Sample<-rownames(pca_data)

# Bind into one data frame
#svm_input <- data.frame(pca_data, esc = esc_data)
#svm_input<-full_join(pca_data, esc_data, "Sample")

#svm_input$Sample=NULL

#intersect(names(esc_data), names(pca_data))
#svm_input<-pca_data
#svm_input<-esc_data
#Mahal-anobis distance
#esc_data$Sample=NULL
center <- colMeans(pca_data)
cov_matrix <- cov(pca_data)
md <- mahalanobis(pca_data, center, cov_matrix)
threshold <- qchisq(0.99, df = ncol(pca_data))  # 99% cutoff
outliers_md <- md > threshold
outliers_md<-data.frame(outliers_md)
outliers_md$Sample<-rownames(outliers_md)

outliers_md$Sample <- sub("^X", "", outliers_md$Sample)
outliers_md$Sample <- gsub("\\.", "-", outliers_md$Sample)       # Replace all periods with dashes

rownames(outliers_md)<-outliers_md$Sample


#seurat_obj@meta.data<-left_join(seurat_obj@meta.data, outliers_md, "Sample")

meta<-seurat_obj@meta.data
meta$outliers_md_pca=NULL
meta<-left_join(meta, outliers_md, "Sample")

seurat_obj@meta.data$outliers_md_pca<-meta$outliers_md

#rownames(seurat_obj@meta.data)<-seurat_obj@meta.data$Sample

DimPlot(seurat_obj, group.by="outliers_md_pca")
DimPlot(seurat_obj, group.by="outliers_md_pca",reduction="pca")

plot(density(md, bw=0.5))
qqplot(qchisq(ppoints(100), df = 3), md, main = expression("Q-Q plot of Mahalanobis" * ~D^2 *
                                                             " vs. quantiles of" * ~ chi[3]^2))
abline(0, 1, col = 'gray')


########mahalanobis on esc
center <- colMeans(esc_data$esc_data)
cov_matrix <- cov(esc_data$esc_data)
md <- mahalanobis(esc_data, center, cov_matrix)
threshold <- qchisq(0.99, df = ncol(esc_data))  # 99% cutoff
outliers_md <- md > threshold
outliers_md<-data.frame(outliers_md)
outliers_md$Sample<-rownames(outliers_md)

outliers_md$Sample <- sub("^X", "", outliers_md$Sample)
outliers_md$Sample <- gsub("\\.", "-", outliers_md$Sample)       # Replace all periods with dashes

rownames(outliers_md)<-outliers_md$Sample


#seurat_obj@meta.data<-left_join(seurat_obj@meta.data, outliers_md, "Sample")

meta<-seurat_obj@meta.data
meta$outliers_md_pca=NULL
meta<-left_join(meta, outliers_md, "Sample")

seurat_obj@meta.data$outliers_md_pca<-meta$outliers_md

#rownames(seurat_obj@meta.data)<-seurat_obj@meta.data$Sample

DimPlot(seurat_obj, group.by="outliers_md_pca")
DimPlot(seurat_obj, group.by="outliers_md_pca",reduction="pca")

plot(density(md, bw=0.5))
qqplot(qchisq(ppoints(100), df = 3), md, main = expression("Q-Q plot of Mahalanobis" * ~D^2 *
                                                             " vs. quantiles of" * ~ chi[3]^2))
abline(0, 1, col = 'gray')

########


valid_cells <- rownames(seurat_obj@meta.data[!is.na(seurat_obj@meta.data$esc), ])
seurat_obj4 <- seurat_obj[, valid_cells]

#Z scores on esc scores
z_scores <- scale(seurat_obj$esc)
outliers_z <- abs(z_scores) > 3  # 3 SDs from the mean

#Z scores on PCs
pc1 <- pca_data[, 1]
z_pc1 <- scale(pc1)
outliers_pc1 <- abs(z_pc1) > 3


#Local outlier factor
library(Rlof)
lof_scores <- lof(svm_input, k = 20)
outliers_lof <- lof_scores > 1.5  # tweakable threshold


#DBSCAN
library(dbscan)
svm_input<-pca_data
kNNdistplot(svm_input, k = 10)  # pcs: your 10-PC matrix
abline(h = 40, col = "red")

db <- dbscan(data.frame(na.omit(seurat_obj@meta.data$nFeature_RNA)), eps = 40, minPts = 10)
outliers_dbscan <- db$cluster
outliers_dbscan<-data.frame(outliers_dbscan)
outliers_dbscan$Sample<-rownames(pca_data)
plot(outliers_dbscan)
seurat_obj2@meta.data<-left_join(seurat_obj2@meta.data, outliers_dbscan, "Sample")
rownames(seurat_obj@meta.data)<-seurat_obj@meta.data$Sample
DimPlot(seurat_obj2, group.by="outliers_dbscan")
# gene_loadings <- data.frame(Loadings(seurat_obj[["pca"]]))
# gene_loadings$gene<-rownames(gene_loadings)

# Extract PCA gene loadings
gene_loadings <- Loadings(seurat_obj[["pca"]])

# How many top genes to fetch?
top_n <- 50

# Get names of principal components (columns)
pcs <- colnames(gene_loadings)

# Initialize a list to hold top genes per PC
top_genes_list <- lapply(pcs, function(pc) {
  # Order genes by absolute loading for this PC
  top_genes <- names(sort(abs(gene_loadings[, pc]), decreasing = TRUE))[1:top_n]
  return(top_genes)
})

# Name the list by PC
names(top_genes_list) <- pcs

# Bind into a tidy dataframe (pad with NA if unequal, but here they'll all be equal)
top_genes_df <- as.data.frame(top_genes_list)


#KNN Distance Outliers
library(FNN)
k <- 10
knn <- get.knn(svm_input, k = k)
distances <- rowMeans(knn$nn.dist)
threshold <- quantile(distances, 0.95)
outliers_knn <- distances > threshold


#Consensus: find cells which are classified as outliers across multiple methods
consensus_outlier <- outliers_svm | outliers_iso | outliers_md
true_eldritch_cells <- outliers_svm & outliers_iso & outliers_md





#Assess PCA gene loadings

# gene_loadings <- data.frame(Loadings(seurat_obj[["pca"]]))
# gene_loadings$gene<-rownames(gene_loadings)

# Extract PCA gene loadings
gene_loadings <- Loadings(merged[["pca"]])

top_n <- 50
pcs <- colnames(gene_loadings)

top_genes_list <- lapply(pcs, function(pc) {
  top_genes <- names(sort(abs(gene_loadings[, pc]), decreasing = TRUE))[1:top_n]
  return(top_genes)
})

names(top_genes_list) <- pcs

top_genes_df <- as.data.frame(top_genes_list)



#cell cycle and more pca plotting
library(Seurat)
#cell cycle score
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

ipsc <- CellCycleScoring(ipsc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Visualize the distribution of cell cycle markers across
RidgePlot(ipsc, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
ipsc <- RunPCA(ipsc, features = c(s.genes, g2m.genes))
DimPlot(ipsc)

DimPlot(ipsc) + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))


#PCA plots
DimHeatmap(ipsc, dims = c(1, 2, 3))

stdev <- ipsc@reductions$pca@stdev

# Compute variance explained
var_explained <- (stdev^2) / sum(stdev^2)

# Create a data frame
pc_data <- data.frame(PC = paste0("PC", 1:length(var_explained)),
                      VarianceExplained = var_explained)

# Keep only the first 20 PCs
pc_data <- pc_data[1:20, ]

# Ensure correct PC order
pc_data$PC <- factor(pc_data$PC, levels = paste0("PC", 1:20))

# Load plotting magic
library(ggplot2)

# Bar plot with numbers on top
ggplot(pc_data, aes(x = PC, y = VarianceExplained)) +
  geom_bar(stat = "identity", fill = "#ADD8E6", color = "black") +  # Baby blue bars
  geom_text(aes(label = round(VarianceExplained, 3)), 
            vjust = -0.3, size = 3.5, color = "black") +           # Baby pink numbers
  theme_minimal() +
  ylab("Proportion of Variance Explained") +
  xlab("Principal Component") +
  ggtitle("Variance Captured by PCs 1–20") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14))

