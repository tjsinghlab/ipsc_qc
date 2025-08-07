library(Biobase)
library(BisqueRNA)
library(Seurat)
library(Matrix)
library(stringr)
library(ggplot2)
library(reshape2)
library(dplyr)
#remotes::install_github("jlaffy/scalop")
library(scalop)

#bulk ipsc object:
bulk_ipsc <- readRDS("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/merged_ipsc_seurat_obj_with_prelim_outlier_res.rds")

#reference scRNA data:
  # setwd("/gpfs/commons/groups/singh_lab/projects/bd2village/dl_datasets/cuomo_2020/scrna/")
  # unzip("raw_counts.csv.zip", exdir = "unzipped_counts")
  # raw_counts <- read.csv("unzipped_counts/raw_counts.csv", header = TRUE, row.names = 1)
  # meta<-read.csv("cell_metadata_cols.csv")
  # cuomo_obj<-CreateSeuratObject(counts=raw_counts, meta.data = meta)
  # saveRDS(cuomo_obj, "cuomo_2020_seurat_object_raw_counts_and_metadata.rds")

ref<-readRDS("/gpfs/commons/groups/singh_lab/projects/bd2village/dl_datasets/cuomo_2020/scrna/cuomo_2020_seurat_object_raw_counts_and_metadata.rds")

#Bisque requires expression data in the ExpressionSet format from the Biobase package:
#convert bulk data from a matrix (columns are samples, rows are genes) to an ExpressionSet:
bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(GetAssayData(bulk_ipsc, assay = "RNA", layer="data")))

sc.counts.matrix<-as.matrix(GetAssayData(ref, assay = "RNA", layer="counts"))
rownames(sc.counts.matrix)<-str_split_i(rownames(sc.counts.matrix), "-", 2)
genes <- str_split_i(rownames(sc.counts.matrix), "-", 1)
sc.counts.matrix <- rowsum(sc.counts.matrix, group=genes)
sample.ids <- colnames(sc.counts.matrix)

sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=sample.ids,
                       SubjectName=ref$donor,
                       cellType=ref$day)
sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc.pdata <- new("AnnotatedDataFrame",
                data=sc.pheno,
                varMetadata=sc.meta)
sc.eset <- Biobase::ExpressionSet(assayData=as.matrix(sc.counts.matrix),
                                  phenoData=sc.pdata)

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, use.overlap=FALSE, verbose=TRUE)

#A list is returned with decomposition estimates in slot bulk.props.
ref.based.estimates <- res$bulk.props
knitr::kable(ref.based.estimates, digits=2)

est<-data.frame(t(ref.based.estimates))
est$Sample<-rownames(est)

#add results to object
meta<-left_join(bulk_ipsc@meta.data, est, "Sample")
rownames(meta)<-meta$Sample

bulk_ipsc@meta.data<-meta

#Plotting:
# Extract metadata
df <- bulk_ipsc@meta.data
df$Sample <- rownames(df)

# Melt into long format
df_melt <- melt(df[, c("Sample", "day0", "day1", "day2", "day3")], id.vars = "Sample")

# Plot
ggplot(df_melt, aes(x = Sample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  #scale_fill_manual(values = c("Success" = "#87CEFA", "Fail" = "#FFB6C1")) +
  theme_minimal() +
  labs(title = "Timepoint Gene Signatures Proportions per Sample",
       x = "Sample", y = "Proportion", fill = "Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 1))

#see "process_puigdevall_kilpinen_data.R" for how below object was created:
puigdevall<-readRDS("bulk_ipsc_dataset_with_puigdevall_d11_differentiationstatus_deconvolution_results.rds")

#add deconvolution results from puigdevall classifications to our object, so we'll now have deconvolution results from both the puigdevall and cuomo publications:
bulk_ipsc$Successful_Puigdevall<-puigdevall$Successful
bulk_ipsc$Failed_Puigdevall<-puigdevall$Failed

saveRDS(bulk_ipsc, "bulk_ipsc_object_with_cuomo_puigdevall_deconvolution_results.rds")


bulk_ipsc<-readRDS("bulk_ipsc_object_with_cuomo_puigdevall_deconvolution_results.rds")

bulk_ipsc$differentiation_status %>% table

#Find DEGs between failed and successful hipsci lines (as defined by puigdevall):
hipsci_classified <- subset(bulk_ipsc, 
                            subset = differentiation_status %in% c("Successful", "Failed"))

Idents(hipsci_classified)<-hipsci_classified$differentiation_status
markers<-FindAllMarkers(hipsci_classified)

#Score all cells (samples) by expression of "successful" geneset
successful_geneset<-subset(markers, 
                           cluster=="Successful" & 
                             avg_log2FC>=1 & 
                             p_val_adj<=0.01)$gene
failed_geneset<-subset(markers, 
                       cluster=="Failed" & 
                         avg_log2FC>=1 & 
                         p_val_adj<=0.01)$gene

#this is just the output of running NormalizeData and FindAllMarkers on the cuomo scRNA dataset:
cuomo_timepoint_markers<-readRDS("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/cuomo_day_findallmarkers.rds")

cuomo_timepoint_markers<-subset(cuomo_timepoint_markers, 
                         avg_log2FC>=1 & 
                         p_val_adj<=0.01)

d0_markers<-subset(cuomo_timepoint_markers, cluster=="day0")$gene
d1_markers<-subset(cuomo_timepoint_markers, cluster=="day1")$gene
d2_markers<-subset(cuomo_timepoint_markers, cluster=="day2")$gene
d3_markers<-subset(cuomo_timepoint_markers, cluster=="day3")$gene

bulk_ipsc$successful_sigscores<-scalop::sigScores(m = as.matrix(GetAssayData(bulk_ipsc, assay = "RNA", layer="counts")), sigs=successful_geneset)

bulk_ipsc$failed_sigscores<-scalop::sigScores(m = as.matrix(GetAssayData(bulk_ipsc, assay = "RNA", layer="counts")), sigs=failed_geneset)

bulk_ipsc$cuomo_d0_sigscores<-scalop::sigScores(m = as.matrix(GetAssayData(bulk_ipsc, assay = "RNA", layer="counts")), sigs=d0_markers)
bulk_ipsc$cuomo_d1_sigscores<-scalop::sigScores(m = as.matrix(GetAssayData(bulk_ipsc, assay = "RNA", layer="counts")), sigs=d1_markers)
bulk_ipsc$cuomo_d2_sigscores<-scalop::sigScores(m = as.matrix(GetAssayData(bulk_ipsc, assay = "RNA", layer="counts")), sigs=d2_markers)
bulk_ipsc$cuomo_d3_sigscores<-scalop::sigScores(m = as.matrix(GetAssayData(bulk_ipsc, assay = "RNA", layer="counts")), sigs=d3_markers)

status_success<-data.frame(
  bulk_ipsc$differentiation_status, 
  bulk_ipsc$successful_sigscores)

status_failure<-data.frame(
  bulk_ipsc$differentiation_status, 
  bulk_ipsc$failed_sigscores)


#Plotting:
# Summarize your data: mean and standard error
colnames(status_success)<-c("Status","Score")
summary_df <- status_success %>%
  group_by(Status) %>%
  summarise(
    mean_score = mean(Score, na.rm = TRUE),
    se_score = sd(Score, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Now plot the magic
ggplot(summary_df, aes(x = Status, y = mean_score, fill = Status)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
                width = 0.2, linewidth = 0.8) +
  labs(x = "Status", y = "Mean Score", title = "Mean Success Score by Group") +
  theme_minimal() +
  theme(legend.position = "none")


# day0_failed_scores<-data.frame(bulk_ipsc$day0, bulk_ipsc$failed_sigscores)
# colnames(day0_failed_scores)<-c("day0","score")

day0_failed_scores<-data.frame(bulk_ipsc$day3, bulk_ipsc$cuomo_d3_sigscores)
colnames(day0_failed_scores)<-c("day0","score")


ggplot(day0_failed_scores, aes(x = day0, y = score)) +
  geom_point(size = 2, alpha = 0.7) +  # baby pink points
  geom_smooth(method = "lm", se = TRUE) +  # baby blue regression
  labs(x = "Day 3 bisque", y = "Day 3 sigScore", title = "bisque vs sigscore: d3") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12)
  )


#Use failed/successful hipsci classifications as validation set for deconvolution
table(hipsci_classified$differentiation_status, hipsci_classified$Successful_Puigdevall)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 35) %>%
  ungroup() -> top10
DoHeatmap(hipsci_classified, features = top10$gene) 

#bMIND for cell type expression
devtools::install_github('randel/MIND')
library(MIND)
bulk_expr <- GetAssayData(bulk_ipsc, assay = "RNA", slot = "data") # log-normalized
bulk_expr <- as.matrix(bulk_expr)  # Ensure it's a matrix

# Extract the columns from metadata
cell_frac <- bulk_ipsc@meta.data[, c("Successful_Puigdevall", "Failed_Puigdevall")]

# Ensure it's a matrix
cell_frac <- as.matrix(cell_frac)

# Set sample names as rownames
# rownames(cell_frac) <- rownames(seurat_object@meta.data)

# Transpose to match bMIND's needs: cell types (rows) × samples (columns)
cell_frac <- t(cell_frac)



bMIND_res <- bMIND(
  bulk = bulk_expr,
  frac = cell_frac,
  sample_id = colnames(bulk_expr)
)



expr_list <- bMIND_res

seurat_ct_list <- lapply(names(expr_list), function(ct) {
  mat <- expr_list[[ct]]
  CreateSeuratObject(counts = mat, project = paste0("bMIND_", ct))
})
names(seurat_ct_list) <- names(expr_list)


