library(dplyr)
library(biomaRt)
library(tidyverse)
library(data.table)
library(stringr)
library(patchwork)
library(ggplot2)
library(Seurat)

base_dir<-"/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/"
setwd(base_dir)
read.file<-function(x){
  x<-fread(x)
  #x<-read.csv(x, header=T, sep=" ")
  return(x)
}
bd2<-lapply(list.files("bd2_village/wdl_pipeline/outs/RSEM_outputs/rsem_genes_mats", full.names = T), read.file)
names(bd2)<-str_split_i(list.files("bd2_village/wdl_pipeline/outs/RSEM_outputs/rsem_genes_mats"),"[.]",1)

hipsci<-lapply(list.files("hipsci/wdl_pipeline/outs/RSEM_outputs/rsem_genes_mats", full.names = T), read.file)
names(hipsci)<-str_split_i(list.files("hipsci/wdl_pipeline/outs/RSEM_outputs/rsem_genes_mats"),"[.]",1)
hipsci <- hipsci[!duplicated(names(hipsci))]

val<-lapply(list.files("cancer_validation/rsem_mats", full.names = T),read.file)
names(val)<-str_split_i(list.files("cancer_validation/rsem_mats"),"[.]",1)

hart1<-lapply(list.files("ron_data_may2025/RSEM_outputs/rsem_genes_mats", full.names=T),read.file)
names(hart1)<-paste0("S",str_split_i(list.files("ron_data_may2025/RSEM_outputs/rsem_genes_mats"),"[.]",1))
hart2<-lapply(list.files("ron_data_may2025/set2/RSEM_outputs/rsem_genes_mats", full.names=T),read.file)
names(hart2)<-paste0("S", str_split_i(list.files("ron_data_may2025/set2/RSEM_outputs/rsem_genes_mats"),"[.]",1))

may_rna<-lapply(list.files("ron_data_may2025_2/RSEM_outputs/rsem_genes_mats", full.names=T),read.file)
names(may_rna)<-paste0("S",str_split_i(list.files("ron_data_may2025_2/RSEM_outputs/rsem_genes_mats"),"[.]",1))

ipsc_mats<-c(bd2, hipsci, val, hart1, hart2, may_rna)
ipsc_mats<-may_rna
genes<-lapply(ipsc_mats, function(x){
  x<-x$gene_id
  return(x)
})

mart<-useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
meta.genes<-lapply(genes, function(x){
  y<-getBM(attributes = c("ensembl_gene_id","external_gene_name", "description"), 
           values = x, mart = mart)
  return(y)
}) 

new_mats<-Map(function(og_map, new_genes){
  og_map<-og_map[,c(1,5)]
  colnames(og_map)<-c("ensembl_gene_id","count")
  og_map$ensembl_gene_id<-str_split_i(og_map$ensembl_gene_id, "[.]", 1)
  out<-full_join(og_map, new_genes, by="ensembl_gene_id")
  return(out)
},ipsc_mats, meta.genes)

new_mats2<-lapply(new_mats, function(df){
  df<-df[,-c(1,4)]
  df[df==''] <- NA
  df<-df[complete.cases(df), ]
  df<-df[!grepl("^[0-9]", df$external_gene_name) &
           !grepl("Metazoa", df$external_gene_name) &
           !grepl("^RNU", df$external_gene_name) &
           !grepl("^SNOR", df$external_gene_name) &
           df$external_gene_name != "Y_RNA" &
           df$external_gene_name != "Vault" &
           !grepl("^U[0-9]+$", df$external_gene_name), ] 
  return(df)
})

new_mats2<-lapply(new_mats2, function(df){
  df[!duplicated(df$external_gene_name) &
       !duplicated(df$external_gene_name, fromLast = TRUE), ]
  return(df)
})
new_mats2<-Map(function(df, namex){
  colnames(df)<-c(namex, "gene")
  return(df)
},new_mats2, names(new_mats2))

saveRDS(new_mats2, "hart_and_other_ipsc_matrices_lists.rds")

#If any sample starts with a number, this will cause issues downstream.
#Make sure any name that starts with a number is preceded by "S"
names(new_mats2) <- sapply(names(new_mats2), function(name) {
  if (grepl("^[0-9]", name)) {
    paste0("S", name)
  } else {
    name
  }
})

#There are duplicate gene names in the matrices? Where one gene has a real value, and the duplicate has expression value of 0
#So we should remove the ones that have a gene value of 0 so we are left with only the real values
process_df <- function(df) {
  # Identify duplicate genes
  duplicates <- duplicated(df$gene) | duplicated(df$gene, fromLast = TRUE)
  # Keep rows where the "value" column (or the first column) is not zero
  df <- df[!(duplicates & df[[1]] == 0), ]
  return(df)
}

#For entries where two different Ensembl Gene IDs mapped onto the same gene, we'll take the average of the two values
aggene<-function(df){
  # Convert to data.table for efficiency
  dt <- as.data.table(df)
  # Calculate the mean for duplicated genes (accessing the first column dynamically)
  dt <- dt[, .(value = mean(.SD[[1]])), by = gene]
  return(dt)
}
ipsc_matrices_lists<-readRDS("hart_and_other_ipsc_matrices_lists.rds")

ipsc2<-lapply(ipsc_matrices_lists, data.frame)
ipsc2<-lapply(ipsc2, process_df)
ipsc2<-lapply(ipsc2, aggene)
ipsc2<-Map(function(df, namex){
  colnames(df)<-c("gene",namex)
  return(df)
},ipsc2, names(ipsc2))
group_batches <- function(prefix, df_list) {
  # Filter data frames by prefix
  group <- df_list[grepl(paste0("^", prefix), names(df_list))]
}
# hart_data <- lapply(ipsc2, function(x){
#   x<-data.frame(x)
#   return(x)
# })

#hart_data2<-Reduce(function(x, y) merge(x, y, all = TRUE), hart_data)
names(ipsc2)[1:24]<-paste0("BD2_", names(ipsc2)[1:24])
# Join data frames by their prefixes
bd2 <- group_batches("BD2", ipsc2)
hipsci <- group_batches("ERR", ipsc2)
val <- group_batches("ipsc", ipsc2)
hart1 <- group_batches("S24246R-09", ipsc2)
hart2 <- group_batches("S24246R-11", ipsc2)
hart3 <- group_batches("S24246R-15", ipsc2)

library(purrr)
bd2<-bd2 %>% purrr::reduce(full_join, by="gene")
val<-val%>% purrr::reduce(full_join, by="gene")
hipsci<-hipsci%>% purrr::reduce(full_join, by="gene")
hart1<-hart1%>% purrr::reduce(full_join, by="gene")
hart2<-hart2%>% purrr::reduce(full_join, by="gene")
hart3<-hart3 %>% purrr::reduce(full_join, by="gene")
ipsc_mats<-list(bd2, val, hipsci, hart1, hart2)

ipsc_mats2<-lapply(ipsc_mats, function(x){
  x<-data.frame(x)
  rownames(x)<-x$gene
  x$gene=NULL
  return(x)
})
# hart_data2<-data.frame(hart_data2)
# rownames(hart_data2)<-hart_data2$gene
# hart_data2$gene<-NULL

names(ipsc_mats2)<-c("bd2","val","hipsci","hart1","hart2")

#saveRDS(ipsc_mats2, "merged_ipsc_mats_by_batch.rds")
# saveRDS(hart_data2, "merged_hart_ipsc_30_samples.rds")
saveRDS(ipsc_mats2, "merged_ipsc_mats_by_batch_plus_30hart_samples.rds")

library(Seurat)
ipsc_mats3<-lapply(ipsc_mats2, function(x){
  x[is.na(x)]<-0
  return(x)
})
#hart_data2[is.na(hart_data2)]<-0
objects<-lapply(ipsc_mats3, function(mat){
  obj<-CreateSeuratObject(counts=mat, assay="RNA")
  obj$Sample<-colnames(obj)
  obj<-PercentageFeatureSet(obj,pattern="^MT-",col.name = "percent.mito")
  return(obj)
})
#names(objects)<-c("bd2","val","hipsci","hart1","hart2","hart3")

#obj<-CreateSeuratObject(counts=hart_data2, assay="RNA")
#obj$Sample<-colnames(obj)
#obj<-PercentageFeatureSet(obj, pattern="MT-", col.name="percent.mito")

# saveRDS(objects, "ipsc_objects_00.rds")
# saveRDS(obj, "hart_30samples_object.rds")
saveRDS(objects, "ipsc_objects_plus30hartsamples_00.rds")

#obj<-readRDS("hart_30samples_object.rds")
samples_key <- read.csv("/gpfs/commons/groups/singh_lab/projects/bd2village/data/hart_data_april2025/samples_key.csv")
samples_key$sample_name<-str_split_i(samples_key$R1_Filename, "_", 1)
samples_key$Sample<-paste0("S",samples_key$Sample)

#objects$hart1@meta.data$sample_name<-str_split_i(objects$hart1@meta.data$Sample, "X",2)
objects$hart1@meta.data$Sample<-str_replace_all(objects$hart1@meta.data$Sample, "[.]", "-")
objects$hart1@meta.data<-left_join(objects$hart1@meta.data, samples_key, "Sample")

#objects$hart2@meta.data$sample_name<-str_split_i(objects$hart2@meta.data$Sample, "X",2)
objects$hart2@meta.data$Sample<-str_replace_all(objects$hart2@meta.data$Sample, "[.]", "-")
objects$hart2@meta.data<-left_join(objects$hart2@meta.data, samples_key, "Sample")

may_meta<-read.csv("/gpfs/commons/groups/singh_lab/projects/bd2village/data/hart_data_may2025/24246-15-final-QC-summary.csv")
may_meta$Sample<-paste0("S", may_meta$Sample)
meta3<-objects$hart3@meta.data
meta3$Sample<-str_replace_all(meta3$Sample, "[.]", "-")
meta3<-left_join(meta3, may_meta, "Sample")
objects$hart3@meta.data<-meta3

hipsci_meta<-read.csv(paste0(base_dir,"hipsci/hipsci-cell-lines.csv"),header=T)
hipsci_files_RNA <- read.csv("hipsci/hipsci-files-RNA.csv")
colnames(hipsci_files_RNA)[6]<-"Sample"
colnames(hipsci_files_RNA)[1]<-"Name"
hipsci_meta<-full_join(hipsci_meta, hipsci_files_RNA, "Name")

hipsci_obj_meta<-objects$hipsci@meta.data
hipsci_obj_meta<-left_join(hipsci_obj_meta, hipsci_meta, "Sample")
colnames(hipsci_obj_meta)[9]<-"Sex"
colnames(hipsci_obj_meta)[15]<-"Culture"
hipsci_obj_meta$Culture.y=NULL
hipsci_obj_meta$Sex.y=NULL
hipsci_obj_meta$Disease.Status<-"Control"
rownames(hipsci_obj_meta)<-rownames(objects$hipsci@meta.data)
objects$hipsci@meta.data<-hipsci_obj_meta
objects$hipsci$Group<-"HIPSCI_BATCH1"

bd2_meta<-read.csv(paste0(base_dir,"bd2_village/bd2_dec2024_meta.csv"),header=T)
colnames(bd2_meta)[3]<-"Name"
colnames(bd2_meta)[1]<-"Sample"
colnames(bd2_meta)[14]<-"Sex"
colnames(bd2_meta)[15]<-"Disease.Status"
bd2_meta$Sex<-ifelse(bd2_meta$Sex=="TRUE", "Female", "Male")
bd2_meta$Group<-"BD2_BATCH1"
bd2_meta<-left_join(objects$bd2@meta.data, bd2_meta, "Sample")
rownames(bd2_meta)<-rownames(objects$bd2@meta.data)
objects$bd2@meta.data<-bd2_meta

objects$val$Group<-"CANCER_VALIDATION"

saveRDS(objects, "ipsc_batch1_and_54hartsamples_objects_with_metadata.rds")

objects<-readRDS("ipsc_batch1_and_54hartsamples_objects_with_metadata.rds")
merged<-merge(objects$bd2, objects$hipsci)
merged<-JoinLayers(merged)
rep(gc(),100)
merged<-merge(merged, objects$val)
merged<-JoinLayers(merged)
hart30 <- readRDS("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/hart_30_ipsc_merged_tsne_no_sct_no_harmony.rds")
merged<-merge(merged, hart30)
merged<-JoinLayers(merged)
rep(gc(),100)
saveRDS(merged, "ipsc_batches_plus_30hart_merged_seurat_raw.rds")
#object<-CreateSeuratObject(counts=merged@assays$RNA@layers$counts, meta.data = merged@meta.data)
merged2<-readRDS("ipsc_batches_plus_30hart_merged_seurat_raw.rds")

merged<-CreateSeuratObject(counts=merged@assays$RNA, meta.data = merged@meta.data)

#merged<-SCTransform(merged, vst.flavor="v2", vars.to.regress=c("Sample","percent.mito","Group"), conserve.memory = T)
DefaultAssay(merged)<-"RNA"
merged@assays$RNA$data<-log(merged@assays$RNA@layers$counts+1)
merged<-FindVariableFeatures(merged)
merged<-ScaleData(merged)
DefaultAssay(merged)<-"RNA"
merged<-RunPCA(merged)
ElbowPlot(merged)
merged<-FindNeighbors(merged,dims=1:12)
merged<-FindClusters(merged,resolution = 0.2)
merged<-RunTSNE(merged,dims = 1:12)
#merged$Group[is.na(merged$Group)]<-"Hart30"
merged$Project[is.na(merged$Project)]<-merged$Group
merged$Sample[is.na(merged$Sample)] <- merged$sample_name[is.na(merged$Sample)]
TSNEPlot(merged,group.by=c("seurat_clusters","Group","Disease.Status","characteristic"))
DimPlot(merged, reduction="pca",group.by=c("seurat_clusters","Group","Disease.Status","characteristic"))
#TSNEPlot(merged,group.by=c("characteristic"))
#saveRDS(merged, "ipsc_merged_tsne_no_sct_no_harmony_increasedresolution0point2.rds")
#saveRDS(merged, "hart_30_ipsc_merged_tsne_no_sct_no_harmony.rds")
# saveRDS(merged, "ipsc_and_hart30_merged_tsne_pca.rds")
# merged<-readRDS("ipsc_and_hart30_merged_tsne_pca.rds")

saveRDS(merged, "ipsc_object_merged_pca_tsne_norm.rds")

merged<-readRDS("ipsc_and_hart30_merged_tsne_pca.rds")
outliers<-c("ERR1243459","ERR1203456","24246R-11-03","24246R-11-04",
            "ERR2039323","ERR2039333","ERR2278299","ERR2278291")
filtered <- subset(merged, cells = setdiff(Cells(merged), outliers))
filtered<-subset(filtered, Sample!="24246R-11-03")
filtered<-subset(filtered, Sample!="24246R-11-04")

DimPlot(filtered)

your_seurat_object<-merged

# Extract PCA embeddings and attach cell names
pca_coords <- Embeddings(your_seurat_object, "pca") %>%
  as.data.frame() %>%
  mutate(cell = rownames(.))

# Attach sample names from meta.data
metadata <- your_seurat_object@meta.data %>%
  mutate(cell = rownames(.)) %>%
  dplyr::select(cell, Sample)

# Merge coords with metadata
pca_coords <- left_join(pca_coords, metadata, by = "cell")

# Compute Mahalanobis distance to identify outliers
center <- colMeans(pca_coords[, c("PC_1", "PC_2")])
cov_matrix <- cov(pca_coords[, c("PC_1", "PC_2")])
pca_coords$mahal_dist <- mahalanobis(pca_coords[, c("PC_1", "PC_2")], center, cov_matrix)

# Flag the top 1% as outliers
threshold <- quantile(pca_coords$mahal_dist, 0.99)
pca_coords$outlier <- pca_coords$mahal_dist > threshold
library(ggrepel)
# Plot with outlier labels
DimPlot(your_seurat_object, reduction = "pca") +
  geom_text_repel(data = pca_coords %>% filter(outlier),
                  aes(x = PC_1, y = PC_2, label = Sample),
                  color = "red", size = 3)

markers<-FindAllMarkers(merged)

expression_matrix <- as.data.frame(GetAssayData(merged, slot = "data"))

write.csv(expression_matrix, file = "hart_30_ipsc_expr_mat.csv", row.names = TRUE)

objs<-SplitObject(merged, "Group")

mat<-expression_matrix
colnames(mat)<-str_split_i(colnames(mat), "X", 2)
colnames(mat)<-str_replace_all(colnames(mat), "[.]", "-")
meta<-data.frame(merged$sample_name, merged$characteristic)
colnames(meta)<-c("sample_name","description1")
write.csv(mat, "hart_30_mat.csv")
write.csv(meta, "hart_30_meta.csv")

#Export matrices for PACnet:
bd2_mat<-as.data.frame(GetAssayData(objs$BD2_BATCH1, slot="data"))
bd2_meta<-data.frame(objs$BD2_BATCH1$Sample, objs$BD2_BATCH1$Group)
colnames(bd2_meta)<-c("sample_name","description1")
write.csv(bd2_mat, "bd2_matrix_from_split.csv")
write.csv(bd2_meta, "bd2_meta_from_split.csv")

val_mat<-as.data.frame(GetAssayData(objs$CANCER_VALIDATION, slot="data"))
val_meta<-data.frame(objs$CANCER_VALIDATION$Sample, objs$CANCER_VALIDATION$Group)
colnames(val_meta)<-c("sample_name","description1")
write.csv(val_mat, "val_matrix_from_split.csv")
write.csv(val_meta, "val_meta_from_split.csv")

objs[["HIPSCI_BATCH1"]]->hipsci
cells <- colnames(hipsci)  # Get all cell names
n <- length(cells)
group_size <- 100

cell_groups <- split(cells, ceiling(seq_along(cells) / group_size))

seurat_chunks <- lapply(cell_groups, function(cell_set) {
  subset(hipsci, cells = cell_set)
})

hipsci_mat<-as.data.frame(GetAssayData(seurat_chunks$`7`, layer="data"))
hipsci_meta<-data.frame(seurat_chunks$`7`$Sample, seurat_chunks$`7`$Group)
colnames(hipsci_meta)<-c("sample_name","description1")
write.csv(hipsci_mat, "hipsci_matrix_from_split_7.csv")
write.csv(hipsci_meta, "hipsci_meta_from_split_7.csv")

#Split hipsci b/c it's too big for pacnet
sample_names <- colnames(hipsci_mat)
split_indices <- split(sample_names, ceiling(seq_along(sample_names) / 100))

split_matrices <- lapply(split_indices, function(samples) hipsci_mat[, samples, drop = FALSE])
split_metadata <- lapply(split_indices, function(samples) hipsci_meta[samples, , drop = FALSE])

DotPlot(merged, features=c(
  #c0 pluripotent, male
  "TMSB4Y","TXLNGY","TTTY14","NLGN4Y","CD24","CCR6","NGB",
  #c1 closer to differentiation to endothelial/circulation cells
  "HOXC8","HOXD8","DCN","FGF7","DPT","NR2E1","IL7R","CCL7","TGFBI",
  #c2 cell death? closer to differentiating into sensory neurons? dna damage and repair? nonsense genes and cell death regulators
  "AGTR2","OR5BE1P","RMRP","H2AC21","ELOAP1",
  #c3 proliferation
  "PRAC1","PRAC2","HOXB13","CALM1P1","FTH1P6","TBX1","SOX14"))+coord_flip()
merged$prediction<-"NA"

merged$prediction[merged$seurat_clusters==0]<-"Pluripotent_Ychr"
merged$prediction[merged$seurat_clusters==1]<-"Slightly_less_pluripotent"
merged$prediction[merged$seurat_clusters==2]<-"cell_death_dna_damage"
merged$prediction[merged$seurat_clusters==3]<-"Proliferative"
TSNEPlot(merged, group.by="prediction",label=T, label.box=T)

ggplot(merged@meta.data, aes(x = as.factor(seurat_clusters), fill = Sex)) +
  geom_bar(position = "dodge") +  # "dodge" places bars side by side
  labs(x = "Seurat Clusters", y = "Count", title = "Sex Distribution in Seurat Clusters") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(stringr)
merged$hipsci_donor<-"NA"
merged$hipsci_donor<-ifelse(merged$Group=="HIPSCI_BATCH1", str_split_i(merged$Name, "-", 2), "NA")
merged$hipsci_donor<-str_split_i(merged$hipsci_donor, "_", 1)
