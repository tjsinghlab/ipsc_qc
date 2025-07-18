pac_res<-read.csv("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ron_data_may2025/classification_results_2025-05-15.csv")
colnames(pac_res)<-str_split_i(colnames(pac_res),"X",2)
colnames(pac_res)<-str_replace_all(colnames(pac_res),"[.]","-")
colnames(pac_res)[1]<-"X"
pac_list<-lapply(list.files("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/files_from_pacnet/", full.names = T), function(x){
  x<-read.csv(x)
  return(x)
})
names(pac_list)<-list.files("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/files_from_pacnet/")
pac_list$new_ron_data<-pac_res

#cell_type <- "neuron"  # your chosen cell type
#cell_type="esc"
#df<-pac_res
pac_list2<-lapply(pac_list, function(df){
  rownames(df)<-df$X
  df$X<-NULL
  df<-t(df)
  df <- df[1:(nrow(df) - 3), ]
  df<-data.frame(df)
  df$sample_name<-rownames(df)
  return(df)
})

obj<-Reduce(rbind, pac_list2)
obj$Sample<-obj$sample_name

obj$Sample <- gsub("[.]", "-", obj$Sample)  # Step 1: replace dots with dashes
obj$Sample <- ifelse(grepl("^[0-9]", obj$Sample),
                    paste0("S", obj$Sample),
                    obj$Sample)             # Step 2: prefix 'S' if starts with a number



#merged<-readRDS("ipsc_and_hart30_merged_tsne_pca.rds")
#merged@meta.data<-left_join(merged@meta.data, obj, "Sample")
#merged2 <- subset(merged, subset = !is.na(Sample))
merged<-readRDS("ipsc_object_merged_pca_tsne_norm.rds")
met<-merged@meta.data
met<-left_join(met, obj, "Sample")
merged$esc<-met$esc
merged$neuron<-met$neuron

#rownames(met)<-met$Sample
#merged@meta.data<-met
#valid_cells <- rownames(merged@meta.data[!is.na(merged@meta.data$esc), ])

# Subset the Seurat object using these valid cells
#merged2 <- merged[, valid_cells]


FeaturePlot(merged, features="esc", reduction="tsne")

merged2$expression<-ifelse(merged2$esc>=0.9, "high_esc", "low_esc")
hipsci<-subset(merged2,Group=="HIPSCI_BATCH1")
bd2<-subset(merged2, Group=="BD2_BATCH1")
val<-subset(merged2, Group=="CANCER_VALIDATION")
hart30<-subset(merged2, Group=="Hart30")
objs<-list(hipsci, bd2, val)
de_objs<-lapply(objs, function(x){
  x<-FindMarkers(
    x,
    ident.1 = "high_esc",
    ident.2 = "low_esc",
    group.by = "expression",
    logfc.threshold = 0,      # Set to 0 to include all genes, adjust as desired
    min.pct = 0.1             # Optional: filter by minimum percent expressed
  )
})
names(de_objs)<-c("hipsci","bd2","val")
# Run differential expression using the 'expression' metadata column
de_results <- FindMarkers(
  merged2,
  ident.1 = "high_esc",
  ident.2 = "low_esc",
  group.by = "expression",
  logfc.threshold = 0,      # Set to 0 to include all genes, adjust as desired
  min.pct = 0.1             # Optional: filter by minimum percent expressed
)


saveRDS(de_results, "de_genes_high_pacnet_esc_vs_low_pacnet_esc_thresh_0point98.rds")

seu<-merged2
# Assume your Seurat object is named 'seu'
# First, extract the 'esc' column into a data frame
esc_data <- seu@meta.data %>%
  dplyr::select(esc) %>%
  dplyr::filter(!is.na(esc))  # Remove NAs if any

# Calculate summary stats
summary_stats <- esc_data %>%
  summarise(
    Mean = mean(esc),
    Median = median(esc),
    SD = sd(esc)
  )

# Print summary for a peek behind the curtain
print(summary_stats)

# Now craft the histogram
ggplot(esc_data, aes(x = esc)) +
  geom_histogram(bins=100, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = mean(esc)), color = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = median(esc)), color = "purple", linetype = "dotted", size = 1) +
  labs(
    title = "Distribution of ESC Values",
    x = "ESC Score",
    y = "Cell Count"
  ) +
  theme_minimal() +
  annotate("text", x = summary_stats$Mean, y = Inf, label = paste0("Mean: ", round(summary_stats$Mean, 2)), vjust = -1, color = "red") +
  annotate("text", x = summary_stats$Median, y = Inf, label = paste0("Median: ", round(summary_stats$Median, 2)), vjust = -1, color = "purple")


de_results2<-subset(de_results, p_val_adj<=0.05)
head(de_results2[order(de_results2$avg_log2FC, decreasing = TRUE), ], 50)->esc_genes


de_results$gene <- rownames(de_results)

padj_thresh <- 0.05

de_results$color <- "not_significant"
de_results$color[de_results$p_val_adj < padj_thresh & de_results$avg_log2FC > 0] <- "up"
de_results$color[de_results$p_val_adj < padj_thresh & de_results$avg_log2FC < 0] <- "down"

de_results$color <- factor(de_results$color, levels = c("down", "not_significant", "up"))

ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("down" = "blue", "not_significant" = "black", "up" = "red")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot: High vs. Low Expression",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value",
    color = "Expression Change"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray")


# filtered3<-subset(filtered2, esc>=0.98)
# FeaturePlot(filtered3, features="esc")
# 
# expr_matrix <- GetAssayData(filtered3, slot = "data")  # Use "counts" for raw counts
# gene_totals <- Matrix::rowSums(expr_matrix)
# top_genes <- sort(gene_totals, decreasing = TRUE)
# 
# head(top_genes, 20)


library(tidyr)

# Assuming your data frame is called `df`
df<-merged@meta.data
# Step 1: Give each sample a unique ID
df <- df %>%
  mutate(sample_id = paste0("Sample_", row_number()))

# Step 2: Pivot longer
df_long <- df %>%
  pivot_longer(cols = c(esc, neuron),
               names_to = "cell_type",
               values_to = "score")

# Step 3: Plot — one bar per sample, two per sample (esc + neuron), colored by cell type
ggplot(df_long, aes(x = sample_id, y = score, fill = cell_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("esc" = "#ADD8E6", "neuron" = "#FFC0CB")) +
  labs(
    x = "Sample",
    y = "Score",
    fill = "Cell Type",
    title = "ESC and Neuron Scores per Sample, Colored by Cell Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(df_long, aes(x = sample_name, y = score, fill = cell_type)) +
  geom_bar(stat = "identity") +  # stacking is default behavior
  scale_fill_manual(values = c("esc" = "#ADD8E6", "neuron" = "#FFC0CB")) +
  labs(
    x = "Sample",
    y = "Score",
    fill = "Cell Type",
    title = "Stacked ESC and Neuron Scores per Sample"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
