ipsc_merged <- readRDS("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ipsc_merged_tsne_no_sct_no_harmony_increasedresolution0point2.rds")

ipsc_merged<-merged
all_features <- rownames(ipsc_merged@assays$RNA@layers$data)  # Get all features (genes)
ipsc_merged <- ScaleData(ipsc_merged, features = all_features)
ipsc_merged<-RunPCA(ipsc_merged, features=all_features)

DimPlot(ipsc_merged, reduction="pca", group.by="Group")
DimPlot(ipsc_merged, reduction="pca", group.by="characteristic")



metadata <- ipsc_merged@meta.data
metadata$Cell <- rownames(metadata)  # Add cell names as a separate column

# Create the bar plot with each cell as a separate bar, grouped by Group
ggplot(metadata, aes(x = reorder(Cell, Group), y = nFeature_RNA, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Cells", y = "nFeature_RNA", fill = "Group") +
  theme(axis.text.x = element_blank(),  # Hide x-axis text for clarity
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  facet_wrap(~ Group, scales = "free_x", nrow = 1)+
  ggtitle("Number of Genes per Line")# Group bars into separate facets


ggplot(metadata, aes(x = reorder(Cell, Group), y = nCount_RNA, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Cells", y = "nCount_RNA", fill = "Group") +
  theme(axis.text.x = element_blank(),  # Hide x-axis text for clarity
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  facet_wrap(~ Group, scales = "free_x", nrow = 1)+
  ggtitle("Number of Cells per Line")# Group bars into separate facets



# Create scatter plot
ggplot(metadata, aes(x = nFeature_RNA, y = nCount_RNA, color = Group)) +
  geom_point(alpha = 0.7) +  # Adjust transparency for visibility
  theme_minimal() +
  labs(x = "nFeature_RNA", y = "nCount_RNA", color = "Group") +
  theme(legend.position = "right") +
  scale_color_brewer(palette = "Set2")+
  ggtitle("Genes x Cells Per Line")



# Create the bar plot with improved colors
ggplot(metadata, aes(x = reorder(Cell, Group), y = percent.mito, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Cells", y = "Percent Mitochondrial", fill = "Group") +
  theme(axis.text.x = element_blank(),  # Hide x-axis text for clarity
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  facet_wrap(~ Group, scales = "free_x", nrow = 1) +  # Group bars into separate facets
  scale_fill_brewer(palette = "Set2")  +
  ggtitle("Percent Mitochondrial Genes by Line")




# Compute IQR and define thresholds for outliers
iqr_nFeature <- IQR(metadata$nFeature_RNA)
iqr_nCount <- IQR(metadata$nCount_RNA)

q1_nFeature <- quantile(metadata$nFeature_RNA, 0.25)
q3_nFeature <- quantile(metadata$nFeature_RNA, 0.75)
q1_nCount <- quantile(metadata$nCount_RNA, 0.25)
q3_nCount <- quantile(metadata$nCount_RNA, 0.75)

lower_threshold_nFeature <- q1_nFeature - 1.5 * iqr_nFeature
upper_threshold_nFeature <- q3_nFeature + 1.5 * iqr_nFeature
lower_threshold_nCount <- q1_nCount - 1.5 * iqr_nCount
upper_threshold_nCount <- q3_nCount + 1.5 * iqr_nCount

# Identify outliers (both low and high)
metadata$Outlier <- (metadata$nFeature_RNA < lower_threshold_nFeature | 
                       metadata$nFeature_RNA > upper_threshold_nFeature |
                       metadata$nCount_RNA < lower_threshold_nCount | 
                       metadata$nCount_RNA > upper_threshold_nCount)

# Scatter plot with labeled outliers
library(ggrepel)
ggplot(metadata, aes(x = nFeature_RNA, y = nCount_RNA, color = Group)) +
  geom_point(alpha = 0.7) +  # Regular points
  geom_text_repel(aes(label = ifelse(Outlier, Sample, "")), size = 3) +  # Label outliers
  theme_minimal() +
  labs(x = "nFeature_RNA", y = "nCount_RNA", color = "Group") +
  theme(legend.position = "right") +
  scale_color_viridis_d(option = "plasma")+
  scale_color_brewer(palette = "Set2")+
  ggtitle("Genes x Cells Per Line")



ipsc<-CellCycleScoring(ipsc_merged, g2m.features = cc.genes$g2m.genes, s.features=cc.genes$s.genes)
metadata <- ipsc@meta.data
metadata$Cell <- rownames(metadata)  # Add cell names as a separate column

ggplot(metadata, aes(x = reorder(Cell, characteristic), y = G2M.Score, fill = characteristic)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Cells", y = "G2M Score", fill = "characteristic") +
  theme(axis.text.x = element_blank(),  # Hide x-axis text for clarity
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  facet_wrap(~ characteristic, scales = "free_x", nrow = 1)+
  ggtitle("G2M Score per Line")# Group bars into separate facets


ggplot(metadata, aes(x = reorder(Cell, characteristic), y = S.Score, fill = characteristic)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Cells", y = "S Score", fill = "characteristic") +
  theme(axis.text.x = element_blank(),  # Hide x-axis text for clarity
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  facet_wrap(~ characteristic, scales = "free_x", nrow = 1)+
  ggtitle("S Score per Line")# Group bars into separate facets


########################################

# library(Seurat)
# library(ggplot2)
# library(dplyr)

# Grab the PCA embeddings
pca_coords <- Embeddings(merged, reduction = "pca") %>%
  as.data.frame() %>%
  mutate(cell = rownames(.))

# Merge with metadata to get sample names
pca_coords <- pca_coords %>%
  left_join(merged@meta.data %>% 
              dplyr::select(Sample) %>%
              mutate(cell = rownames(.)), 
            by = "cell")

# Identify outliers - here's a basic method using Mahalanobis distance (feel free to replace!)
# For dramatic flair, let's say the top 1% furthest from the center are our outliers
center <- colMeans(pca_coords[, c("PC_1", "PC_2")])
cov_matrix <- cov(pca_coords[, c("PC_1", "PC_2")])
pca_coords$mahal_dist <- mahalanobis(pca_coords[, c("PC_1", "PC_2")], center, cov_matrix)

threshold <- quantile(pca_coords$mahal_dist, 0.99)  # Top 1% most distant
pca_coords$outlier <- pca_coords$mahal_dist > threshold

# Plot with DimPlot first for background
DimPlot(merged, reduction = "pca") +
  geom_text(data = pca_coords %>% filter(outlier),
            aes(x = PC_1, y = PC_2, label = Sample),
            color = "red", size = 3, vjust = -0.5)

#Mature neuron markers
FeaturePlot(merged, features=c("MAP2","TUBB3","NEUN","RBFOX3","SYN1","SYP","DCX","GAP43","VGLUT1"),cols=c("gray","red"),reduction="pca")
#Pluripotent markers
FeaturePlot(merged, features=c("POU5F1","OCT4","NANOG","SOX2","LIN28A","LIN28B","DPPA4","TDGF1","ZFP42","REX1"),cols=c("gray","red"),reduction="pca")

###############################
###############################


###############
merged<-readRDS("ipsc_object_merged_pca_tsne_norm.rds")

library(ggplot2)
library(dplyr)

# Combine into a data frame
plot_df <- data.frame(
  esc = esc,
  Project = seurat_object$Project
)

# Compute IQR per project to define outliers
outlier_df <- merged@meta.data %>%
  group_by(Project) %>%
  mutate(Q1 = quantile(esc, 0.25),
         Q3 = quantile(esc, 0.75),
         IQR = Q3 - Q1,
         is_outlier = esc < Q1 - 1.5*IQR | esc > Q3 + 1.5*IQR)

# Plot
ggplot(outlier_df, aes(x = Project, y = esc)) +
  geom_violin(fill = "skyblue", alpha = 0.5) +
  geom_jitter(aes(color = is_outlier), width = 0.2, alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Tilt labels
  labs(title = "ESC Value Distribution by Project",
       subtitle = "is_outlier = esc < Q1 - 1.5*IQR | esc > Q3 + 1.5*IQR",
       y = "ESC",
       x = "Project")



library(ggridges)

ggplot(merged@meta.data, aes(x = esc, y = Project, fill = Project)) +
  geom_density_ridges(scale = 3, alpha = 0.6) +
  theme_minimal() +
  labs(title = "ESC Distributions Across Projects",
       x = "ESC", y = "Project")


plot_df<-merged@meta.data
plot_df$rank <- rank(-plot_df$esc)  # descending rank

ggplot(plot_df, aes(x = rank, y = esc, color = Project)) +
  geom_point(alpha = 0.7, size = 1) +
  theme_minimal() +
  labs(title = "Ranked ESC Values Across Projects",
       x = "ESC Rank",
       y = "ESC Value")


library(ggbeeswarm)

ggplot(outlier_df, aes(x = Project, y = esc)) +
  geom_boxplot(outlier.shape = NA, fill = "lightpink", alpha = 0.4) +
  geom_beeswarm(aes(color = is_outlier), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "firebrick")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "ESC Boxplots with Swarming Outliers",
       x = "Project", y = "ESC")


top_outliers <- outlier_df %>%
  filter(is_outlier) %>%
  top_n(400, esc)

ggplot(top_outliers, aes(x = reorder(Sample, esc), y = esc, color = Project)) +
  geom_segment(aes(xend = Sample, y = 0, yend = esc), size = 0.8) +
  geom_point(size = 3) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 400 ESC Outliers by Sample",
       x = "Sample", y = "ESC")



esc_threshold <- mean(plot_df$esc, na.rm = TRUE) + 2 * sd(plot_df$esc, na.rm = TRUE)

plot_df$outlier_label <- ifelse(plot_df$esc > esc_threshold,
                                as.character(seurat_object$Sample),
                                NA)

ggplot(plot_df, aes(x = Project, y = esc)) +
  geom_jitter(aes(color = Project), width = 0.3, alpha = 0.5) +
  geom_text_repel(aes(label = outlier_label), size = 2.5, max.overlaps = 10) +
  geom_hline(yintercept = esc_threshold, linetype = "dashed", color = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "ESC Outliers Above Threshold",
       subtitle = paste("Threshold =", round(esc_threshold, 2)),
       y = "ESC Value")



library(pheatmap)

top_samples <- plot_df %>%
  top_n(50, esc) %>%
  arrange(desc(esc))

mat <- matrix(top_samples$esc, nrow = 1)
colnames(mat) <- top_samples$Sample

pheatmap(mat,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = colorRampPalette(c("skyblue", "white", "pink"))(100),
         main = "Top 50 ESC Samples")

