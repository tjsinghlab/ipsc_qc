library(Seurat)

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

