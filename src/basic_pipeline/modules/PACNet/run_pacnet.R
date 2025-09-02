#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(data.table)
  library(stringr)
  library(biomaRt)
  library(purrr)
  library(CellNet)
  library(cancerCellNet)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

# ---------------- CLI options ----------------
option_list <- list(
  make_option(c("-p", "--project"), type="character", default="default_project",
              help="Project name"),
  make_option(c("-o", "--output_dir"), type="character", default="outputs",
              help="Output directory"),
  make_option(c("-i", "--input_dir"), type="character", default="rsem_outputs",
              help="Directory with RSEM outputs"),
  make_option(c("-r", "--ref"), type="character", default="ref",
              help="Reference directory for PACNet")
)

opt <- parse_args(OptionParser(option_list=option_list))
project_name <- opt$project
output_dir <- opt$output_dir
input_dir <- opt$input_dir
ref_dir <- opt$ref

# ---------------- Prepare folders ----------------
pacnet_dir <- file.path(output_dir, "pacnet")
if (!dir.exists(pacnet_dir)) dir.create(pacnet_dir, recursive = TRUE)

# ---------------- Step 1: Build Expression Matrix ----------------
message(">>> Reading RSEM outputs")

mats <- lapply(list.files(input_dir, pattern = "genes", full.names = TRUE), fread)
names(mats) <- str_split_i(list.files(input_dir, pattern = "genes"), "[.]", 1)

genes <- lapply(mats, function(x) x$gene_id)
mart  <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

meta.genes <- lapply(genes, function(x) {
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),
        values = x, mart = mart)
})

mats <- Map(function(og_map, new_genes) {
  og_map <- data.frame(og_map[["gene_id"]], og_map[["expected_count"]])
  colnames(og_map) <- c("ensembl_gene_id", "count")
  og_map$ensembl_gene_id <- str_split_i(og_map$ensembl_gene_id, "[.]", 1)
  full_join(og_map, new_genes, by = "ensembl_gene_id")
}, mats, meta.genes)

mats <- lapply(mats, function(df) {
  df <- df[,-c(1,4)] # drop ensembl + description
  df[df == ""] <- NA
  df <- df[complete.cases(df), ]
  df <- df[!grepl("^[0-9]", df$external_gene_name) &
             !grepl("Metazoa", df$external_gene_name) &
             !grepl("^RNU", df$external_gene_name) &
             !grepl("^SNOR", df$external_gene_name) &
             df$external_gene_name != "Y_RNA" &
             df$external_gene_name != "Vault" &
             !grepl("^U[0-9]+$", df$external_gene_name), ]
  df
})

# Handle sample names starting with numbers
names(mats) <- sapply(names(mats), function(name) {
  if (grepl("^[0-9]", name)) paste0("S", name) else name
})

# Process duplicates
process_df <- function(df) {
  duplicates <- duplicated(df$external_gene_name) | duplicated(df$external_gene_name, fromLast = TRUE)
  df <- df[!(duplicates & df[[1]] == 0), ]
  df
}
aggene <- function(df) {
  dt <- as.data.table(df)
  dt <- dt[, .(value = max(.SD[[1]])), by = external_gene_name]
  setnames(dt, c("gene", names(df)[1]))
  dt
}

mats <- lapply(mats, data.frame)
mats <- lapply(mats, process_df)
mats <- lapply(mats, aggene)
mats <- Map(function(df, namex) {
  colnames(df) <- c("gene", namex)
  df
}, mats, names(mats))

full_matrix <- mats %>% purrr::reduce(full_join, by = "gene")
genes <- full_matrix$gene
rownames(full_matrix) <- genes
full_matrix$gene <- NULL

metadata <- setNames(
  data.frame(colnames(full_matrix), rep(project_name, ncol(full_matrix))),
  c("sample_name", "description1")
)

write.csv(full_matrix, file.path(output_dir, "query_matrix.csv"))
write.csv(metadata,    file.path(output_dir, "query_meta.csv"))

# ---------------- Step 2: PACNet Classifier ----------------
message(">>> Running PACNet classifier")

expTrain <- utils_loadObject(file.path(ref_dir, "Hs_expTrain_Jun-20-2017.rda"))
stTrain  <- utils_loadObject(file.path(ref_dir, "Hs_stTrain_Jun-20-2017.rda"))

queryExpDat  <- read.csv(file.path(output_dir, "query_matrix.csv"), row.names = 1)
querySampTab <- read.csv(file.path(output_dir, "query_meta.csv"), row.names = 1)

iGenes <- intersect(rownames(expTrain), rownames(queryExpDat))
expTrain <- expTrain[iGenes, ]

set.seed(99)
stList <- splitCommon_proportion(stTrain, proportion = 0.66, dLevel = "description1")
stTrainSubset <- stList$trainingSet
expTrainSubset <- expTrain[, rownames(stTrainSubset)]
stValSubset <- stList$validationSet
expValSubset <- expTrain[, rownames(stValSubset)]

my_classifier <- broadClass_train(
  stTrain = stTrainSubset,
  expTrain = expTrainSubset,
  colName_cat = "description1",
  colName_samp = "sra_id",
  nRand = 70,
  nTopGenes = 100,
  nTopGenePairs = 100,
  nTrees = 2000,
  stratify = TRUE,
  sampsize = 25,
  quickPairs = TRUE
)

save(my_classifier, file = file.path(pacnet_dir, "classifier.rda"))

# Validation heatmap
classMatrix <- broadClass_predict(my_classifier$cnProc, expValSubset, nrand = 60)
stValRand <- addRandToSampTab(classMatrix, stValSubset, desc = "description1", id = "sra_id")
grps <- as.vector(stValRand$description1)
names(grps) <- rownames(stValRand)

png(file.path(pacnet_dir, "classification_validation_hm.png"), height=6, width=10, units="in", res=300)
ccn_hmClass(classMatrix, grps=grps, fontsize_row=10)
dev.off()

# Query predictions
classMatrixEx <- broadClass_predict(my_classifier$cnProc, queryExpDat, nrand = 10)
grp_names1 <- c(as.character(querySampTab$description1), rep("random", 10))
names(grp_names1) <- c(as.character(rownames(querySampTab)), paste0("rand_", c(1:10)))
classMatrixEx <- classMatrixEx[, names(grp_names1)]
write.csv(classMatrixEx, file.path(pacnet_dir, "classification_scores.csv"))

######### Generate Heatmap #########
#Heatmap code adapted from pacnet_utils heatmapRef function)
heatmapRef_split <- function(c_scoresMatrix, sampTab) {
  melted_tab = data.frame(
    "classificationScore" = numeric(),
    "sampleName" = character(),
    "tissueType" = character(),
    "description" = character(),
    stringsAsFactors = FALSE
  )
  for (sampleName in colnames(c_scoresMatrix)) {
    temp_cscore = c_scoresMatrix[, sampleName]
    # Handle metadata
    if (sampleName %in% rownames(sampTab)) {
      samp_descrip <- sampTab[sampleName, "description1"]
    } else {
      samp_descrip <- "random"
    }
    tempTab = data.frame(
      "classificationScore" = temp_cscore,
      "sampleName" = sampleName,
      "tissueType" = rownames(c_scoresMatrix),
      "description" = samp_descrip,
      stringsAsFactors = FALSE
    )
    melted_tab = rbind(melted_tab, tempTab)
  }
  #Label added randoms
  melted_tab$isRandom <- melted_tab$description == "random" | melted_tab$sampleName == "random"
  melted_tab$tissueType <- factor(
    melted_tab$tissueType,
    levels = sort(unique(melted_tab$tissueType), decreasing = TRUE)
  )
  # color palette
  cools <- colorRampPalette(c("black", "limegreen", "yellow"))(100)
  #Plot for query data
  p1 <- ggplot(subset(melted_tab, !isRandom), aes(x = sampleName, y = tissueType, fill = classificationScore)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = cools[1],
      high = cools[length(cools)],
      mid = cools[length(cools)/2],
      midpoint = 0.5,
      limit = c(0, 1),
      space = "Lab",
      name = "Classification Score"
    ) +
    xlab("Samples") + ylab("Tissue Types") +
    ggtitle("Query Samples") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
  #Plot for random controls
  p2 <- ggplot(subset(melted_tab, isRandom), aes(x = sampleName, y = tissueType, fill = classificationScore)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = cools[1],
      high = cools[length(cools)],
      mid = cools[length(cools)/2],
      midpoint = 0.5,
      limit = c(0, 1),
      space = "Lab",
      name = "Classification Score"
    ) +
    xlab("Random Samples") + ylab("Tissue Types") +
    ggtitle("Random Controls") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
  p <- p1 + p2 + plot_layout(ncol = 2, widths = c(2, 1))
  return(p)
}

png(file.path(pacnet_dir, "heatmapEx.png"), height=12, width=20, units="in", res=300)
heatmapRef_split(classMatrixEx, querySampTab)
dev.off()

# ---------------- Step 3: Metrics JSON ----------------
metrics <- list(
  project = project_name,
  n_samples = ncol(queryExpDat),
  n_genes   = nrow(queryExpDat),
  classifier_trees = 2000,
  top_genes = 100
)
write_json(metrics, file.path(pacnet_dir, "metrics.json"), pretty = TRUE)

message(">>> PACNet pipeline completed. Outputs written to: ", pacnet_dir)