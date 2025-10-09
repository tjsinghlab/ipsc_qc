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
  make_option(c("-s", "--sample"), type="character",
              help="Sample name (basename before underscore in FASTQ)"),
  make_option(c("-p", "--project"), type="character", default="default_project",
              help="Project name"),
  make_option(c("-o", "--output_dir"), type="character", default="outputs",
              help="Top-level output directory"),
  make_option("--ref_dir", type="character", default="/ref",
              help="Reference directory inside Docker image (default: /ref)")
)

opt <- parse_args(OptionParser(option_list=option_list))

sample_name <- opt$sample
project_name <- opt$project
output_dir <- opt$output_dir
ref_dir <- opt$ref

# ---------------- Locate RSEM outputs automatically ----------------
message("Searching for RSEM outputs under: ", output_dir)

# Find all subdirectories inside output_dir (each should correspond to a sample)
sample_dirs <- list.dirs(output_dir, recursive = FALSE, full.names = TRUE)

# Collect RSEM genes.results files from each sample directory
rsem_files <- unlist(lapply(sample_dirs, function(sdir) {
  list.files(
    file.path(sdir, "RSEM_outputs"),
    pattern = "genes.results(\\.gz)?$",
    full.names = TRUE
  )
}))

if (length(rsem_files) == 0) {
  stop("No RSEM genes.results files found under ", output_dir)
}

message("Found ", length(rsem_files), " RSEM file(s).")

# ---------------- Prepare folders ----------------
pacnet_outdir <- file.path(output_dir, "pacnet")
if (!dir.exists(pacnet_outdir)) dir.create(pacnet_outdir, recursive = TRUE)

# ---------------- Step 1: Build Expression Matrix ----------------
mats <- lapply(rsem_files, fread)
names(mats) <- str_split_i(basename(rsem_files), "[.]", 1)

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

# ---------------- Clean up repetitive genes ----------------
mats <- lapply(mats, function(df) {
  df <- df[,-c(1,4)]
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

# Collapse duplicates
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

write.csv(full_matrix, file.path(pacnet_outdir, "query_matrix.csv"), row.names = TRUE)
write.csv(metadata,    file.path(pacnet_outdir, "query_meta.csv"))

#Load in training data from PACnet repo 
message(">>> Running PACNet classifier")

expTrain <- utils_loadObject(file.path(ref_dir, "Hs_expTrain_Jun-20-2017.rda"))
stTrain  <- utils_loadObject(file.path(ref_dir, "Hs_stTrain_Jun-20-2017.rda"))

#Load in query matrix (rownames are genes, col names are samples)
queryExpDat<-read.csv(file.path(pacnet_outdir, "query_matrix.csv"),row.names=1, check.names=FALSE)

#Load in query sample data (cols sample_name and description1)
querySampTab<-read.csv(file.path(pacnet_outdir, "query_meta.csv"), row.names=1, check.names=FALSE)

#Ensure overlapping genes
iGenes <- intersect(rownames(expTrain), rownames(queryExpDat))

#Subset training expression matrix based on iGenes
expTrain <- expTrain[iGenes,]

set.seed(99) #Set seed for reproducibility with classifier
stList <- splitCommon_proportion(sampTab = stTrain, proportion = 0.66, dLevel = "description1") # Use 2/3 of training data for training and 1/3 for validation
stTrainSubset <- stList$trainingSet
expTrainSubset <- expTrain[,rownames(stTrainSubset)]

stValSubset <- stList$validationSet
expValSubset <- expTrain[,rownames(stValSubset)]

#####################################
######### Create Classifier #########
#####################################

system.time(my_classifier <- broadClass_train(stTrain = stTrainSubset, 
                                              expTrain = expTrainSubset, 
                                              colName_cat = "description1", 
                                              colName_samp = "sra_id", 
                                              nRand = 70, 
                                              nTopGenes = 100, 
                                              nTopGenePairs = 100, 
                                              nTrees = 2000, 
                                              stratify=TRUE, 
                                              sampsize=25, # Must be less than the smallest number in table(stTrainSubset$description1)
                                              quickPairs=TRUE)) # Increasing the number of top genes and top gene pairs increases the resolution of the classifier but increases the computing time
save(my_classifier, file = file.path(pacnet_outdir, "classifier.rda"))

stValSubsetOrdered <- stValSubset[order(stValSubset$description1), ] #order samples by classification name
expValSubset <- expValSubset[,rownames(stValSubsetOrdered)]
cnProc <- my_classifier$cnProc #select the cnProc from the earlier class training

classMatrix <- broadClass_predict(cnProc, expValSubset, nrand = 60) 
stValRand <- addRandToSampTab(classMatrix, stValSubsetOrdered, desc="description1", id="sra_id")

grps <- as.vector(stValRand$description1)
names(grps)<-rownames(stValRand)

# Create validation heatmap
png(file="classification_validation_hm.png", height=6, width=10, units="in", res=300)
ccn_hmClass(classMatrix, grps=grps, fontsize_row=10)
dev.off()

genePairs <- cnProc$xpairs

message(">>> Query expression matrix dims: ", paste(dim(queryExpDat), collapse=" x "))
message(">>> First few query sample names: ", paste(colnames(queryExpDat)[1:5], collapse=", "))

#Gene to gene comparison of each gene pair in the expression table
expTransform <- query_transform(expTrainSubset, genePairs)
avgGenePair_train <- avgGeneCat(expDat = expTransform, sampTab = stTrainSubset, 
                                dLevel = "description1", sampID = "sra_id")

genePairs_val <- query_transform(expValSubset, genePairs)
geneCompareMatrix <- makeGeneCompareTab(queryExpTab = genePairs_val,
                                        avgGeneTab = avgGenePair_train, geneSamples = genePairs)
val_grps <- stValSubset[,"description1"]
val_grps <- c(val_grps, colnames(avgGenePair_train))
names(val_grps) <- c(rownames(stValSubset), colnames(avgGenePair_train))

xpairs_list <- vector("list", 14) 
for (pair in rownames(avgGenePair_train)) {
  for (j in 1:ncol(avgGenePair_train)) {
    if (avgGenePair_train[pair,j] >= 0.5) {
      if (is.null(xpairs_list[[j]])) {
        xpairs_list[[j]] <- c(pair)
      } else { 
        xpairs_list[[j]] <- c(xpairs_list[[j]], pair)
      }
    }  
  }
}

xpair_names <- colnames(avgGenePair_train)
xpair_names <- sub(pattern="_Avg", replacement="", x=xpair_names)
names(xpairs_list) <- xpair_names

for (type in names(xpairs_list)) {
  names(xpairs_list[[type]]) <- xpairs_list[[type]]
}
save(xpairs_list, file="Hs_xpairs_list.rda")

#####################################
###### Querying the Classifier ######
#####################################
# Query predictions
nrand_val <- min(10, ncol(queryExpDat))
queryExpDat[is.na(queryExpDat)] <- 0

rownames(querySampTab)<-querySampTab$sample_name

classMatrixEx <- broadClass_predict(cnProc = cnProc, expDat = queryExpDat, nrand = nrand_val) 
grp_names1 <- c(as.character(querySampTab$description1), rep("random", nrand_val))
names(grp_names1) <- c(as.character(rownames(querySampTab)), paste0("rand_", c(1:nrand_val)))

cat("colnames(classMatrixEx):", colnames(classMatrixEx), "\n")
cat("names(grp_names1):", names(grp_names1), "\n")

common_names <- intersect(colnames(classMatrixEx), names(grp_names1))

cat("common_names:", common_names, "\n")

classMatrixEx <- classMatrixEx[, common_names, drop = FALSE]

cat("colnames(classMatrixEx):", colnames(classMatrixEx), "\n")

grp_names1 <- grp_names1[common_names]

cat("names(grp_names1):", names(grp_names1), "\n")


# Re-order classMatrixQuery to match order of rows in querySampTab
#classMatrixEx <- classMatrixEx[,names(grp_names1)]

write.csv(classMatrixEx, file=file.path(pacnet_outdir, "classification_scores.csv"))
write.csv(querySampTab, file="classification_sample_table.csv")

####################################
######### Generate Heatmap #########
####################################

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

png(file = file.path(pacnet_outdir, "PACNet_heatmap.png"), height = 12, width = 20, units = "in", res = 300)
print(heatmapRef_split(classMatrixEx, querySampTab))  # <-- must explicitly print
dev.off()

message(">>> PACNet completed. Outputs in: ", pacnet_outdir)

