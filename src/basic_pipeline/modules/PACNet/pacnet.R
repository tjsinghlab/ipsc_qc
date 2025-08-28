library(CellNet)
library(cancerCellNet)
library(ggplot2)
library(patchwork)

####################################
########### Load in Data ########### 
####################################

#Load in training data from PACnet repo 
expTrain <- utils_loadObject("Hs_expTrain_Jun-20-2017.rda") #run aws s3 cp s3://cellnet-rnaseq/ref/cnproc/HS/Hs_stTrain_Jun-20-2017.rda . --no-sign-request
stTrain <- utils_loadObject("Hs_stTrain_Jun-20-2017.rda") # run aws s3 cp s3://cellnet-rnaseq/ref/cnproc/HS/Hs_expTrain_Jun-20-2017.rda . --no-sign-request

#Load in query matrix (rownames are genes, col names are samples)
queryExpDat<-read.csv("query_matrix.csv",row.names=1)

#Load in query sample data (cols sample_name and description1)
querySampTab<-read.csv("query_meta.csv", row.names=1)

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
save(my_classifier, file="/test_cellnet_classifier_100topGenes_100genePairs.rda")

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

classMatrixEx <- broadClass_predict(cnProc = cnProc, expDat = queryExpDat, nrand = 10) 
grp_names1 <- c(as.character(querySampTab$description1), rep("random", 10))
names(grp_names1) <- c(as.character(rownames(querySampTab)), paste0("rand_", c(1:10)))

# Re-order classMatrixQuery to match order of rows in querySampTab
classMatrixEx <- classMatrixEx[,names(grp_names1)]

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

png(file="heatmapEx.png", height=12, width=20, units="in", res=300)
heatmapRef_split(classMatrixEx, querySampTab)
dev.off()

