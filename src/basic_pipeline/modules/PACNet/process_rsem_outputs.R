library(dplyr)
library(data.table)
library(stringr)
library(biomaRt)
library(purrr)

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-p", "--project"), type="character", default="default_project",
              help="Project name", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=".",
              help="Output directory", metavar="character"),
  make_option(c("-i", "--input_dir"), type="character", default=".",
              help="Input directory", metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

project_name <- opt$project
output_dir <- opt$output_dir
input_dir  <- opt$input_dir

# Example use
metadata <- setNames(
  data.frame(colnames(full_matrix2), rep(project_name, ncol(full_matrix2))),
  c("sample_name", "description1")
)


setwd(output_dir)

#Specify RSEM output directory. File names will be used as sample names.
rsem_output_directory<-input_dir

#Read in and name files
mats<-lapply(list.files(rsem_output_directory, pattern = "genes", full.names = T), fread)
names(mats)<-str_split_i(list.files(rsem_output_directory, pattern = "genes"),"[.]",1)

genes<-lapply(mats, function(x){
  x<-x$gene_id
  return(x)
})

mart<-useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
meta.genes<-lapply(genes, function(x){
  y<-getBM(attributes = c("ensembl_gene_id","external_gene_name", "description"), 
           values = x, mart = mart)
  return(y)
}) 

mats<-Map(function(og_map, new_genes){
  #og_map<-data.frame(og_map[["gene_id"]], og_map[["TPM"]]) #If you want to use TPMs
  og_map<-data.frame(og_map[["gene_id"]], og_map[["expected_count"]]) #If you want to use RSEM's expected counts
  colnames(og_map)<-c("ensembl_gene_id","count")
  og_map$ensembl_gene_id<-str_split_i(og_map$ensembl_gene_id, "[.]", 1)
  out<-full_join(og_map, new_genes, by="ensembl_gene_id")
  return(out)
},mats, meta.genes)

mats<-lapply(mats, function(df){
  df<-df[,-c(1,4)] #Removing emsembl ID and description columns
  df[df==''] <- NA #Replace empty entries with NA
  df<-df[complete.cases(df), ]
  #Removing repetitive and/or irrelevant genes:
  df<-df[!grepl("^[0-9]", df$external_gene_name) &
           !grepl("Metazoa", df$external_gene_name) &
           !grepl("^RNU", df$external_gene_name) &
           !grepl("^SNOR", df$external_gene_name) &
           df$external_gene_name != "Y_RNA" &
           df$external_gene_name != "Vault" &
           !grepl("^U[0-9]+$", df$external_gene_name), ] 
  return(df)
})

mats<-Map(function(df, namex){
  colnames(df)<-c(namex, "gene")
  return(df)
},mats, names(mats))

#If any sample starts with a number, this will cause issues downstream.
#Make sure any name that starts with a number is preceded by "S"
names(mats) <- sapply(names(mats), function(name) {
  if (grepl("^[0-9]", name)) {
    paste0("S", name)
  } else {
    name
  }
})

process_df <- function(df) {
  # Identify duplicate genes
  duplicates <- duplicated(df$gene) | duplicated(df$gene, fromLast = TRUE)
  # Keep rows where the "value" column (or the first column) is not zero
  df <- df[!(duplicates & df[[1]] == 0), ]
  return(df)
}

#To really make sure all duplicates are dealt with, if duplicated genes persist, take the max of the values:
#For entries where two different Ensembl Gene IDs mapped onto the same gene, we'll take the max of the two values (this is not common)
aggene<-function(df){
  # Convert to data.table for efficiency
  dt <- as.data.table(df)
  # Calculate the maximum for duplicated genes (the duplicated value is generally 0)
  dt <- dt[, .(value = max(.SD[[1]])), by = gene]
  return(dt)
}

mats<-lapply(mats, data.frame)
mats<-lapply(mats, process_df)
mats<-lapply(mats, aggene)
mats<-Map(function(df, namex){
  colnames(df)<-c("gene",namex)
  return(df)
},mats, names(mats))

full_matrix<-mats %>% purrr::reduce(full_join, by="gene")

full_matrix<-data.frame(full_matrix)
genes<-full_matrix$gene
full_matrix$gene=NULL
rownames(full_matrix)<-genes

metadata <- setNames(
  data.frame(colnames(full_matrix), rep(project_name, ncol(full_matrix))),
  c("sample_name", "description1")
)

write.csv(full_matrix, "query_matrix.csv")
write.csv(metadata, "query_meta.csv")