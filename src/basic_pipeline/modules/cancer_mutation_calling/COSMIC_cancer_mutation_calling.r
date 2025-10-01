#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(data.table)
  library(stringr)
})

option_list <- list(
  make_option("--sample", type = "character", default = NULL,
              help = "Sample name (used for labeling outputs)", metavar = "character"),
  make_option("--ref_dir", type = "character", default = "/ref",
            help = "Reference directory inside Docker image (default: /ref)"),
  make_option("--output_dir", type = "character", default = ".",
              help = "Output directory for results", metavar = "character"),
  make_option("--cosmic_dir", type = "character", default = ".",
              help = "Path to cosmic database downloaded by user", metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

sample_id <- opt$sample
ref_dir <- opt$ref_dir
output_dir <- opt$output_dir
cosmic_dir <- opt$cosmic_dir
dir.create(file.path(output_dir, "cosmic_calling"), recursive = TRUE, showWarnings = FALSE)

#extract COSMIC database
untar(file.path(cosmic_dir,"Cosmic_CancerGeneCensus_Tsv_v101_GRCh38.tar"), exdir = cosmic_dir)

Organism="Human"
print("Reading VCF file...")
#Directory  
vcf_file = file.path(output_dir, "variant_calling", paste0(sample_id, ".variant_filtered.vcf.gz"))
output_dir <- file.path(output_dir, "cosmic_calling")
readData = read.delim(vcf_file,as.is=T)
readData=as.character(readData[-c(1:which(readData=="#CHROM")-1),1])

print("Data processing...")
jump = 10
startChr = 1+jump
startPos = 2+jump
startInfo = 10+jump

print("Extracting VCF data...")
len = length(readData)
chrRegex ="^chr(\\w+)$"
infoRegex = "^([01])\\/([01]):(\\d+)\\,(\\d+):(\\d+):\\d+:\\d+\\,\\d+\\,\\d+$"

chrVector =  readData[startChr]
posVector =  readData[startPos]
infoVector = readData[startInfo]

print("Looping through VCF entries...")
while (startInfo + jump < len) {
  startChr = startChr + jump
  startPos = startPos + jump
  startInfo = startInfo + jump
  chrVector = append(chrVector, readData[startChr])
  posVector = append(posVector, readData[startPos])
  infoVector = append(infoVector, readData[startInfo])
}

print("Defining X and Y chr numbers...")
chrNum = gsub(chrRegex, "\\1", chrVector)
if (Organism=="Human"){
  chrNum[chrNum=="X"]="23"
  chrNum[chrNum=="Y"]="24"}

print("Extracting allele info...")
chrNum =  as.numeric(chrNum)
Karyotape = 10*abs(as.numeric(gsub(infoRegex, "\\1", infoVector))-as.numeric(gsub(infoRegex, "\\2", infoVector)))
AD1 = as.numeric(gsub(infoRegex, "\\3", infoVector))
AD2 = as.numeric(gsub(infoRegex, "\\4", infoVector))
DP = as.numeric(gsub(infoRegex, "\\5", infoVector))

posVector = as.numeric(posVector)

print("Constructing dataframe from vcf...")
table = data.frame("chr" = chrNum, "position" = posVector, "AD1" = AD1, "AD2" = AD2, "DP" = DP, "Karyotape" = Karyotape)

print("Ordering dataframe...")
table$chr=as.numeric(table$chr)
table=table[order(table$chr,table$position),]
table=table[table$chr>0,]
colnames(table)<-c("CHROMOSOME","POS","AD1","AD2","DP","KARYOTYPE")

print("Reading COSMIC file...")
cosmic = fread(file.path(cosmic_dir, "Cosmic_MutantCensus_v101_GRCh38.tsv.gz"))
db<-cosmic[,-c(2:8,19:25)]
colnames(table)<-c("CHROMOSOME","POS","AD1","AD2","DP","Karyotype")

print("Reading gene list...")
genes<-read.table(file.path(ref_dir, "genes.txt"))

print("Subsetting on genes of interest...")
db<-db[db$GENE_SYMBOL%in%genes$V1]
db<-subset(db, GENE_SYMBOL%in%genes$V1)
# Function to match POS with genome regions and retain all columns
match_genes_full <- function(df1_chr, df2_chr) {
  result <- list()  # Store results
  
  for (i in seq_len(nrow(df1_chr))) {
    pos <- df1_chr$POS[i]
    
    # Find matches in df2
    matches <- df2_chr %>%
      filter(pos >= GENOME_START & pos <= GENOME_STOP)
    
    if (nrow(matches) > 0) {
      # Merge matched rows with the current row of df1
      matched_rows <- cbind(df1_chr[rep(i, nrow(matches)), ], matches)
      result[[i]] <- matched_rows
    }
  }
  
  if (length(result) > 0) {
    return(do.call(rbind, result))  # Combine list into a dataframe
  } else {
    return(NULL)
  }
}

#########################################
VCF_FILE=vcf_file
tmp_vcf<-readLines(VCF_FILE)
tmp_vcf_data<-read.table(VCF_FILE, stringsAsFactors = FALSE)

# filter for the columns names
tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
names(tmp_vcf_data)<-vcf_names
vcf<-tmp_vcf_data
colnames(vcf)[1]<-"CHROMOSOME"

vcf$CHROMOSOME<-str_split_i(vcf$CHROMOSOME, "r", 2)
vcf$CHROMOSOME<-as.numeric(vcf$CHROMOSOME)
#########################################

table1<-vcf
# Split df1 and df2 by chromosome
df1_split <- split(table1, table$CHROMOSOME)
df2_split <- split(db, db$CHROMOSOME)

print("Constructing final dataframe...")
# Iterate over chromosomes and find matches, keeping all columns
final_results <- lapply(names(df1_split), function(chr) {
  if (chr %in% names(df2_split)) {
    match_genes_full(df1_split[[chr]], df2_split[[chr]])
  }
})

# Combine results into a final dataframe
final_df <- do.call(rbind, final_results)
summary_genes<-final_df$GENE_SYMBOL %>% unique()

print("Creating plot...")
df_plot<-final_df
df_plot$CHROMOSOME=NULL
df_plot <- df_plot %>%
  group_by(GENE_SYMBOL, MUTATION_DESCRIPTION) %>%
  summarise(unique_pos_count = n_distinct(POS), .groups = "drop")

plotty<-ggplot(df_plot, aes(x = GENE_SYMBOL, y = unique_pos_count, fill = MUTATION_DESCRIPTION)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) + 
  labs(
    x = "Gene Symbol", 
    y = "Unique Position Count", 
    fill = "Mutation Type",
    title = "SNP Hits in Cancer Mutation Database Metrics",
    subtitle = "Cosmic_MutantCensus_Tsv_v101_GRCh38"
  ) +
  theme_minimal(base_size = 14) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title = element_text(face = "bold"), 
    axis.text = element_text(size = 10), 
    legend.title = element_text(size = 12, face = "bold"), 
    legend.text = element_text(size = 10), 
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic")
  ) 

print("Saving results...")
ggsave(
  filename="CancerMutationPlot.png",
  plot = plotty,
  device = "png",
  path = output_dir,
  width=16,
  height=8)

write.csv(final_df, file.path(output_dir, paste0(sample_id, "_CancerMutations.tsv")), sep = "\t", row.names = FALSE)
