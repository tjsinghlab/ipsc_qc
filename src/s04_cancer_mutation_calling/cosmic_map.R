library(dplyr)
library(patchwork)
library(ggplot2)
library(data.table)
setwd("/path/to/wd/")

#loop over all files in directory containing "variant_filtered.vcf.gz" but not containing the string "tbi"
#extract COSMIC database
untar("/gpfs/commons/groups/singh_lab/users/kjakubiak/COSMIC_DB/Cosmic_CancerGeneCensus_Tsv_v101_GRCh37.tar", exdir = "COSMIC_DB")
#Directory="/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ekaryo/"
Organism="Human"
#print("Editing VCF File")
#Directory  
Dir="/directory/to/vcf/files/"
file = "FILENAME.variant_filtered.vcf.gz"
path = paste(Dir, file, sep="")
readData = read.delim(VCF_FILE,as.is=T)
readData=as.character(readData[-c(1:which(readData=="#CHROM")-1),1])

jump = 10
startChr = 1+jump
startPos = 2+jump
startInfo = 10+jump

len = length(readData)
chrRegex ="^chr(\\w+)$"
infoRegex = "^([01])\\/([01]):(\\d+)\\,(\\d+):(\\d+):\\d+:\\d+\\,\\d+\\,\\d+$"

chrVector =  readData[startChr]
posVector =  readData[startPos]
infoVector = readData[startInfo]

while (startInfo + jump < len) {
  startChr = startChr + jump
  startPos = startPos + jump
  startInfo = startInfo + jump
  chrVector = append(chrVector, readData[startChr])
  posVector = append(posVector, readData[startPos])
  infoVector = append(infoVector, readData[startInfo])
}

chrNum = gsub(chrRegex, "\\1", chrVector)
if (Organism=="Human"){
  chrNum[chrNum=="X"]="23"
  chrNum[chrNum=="Y"]="24"}

if (Organism=="Mouse"){
  chrNum[chrNum=="X"]="20"
  chrNum[chrNum=="Y"]="21"}

chrNum =  as.numeric(chrNum)
Karyotape = 10*abs(as.numeric(gsub(infoRegex, "\\1", infoVector))-as.numeric(gsub(infoRegex, "\\2", infoVector)))
AD1 = as.numeric(gsub(infoRegex, "\\3", infoVector))
AD2 = as.numeric(gsub(infoRegex, "\\4", infoVector))
DP = as.numeric(gsub(infoRegex, "\\5", infoVector))

posVector = as.numeric(posVector)

table = data.frame("chr" = chrNum, "position" = posVector, "AD1" = AD1, "AD2" = AD2, "DP" = DP, "Karyotape" = Karyotape)

table$chr=as.numeric(table$chr)
table=table[order(table$chr,table$position),]
table=table[table$chr>0,]
colnames(table)<-c("CHROMOSOME","POS","AD1","AD2","DP","KARYOTYPE")
cosmic = fread("/gpfs/commons/groups/singh_lab/users/kjakubiak/COSMIC_DB/Cosmic_MutantCensus_v101_GRCh38.tsv.gz")
#db<-cosmic[,-c(2:8,19:25)]
colnames(table)<-c("CHROMOSOME","POS","AD1","AD2","DP","Karyotype")

genes<-read.table("/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/gene_list.txt")
db<-cosmic_tier1_10tumors_confirmed
#db<-db[db$GENE_SYMBOL%in%genes$V1]
#db<-subset(db, GENE_SYMBOL%in%genes$V1)
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
table1<-vcf
# Split df1 and df2 by chromosome
df1_split <- split(table1, table$CHROMOSOME)
df2_split <- split(db, db$CHROMOSOME)

# Iterate over chromosomes and find matches, keeping all columns
final_results <- lapply(names(df1_split), function(chr) {
  if (chr %in% names(df2_split)) {
    match_genes_full(df1_split[[chr]], df2_split[[chr]])
  }
})

# Combine results into a final dataframe
final_df <- do.call(rbind, final_results)
summary_genes<-final_df$GENE_SYMBOL %>% unique()

df_plot<-final_df
df_plot$CHROMOSOME=NULL
df_plot <- df_plot %>%
  group_by(GENE_SYMBOL, MUTATION_DESCRIPTION) %>%
  summarise(unique_pos_count = n_distinct(POS), .groups = "drop")
ggplot(df_plot, aes(x = GENE_SYMBOL, y = unique_pos_count, fill = MUTATION_DESCRIPTION)) +
  geom_bar(stat = "identity", position = "dodge", width=0.6) +
  labs(x = "Gene Symbol", y = "Unique Position Count", fill = "Mutation Type") +
  theme_minimal() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("SNP Hits in Cancer Mutation Database Metrics")

library(ggplot2)

ggplot(df_plot, aes(x = GENE_SYMBOL, y = unique_pos_count, fill = MUTATION_DESCRIPTION)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Adjusted bar width for a cleaner look
  #scale_fill_brewer(palette = "Set3") +  # Use a nice color palette
  labs(
    x = "Gene Symbol", 
    y = "Unique Position Count", 
    fill = "Mutation Type",
    title = "SNP Hits in Cancer Mutation Database Metrics",
    subtitle = "Cosmic_MutantCensus_Tsv_v101_GRCh38"
  ) +
  theme_minimal(base_size = 14) +  # Clean minimal theme with larger text
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle for readability
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(size = 10),  # Font size for axis text
    legend.title = element_text(size = 12, face = "bold"),  # Refined legend title
    legend.text = element_text(size = 10),  # Refined legend text
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centered and bold title
    plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic")  # Centered subtitle
  ) 



ggsave(
  filename=paste(namex, "_CancerMutationPlot.pdf", sep = ""),
  plot = plotty,
  device = "pdf",
  path = out_dir,
  width=10,
  height=8)


##########################################################################
############################ --------------- #############################
##########################################################################

setwd("/gpfs/commons/groups/singh_lab/users/kjakubiak")
fread("COSMIC_DB/Cosmic_MutantCensus_v101_GRCh38.tsv.gz")->cosmic
VCF_FILE="/gpfs/commons/groups/singh_lab/users/dongwang/RNA-variant-calling/24246R-06-15.variant_filtered.vcf.gz"
GENE_LIST="/gpfs/commons/groups/singh_lab/users/dongwang/RNA-variant-calling/"
genes[26,]<-"BRCA2"
readData = read.delim(VCF_FILE,as.is=T)
readData=as.character(readData[-c(1:which(readData=="#CHROM")-1),1])
jump = 10
startChr = 1+jump
startPos = 2+jump
startInfo = 10+jump
len = length(readData)
chrRegex ="^chr(\\w+)$"
infoRegex = "^([01])\\/([01]):(\\d+)\\,(\\d+):(\\d+):\\d+:\\d+\\,\\d+\\,\\d+$"
chrVector =  readData[startChr]
posVector =  readData[startPos]
infoVector = readData[startInfo]
while (startInfo + jump < len) {
  startChr = startChr + jump
  startPos = startPos + jump
  startInfo = startInfo + jump
  chrVector = append(chrVector, readData[startChr])
  posVector = append(posVector, readData[startPos])
  infoVector = append(infoVector, readData[startInfo])
}
chrNum = gsub(chrRegex, "\\1", chrVector)

Organism="Human"
if (Organism=="Human"){
  chrNum[chrNum=="X"]="23"
  chrNum[chrNum=="Y"]="24"}
chrNum =  as.numeric(chrNum)

Karyotape = 10*abs(as.numeric(gsub(infoRegex, "\\1", infoVector))-as.numeric(gsub(infoRegex, "\\2", infoVector)))
AD1 = as.numeric(gsub(infoRegex, "\\3", infoVector))
AD2 = as.numeric(gsub(infoRegex, "\\4", infoVector))
DP = as.numeric(gsub(infoRegex, "\\5", infoVector))
posVector = as.numeric(posVector)
table = data.frame("chr" = chrNum, "position" = posVector, "AD1" = AD1, "AD2" = AD2, "DP" = DP, "Karyotape" = Karyotape)
table$chr=as.numeric(table$chr)
table=table[order(table$chr,table$position),]
table=table[table$chr>0,]
db<-cosmic
colnames(table)<-c("CHROMOSOME","POS","AD1","AD2","DP","KARYOTYPE")

# Load required libraries
library(dplyr)
library(ggplot2)
library(data.table)
library(fuzzyjoin)
library(tidyr)
# Step 1: Load Data
#vcf <- fread(VCF_FILE,  header = FALSE)
#colnames(vcf)[1:2] <- c("CHROMOSOME", "POS")  # Assuming first two columns are CHR and POS

cancer_db<-cosmic_tier1_10tumors
cancer_db <- cancer_db %>% select(CHROMOSOME, GENOME_START, GENOME_STOP, GENE_SYMBOL, MUTATION_DESCRIPTION)

# Group by GENE_SYMBOL and CHROMOSOME to summarize the genomic range
cancer_gene_ranges <- cancer_db %>%
  group_by(GENE_SYMBOL, CHROMOSOME) %>%
  summarise(GENOME_START = min(GENOME_START), 
            GENOME_STOP = max(GENOME_STOP), .groups = "drop")

cancer_gene_ranges<-na.omit(cancer_gene_ranges)
vcf$CHROMOSOME %>% typeof
cancer_gene_ranges$CHROMOSOME %>% typeof
cancer_gene_ranges$CHROMOSOME<-as.numeric(cancer_gene_ranges$CHROMOSOME)
# Step 3: Match SNPs to Cancer Gene Ranges (Ensuring Chromosomes Match)
matched_snps <- genome_inner_join(
  vcf, cancer_gene_ranges,
  by = c("CHROMOSOME" = "CHROMOSOME", "POS" = "GENOME_START", "POS" = "GENOME_STOP")
)

# Step 4: Fix Duplicate Chromosome Columns
matched_snps <- matched_snps %>%
  rename(CHROMOSOME = CHROMOSOME.x) %>%  # Keep CHROMOSOME from vcf
  select(-CHROMOSOME.y)  # Remove duplicate column from cancer_gene_ranges

# Step 5: Count Unique Gene Matches Per Chromosome
gene_match_counts <- matched_snps %>%
  group_by(CHROMOSOME) %>%
  summarise(unique_genes_matched = n_distinct(GENE_SYMBOL), .groups = "drop")
gene_match_counts$CHROMOSOME<-as.factor(gene_match_counts$CHROMOSOME)
# Step 6: Plot the Number of Cancer Genes with SNP Matches Per Chromosome
library(viridis)

ggplot(gene_match_counts, aes(x = CHROMOSOME, y = unique_genes_matched, fill = CHROMOSOME)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_viridis_d(option = "turbo") +  # "turbo" is colorful and distinguishable
  theme_minimal(base_size = 15) +
  labs(title = "Number of Cancer Genes with SNP Matches Per Chromosome",
       x = "Chromosome", y = "Unique Genes Matched") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Step 1: Filter VCF file for matching chromosomes
vcf_filtered <- vcf %>%
  filter(CHROMOSOME %in% cancer_db$CHROMOSOME)

# Step 2: For each SNP in the VCF file, check if POS is within the range of GENOME_START and GENOME_STOP
matched_snps <- vcf_filtered %>%
  rowwise() %>%  # Iterate row-by-row
  mutate(
    matched_gene = list(cancer_db %>%
                          filter(
                            CHROMOSOME == CHROMOSOME & POS >= GENOME_START & POS <= GENOME_STOP  # Check for positional overlap
                          ))
  ) %>%
  unnest(matched_gene)  # Expand matched genes into separate rows

# Step 3: Group by CHROMOSOME, GENE_SYMBOL, and MUTATION_DESCRIPTION, and count unique gene symbols
gene_match_counts <- matched_snps %>%
  group_by(CHROMOSOME, GENE_SYMBOL) %>%
  summarise(unique_genes_matched = n_distinct(GENE_SYMBOL), .groups = "drop")

# Step 4: Plot the counts of unique gene matches per chromosome and mutation type
ggplot(gene_match_counts, aes(x = factor(CHROMOSOME), y = unique_genes_matched, fill = factor(CHROMOSOME))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Separate bars by mutation type
  scale_fill_viridis_d(option = "turbo") +  # Color palette for mutation types
  theme_minimal(base_size = 15) +
  labs(title = "Number of Unique Genes Per Mutation Type on Each Chromosome",
       x = "Chromosome", y = "Unique Genes Matched",
       fill = "Mutation Type") +  # Label for the mutation type legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability



# # Load required libraries
# library(dplyr)
# library(ggplot2)
# library(data.table)  # For fast file reading
# library(fuzzyjoin)   # To handle genomic position-based joins
# 
# # Step 1: Load Data
# # Read VCF File (assuming tab-separated without headers)
vcf2 <- fread(VCF_FILE, header = FALSE)
# colnames(vcf)[1:2] <- c("CHROMOSOME", "POS")  # Assuming first two columns are CHR and POS
# 
# # Read Cancer Mutation Database (TSV)
# cancer_db <- fread("your_cancer_db.tsv")
# 
# # Ensure column names match
# colnames(cancer_db) <- toupper(colnames(cancer_db))  # Normalize column names
# cancer_db <- cancer_db %>% select(CHROMOSOME, GENOME_START, GENOME_STOP, GENE_SYMBOL, MUTATION_DESCRIPTION)
# 
# # Convert Chromosome Names to Character (for consistency)
# vcf$CHROMOSOME <- as.character(vcf$CHROMOSOME)
# cancer_db$CHROMOSOME <- as.character(cancer_db$CHROMOSOME)
# 
# # Step 2: Perform Position-Based Join Efficiently
# matched_snps <- genome_inner_join(
#   vcf, cancer_db,
#   by = c("CHROMOSOME" = "CHROMOSOME", "POS" = "GENOME_START", "POS" = "GENOME_STOP")
# )
# 
# # Step 3: Summarize Matches
# match_counts <- matched_snps %>%
#   group_by(CHROMOSOME.y) %>%
#   summarise(total_matches = n())
# 
# # Step 4: Plot - SNPs Matching Cancer Mutations Per Chromosome
ggplot(gene_match_counts, aes(x = factor(CHROMOSOME), y = unique_genes_matched, fill = CHROMOSOME)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  #scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 15) +
  labs(title = "Total Cancerous SNPs per Chromosome", x = "Chromosome", y = "SNP Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # Step 5: Print total SNP matches
# cat("Total SNPs matching cancer mutations:", nrow(matched_snps), "\n")


