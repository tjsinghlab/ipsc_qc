#!/bin/bash

# Set default values for the variables
echo "Ready..."
echo "Set..."
echo "GO!!"
d=""
v=""
g=""
o=""
r=""

while getopts "d:v:g:o:r:" opt; do
  echo "Processing option: $opt with argument: $OPTARG"
  case ${opt} in
    d)
      db_dir="$OPTARG"
      ;;
    v)
      directory_containing_vcf="$OPTARG"
      ;;
    g)
      txt_file="$OPTARG"
      ;;
    o)
      out_dir="$OPTARG"
      ;;
    r)
      Organism="$OPTARG"
      ;;
    \?)
      echo "Usage: $0 [-d db_dir] [-v directory_containing_vcf] [-g txt_file] [-o out_dir] [-r Organism]"
      exit 1
      ;;
  esac
done


# Validate inputs
if [[ -z "$db_dir" || -z "$directory_containing_vcf" || -z "$txt_file" || -z "$out_dir" || -z "$Organism" ]]; then
  echo "All arguments (-d, -v, -g, -o, -r) must be provided"
  exit 1
fi

echo "Loading module R 4.3.3"
module load R/4.3.3

# Now we will run R code embedded in this script
Rscript - <<EOF

# Load necessary libraries
print("Importing R libraries...")
library(dplyr)
library(patchwork)
library(ggplot2)
library(data.table)
library(stringr)

# Get command line arguments from shell script
print("Passing command line entries into R script...")
args <- commandArgs(trailingOnly = TRUE)
db_dir <- "$db_dir"
print(paste("db_dir: ", db_dir))
directory_containing_vcf <- "$directory_containing_vcf"
print(paste("directory_containing_vcf: ", directory_containing_vcf))
txt_file <- "$txt_file"
print(paste("txt_file: ", txt_file))
out_dir <- "$out_dir"
print(paste("out_dir: ", out_dir))
Organism <- "$Organism"
print(paste("Organism: ", Organism))
print("Reading database...")
cosmic = fread(db_dir)

db<-cosmic#[,-c(2:8,19:25)]


print("Reading genes list...")
genes<-read.table(txt_file)

print("Filtering database for genes in list...")
db<-db[db[["GENE_SYMBOL"]]%in%genes[["V1"]]]
print("Defining function.")
# Function to process each file; comes directly from eSNPKaryotyping package
process_vcf_file <- function(db_dir, Organism, directory_containing_vcf, name_vcf, out_dir) {

  # Edit VCF File
  print("Editing VCF File")
  Dir <- directory_containing_vcf
  file <- name_vcf
  path <- paste(Dir, file, sep="")
  readData <- read.delim(path, as.is = TRUE)
  readData <- as.character(readData[-c(1:which(readData == "#CHROM")-1), 1])
  print(paste("Reading file", name_vcf))
  
  jump <- 10
  startChr <- 1 + jump
  startPos <- 2 + jump
  startInfo <- 10 + jump
  len <- length(readData)
  
  print("Defining chrRegex and infoRegex")
  chrRegex <- "^chr(\\\\w+)$"
  infoRegex <- "^([01])\\\\/([01]):(\\\\d+)\\\\,(\\\\d+):(\\\\d+):\\\\d+:\\\\d+\\\\,\\\\d+\\\\,\\\\d+$"
  
  print("Reading data...")
  chrVector <- readData[startChr]
  posVector <- readData[startPos]
  infoVector <- readData[startInfo]
  
  print(paste0("Processing file ",name_vcf))
  while (startInfo + jump < len) {
    startChr <- startChr + jump
    startPos <- startPos + jump
    startInfo <- startInfo + jump
    chrVector <- append(chrVector, readData[startChr])
    posVector <- append(posVector, readData[startPos])
    infoVector <- append(infoVector, readData[startInfo])
  }
  
  print("Defining chrNum")
  chrNum <- gsub(chrRegex, "\\\\1", chrVector)
  
  print("Defining chr for human...")
  if (Organism == "Human") {
    chrNum[chrNum == "X"] <- "23"
    chrNum[chrNum == "Y"] <- "24"
  }
  
   if (Organism == "Mouse") {
     chrNum[chrNum == "X"] <- "20"
     chrNum[chrNum == "Y"] <- "21"
   }
  
  print("Setting chr as numeric...")
  chrNum <- as.numeric(chrNum)

  print("Defining table columns...")
  Karyotape = 10*abs(as.numeric(gsub(infoRegex, "\\\\1", infoVector))-as.numeric(gsub(infoRegex, "\\\\2", infoVector)))
  AD1 = as.numeric(gsub(infoRegex, "\\\\3", infoVector))
  AD2 = as.numeric(gsub(infoRegex, "\\\\4", infoVector))
  DP = as.numeric(gsub(infoRegex, "\\\\5", infoVector))
  
  posVector <- as.numeric(posVector)

  print("Unlisting.")
  chrNum <- as.numeric(unlist(chrNum))
  posVector <- as.numeric(unlist(posVector))
  AD1 <- as.numeric(unlist(AD1))
  AD2 <- as.numeric(unlist(AD2))
  DP <- as.numeric(unlist(DP))
  Karyotape <- as.numeric(unlist(Karyotape))

  print("Constructing table...")
  table <- data.frame("chr" = chrNum, "position" = posVector, "AD1" = AD1, "AD2" = AD2, "DP" = DP, "Karyotape" = Karyotape)
  
  # Further processing
  table[["chr"]] <- as.numeric(table[["chr"]])
  table <- table[order(table[["chr"]], table[["position"]]), ]
  table <- table[table[["chr"]] > 0, ]
  
  

# Function to match POS with genome regions and retain all columns

match_genes_full <- function(df1_chr, df2_chr) {
  result <- list()  # Store results
  
  for (i in seq_len(nrow(df1_chr))) {
    pos <- df1_chr[["POS"]][i]
    
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

# Split df1 and df2 by chromosome
colnames(table)<-c("CHROMOSOME","POS","AD1","AD2","DP","Karyotype")
print("Splitting by chromosome...")
df1_split <- split(table, table[["CHROMOSOME"]])
df2_split <- split(db, db[["CHROMOSOME"]])

df2_split[[1]]=NULL

print("Finding matching positions...")
# Iterate over chromosomes and find matches, keeping all columns
final_results <- lapply(names(df1_split), function(chr) {
  if (chr %in% names(df2_split)) {
    match_genes_full(df1_split[[chr]], df2_split[[chr]])
  }
})
namex<-str_split_i(name_vcf, "[.]", 1)
# Combine results into a final dataframe
print("Combining results into final df...")
final_df <- do.call(rbind, final_results)
saveRDS(final_df, paste0(out_dir,namex,"_cancer_mutation_screen_results.rds"))
summary_genes<-final_df[["GENE_SYMBOL"]] %>% unique()

df_plot<-final_df
df_plot[["CHROMOSOME"]]=NULL
df_plot <- df_plot %>%
  group_by(GENE_SYMBOL, MUTATION_DESCRIPTION) %>%
  summarise(unique_pos_count = n_distinct(POS), .groups = "drop")

#pdf(paste(out_dir,"/", namex, "_CancerMutationPlot.pdf", sep = ""), width = 10, height = 8)

plotty<-ggplot(df_plot, aes(x = GENE_SYMBOL, y = unique_pos_count, fill = MUTATION_DESCRIPTION)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Adjusted bar width for a cleaner look
  scale_fill_brewer(palette = "Set3") +  # Use a nice color palette
  labs(
    x = "Gene Symbol", 
    y = "Unique Position Count", 
    fill = "Mutation Type",
    title = "SNP Cancer Mutation Hits",
    subtitle = namex
  ) +
  theme_minimal(base_size = 14) +  # Clean minimal theme with larger text
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle for readability
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(size = 12),  # Font size for axis text
    legend.title = element_text(size = 12, face = "bold"),  # Refined legend title
    legend.text = element_text(size = 11),  # Refined legend text
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centered and bold title
    plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic")  # Centered subtitle
  ) 

 ggsave(
  filename=paste(namex, "_CancerMutationPlot.pdf", sep = ""),
  plot = plotty,
  device = "pdf",
  path = out_dir,
  width=12,
  height=8)
  

  #dev.off()  # Close the PDF device
}

# Apply function to all .variant_filtered.vcf.gz files in the directory
files <- list.files(directory_containing_vcf, pattern = "variant_filtered.vcf.gz", full.names = F)
files <- files[!grepl("tbi", files)]
print("Running loop for the following files:")
print(files)
for (file in files) {
  process_vcf_file(db_dir, Organism, directory_containing_vcf, basename(file), out_dir)
}
print("Finished screening all samples for cancerous mutations. Huzzah!")
EOF