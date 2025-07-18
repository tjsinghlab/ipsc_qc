#!/bin/bash
#Usage: esnpkaryo_loop.sh -d /path/to/dbsnp -o Human -p /path/to/directory_containing_vcf -n name.vcf.gz -c /path/to/directory_for_csv
# Set default values for the variables
echo "Ready..."
echo "Set..."
echo "GO!!"
dbsnp_dir=""
Organism=""
directory_containing_vcf=""
directory_for_csv=""

# Parse command-line arguments
while getopts "d:o:p:c:" opt; do
  case ${opt} in
    d)
      dbsnp_dir="$OPTARG"
      ;;
    o)
      Organism="$OPTARG"
      ;;
    p)
      directory_containing_vcf="$OPTARG"
      ;;
    c)
      directory_for_csv="$OPTARG"
      ;;
    \?)
      echo "Usage: $0 [-d dbsnp_dir] [-o Organism] [-p directory_containing_vcf] [-c directory_for_csv]"
      exit 1
      ;;
  esac
done

# Validate inputs
if [[ -z "$dbsnp_dir" || -z "$Organism" || -z "$directory_containing_vcf" || -z "$directory_for_csv" ]]; then
  echo "All arguments (-d, -o, -p, -c) must be provided"
  exit 1
fi

echo "Loading module R 4.3.3"
module load R/4.3.3

# Now we will run R code embedded in this script
Rscript - <<EOF

# Load necessary libraries
print("Importing R libraries...")
library(devtools)
library(eSNPKaryotyping)
library(zoo)
library(gplots)
library(patchwork)
library(ggplot2)
library(stringr)

# Get command line arguments from shell script
print("Passing command line entries into R script...")
args <- commandArgs(trailingOnly = TRUE)
dbsnp_dir <- "$dbsnp_dir"
print(paste("dbsnp_dir: ", dbsnp_dir))
Organism <- "$Organism"
print(paste("Organism: ", Organism))
directory_containing_vcf <- "$directory_containing_vcf"
print(paste("directory_containing_vcf: ", directory_containing_vcf))
directory_for_csv <- "$directory_for_csv"
print(paste("directory_for_csv: ", directory_for_csv))

# Edit dbSNP files (this only occurs once)
print("Editing dbSNP files...")
#Edit_dbSNP_Files(Directory = dbsnp_dir, File_Name = "chr", Organism = Organism)

print("Defining function.")
# Function to process each file; comes directly from eSNPKaryotyping package
process_vcf_file <- function(dbsnp_dir, Organism, directory_containing_vcf, name_vcf, directory_for_csv) {

  # Edit VCF File
  print("Editing VCF File")
  Dir <- directory_containing_vcf
  file <- name_vcf
  path <- paste(Dir, file, sep="")
  readData <- read.delim(path, as.is = TRUE)
  readData <- as.character(readData[-c(1:which(readData == "#CHROM")-1), 1])
  print(paste("ReadVCF file for ", name_vcf))
  
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
  
  print("Assessing VCF...")
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
  
  # Save table as CSV
  fileName = paste0(name_vcf, "_variantTable.csv")
  pathToSave <- paste(directory_for_csv, fileName, sep="")
  table[is.na(table)] <- 0
  write.table(table, pathToSave, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Further processing
  table[["chr"]] <- as.numeric(table[["chr"]])
  table <- table[order(table[["chr"]], table[["position"]]), ]
  table <- table[table[["chr"]] > 0, ]
  
  # Run the MajorMinorCalc function
  print("Running MajorMinorCalc...")
  table2 <- MajorMinorCalc(Table = table, minDP = 20, maxDP = 100000, minAF = 0.2)
  saveRDS(table2, paste0(directory_for_csv, name_vcf, "major_minor_table.rds"))
  # Save the plot to a PDF file in the specified directory
  namex<-str_split_i(name_vcf, "[.]", 1)
  pdf(paste(directory_for_csv,"/", namex, "_plotGenome_output.pdf", sep = ""), width = 10, height = 8)
  PlotGenome(table2, Window = 151, Organism = Organism, Ylim = 3, PValue = TRUE)
  title(main = paste("Allelic Frequency for ", namex))
  dev.off()  # Close the PDF device
}

# Apply function to all .variant_filtered.vcf.gz files in the directory
files <- list.files(directory_containing_vcf, pattern = "variant_filtered.vcf.gz", full.names = F)
files <- files[!grepl("tbi", files)]
print("Running loop for the following files:")
print(files)
for (file in files) {
  process_vcf_file(dbsnp_dir, Organism, directory_containing_vcf, basename(file), directory_for_csv)
  #saveRDS(t2, paste0(file, "major_minor_table.rds))
}
print("Finished running eSNPKaryotyping for all samples. Hooray!")
EOF
