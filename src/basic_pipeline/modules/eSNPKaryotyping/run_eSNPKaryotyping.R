#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(devtools)
  library(eSNPKaryotyping)
  library(zoo)
  library(gplots)
  library(patchwork)
  library(ggplot2)
  library(optparse)
})

#Define variables from arguments
cwd <- getwd()

option_list <- list(
  make_option("--output_dir", type="character", default=file.path(cwd, "outputs"),
              help="Path to desired output directory", metavar="character"),
  make_option("--ref_dir", type = "character", default = "/ref",
              help = "Reference directory inside Docker image (default: /ref)"),
  make_option("--sample", type="character", default=NULL,
              help="Sample name", metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

output_dir <- opt$output_dir
output_dir <- file.path(output_dir, "eSNPKaryotyping")
ref_dir <- opt$ref_dir
sample <- opt$sample

if (is.null(sample)) stop("Error: --sample must be provided!")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
bam_dir <- file.path(opt$output_dir, "Mark_duplicates_outputs")
# -----------------------------
# Edit dbSNP files (only needs to happen once)
# If you've run this before, you may comment out this chunk to save time.
# -----------------------------
Edit_dbSNP_Files(
  Directory = file.path(ref_dir, "chr/"), 
  File_Name = "chr", 
  Organism = "Human"
)

# -----------------------------
# Edit VCF & generate variant table
# -----------------------------
vcf_file <- file.path(opt$output_dir, "variant_calling", paste0(sample, ".variant_filtered.vcf.gz"))
if (!file.exists(vcf_file)) stop(paste("VCF file not found:", vcf_file))

readData <- read.delim(vcf_file, as.is = TRUE)
header_index <- which(readData == "#CHROM")
readData <- as.character(readData[-c(1:(header_index - 1)), 1])

jump <- 10
startChr <- 1 + jump
startPos <- 2 + jump
startInfo <- 10 + jump
len <- length(readData)

chrRegex <- "^chr(\\w+)$"
infoRegex <- "^([01])/([01]):(\\d+),(\\d+):(\\d+):\\d+:\\d+,(\\d+),(\\d+)$"

chrVector <- readData[startChr]
posVector <- readData[startPos]
infoVector <- readData[startInfo]

while (startInfo + jump < len) {
  startChr <- startChr + jump
  startPos <- startPos + jump
  startInfo <- startInfo + jump
  
  chrVector <- append(chrVector, readData[startChr])
  posVector <- append(posVector, readData[startPos])
  infoVector <- append(infoVector, readData[startInfo])
}

chrNum <- gsub(chrRegex, "\\1", chrVector)
chrNum[chrNum == "X"] <- "23"
chrNum[chrNum == "Y"] <- "24"
chrNum <- as.numeric(chrNum)

Karyotape <- 10 * abs(as.numeric(gsub(infoRegex, "\\1", infoVector)) - as.numeric(gsub(infoRegex, "\\2", infoVector)))
AD1 <- as.numeric(gsub(infoRegex, "\\3", infoVector))
AD2 <- as.numeric(gsub(infoRegex, "\\4", infoVector))
DP <- as.numeric(gsub(infoRegex, "\\5", infoVector))
posVector <- as.numeric(posVector)

variant_table <- data.frame(
  "chr" = chrNum,
  "position" = posVector,
  "AD1" = AD1,
  "AD2" = AD2,
  "DP" = DP,
  "Karyotape" = Karyotape
)

# Save variant table
variant_file <- file.path(output_dir, paste0(sample, "_variantTable.csv"))
variant_table[is.na(variant_table)] <- 0
write.table(variant_table, variant_file, sep = "\t", row.names = FALSE, quote = FALSE)

# -----------------------------
# Run MajorMinorCalc and plot
# -----------------------------
variant_table$chr <- as.numeric(variant_table$chr)
variant_table <- variant_table[order(variant_table$chr, variant_table$position), ]
variant_table <- variant_table[variant_table$chr > 0, ]

table2 <- MajorMinorCalc(Table = variant_table, minDP = 20, maxDP = 1000000, minAF = 0.2)

#export results
png(file.path(output_dir, paste0(sample, "_PlotGenome.png")), width = 2000, height = 1200, res = 150)
PlotGenome(table2, Window = 151, Organism = "Human", Ylim = 3, PValue = TRUE)
dev.off()
