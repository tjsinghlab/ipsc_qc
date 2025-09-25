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

# -----------------------
# Options
# -----------------------
cwd <- getwd()

option_list <- list(
  make_option("--output_dir", type="character", default=file.path(cwd, "outputs"),
              help="Path to desired output directory", metavar="character"),
  make_option("--bam", type="character", default=file.path(cwd, "bams"),
              help="Path to BAM file for sample", metavar="character"),
  make_option("--ref_dir", type = "character", default = "/ref",
              help = "Reference directory inside Docker image (default: /ref)"),
  make_option("--vcf", type="character", default=file.path(cwd, "vcfs"),
              help="Path to VCF file for sample", metavar="character"),
  make_option("--sample", type="character", default=NULL,
              help="Sample name", metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

output_dir <- opt$output_dir
output_dir <- file.path(output_dir, "eSNPKaryotyping")
ref_dir <- opt$ref_dir
bam_dir <- opt$bam
vcf_dir <- opt$vcf
sample <- opt$sample

if (is.null(sample)) stop("Error: --sample must be provided!")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Edit dbSNP files (once only)
# -----------------------------
Edit_dbSNP_Files(
  Directory = file.path(ref_dir, "chr/"), 
  File_Name = "chr", 
  Organism = "Human"
)

# -----------------------------
# Edit VCF & generate variant table
# -----------------------------
vcf_file <- file.path(vcf_dir, paste0(sample, ".variant_filtered.vcf.gz"))
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

png(file.path(output_dir, paste0(sample, "_PlotGenome.png")), width = 2000, height = 1200, res = 150)
PlotGenome(table2, Window = 151, Organism = "Human", Ylim = 3, PValue = TRUE)
dev.off()

# -----------------------------
# Run DeletionTable
# -----------------------------
bam_file <- file.path(bam_dir, paste0(sample, ".bam"))
if (!file.exists(bam_file)) stop(paste("BAM file not found:", bam_file))

DeletionTable <- function(Directory, Table, dbSNP_Data_Directory, dbSNP_File_Name, Genome_Fa_dict, Organism, bam_file) {
  message("Reading SNPs table...")
  i=1
  
  if (Organism=="Human"){mx=25} else if (Organism=="Mouse"){mx=22}
  
  snpTable <- NULL
  while (i < mx){
    chrTable <- read.delim(paste0(dbSNP_Data_Directory, dbSNP_File_Name, i))
    snpTable <- if (is.null(snpTable)) chrTable else rbind(snpTable, chrTable)
    i=i+1
  }
  
  x <- merge(snpTable, Table, by = c("chr","position"), all.x = TRUE)
  x <- x[order(x$chr,x$position),]
  
  dir.create(Directory, recursive = TRUE, showWarnings = FALSE)
  setwd(Directory)
  system(paste("samtools index", shQuote(bam_file)))
  
  dict <- read.csv(Genome_Fa_dict, as.is=TRUE)
  dict_type <- ifelse(length(grep("chr", dict[1,1]))==0, 0, 1)
  
  tbl <- NULL
  for (i in 1:(mx-1)) {
    message(paste("Chromosome ", i, " | ", Sys.time()))
    x1 <- x[x$chr==i,]
    if (nrow(x1) == 0) next
    
    if(dict_type==0){
      loc <- paste0("chr", i, ":", x1$position[1], "-", x1$position[nrow(x1)])
      if (i==23) loc <- paste0("chrX:", x1$position[1], "-", x1$position[nrow(x1)])
      if (i==24) loc <- paste0("chrY:", x1$position[1], "-", x1$position[nrow(x1)])
    } else {
      loc <- paste0(i, ":", x1$position[1], "-", x1$position[nrow(x1)])
    }
    
    command <- paste("samtools depth -r", shQuote(loc), bam_file, "> reads-per-position.txt")
    system(command)
    system("awk -F \" \" '($3 >20){print $0}' reads-per-position.txt > reads-per-position2.txt")
    
    chr <- read.delim("reads-per-position2.txt", header = FALSE)
    colnames(chr) <- c("chr","position","Depth")
    x2 <- merge(x1, chr, by="position")
    x2$Depth_group <- cut(x2$Depth,
                          breaks = c(20,50,100,200,500,Inf),
                          labels = c("20-50","50-100","100-200","200-500",">500"),
                          include.lowest = TRUE)
    tbl <- if (is.null(tbl)) x2 else rbind(tbl,x2)
  }
  
  tbl[is.na(tbl)] <- 0
  write.table(tbl, "Deletions.txt", sep="\t", row.names=FALSE, quote=FALSE)
  return(tbl)
}

tbl <- DeletionTable(
  Directory = file.path(output_dir, "tmp"),
  Table = table2,
  dbSNP_Data_Directory = file.path(ref_dir, "chr"),
  dbSNP_File_Name = "Edited_Common_chr",
  Genome_Fa_dict = file.path(ref_dir, "ncbi_dataset/data/GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna"),
  Organism = "Human",
  bam_file = bam_file
)

png(file.path(output_dir, paste0(sample, "_Zygosity_Single.png")), width = 2000, height = 1200, res = 150)
Plot_Zygosity_Sinle(Table = tbl, Organism = "Human")
dev.off()

png(file.path(output_dir, paste0(sample, "_Zygosity_Blocks.png")), width = 2000, height = 1200, res = 150)
Plot_Zygosity_Blocks(Table = tbl, Window = 1500000, Max = 6, Max2 = 60, Organism = "Human")
dev.off()
