#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(tidyr)
  library(fuzzyjoin)
})

# -----------------------
# Options
# -----------------------
option_list <- list(
  make_option("--ref_dir", type = "character", default = "./ref",
              help = "Reference directory (default: ./ref)", metavar = "character"),
  make_option("--vcf_dir", type = "character", default = ".",
              help = "Directory containing VCF files", metavar = "character"),
  make_option("--output_dir", type = "character", default = ".",
              help = "Output directory", metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))
ref_dir <- opt$ref_dir
vcf_dir <- opt$vcf_dir
output_dir <- opt$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------
# Load COSMIC data
# -----------------------
cosmic_tar <- file.path(ref_dir, "Cosmic_CancerGeneCensus_Tsv_v101_GRCh38.tar")
untar(cosmic_tar, exdir = ref_dir)

cosmic_file <- file.path(ref_dir, "Cosmic_MutantCensus_v101_GRCh38.tsv.gz")
if (!file.exists(cosmic_file)) {
  stop("COSMIC file not found in reference directory.")
}
cosmic <- fread(cosmic_file)

cancer_db <- cosmic %>%
  select(CHROMOSOME, GENOME_START, GENOME_STOP, GENE_SYMBOL, MUTATION_DESCRIPTION) %>%
  na.omit()

# -----------------------
# Process VCFs
# -----------------------
vcf_files <- list.files(vcf_dir, pattern = "variant_filtered.vcf.gz$", full.names = TRUE)
vcf_files <- vcf_files[!grepl(".tbi$", vcf_files)]  # skip index files

if (length(vcf_files) == 0) stop("No VCF files found in ", vcf_dir)

for (vcf_file in vcf_files) {
  sample_id <- sub("\\.variant_filtered\\.vcf\\.gz$", "", basename(vcf_file))
  message("[INFO] Processing sample: ", sample_id)

  # Load VCF (assumes standard tab-delimited format)
  vcf <- fread(vcf_file, data.table = FALSE)
  colnames(vcf)[1:2] <- c("CHROMOSOME", "POS")
  vcf$CHROMOSOME <- as.numeric(gsub("^chr", "", vcf$CHROMOSOME))

  # Match SNPs to cancer genes
  matched_snps <- genome_inner_join(
    vcf, cancer_db,
    by = c("CHROMOSOME" = "CHROMOSOME", "POS" = "GENOME_START", "POS" = "GENOME_STOP")
  )

  if (nrow(matched_snps) == 0) {
    message("[INFO] No COSMIC matches found for ", sample_id)
    next
  }

  # Summarize matches
  gene_match_counts <- matched_snps %>%
    group_by(CHROMOSOME, GENE_SYMBOL, MUTATION_DESCRIPTION) %>%
    summarise(unique_pos_count = n_distinct(POS), .groups = "drop")

  # -----------------------
  # Plots
  # -----------------------
  p1 <- ggplot(gene_match_counts, aes(x = GENE_SYMBOL, y = unique_pos_count, fill = MUTATION_DESCRIPTION)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.6) +
    labs(x = "Gene Symbol", y = "Unique Position Count", fill = "Mutation Type",
         title = "SNP Hits in Cancer Mutation Database Metrics") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p2 <- ggplot(gene_match_counts, aes(x = factor(CHROMOSOME), y = unique_pos_count, fill = MUTATION_DESCRIPTION)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    labs(x = "Chromosome", y = "Unique Genes Matched", fill = "Mutation Type",
         title = "Number of Cancer Genes with SNP Matches Per Chromosome") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  pdf_out <- file.path(output_dir, paste0(sample_id, "_CancerMutationPlots.pdf"))
  ggsave(pdf_out, p1 / p2, width = 12, height = 10)
  message("[INFO] Plots saved to ", pdf_out)

  # Save matched table
  out_tsv <- file.path(output_dir, paste0(sample_id, "_CancerMutations.tsv"))
  fwrite(matched_snps, out_tsv, sep = "\t")
  message("[INFO] Results saved to ", out_tsv)
}
