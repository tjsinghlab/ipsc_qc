#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(tidyr)
  library(fuzzyjoin)
  library(patchwork)
})

# -----------------------
# Options
# -----------------------
option_list <- list(
  make_option("--sample", type = "character", default = NULL,
              help = "Sample name (used for labeling outputs)", metavar = "character"),
  make_option("--ref_dir", type = "character", default = "./ref",
              help = "Reference directory (default: ./ref)", metavar = "character"),
  make_option("--vcf", type = "character", default = NULL,
              help = "VCF file for this sample", metavar = "character"),
  make_option("--output_dir", type = "character", default = ".",
              help = "Output directory for results", metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$sample) || is.null(opt$vcf)) {
  stop("Both --sample and --vcf arguments are required.")
}

sample_id <- opt$sample
ref_dir <- opt$ref_dir
vcf_file <- opt$vcf
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
# Load VCF
# -----------------------
message("[INFO][", sample_id, "] Processing VCF: ", vcf_file)

vcf <- fread(vcf_file, data.table = FALSE)
colnames(vcf)[1:2] <- c("CHROMOSOME", "POS")
vcf$CHROMOSOME <- as.numeric(gsub("^chr", "", vcf$CHROMOSOME))

# -----------------------
# Match SNPs to cancer genes
# -----------------------
matched_snps <- genome_inner_join(
  vcf, cancer_db,
  by = c("CHROMOSOME" = "CHROMOSOME", "POS" = "GENOME_START", "POS" = "GENOME_STOP")
)

if (nrow(matched_snps) == 0) {
  message("[INFO][", sample_id, "] No COSMIC matches found.")
  quit(status = 0)
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
       title = paste("SNP Hits in COSMIC -", sample_id)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(gene_match_counts, aes(x = factor(CHROMOSOME), y = unique_pos_count, fill = MUTATION_DESCRIPTION)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(x = "Chromosome", y = "Unique Genes Matched", fill = "Mutation Type",
       title = paste("COSMIC Matches Per Chromosome -", sample_id)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf_out <- file.path(output_dir, paste0(sample_id, "_CancerMutationPlots.pdf"))
ggsave(pdf_out, p1 / p2, width = 12, height = 10)
message("[INFO][", sample_id, "] Plots saved to ", pdf_out)

# -----------------------
# Save matched table
# -----------------------
out_tsv <- file.path(output_dir, paste0(sample_id, "_CancerMutations.tsv"))
fwrite(matched_snps, out_tsv, sep = "\t")
message("[INFO][", sample_id, "] Results saved to ", out_tsv)
