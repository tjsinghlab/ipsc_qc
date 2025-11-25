#!/usr/bin/env Rscript

# ===========================================
# Construct final summary table (generalized genes)
# ===========================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

# -------------------------------------------
# Command-line options
# -------------------------------------------
option_list <- list(
  make_option(c("--output_dir"), type="character", help="Output directory containing per-sample results"),
  make_option(c("--pacnet_file"), type="character", help="Path to PACNet classification_scores.csv"),
  make_option(c("--samples"), type="character", help="Comma-separated list of sample names"),
  make_option(c("--project"), type="character", default="Project", help="Project name (optional)")
)

opt <- parse_args(OptionParser(option_list=option_list))
OUTPUT_DIR <- opt$output_dir
PACNET_FILE <- opt$pacnet_file
SAMPLES <- strsplit(opt$samples, ",")[[1]]
#SAMPLES<-c("24246R-06-01_S0_L001", "24246R-06-08_S415_L008", "24246R-06-11_S198_L008")

# -------------------------------------------
# Discover all genes across all samples
# -------------------------------------------
all_genes <- character()
for (sample in SAMPLES) {
  cosmic_file <- file.path(OUTPUT_DIR, sample, "cosmic_calling", paste0(sample, "_CancerMutations_summary.tsv"))
  if (file.exists(cosmic_file)) {
    df <- fread(cosmic_file)
    if (ncol(df) >= 3) {
      setnames(df, c("GENE_SYMBOL", "unique_pos_count", "MUTATION_DESCRIPTION"))
      all_genes <- union(all_genes, unique(df$GENE_SYMBOL))
    }
  }
}

# -------------------------------------------
# Prepare summary table
# -------------------------------------------
summary_list <- list()

for (sample in SAMPLES) {
  sample_dir <- file.path(OUTPUT_DIR, sample)
  
  # --------------------------
  # Mycoplasma %
  # --------------------------
  myco_file <- file.path(sample_dir, "mycoplasma/mycoplasma_alignment_stats.tsv")
  myco_pct <- NA
  if (file.exists(myco_file)) {
    df <- fread(myco_file, skip=1)
    if (nrow(df) > 0) {
      myco_pct <- sum(as.numeric(df[[5]]))
    }
  }
  
  # --------------------------
  # Fraction abnormal karyotype
  # --------------------------
  karyo_file <- file.path(sample_dir, "eSNPKaryotyping", paste0(sample, "_variantTable.csv"))
  frac_karyo <- NA
  if (file.exists(karyo_file)) {
    df <- fread(karyo_file)
    total <- nrow(df)
    abnormal <- sum(df[[6]] != 10)
    if (total > 0) frac_karyo <- round(abnormal/total, 6)
  }
  
  # --------------------------
  # PACNet esc score
  # --------------------------
  pacnet_scores <- fread(PACNET_FILE)
  esc_scores <- subset(pacnet_scores, V1=="esc")

  print(names(esc_scores))
  print(sample)

  esc_score <- NA
  col <- grep(sample, names(esc_scores), value = TRUE)

  if (length(col) == 1) {
      esc_score <- esc_scores[[col]]
  } else if (length(col) > 1) {
      warning(paste("Multiple PACNet columns matched sample:", sample))
      esc_score <- esc_scores[[col[1]]]
  } else {
      warning(paste("No PACNet ESC score column matched sample:", sample))
      esc_score <- NA
  }
  
  # --------------------------
  # COSMIC mutations (dynamic genes)
  # --------------------------
  cosmic_file <- file.path(sample_dir, "cosmic_calling", paste0(sample, "_CancerMutations_summary.tsv"))
  
  frame_counts <- setNames(rep(0, length(all_genes)), all_genes)
  total_counts <- setNames(rep(0, length(all_genes)), all_genes)
  
  if (file.exists(cosmic_file)) {
    df <- fread(cosmic_file)
    if (ncol(df) >= 3) {
      setnames(df, c("GENE_SYMBOL", "MUTATION_DESCRIPTION", "unique_pos_count" ))
      for (gene in all_genes) {
        sub <- df[GENE_SYMBOL == gene]
        total_counts[gene] <- sum(as.numeric(sub$unique_pos_count))
        frame_counts[gene] <- sum(as.numeric(sub$unique_pos_count[grep("frame", sub$MUTATION_DESCRIPTION, ignore.case=TRUE)]))
      }
    }
  }
  
  # --------------------------
  # Construct row
  # --------------------------
  row <- c(
    Sample = sample,
    Total_Mycoplasma_Percent_Alignment = myco_pct,
    Fraction_Abnormal_Karyotype = frac_karyo,
    esc_score = esc_score
  )
  
  # Append dynamic genes
  for (gene in all_genes) {
    row[paste0(gene, "_Frame_Mutations")] <- frame_counts[gene]
    row[paste0(gene, "_Total_Mutations")] <- total_counts[gene]
  }
  
  summary_list[[sample]] <- row
}

# -------------------------------------------
# Write summary
# -------------------------------------------
summary_df <- rbindlist(lapply(summary_list, as.list), fill=TRUE)
output_file <- file.path(OUTPUT_DIR, "final_summary_table.tsv")
fwrite(summary_df, output_file, sep="\t", na="NA", quote=FALSE)

# Construct and save summary heatmap
require(lattice)
df<-summary_df
rnamx<-df$Sample
df$Sample<-NULL
df <- df %>% mutate_if(is.character,as.numeric)

mat<-as.matrix(df)
rownames(mat)<-rnamx

png(
  file.path(OUTPUT_DIR, "final_summary_heatmap.png"),
  width = 1000 * length(SAMPLES),
  height = 5000,
  #units = "in",
  res = 300,
  type = "cairo"
)

levelplot(
  mat,
  scale = list(
    x = list(rot = 45, cex = 1.2),
    y = list(cex = 1.2)
  ),
  par.settings = list(
    fontsize = list(
      text = 14,
      points = 12
    ),
    axis.text = list(cex = 1.2),
    axis.line = list(lwd = 1.2)
  ),
  colorkey = list(
    labels = list(cex = 1.2)
  )
)

dev.off()

message("[INFO] Summary files written to output path.")


