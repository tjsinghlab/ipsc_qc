#!/usr/bin/env Rscript

library(pdftools)

message("Setting directories...")
# Set directories
output_dir <- "my_outputs"
summary_pdf <- file.path(output_dir, "Analysis_Summary.pdf")

message("Creating empty list to store PDF pages...")
# Create an empty list to store PDF pages
pdf_pages <- list()
library(magick)
# ----------------------------
# First page: date/time
# ----------------------------
message("Writing date/time page...")
date_page <- image_blank(width = 800, height = 600, color = "white") %>%
  image_annotate(text = paste("Summary PDF\nCreated on:", Sys.time()),
                 size = 30, gravity = "center", color = "black")
tmp_file <- tempfile(fileext = ".pdf")
image_write(date_page, tmp_file)
pdf_pages <- c(pdf_pages, tmp_file)

# ----------------------------
# All-sample analyses
# ----------------------------
message("Finding PACNet_heatmap.png...")
all_sample_files <- c(
  file.path(output_dir, "pacnet", "PACNet_heatmap.png")
)

message("Looking for outlier analysis plot...")
# Include outlier_analysis if present
outlier_dir <- file.path(output_dir, "outlier_analysis")
if (dir.exists(outlier_dir)) {
  outlier_files <- list.files(outlier_dir, pattern = "\\.png$", full.names = TRUE)
  all_sample_files <- c(all_sample_files, outlier_files)
}

message("Converting all-sample images to PDF pages...")
# Convert each image to PDF page
for (f in all_sample_files) {
  img <- image_read(f)
  tmp_file <- tempfile(fileext = ".pdf")
  image_write(img, tmp_file)
  pdf_pages <- c(pdf_pages, tmp_file)
}

# ----------------------------
# Per-sample pages
# ----------------------------
message("Finding per-sample directories...")
sample_dirs <- list.dirs(output_dir, recursive = FALSE)
for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  
  # Skip non-sample dirs like pacnet and outlier_analysis
  if (sample_name %in% c("pacnet", "outlier_analysis")) next
  
  # Header page for sample
  header_page <- image_blank(width = 800, height = 200, color = "white") %>%
    image_annotate(text = sample_name, size = 40, gravity = "center", color = "black")
  tmp_file <- tempfile(fileext = ".pdf")
  image_write(header_page, tmp_file)
  pdf_pages <- c(pdf_pages, tmp_file)
  
  # CancerMutationPlots.pdf
  cancer_plot <- list.files(file.path(sample_dir, "cosmic_calling"), 
                           pattern = "\\.png$", full.names = TRUE)
  if (length(cancer_plot) > 0) pdf_pages <- c(pdf_pages, cancer_plot)
  
  # eSNPKaryotyping plots
  esnpk_files <- list.files(file.path(sample_dir, "eSNPKaryotyping"), 
                            pattern = "_PlotGenome\\.png$", full.names = TRUE)
  for (f in esnpk_files) {
    img <- image_read(f)
    tmp_file <- tempfile(fileext = ".pdf")
    image_write(img, tmp_file)
    pdf_pages <- c(pdf_pages, tmp_file)
  }
  
  # mycoplasma alignment
  myco_file <- file.path(sample_dir, "mycoplasma", "mycoplasma_alignment_summary.png")
  if (file.exists(myco_file)) {
    img <- image_read(myco_file)
    tmp_file <- tempfile(fileext = ".png")
    image_write(img, tmp_file)
    pdf_pages <- c(pdf_pages, tmp_file)
  }
}

# ----------------------------
# Combine all pages into a single PDF
# ----------------------------
message("Combining all pages into summary PDF...")
pdf_combine(input = pdf_pages, output = summary_pdf)
cat("Summary PDF created at:", summary_pdf, "\n")
