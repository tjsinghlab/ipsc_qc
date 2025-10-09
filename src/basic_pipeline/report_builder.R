#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(magick)
  library(pdftools)
  library(grid)
  library(fs)
})

# ----------------------------
# Parse arguments
# ----------------------------
option_list <- list(
  make_option(c("--output_dir"), type = "character", help = "Main output directory containing per-sample folders"),
  make_option(c("--project"), type = "character", help = "Project name")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.na(opt$output_dir) || is.na(opt$project)) {
  stop("Usage: Rscript report_builder.R --output_dir <dir> --project <name>")
}

output_dir <- normalizePath(opt$output_dir, mustWork = TRUE)
project <- opt$project
pdf_out <- file.path(output_dir, "final_summary.pdf")

cat("[INFO] Building project summary for:", project, "\n")

# ----------------------------
# Helper functions
# ----------------------------

# Create a header image with project or sample name
header_img <- function(text, width = 8, height = 1) {
  tf <- tempfile(fileext = ".png")
  png(tf, width = width, height = height, units = "in", res = 150)
  grid.newpage()
  grid.text(text, gp = gpar(fontsize = 20, fontface = "bold"))
  dev.off()
  image_read(tf)
}

# Combine multiple images nicely on one page
arrange_images_on_page <- function(imgs, ncol = 2) {
  # resize all to the same width
  imgs <- lapply(imgs, function(im) image_scale(im, "1000x"))
  # tile them in grid
  page <- image_montage(do.call(c, imgs),
                        tile = paste0(ncol, "x"),
                        geometry = "+10+10")
  return(page)
}

# Build a sample page with header + plots
make_sample_page <- function(sample_name, png_files) {
  imgs <- lapply(png_files, image_read)
  head <- header_img(paste("Sample:", sample_name))
  body <- arrange_images_on_page(imgs)
  image_append(c(head, body), stack = TRUE)
}

# ----------------------------
# Gather samples
# ----------------------------
sample_dirs <- dir_ls(output_dir, type = "directory", recurse = FALSE)
sample_dirs <- sample_dirs[!basename(sample_dirs) %in% c("logs", "tmp", "qc", "figures")]

if (length(sample_dirs) == 0) {
  stop("[ERROR] No sample directories found in ", output_dir)
}

cat("[INFO] Found", length(sample_dirs), "samples\n")

# ----------------------------
# Make a project title page
# ----------------------------
cat("[INFO] Generating project title page...\n")
title_img <- header_img(paste("Project Summary:", project), height = 2)
tf_title <- tempfile(fileext = ".pdf")
image_write(title_img, tf_title, format = "pdf")

# ----------------------------
# Make per-sample pages
# ----------------------------
pdf_pages <- c(tf_title)

for (sample_dir in sample_dirs) {
  sample <- basename(sample_dir)
  cat("[INFO] Processing sample:", sample, "\n")

  png_files <- c(
    list.files(file.path(sample_dir, "cosmic_calling"), pattern = "\\.png$", full.names = TRUE),
    list.files(file.path(sample_dir, "esnp_karyotyping"), pattern = "\\.png$", full.names = TRUE),
    file.path(sample_dir, "mycoplasma", "mycoplasma_alignment_summary.png")
  )
  png_files <- png_files[file.exists(png_files)]

  if (length(png_files) == 0) {
    cat("[WARN] No plots found for sample:", sample, "\n")
    next
  }

  sample_page <- make_sample_page(sample, png_files)
  tf <- tempfile(fileext = ".pdf")
  image_write(sample_page, tf, format = "pdf")
  pdf_pages <- c(pdf_pages, tf)
}

# ----------------------------
# Combine into final PDF
# ----------------------------
if (length(pdf_pages) > 0) {
  pdf_combine(pdf_pages, output = pdf_out)
  cat("[SUCCESS] Final report written to:", pdf_out, "\n")
} else {
  cat("[WARN] No pages generated.\n")
}
