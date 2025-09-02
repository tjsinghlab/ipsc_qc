#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(rmarkdown)
  library(knitr)
  library(ggplot2)
})

# -------------------------------
# Command-line arguments (named)
# -------------------------------
option_list <- list(
  make_option(c("-p", "--project"), type="character", default="default_project",
              help="Project name", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL,
              help="Path to output directory", metavar="character")
)

opt <- parse_args(OptionParser(option_list = option_list))
project_name <- opt$project
output_dir <- ifelse(is.null(opt$output_dir), file.path(getwd(), "outputs"), opt$output_dir)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# Define paths to previously generated plots
# Adjust filenames as needed
# -------------------------------
plots <- list(
  PCA_counts      = file.path(output_dir, paste0(project_name, "_PCA_counts.pdf")),
  PCA_scores      = file.path(output_dir, paste0(project_name, "_PCA_pacnet_scores.pdf")),
  Mycoplasma      = file.path(output_dir, "mycoplasma", "mycoplasma_alignment_summary.pdf")
  # Add more plots here as needed
)

# -------------------------------
# Generate a single RMarkdown file
# -------------------------------
report_rmd <- file.path(output_dir, paste0(project_name, "_summary_report.Rmd"))

cat(file = report_rmd, "---
title: 'Project Summary Report'
output: 
  html_document:
    toc: true
    toc_depth: 2
    theme: united
    highlight: tango
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE)
library(ggplot2)
```

# Project: `r project_name`

## PCA Analyses

### Gene Expression PCA
```{r}
if (file.exists('", plots$PCA_counts, "')) {
  knitr::include_graphics('", plots$PCA_counts, "')
} else {
  cat('PCA counts plot not found.')
}
```

### PACNet ESC Scores PCA
```{r}
if (file.exists('", plots$PCA_scores, "')) {
  knitr::include_graphics('", plots$PCA_scores, "')
} else {
  cat('PCA PACNet scores plot not found.')
}
```

## Mycoplasma Contamination
```{r}
if (file.exists('", plots$Mycoplasma, "')) {
  knitr::include_graphics('", plots$Mycoplasma, "')
} else {
  cat('Mycoplasma plot not found.')
}
```

## Notes
- PCA plots highlight outliers using Mahalanobis distance.
- Mycoplasma plot shows % reads aligned per sample and correlation with sequencing depth.
- Additional plots from pipeline modules can be added to this summary by appending them to the plots list.
")

# -------------------------------
# Render HTML report
# -------------------------------
rmarkdown::render(
  report_rmd,
  output_file = paste0(project_name, "_summary_report.html"),
  output_dir = output_dir
)

message("[INFO] Summary report created at: ", file.path(output_dir, paste0(project_name, "_summary_report.html")))