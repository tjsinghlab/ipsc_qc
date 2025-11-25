#!/usr/bin/env bash

suppressPackageStartupMessages({
  library(optparse)
  library(fs)
  library(tidyverse)
  library(rmarkdown)
})

# -------------------------------
# Command-line options
# -------------------------------
option_list <- list(
  make_option(c("-o", "--output_dir"), type="character", default=".",
              help="Directory containing sample outputs [default: %default]"),
  make_option(c("-p", "--project"), type="character", default="MyProject",
              help="Project name for report [default: %default]"),
  make_option(c("-s", "--samples"), type="character", default=NULL,
              help="Comma-separated list of sample names to include [required]"),
  make_option(c("--html"), action="store_true", default=TRUE,
              help="Render HTML report [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

output_dir <- normalizePath(opt$output_dir)
project <- opt$project

# -------------------------------
# Parse samples
# -------------------------------
if (is.null(opt$samples)) {
  stop("Please provide samples via --samples argument (comma-separated).")
}

samples_to_include <- str_split(opt$samples, ",")[[1]] %>% str_trim()

# -------------------------------
# Build temporary Rmd
# -------------------------------
tmp_rmd <- tempfile(fileext = ".Rmd")

cat(
paste0(
'---
title: "iPSC QC Summary Report"
output:
  html_document:
    toc: true
    toc_depth: 2
    theme: flatly
    df_print: paged
date: "', format(Sys.Date()), '"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, results="asis", fig.align="center")
library(fs)
library(tidyverse)

cat("## Project Summary: ', project, ' \n")
```

'
), file = tmp_rmd
)

# -------------------------------
# PACNet + Final heatmap block
# -------------------------------
pacnet_png <- file.path(output_dir, "pacnet", "PACNet_heatmap.png")
final_pdf  <- file.path(output_dir, "final_summary_heatmap.pdf")

cat(
  paste0(
'```{r pacnet, echo=FALSE, results="asis"}\n',
'cat("\\n\\n# PACNet Overview\\n\\n")\n',
'if (file.exists("', pacnet_png, '")) {\n',
'  cat("## PACNet Heatmap\\n\\n")\n',
'  cat("![](", normalizePath("', pacnet_png, '"), ")\\n\\n", sep="")\n',
'} else {\n',
'  cat("*(PACNet heatmap not found.)*\\n")\n',
'}\n',
'\n',
'cat("\\n\\n# Final Summary Heatmap\\n\\n")\n',
'if (file.exists("', final_pdf, '")) {\n',
'  cat("## Combined Heatmap Overview\\n\\n")\n',
'  cat("<iframe src=\'", normalizePath("', final_pdf, '"), "\' width=\'100%\' height=\'800px\'></iframe>\\n", sep="")\n',
'} else {\n',
'  cat("*(File not found.)*\\n")\n',
'}\n',
'```\n'
  ),
  file = tmp_rmd,
  append = TRUE
)

# -------------------------------
# Loop over samples
# -------------------------------
for (sample in samples_to_include) {

  cat(
    paste0(
'```{r sample_', sample, ', echo=FALSE, results="asis"}\n',
'cat("\\n\\n# Sample: ', sample, '\\n\\n")\n',
'pngs <- c(\n',
'  dir_ls(file.path("', output_dir, '", "', sample, '", "cosmic_calling"), glob="*.png", fail=FALSE),\n',
'  dir_ls(file.path("', output_dir, '", "', sample, '", "eSNPKaryotyping"), glob="*.png", fail=FALSE),\n',
'  file.path("', output_dir, '", "', sample, '", "mycoplasma", "mycoplasma_alignment_summary.png")\n',
')\n',
'pngs <- pngs[file.exists(pngs)]\n',
'if (length(pngs) == 0) {\n',
'  cat("*(No plots found for this sample.)*\\\\n")\n',
'} else {\n',
'  for (img in pngs) {\n',
'    cat("\\\\n\\\\n## ", basename(img), "\\\\n\\\\n", sep="")\n',
'    cat("![](", normalizePath(img), ")\\\\n\\\\n", sep="")\n',
'  }\n',
'}\n',
'```\n'
    ),
    file = tmp_rmd,
    append = TRUE
  )
}

# -------------------------------
# Render HTML
# -------------------------------
out_file <- file.path(output_dir, "Pipeline_Summary_Report.html")
rmarkdown::render(tmp_rmd, output_file = out_file, quiet = FALSE)

cat("\nReport written to: ", out_file, "\n")