# ipsc_qc
QC pipeline for iPSC cultures

# Overview

Performing quality control from multiple angles on iPSC lines to ensure viability, reproducibility, and scability in cultures.

# Directory Structure and Organization

```
.
├── README.md
└── src
    ├── s01_bulk_rna_preprocessing
    │   ├── bulk_rna_processing.py
    │   ├── cram_to_fastq.sh
    │   ├── create_variant_tables.sh
    │   ├── plotting.R
    │   └── process_rsem_outputs.R
    ├── s02_germline_calling
    │   └── gatk4_rna_germline_calling_run.py
    ├── s03_pluripotency
    │   ├── pacnet.R
    │   └── process_pacnet_outputs.R
    ├── s04_cancer_mutation_calling
    │   ├── cosmic_map.R
    │   ├── cosmic_mapping.sh
    │   ├── filter_vcf_cosmic.sh
    │   └── screen_by_mutation_tables.R
    ├── s05_demultiplexing_donors
    │   ├── demuxlet_and_reorder_chromosomes.sh
    │   ├── demuxlet.sh
    │   └── souporcell.sh
    ├── s06_mycoplasma_detection
    │   └── mycoplasma_bowtie.sh
    ├── s07_eKaryotyping
    │   ├── delTabl.R
    │   ├── eSNPKaryo_loop.sh
    │   └── eSNPKaryo.R
    └── s08_outlier_analysis
        ├── outlier_detection.R
        └── process_vcf_in_R.R
```

# Input Files

- Bulk RNAseq fastq.gz files
- WGS data (for multiplexed samples)

# Script Descriptions

- `s01_bulk_rna_preprocessing` contains scripts to run fastqc, STAR alignment, RSEM, and doublet calling on bulk RNAseq fastq files
- `s02_germline_calling` takes STAR alignment outputs to call germline variants using GATK pipeline
- `s03_pluripotency` uses the CellNet script to assess pluripotency (alignment to embryonic stem cell geneset)
- `s04_cancer_mutation_calling` uses the COSMIC database to asses bulkRNAseq data for presence of cancerous mutations
- `s05_demultiplexing_donors` contains scripts for both demuxlet and souporcell used on other RNAseq datasets which can be applied to this data
- `s06_mycoplasma_detection` contains scripts using bowtie to assess alignment of the bulk RNA data to the mycoplasma genome
- `s07_eKaryotyping` contains scripts to run eSNPKaryo in R, and create deletion tables to generate useful eSNPKaryo plots
- `s06_outlier_analysis` contains scripts for running outlier analysis algorithms on the processed data

# Main Features

# Data Sources

# Installation

# Usage

# Cite

# Maintainer

# Acknowledgements

# Release Notes
