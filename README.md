# ipsc_qc
QC pipeline for iPSC cultures

# Overview

Performing quality control from multiple angles on iPSC lines to ensure viability, reproducibility, and scability in cultures.

# Directory Structure and Organization

```
.
├── README.md
└── src
    ├── basic_pipeline
    │   ├── Dockerfile
    │   ├── modules
    │   │   ├── cancer_mutation_calling
    │   │   │   ├── cancer_mutation_mapping.R
    │   │   │   ├── cosmic_loop.sh
    │   │   │   └── filter_vcf_on_cosmic.sh
    │   │   ├── eSNPKaryo
    │   │   │   └── eSNPKaryotyping.R
    │   │   ├── mycoplasma_detection
    │   │   │   └── mycoplasma_detection.sh
    │   │   ├── outlier_detection
    │   │   │   └── outlier_detection.R
    │   │   └── PACNet
    │   │       ├── pacnet.R
    │   │       └── process_rsem_outputs.R
    │   ├── pipeline_runner.sh
    │   ├── README.md
    │   ├── requirements.txt
    │   ├── scripts
    │   │   └── reference_check.sh
    │   └── tarballs
    │       └── CellNet_master.tar.gz
    ├── s01.1_bulk_rna_preprocessing
    │   ├── bulk_rna_processing.py
    │   └── cram_to_fastq.sh
    ├── s01.2_analysis
    │   ├── batch_correction
    │   │   ├── deseq.R
    │   │   ├── limma_voom.R
    │   │   ├── peer_prep_and_process.R
    │   │   └── run_peer.py
    │   ├── create_variant_tables.sh
    │   ├── deconvolution
    │   │   ├── bisque.R
    │   │   ├── bMIND.R
    │   │   ├── bMIND.sh
    │   │   ├── CIBERSORT.R
    │   │   └── plot_deconvolution.R
    │   ├── outlier_detection
    │   │   └── outlier_detection.R
    │   ├── plotting.R
    │   ├── process_puigdevall_data.R
    │   ├── process_rsem_outputs.R
    │   ├── process_vcf_in_R.R
    │   └── read_digital_gene_expression.R
    ├── s02_germline_calling
    │   ├── gatk4_rna_germline_calling_run.py
    │   ├── gatk4-rna-best-practices.wdl
    │   └── gatk4-rna-germline-variant-calling.inputs.json
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
    └── s07_eKaryotyping
        ├── delTabl.R
        ├── eSNPKaryo_loop.sh
        └── eSNPKaryo.R
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
