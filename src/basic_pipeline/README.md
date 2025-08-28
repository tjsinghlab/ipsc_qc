# Pipeline Runner

This pipeline will call cancer mutations in select cancerous genes, run eSNPKaryotyping, determine mycoplasma detection, run PACNet, and identify outliers in bulk RNA sequencing data.
This pipeline expects outputs from STAR alignment, RSEM, and GATK variant calling pipelines.

---

## Table of Contents

- [Requirements](#requirements)
- [Directory Structure](#directory-structure)
- [Inputs](#inputs)
- [Pipeline Runner Arguments](#pipeline-runner-arguments)
- [Reference Files](#reference-files)
- [Running the Pipeline](#running-the-pipeline)
- [Troubleshooting](#troubleshooting)

---

## Table of Contents

- [Requirements](#requirements)
- [Directory Structure](#directory-structure)
- [Inputs](#inputs)
- [Pipeline Runner Arguments](#pipeline-runner-arguments)
- [Reference Files](#reference-files)
- [Running the Pipeline](#running-the-pipeline)
- [Troubleshooting](#troubleshooting)

---

## Requirements

### Tools

The pipeline depends on the following command-line tools:

- `bowtie2`
- `samtools`
- `bcftools`
- `vcftools`
- `aws` CLI (only needed if PACNet training files are missing)

Ensure these tools are installed or available via a module system. For AWS CLI, make sure v2 is installed.

### R Packages

The pipeline uses R scripts in several modules. Install the necessary R packages as listed in the module scripts.

---

## Directory Structure

The pipeline expects separate **reference** and **input** directories:

project/
├── refs/ # Reference files (unchanging across runs)
├── inputs/ # Input files for a specific sample run
├── outputs/ # Output directory
├── modules/ # Contains all pipeline modules
└── pipeline_runner.sh

---

## Inputs

### Required Directories

| Argument         | Description                                                                 | Default                        | Required |
|-----------------|-----------------------------------------------------------------------------|--------------------------------|----------|
| `--ref_dir`      | Directory containing reference files                                        | `/refs`                        | Yes      |
| `--inputs_dir`   | Directory containing input files for a run                                  | `./inputs`                     | Yes      |
| `--output_dir`   | Output directory                                                           | `./outputs`                    | No       |
| `--project`      | Project name                                                               | `default_project`              | No       |

### Module-specific inputs

#### Cancer Mutation Calling

- **Refs**:  
  - `ref/Cosmic_CancerGeneCensus_Tsc_v101_GRCh38.tar` (COSMIC Database)  
  - `ref/genes.txt` — plain text file with genes of interest; default list if missing:  
    ```
    BRCA1
    BRCA2
    TP53
    BCOR
    EGFR
    ```  
- **Inputs**:  
  - `inputs/vcfs/` — folder of VCF files from GATK variant calling  

#### Mycoplasma Detection

- **Refs**:  
  - `ref/GCF_000027235.1_ASM273v1_genomic.fna.gz` — Mycoplasma genome  
    > If missing, automatically downloaded from NCBI:  
    ```bash
    wget -O GCF_000027235.1_ASM273v1_genomic.fna.gz \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/235/GCF_000027235.1_ASM273v1/GCF_000027235.1_ASM273v1_genomic.fna.gz
    ```  

#### PACNet

- **Refs**:  
  - `ref/Hs_expTrain_Jun-20-2017.rda`  
  - `ref/Hs_stTrain-Jun-20-2017.rda`  
    > Can be downloaded automatically using AWS CLI (public S3, no credentials):  
    ```bash
    aws s3 cp s3://cellnet-rnaseq/ref/cnproc/HS/Hs_expTrain_Jun-20-2017.rda . --no-sign-request
    aws s3 cp s3://cellnet-rnaseq/ref/cnproc/HS/Hs_stTrain-Jun-20-2017.rda . --no-sign-request
    ```  
- **Inputs**:  
  - `inputs/RSEM/` — folder of RSEM outputs from preprocessing (pipeline will create gene expression matrix and metadata)  

#### eSNPKaryo

- **Refs**:  
  - `ref/chr/` — directory containing dbSNP build 142 common SNP files (one per chromosome)  
    - File names: `chr1` … `chr24` (human chromosomes 1–22, X=23, Y=24)  
- **Inputs**:  
  - `inputs/bams/` — BAM files from STAR alignment  
  - `inputs/vcfs/` — same VCFs used for cancer mutation calling  

#### Outlier Detection / Mycoplasma Detection (fastq)

- **Inputs**:  
  - `inputs/fastqs/` — folder of fastq.gz files from sequencing run (bulk RNA-seq)

---

## Running the Pipeline

Example command:

```bash
bash pipeline_runner.sh \
    --project my_project \
    --ref_dir /path/to/refs \
    --vcf_dir /path/to/vcfs \
    --rsem_dir /path/to/RSEM \
    --bam_dir /path/to/bams \
    --fastq_dir /path/to/fastqs \
    --output_dir /path/to/outputs
