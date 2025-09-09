# iPSC QC Pipeline

This pipeline will screen bulk RNA sequencing data for several QC metrics, including cancer mutation calling (in select oncogenes), eSNPKaryotyping, mycoplasma detection, PACNet classification, and outlier assessment.
This pipeline expects a directory of fastq.gz (paired end: R1 and R2) files from bulk RNA sequencing run.

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
### Reference Files (place in `ref/` directory; default is ./ref)
#### Must be provided by user
- `Cosmic_CancerGeneCensus_Tsv_v101_GRCh37.tar`  
  *Cosmic Database; MUST BE DOWNLOADED BY USER AND PROVIDED IN COSMIC_DIR ARGUMENT. If not, cancer mutation calling will not be performed.*
#### Must be downloaded by user from *insert cloud link here*
- `star_index_oh75`
  *reference files for STAR; included in docker image*
- `rsem_reference`
  *reference files for RSEM; included in docker image*
- `gencode.v39.GRCh38.genes.collapsed_only.gtf`
  *reference GTF; included in docker image*
- `Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta`
  *refFasta; included in docker image*
- `Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta.fai`
  *refFastaIndex; included in docker image*
- `Homo_sapiens_assembly38_noALT_noHLA_noDecoy.dict`
  *refDict; included in docker image*
- `Homo_sapiens_assembly38.dbsnp138.vcf`
  *dbSnpVcf; included in docker image*
- `Homo_sapiens_assembly38.dbsnp138.vcf.idx` 
  *dbSnpVcfIndex; included in docker image*
- `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`
  *knownVcf 1; included in docker image*
- `Homo_sapiens_assembly38.known_indels.vcf.gz`
  *knownVcf 2; included in docker image*
- `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi`
  *knownVcfIndex 1; included in docker image*
- `Homo_sapiens_assembly38.known_indels.vcf.gz.tbi`
  *knownVcfIndex 2; included in docker image*
- `gencode.v39.GRCh38.genes.collapsed_only.gtf`
  *annotationsGTF; included in docker image*
- `images/gatk_4.6.1.0.sif`
  *gatk4_docker; included in docker image*
- `genes.txt`  
  *List of genes of interest (default: BRCA1, BRCA2, TP53, BCOR, EGFR); included in docker image*
- `GCF_000027325.1_ASM2732v1_genomic.fna.gz`  
  *Mycoplasma genome; included in docker image but can also be downloaded with `wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz`*
- `Hs_expTrain_Jun-20-2017.rda`  
  *PACNet training data; included in docker image but can also be downloaded with `aws s3 cp s3://cellnet-rnaseq/ref/cnproc/HS/Hs_expTrain_Jun-20-2017.rda . --no-sign-request`*
- `Hs_stTrain_Jun-20-2017.rda`  
  *PACNet training data; included in docker image but can also be downloaded with `aws s3 cp s3://cellnet-rnaseq/ref/cnproc/HS/Hs_stTrain_Jun-20-2017.rda . --no-sign-request`*
- `chr/`  
  *Directory with dbSNP build 142 common SNP files (GTF format), one per chromosome. Human: X = 23, Y = 24. Source: UCSC Table Browser, snp142common table. Included in docker image.*
- `GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna`  
  *NCBI reference genome FASTA file; included in docker image, but can also be downloaded from `ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38.p13`*

### Tools

The pipeline depends on the following command line tools:

- `bowtie2`
- `samtools`
- `bcftools`
- `vcftools`

Ensure these tools are installed or available via a module system (ex. HPC users). For AWS CLI, make sure v2 is installed.
R packages will be included as part of the docker image.

R packages included in the docker image:
- devtools
- eSNPKaryotyping
- zoo
- gplots
- patchwork
- ggplot2
- optparse
- dplyr
- data.table
- tidyr
- fuzzyjoin
- stringr
- biomaRt
- purrr
- CellNet (tarball included in docker image)
- cancerCellNet
- jsonlite
- edgeR

Python packages included in the docker image:
- os
- json
- logzero (logger)

---

## Directory Structure

The pipeline expects separate **reference** and **input** directories:

project/
- ├── refs/ # Reference files (unchanging across runs; only need to provide cosmic ref file)
- ├── inputs/ # Input files for a run (pipeline will loop over samples)
- ├── outputs/ # Output directory
- ├── modules/ # Contains all pipeline modules
- └── pipeline_runner.sh

---

## Inputs

### Input Directories
- `fastq_dir/`  
  *FASTQ files from sequencing run (bulk RNAseq)*
- `--output_dir`  
  *Path to desired output directory (created if it doesn't exist; default: `./outputs`)*

### Additional Arguments
- `--project`  
  *Name of project for this batch; appears in plots*

**Note:**  
File base names should be unique for each sample and consistent across file types  
(e.g., `sample1.vcf`, `sample1.bam`, `sample1.genes.results`, `sample1_R1.fastq.gz`, `sample1_R2.fastq.gz`)


### Required Directories

| Argument         | Description                                                                 | Default                        | Required |
|-----------------|-----------------------------------------------------------------------------|--------------------------------|----------|
| `--ref_dir`      | Directory containing reference files                                        | `/refs`                        | Yes      |
| `--fastq_dir`   | Directory containing input files (fastq.gz sequencing files) for a run                                  | `./fastq`                     | Yes      |
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
    > Can be downloaded automatically using AWS CLI (public S3, no credentials needed) in accordance with PACNet repo:  
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

Example command (replace all "my_***" directories with your actual directories)

```
docker run -it --rm \
    -v ./my_fastqs:/data \
    -v ./my_refs:/ref \
    -v ./my_cosmic:/cosmic \
    -v ./my_outputs:/output \
    ipsc_qc_image
```


