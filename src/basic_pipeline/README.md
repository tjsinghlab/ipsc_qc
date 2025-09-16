# iPSC QC Pipeline

This pipeline will first process raw fastq files (bulk RNA sequencing) using fastqc, STAR alignment, and RSEM. GATK germline variant calling will be performed, in addition to several QC metrics, including cancer mutation calling (in select oncogenes), eSNPKaryotyping, mycoplasma detection, PACNet classification, and outlier assessment.
This pipeline expects a directory of fastq.gz (paired end: R1 and R2) files from a bulk RNA sequencing run.

---

## Table of Contents

- [Requirements](#requirements)
- [Inputs](#inputs)
- [Running the Pipeline](#running-the-pipeline)
- [Outputs](#outputs)
- [Directory Structure](#directory-structure)

---

## Requirements
### Reference COSMIC Database Files (provide in --cosmic_dir argument; default is ./cosmic)
#### Must be provided by user
- `Cosmic_CancerGeneCensus_Tsv_v101_GRCh37.tar`  
  *Cosmic Database; MUST BE DOWNLOADED BY USER AND PROVIDED IN COSMIC_DIR ARGUMENT. If not, cancer mutation calling will not be performed.*
### Reference Files (place in `ref/` directory; default is ./ref)
- `genes.txt`
  *File containing a list of oncogenes for assessment. Lives in this repo; download into your /ref dir (the one you provide to the --ref_dir argument), and add whichever genes you want analyzed against COSMIC.*
#### If not included in user-provided reference directory (ref_dir), the following files will be automatically downloaded to ref_dir:
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
#### Must be installed by user (Ensure these tools are installed or available via a module system (ex. HPC users)):
- `singularity`
- `parallel`
- `squashfuse`
#### Command line tools included in Docker image:
- `bowtie2`
- `samtools`
- `bcftools`
- `vcftools`
#### R packages included in the docker image:
- `devtools`
- `eSNPKaryotyping`
- `zoo`
- `gplots`
- `patchwork`
- `ggplot2`
- `optparse`
- `dplyr`
- `data.table`
- `tidyr`
- `fuzzyjoin`
- `stringr`
- `biomaRt`
- `purrr`
- `CellNet` (tarball included in docker image)
- `cancerCellNet`
- `jsonlite`
- `edgeR`
#### Python packages included in the docker image:
- `os`
- `json`
- `logzero` (logger)

## Inputs

### Input Directories
- `fastq_dir/`  
  *FASTQ files from sequencing run (bulk RNAseq)*
- `--output_dir`  
  *Path to desired output directory (created if it doesn't exist; default: `./outputs`)*
- `--ref_dir`
  *Path to directory where reference files will live. Absent ref files will be downloaded within the pipeline.*

### Additional Arguments
- `--project`  
  *Name of project for this batch; appears in plots*

**Note:**  
File base names should be unique for each sample and consistent across file types  
(e.g., `sample1.vcf`, `sample1.bam`, `sample1.genes.results`, `sample1_R1.fastq.gz`, `sample1_R2.fastq.gz`)


### Required Directories

| Argument         | Description                                                                | Default                        | Required |
|------------------|----------------------------------------------------------------------------|--------------------------------|----------|
| `--ref_dir`      | Directory for reference files                                       | `/refs`                        | Yes      |
| `--fastq_dir`    | Directory containing input files (fastq.gz sequencing files) for a run     | `./fastq`                      | Yes      |
| `--output_dir`   | Output directory                                                           | `./outputs`                    | No       |
| `--cosmic_dir`   | Directory containing downloaded COSMIC database files.                     | `cosmic`                       | No       |
| `--project`      | Project name                                                               | `default_project`              | No       |

---

## Running the Pipeline

Example command (replace all "my_***" directories with your actual directories)

```
#Download docker image as singularity image (for HPC compatibility):
singularity pull ipsc_image_slim_withmake.sif docker://kjakubiak16/ipsc_image_slim_withmake:latest

#Run the singularity command. First bind the directories, then supply bound directories as arguments.
singularity exec \
  -B ./test_fastqs:/data \
  -B ./ref:/ref \
  -B ./COSMIC_DB:/cosmic \
  -B ./my_outputs:/output \
  ipsc_image_slim_withmake.sif \
  /pipeline/pipeline_runner.sh \
  --fastq_dir /data \
  --output_dir /output \
  --ref_dir /ref \
  --cosmic_dir /cosmic
```

## Outputs
Expected outputs include:

## Directory Structure

### Basic pipeline repo structure
```
basic_pipeline
       ├── docker_image_structure
       ├── docker-compose.yaml
       ├── Dockerfile_slim
       ├── modules
       │   ├── cancer_mutation_calling
       │   │   ├── cancer_mutation_mapping.R
       │   │   └── filter_vcf_on_cosmic.sh
       │   ├── eSNPKaryotyping
       │   │   └── run_eSNPKaryotyping.R
       │   ├── mycoplasma_detection
       │   │   └── detect_mycoplasma.sh
       │   ├── outlier_detection
       │   │   └── outlier_detection.R
       │   ├── PACNet
       │   │   └── run_pacnet.R
       │   └── preprocessing
       │       └── wdlplay
       │           ├── LICENSE.md
       │           ├── Makefile
       │           ├── MANIFEST.in
       │           ├── README.md
       │           ├── requirements.txt
       │           ├── settings.json
       │           ├── setup.cfg
       │           ├── setup.py
       │           ├── versioneer.py
       │           ├── warp-pipelines
       │           │   ├── bulk_RNAseq_preprocess
       │           │   │   ├── RNAseq_pipeline_fastq.wdl
       │           │   │   ├── run_wdl.py
       │           │   │   └── wdl
       │           │   │       ├── fastqc.wdl
       │           │   │       ├── markduplicates.wdl
       │           │   │       ├── rnaseqc2.wdl
       │           │   │       ├── rsem.wdl
       │           │   │       └── star.wdl
       │           │   └── GATK_variant_calling
       │           │       ├── gatk4-rna-best-practices.wdl
       │           │       ├── gatk4-rna-germline-calling_run.py
       │           │       └── gatk4-rna-germline-variant-calling.inputs.json
       │           └── wdlplay
       │               ├── __init__.py
       │               ├── _version.py
       │               ├── db.json
       │               └── wdlplayer.py
       ├── pipeline_runner.sh
       ├── README.md
       ├── report_builder.R
       └── tarballs
           └── CellNet_master.tar.gz
```
---
