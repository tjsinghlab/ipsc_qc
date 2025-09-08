# ----------------------
# Base image: R with version pinned
# ----------------------
FROM rocker/r-ver:4.3.1

LABEL maintainer="Kelly Jakubiak <kjakubiak@nygenome.org>"
LABEL description="Bulk RNA-seq pipeline with pre-processing, GATK germline variant calling, and QC modules."

# ----------------------
# System dependencies
# ----------------------
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    wget curl unzip tar gzip bzip2 git \
    build-essential libcurl4-openssl-dev libssl-dev libxml2-dev \
    zlib1g-dev libbz2-dev liblzma-dev libpcre3-dev \
    bowtie2 samtools bcftools vcftools default-jdk python3 python3-pip && \
    rm -rf /var/lib/apt/lists/*

# ----------------------
# R packages
# ----------------------
ENV CRAN_REPO=https://packagemanager.posit.co/cran/2024-01-01
RUN R -e "install.packages(c('devtools','zoo','gplots','patchwork','ggplot2','optparse','dplyr','data.table','tidyr','fuzzyjoin','stringr','biomaRt','purrr','jsonlite','edgeR','rmarkdown','knitr'), repos='${CRAN_REPO}')"

# ----------------------
# Install CellNet from tarball
# ----------------------
COPY CellNet_*.tar.gz /opt/CellNet/
RUN R -e "install.packages('/opt/CellNet/CellNet_*.tar.gz', repos=NULL, type='source')"

# ----------------------
# Cromwell JAR
# ----------------------
ADD https://github.com/broadinstitute/cromwell/releases/download/86/cromwell-86.jar /pipeline/cromwell.jar

# ----------------------
# Set up directories
# ----------------------
WORKDIR /workspace
ENV REF_DIR=/ref

# -------------------------------
# Copy reference files into container
# -------------------------------
COPY ref/ /ref/

# The /ref folder in the above command contains the following:
# star_index_oh75
# rsem_reference
# gencode.v39.GRCh38.genes.collapsed_only.gtf
# Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fastaREF_DIR
# Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta.fai
# Homo_sapiens_assembly38_noALT_noHLA_noDecoy.dict
# Homo_sapiens_assembly38.dbsnp138.vcf
# Homo_sapiens_assembly38.dbsnp138.vcf.idx
# Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
# Homo_sapiens_assembly38.known_indels.vcf.gz
# Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
# Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
# gatk_4.6.1.0.sif
# genes.txt
# chr
# GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna
# Hs_expTrain_Jun-20-2017.rda
# Hs_stTrain_Jun-20-2017.rda
# GCF_000027325.1_ASM2732v1_genomic.fna.gz
# gtex_rnaseq_V10.sif

# Cosmic database is **not included** due to licensing, user must mount/download
# Cosmic_CancerGeneCensus_Tsv_v101_GRCh37.tar

# -------------------------------
# Copy in all scripts from local /scripts into container /pipeline
# -------------------------------
COPY scripts/ /pipeline/

# -------------------------------
# Make scripts executable
# -------------------------------
RUN chmod +x /pipeline/*.sh && chmod +x /pipeline/modules/*/*

# -------------------------------
# Install wdlplay requirements
# -------------------------------
RUN make -C /pipeline/preprocessing/wdlplay/Makefile install

# ----------------------
# Environment variables (user MUST override these when running container)
# ----------------------
# User must supply their own output and cosmic directories when running the container
ENV OUTPUT_DIR=""
ENV COSMIC_DIR=""

# ----------------------
# Default entrypoint
# ----------------------
ENTRYPOINT ["bash", "/pipeline/pipeline_runner.sh"]
