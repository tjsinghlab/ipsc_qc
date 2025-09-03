# ----------------------
# Base image: pinned R version for reproducibility
# ----------------------
FROM rocker/r-ver:4.3.1

LABEL maintainer="Kelly Jakubiak <kjakubiak@nygenome.org>"
LABEL description="RNA-seq pipeline image with STAR, RSEM, GATK, PACNet, eSNPKaryotyping, mycoplasma detection, and per-sample analysis"

# ----------------------
# System dependencies for Bash, bioinformatics tools, Python, Java
# ----------------------
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    wget curl unzip tar gzip bzip2 git \
    build-essential libcurl4-openssl-dev libssl-dev libxml2-dev \
    zlib1g-dev libbz2-dev liblzma-dev libpcre3-dev \
    bowtie2 samtools bcftools vcftools default-jdk python3 python3-pip && \
    rm -rf /var/lib/apt/lists/*

# ----------------------
# CRAN snapshot for reproducible R packages
# ----------------------
ENV CRAN_REPO=https://packagemanager.posit.co/cran/2024-01-01
RUN R -e "install.packages(c('devtools','zoo','gplots','patchwork','ggplot2','optparse','dplyr','data.table','tidyr','fuzzyjoin','stringr','biomaRt','purrr','jsonlite','edgeR'), repos='https://cloud.r-project.org/')"

# Install CellNet and cancerCellNet from tarball if included
COPY CellNet_*.tar.gz /opt/CellNet/
RUN R -e "install.packages('/opt/CellNet/CellNet_*.tar.gz', repos=NULL, type='source')"

# ----------------------
# Add Cromwell
# ----------------------
ADD https://github.com/broadinstitute/cromwell/releases/download/86/cromwell-86.jar /pipeline/cromwell.jar

# -------------------------------
# Set up working directories
# -------------------------------
WORKDIR /workspace
ENV REF_DIR=/opt/ref
ENV OUTPUT_DIR=/workspace/outputs
ENV FASTQ_DIR=/workspace/fastq_dir

RUN mkdir -p $REF_DIR $OUTPUT_DIR $FASTQ_DIR

# -------------------------------
# Copy reference files into container
# -------------------------------
COPY ref/star_index_oh75 $REF_DIR/star_index_oh75
COPY ref/rsem_reference $REF_DIR/rsem_reference
COPY ref/gencode.v39.GRCh38.genes.collapsed_only.gtf $REF_DIR/
COPY ref/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta $REF_DIR/
COPY ref/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta.fai $REF_DIR/
COPY ref/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.dict $REF_DIR/
COPY ref/Homo_sapiens_assembly38.dbsnp138.vcf $REF_DIR/
COPY ref/Homo_sapiens_assembly38.dbsnp138.vcf.idx $REF_DIR/
COPY ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz $REF_DIR/
COPY ref/Homo_sapiens_assembly38.known_indels.vcf.gz $REF_DIR/
COPY ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi $REF_DIR/
COPY ref/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi $REF_DIR/
COPY ref/images/gatk_4.6.1.0.sif $REF_DIR/images/
COPY ref/genes.txt $REF_DIR/
COPY ref/chr $REF_DIR/chr
COPY ref/GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna $REF_DIR/
COPY ref/Hs_expTrain_Jun-20-2017.rda $REF_DIR/
COPY ref/Hs_stTrain_Jun-20-2017.rda $REF_DIR/
# Cosmic database is **not included** due to licensing, user must mount/download
# COPY ref/Cosmic_CancerGeneCensus_Tsv_v101_GRCh37.tar $REF_DIR/

# Mycoplasma genome
COPY ref/GCF_000027325.1_ASM2732v1_genomic.fna.gz $REF_DIR/

# -------------------------------
# Set PATH for tools
# -------------------------------
ENV PATH=/usr/bin:$PATH

# ----------------------
# Make shell scripts executable
# ----------------------
RUN chmod +x /pipeline/*.sh

# ----------------------
# Default entrypoint
# ----------------------
ENTRYPOINT ["bash", "/pipeline/pipeline_runner.sh"]
