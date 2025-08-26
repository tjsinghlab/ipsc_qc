# Base image: R 4.3.1 pinned for reproducibility
FROM rocker/r-ver:4.3.1

# ----------------------
# System dependencies for Bash and bioinformatics tools
# ----------------------
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    curl git make bash wget gzip unzip \
    bowtie2=2.5.1-1build1 \
    samtools=1.17-1build1 \
    bcftools=1.17-1build1 \
    vcftools \
    && rm -rf /var/lib/apt/lists/*

# ----------------------
# CRAN snapshot for reproducible R packages
# ----------------------
ENV CRAN_REPO=https://packagemanager.posit.co/cran/2024-01-01

RUN R -e "install.packages(c('dplyr','ggplot2','data.table','devtools','gplots','zoo','plyr','plotly','igraph','pheatmap','RColorBrewer','tidyverse','stringr','e1071','Seurat'), repos='${CRAN_REPO}')"
RUN R -e "install.packages('BiocManager', repos='${CRAN_REPO}')"
RUN R -e "BiocManager::install(c('GenomicFeatures','GenomicRanges'), ask=FALSE)"

# ----------------------
# R packages from tarballs
# ----------------------
COPY tarballs/ /tarballs/
RUN R CMD INSTALL /tarballs/CellNet_c7f1d6b.tar.gz
RUN R CMD INSTALL /tarballs/eSNPKaryotyping_vX.Y.Z.tar.gz

# ----------------------
# Set working directory
# ----------------------
WORKDIR /pipeline

# Copy pipeline scripts and modules
COPY pipeline_runner.sh /pipeline/
COPY modules/ /pipeline/modules/
COPY scripts/bash/ /pipeline/scripts/bash/

# ----------------------
# Default entrypoint
# ----------------------
ENTRYPOINT ["bash", "pipeline_runner.sh"]
