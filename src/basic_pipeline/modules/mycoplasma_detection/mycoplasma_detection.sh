#!/bin/bash
#MYCOPLASMA DETECTION
#NEEDS FASTQ FILES AND MYCOPLASMA REFERENCE GENOME
#Use bowtie2 to align reads to mycoplasma genome, and assess alignment
#The alignment will be a proxy for level of contamination

#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz
#--2025-01-24 16:08:05--  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz
#/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/GCF_000027325.1_ASM2732v1_genomic.fna.gz

#!/usr/bin/env bash
set -euo pipefail

# Directory with user-supplied references
REF_DIR="${1:-/refs}"  # default to /refs if not provided

# Make sure the directory exists
mkdir -p "$REF_DIR"

# File name and path
GENOME_FILE="GCF_000027325.1_ASM2732v1_genomic.fna.gz"
GENOME_PATH="$REF_DIR/$GENOME_FILE"

# URL to download if file is missing
GENOME_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/$GENOME_FILE"

# Check if file exists
if [ -f "$GENOME_PATH" ]; then
    echo "[INFO] Using existing genome file: $GENOME_PATH"
else
    echo "[INFO] Genome file not found, downloading..."
    wget -O "$GENOME_PATH" "$GENOME_URL"
    echo "[INFO] Download complete: $GENOME_PATH"
fi





#module load bowtie2
# Check if bowtie2 is available
if command -v bowtie2 >/dev/null 2>&1; then
    echo "[INFO] bowtie2 found in PATH"
else
    # Try HPC module
    if command -v module >/dev/null 2>&1; then
        echo "[INFO] Loading bowtie2 module..."
        module load bowtie2
    else
        echo "[ERROR] bowtie2 not found. Please install it."
        exit 1
    fi
fi


#Build mycoplasma reference
bowtie2-build mycoplasma_genome.fna mycoplasma_index

#Align reads to mycoplasma reference
#Single-end:
bowtie2 -x mycoplasma_index -U sample_reads.fastq -S output.sam
#Paired-end:
bowtie2 -x mycoplasma_index -1 sample_R1.fastq -2 sample_R2.fastq -S output.sam

#Get alignment summary
bowtie2 -x mycoplasma_index -U sample_reads.fastq -S output.sam --no-unal 2> bowtie2_log.txt
cat bowtie2_log.txt
#Check alignment rate in the log to assess mycoplasma contamination

#Convert SAM to BAM to make our lives easier later
samtools view -bS output.sam > output.bam
samtools sort output.bam -o output_sorted.bam
samtools index output_sorted.bam