#!/bin/bash

# Input files
VCF_FILE="/gpfs/commons/groups/singh_lab/users/dongwang/RNA-variant-calling/24246R-06-01.variant_filtered.vcf.gz"
TAR_FILE="/gpfs/commons/groups/singh_lab/users/kjakubiak/COSMIC_DB/Cosmic_MutantCensus_Tsv_v101_GRCh38.tar"
GENE_LIST="/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/gene_list.txt"
OUTPUT_FILE="results_1-31.tsv"


# Extract TSV file from tar
TSV_FILE=$(tar -tf "$TAR_FILE" | grep '\.tsv.gz$')

# Check if the file was found in the archive
if [[ -z "$TSV_FILE" ]]; then
    echo "Error: No .tsv.gz file found in $TAR_FILE"
    exit 1
fi

echo "Extracting $TSV_FILE..."
tar -xvf "$TAR_FILE" "$TSV_FILE"

# Unzip the TSV file
gunzip -c "$TSV_FILE" > extracted_data.tsv

# Check if extracted data is not empty
if [[ ! -s extracted_data.tsv ]]; then
    echo "Error: Extracted TSV file is empty!"
    exit 1
fi

# Filter the TSV file based on gene symbols
grep -F -f "$GENE_LIST" extracted_data.tsv > filtered_data.tsv

# Check if any genes were found
if [[ ! -s filtered_data.tsv ]]; then
    echo "Warning: No matching genes found in extracted data."
    exit 0
fi

# Convert filtered data into BED format for VCF intersection
awk 'BEGIN {OFS="\t"} {print $2, $3, $4}' filtered_data.tsv > regions.bed

# Check if regions file is created
if [[ ! -s regions.bed ]]; then
    echo "Error: No regions found to intersect with VCF."
    exit 1
fi

# Determine if the VCF is compressed
if [[ "$VCF_FILE" == *.gz ]]; then
    VCF_FLAG="--gzvcf"
else
    VCF_FLAG="--vcf"
fi

# Use vcftools to extract VCF records that fall within the BED regions
vcftools $VCF_FLAG "$VCF_FILE" --bed regions.bed --out filtered_vcf --recode --keep-INFO-all

# Check if vcftools produced output
if [[ ! -s filtered_vcf.recode.vcf ]]; then
    echo "No VCF records matched the specified regions."
    exit 0
fi

# Extract relevant VCF data and merge with filtered TSV info
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO\n' filtered_vcf.recode.vcf > matched_vcf.tsv

# Merge matched VCF entries with filtered_data.tsv based
