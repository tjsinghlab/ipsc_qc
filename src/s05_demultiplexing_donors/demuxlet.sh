#!/bin/bash
#SBATCH --job-name=demuxlet
#SBATCH --partition=gpu
#SBATCH --mem=16G
#SBATCH --output=stdout_%j.log
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kjakubiak@nygenome.org

module load demuxlet

# Define directories
VCF_DIR="/path/to/vcf/files"  # Directory where VCF files are stored
#/gpfs/commons/groups/singh_lab/users/dongwang/hbccseq/data/snRNA
BAM_DIR="/path/to/bam/files"  # Directory where BAM files are stored
#/gpfs/commons/groups/singh_lab/users/kjakubiak/hbcc_snRNAseq/sorted_bams
BARCODE_DIR="/path/to/barcode/files"  # Directory where barcode files are stored
#/gpfs/commons/groups/singh_lab/users/kjakubiak/hbcc_snR`NAseq/[SAMPLE_NAME]/cellbender_output_cell_barcodes.csv
OUTPUT_DIR="/path/to/demuxlet_outputs"  # Directory where outputs will be stored
MAP_FILE="/path/to/donor_pool_map.csv"  # Path to donor_pool_map.csv

# Loop through the donor pool map CSV file
while IFS=, read -r pooled_id VCF_file donor1 donor2 donor3 donor4; do
    # Skip header line
    [[ $pooled_id == "pooled_id" ]] && continue

    # Full path to the VCF file for this pooled sample
    VCF_PATH="${VCF_DIR}/${VCF_file}"

    # BAM file for this pooled sample
    BAM_FILE="${BAM_DIR}/${pooled_id}/outs/possorted_genome_bam.bam"

    # Barcode file for this pooled sample
    BARCODE_FILE="${BARCODE_DIR}/${pooled_id}/cellbender_output_cell_barcodes.csv"

    # Output prefix for this pooled sample
    OUTPUT_PREFIX="${OUTPUT_DIR}/${pooled_id}_demuxlet"

    # Check if BAM file exists
    if [[ ! -f "$BAM_FILE" ]]; then
        echo "Missing BAM file for pooled sample $pooled_id, skipping."
        continue
    fi

    # Check if Barcode file exists
    if [[ ! -f "$BARCODE_FILE" ]]; then
        echo "Missing Barcode file for pooled sample $pooled_id, skipping."
        continue
    fi

    # Check if VCF file exists
    if [[ ! -f "$VCF_PATH" ]]; then
        echo "Missing VCF file for pooled sample $pooled_id, skipping."
        continue
    fi

    # Run demuxlet for this pooled sample
    echo "Running demuxlet for pooled sample $pooled_id with VCF file $VCF_PATH."

    demuxlet \
        --sam "$BAM_FILE" \
        --vcf "$VCF_PATH" \
        --field GT \
        --out "$OUTPUT_PREFIX" \
        --group-list "$BARCODE_FILE"
    
done < "$MAP_FILE"


