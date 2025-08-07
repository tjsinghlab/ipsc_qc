#!/bin/bash
#SBATCH --job-name=cram_to_fastq
#SBATCH --partition=pe2
#SBATCH --mem=12G
#SBATCH --output=stdout_%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kjakubiak@nygenome.org

set -e  # Exit immediately if a command exits with a non-zero status

# Define directories
CRAM_DIR="/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/differentiated_datasets/PRJEB55440/cram_files"
FASTQ_DIR="/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/differentiated_datasets/PRJEB55440/fastq_files"

#mkdir -p "$FASTQ_DIR"
module load samtools
# Loop through all .cram files
for cram in "$CRAM_DIR"/*.cram; do
    base=$(basename "$cram" .cram)
    echo "Processing $cram..."

    # Convert CRAM to FASTQ
    samtools fastq -1 "$FASTQ_DIR/${base}_1.fastq.gz" \
                   -2 "$FASTQ_DIR/${base}_2.fastq.gz" \
                   -0 /dev/null -s /dev/null -n "$cram"

    echo "Finished $cram."
done

echo "All CRAM files processed!"