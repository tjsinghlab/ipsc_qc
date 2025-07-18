#!/bin/bash
#SBATCH --job-name=demuxlet
#SBATCH --partition=gpu
#SBATCH --mem=32G
#SBATCH --output=stdout_%j.log
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kjakubiak@nygenome.org

module load samtools
module load bcftools

/gpfs/commons/groups/singh_lab/users/kjakubiak/hbcc_snRNAseq/pool_cellranger_outs_from_hbcc

# Define the base directory containing the sample folders
BASE_DIR="/gpfs/commons/groups/singh_lab/users/kjakubiak/hbcc_snRNAseq/"

# Define the temporary directory for intermediate files
TMPDIR="/gpfs/commons/groups/singh_lab/users/kjakubiak/hbcc_snRNAseq/tmp"

# Define the output directory for the final BAM files
OUTDIR="/gpfs/commons/groups/singh_lab/users/kjakubiak/hbcc_snRNAseq/pooled_sorted_bams"

# Get the list of sample-specific folders (the directories in the base directory)
SAMPLE_FOLDERS=("$BASE_DIR"/*/)

# Get the specific sample folder for this job in the array
SAMPLE_FOLDER="${SAMPLE_FOLDERS[$SLURM_ARRAY_TASK_ID]}"

# Define the BAM file path
BAM_FILE="${SAMPLE_FOLDER}/outs/possorted_genome_bam.bam"

# Define output names
OUTPUT_BAM="${OUTDIR}/$(basename "$SAMPLE_FOLDER")_reordered.bam"
OUTPUT_HEADER="${TMPDIR}/$(basename "$SAMPLE_FOLDER")/new_header.sam"
OUTPUT_TMPDIR="${TMPDIR}/$(basename "$SAMPLE_FOLDER")"

# Create a sample-specific temporary directory
mkdir -p "$OUTPUT_TMPDIR"

# Step 1: Extract the BAM header
samtools view -H "$BAM_FILE" > "$OUTPUT_HEADER"

# Step 2: Reorder the @SQ lines (using a bash script or Python as discussed earlier)
grep "^@SQ" "$OUTPUT_HEADER" > "$OUTPUT_TMPDIR/sq_lines.txt"
grep -v "^@SQ" "$OUTPUT_HEADER" > "$OUTPUT_TMPDIR/other_lines.txt"

# Sort chromosomes numerically and place chrX, chrY, chrM in the correct order
grep -E "SN:chr[0-9]+" "$OUTPUT_TMPDIR/sq_lines.txt" | sort -V -k2,2 -t: > "$OUTPUT_TMPDIR/sorted_numeric_sq.txt"
grep -E "SN:chrX" "$OUTPUT_TMPDIR/sq_lines.txt" >> "$OUTPUT_TMPDIR/sorted_numeric_sq.txt"
grep -E "SN:chrY" "$OUTPUT_TMPDIR/sq_lines.txt" >> "$OUTPUT_TMPDIR/sorted_numeric_sq.txt"
grep -E "SN:chrM" "$OUTPUT_TMPDIR/sq_lines.txt" >> "$OUTPUT_TMPDIR/sorted_numeric_sq.txt"
grep -vE "SN:chr[0-9]+|SN:chrX|SN:chrY|SN:chrM" "$OUTPUT_TMPDIR/sq_lines.txt" >> "$OUTPUT_TMPDIR/sorted_numeric_sq.txt"

# Combine the sorted header
cat "$OUTPUT_TMPDIR/other_lines.txt" "$OUTPUT_TMPDIR/sorted_numeric_sq.txt" > "$OUTPUT_HEADER"

# Step 3: Sort the BAM file by chromosome and position
samtools sort -o "$OUTPUT_TMPDIR/$(basename "$SAMPLE_FOLDER")_sorted.bam" "$BAM_FILE"

# Step 4: Replace the header in the sorted BAM file
samtools reheader "$OUTPUT_HEADER" "$OUTPUT_TMPDIR/$(basename "$SAMPLE_FOLDER")_sorted.bam" > "$OUTPUT_BAM"

# Step 5: Index the reordered BAM file
samtools index "$OUTPUT_BAM"

# Clean up intermediate files
rm -r "$OUTPUT_TMPDIR"

base_dir="/gpfs/commons/groups/singh_lab/users/kjakubiak/hbcc_snRNAseq/"
BAMBAM=("L3-21309" "L3-21310" "L3-21311" "L3-21312")

module load demuxlet

# Define directories
VCF_FILE="/gpfs/commons/groups/singh_lab/users/dongwang/hbccseq/data/snRNA/hbcc_2257_2298_2351_2470.vcf.bgz"  # Directory where VCF files are stored

OUTPUT_DIR="${base_dir}/case_con_demuxlet"  # Directory where outputs will be stored

for sample in "${BAMBAM[@]}"; do

#Specify location of cellbender barcodes
BARCODES="${base_dir}/case_con_cellbender/pooled/${sample}/cellbender_output_cell_barcodes.csv"  # Directory where barcode files are stored

# BAM file for this pooled sample
BAM_FILE="${base_dir}/sorted_bams/${sample}_cellranger_output_reordered.bam"  # Directory where BAM files are stored

# Output prefix for this pooled sample
SAMPLE_OUTDIR="${OUTPUT_DIR}/${sample}_demuxlet"
mkdir $SAMPLE_OUTDIR

# Check if BAM file exists
if [[ ! -f "$BAM_FILE" ]]; then
    echo "WHERE IS THE BAM FILE FOR $sample????"
    continue
fi

# Check if Barcode file exists
if [[ ! -f "$BARCODES" ]]; then
    echo "WHERE IS THE BARCODE FILE FOR $sample????"
    continue
fi

# Check if VCF file exists
if [[ ! -f "$VCF_FILE" ]]; then
    echo "WHERE IS THE VCF FOR $sample????"
    continue
fi

# Run demuxlet for this pooled sample
echo "Running demuxlet for pooled sample $sample with VCF file $VCF_FILE!!!! :3"

demuxlet \
    --sam "$BAM_FILE" \
    --vcf "$VCF_FILE" \
    --field GT \
    --out "$SAMPLE_OUTDIR" \
    --group-list "$BARCODES"
    
done