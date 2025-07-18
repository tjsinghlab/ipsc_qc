#!/bin/bash
#SBATCH --job-name=soup
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --output=stdout_%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kjakubiak@nygenome.org
#SBATCH --time=168:00:00
#SBATCH --array=1-16

conda activate soup
module load singularity
module load bcftools
module load htslib

# Define paths
WGS_VCF="/gpfs/commons/groups/singh_lab/users/dongwang/hbccseq/data/analysis/2.0_eQTL/inputs/hbcc_common_variant_for_eqtl.vcf.bgz"
CELLRANGER_DIR1="/gpfs/commons/groups/singh_lab/projects/hbccseq/data/raw/2025_01_06_pooled_data/CellRangerOuts/pooledLibrary_CellrangerCounts/projectL3_121621yp_cellrangerCounts"  # First directory containing some samples
CELLRANGER_DIR2="/gpfs/commons/groups/singh_lab/projects/hbccseq/data/raw/2025_01_06_pooled_data/CellRangerOuts/pooledLibrary_CellrangerCounts/projectL3_031722yp_cellrangerCounts"  # Second directory containing remaining samples
CELLBENDER_DIR="/gpfs/commons/groups/singh_lab/users/kjakubiak/hbcc_snRNAseq/pool_analysis_2/cellbender"
OUTPUT_DIR="/gpfs/commons/groups/singh_lab/users/kjakubiak/hbcc_snRNAseq/pool_analysis_2/souporcell_with_WGS"
SINGULARITY_IMAGE="/gpfs/commons/groups/singh_lab/users/kjakubiak/bin/souporcell_release.sif"
REF_FASTA="/gpfs/commons/groups/singh_lab/users/kjakubiak/ref_fasta/fasta"
BATCH_SAMPLE_CSV="/gpfs/commons/groups/singh_lab/users/kjakubiak/hbcc_snRNAseq/pool_analysis_2/batch_sample.csv"


# Ensure necessary environment variables are set
if [[ -z "$BATCH_SAMPLE_CSV" || -z "$WGS_VCF" || -z "$OUTPUT_DIR" || -z "$CELLBENDER_DIR" ]]; then
    echo "Error: Required environment variables BATCH_SAMPLE_CSV, WGS_VCF, OUTPUT_DIR, or CELLBENDER_DIR not set."
    exit 1
fi

# Step 1: Precompute Subset VCFs (Only for VCF Versions 1, 2, 3)
echo "Precomputing subset VCFs..."

# Only process versions 1, 2, and 3
for VCF_VERSION in 1 2 3; do
    VCF_FILE="$OUTPUT_DIR/subset_vcf_${VCF_VERSION}.vcf.gz"

    # Skip if already created
    if [[ -f "$VCF_FILE" ]]; then
        echo "Subset VCF already exists: $VCF_FILE"
        continue
    fi

    echo "Creating subset VCF for version: $VCF_VERSION"

    # Extract unique sample names for this VCF version
    SAMPLES=$(awk -F',' -v vcf_version="$VCF_VERSION" 'NR>1 && $6 == vcf_version {print $2"\n"$3"\n"$4"\n"$5}' "$BATCH_SAMPLE_CSV" | grep -v '^$' | sort -u | tr '\n' ',' | sed 's/,$//')

    if [[ -z "$SAMPLES" ]]; then
        echo "Error: No samples found for VCF version $VCF_VERSION in $BATCH_SAMPLE_CSV"
        exit 1
    fi

    echo "Creating BGZF-compressed subset VCF: $VCF_FILE for samples: $SAMPLES"

    # Generate subset VCF and ensure BGZF compression
    bcftools view -s "$SAMPLES" "$WGS_VCF" -Oz -o "$VCF_FILE"
    
    # Index the VCF
    tabix -p vcf "$VCF_FILE"

    echo "Subset VCF created and indexed: $VCF_FILE"
done

echo "All required subset VCFs (1, 2, 3) are ready. Moving to Souporcell execution."

# Step 2: Run Souporcell for each pool
POOL_LIST=($(ls "$CELLBENDER_DIR"))
POOL=${POOL_LIST[$SLURM_ARRAY_TASK_ID - 1]}

echo "Processing pool: $POOL"

# Check if consensus step is already completed
CONSENSUS_DONE="$OUTPUT_DIR/$POOL/consensus.done"
if [[ -f "$CONSENSUS_DONE" ]]; then
    echo "Skipping $POOL: consensus.done file found."
    exit 0
fi

# Define the path to the filtered barcodes file
BARCODES_FILE="$CELLBENDER_DIR/$POOL/cellbender_output_cell_barcodes.csv"

# Find matching BAM file
if [[ -d "$CELLRANGER_DIR1/L3-$POOL" ]]; then
    BAM_FILE="$CELLRANGER_DIR1/L3-$POOL/outs/possorted_genome_bam.bam"
elif [[ -d "$CELLRANGER_DIR2/counts_bothIndex_L3-$POOL" ]]; then
    BAM_FILE="$CELLRANGER_DIR2/counts_bothIndex_L3-$POOL/outs/possorted_genome_bam.bam"
else
    echo "Skipping $POOL: No matching BAM file found."
    exit 1
fi

echo "Found BAM file: $BAM_FILE"

# Determine which VCF file (if any) should be used
VCF_FILE=""
VCF_VERSION=$(awk -F',' -v pool="$POOL" 'NR>1 && $1 == pool {print $6}' "$BATCH_SAMPLE_CSV" | xargs)  # Strip any leading/trailing spaces

# Skip Souporcell if VCF_version is NA
if [[ "$VCF_VERSION" == "NA" ]]; then
    echo "Skipping Souporcell for $POOL (VCF_version is NA)."
    exit 0
else
    # Ensure that the VCF version is valid
    VCF_VERSION=$(echo "$VCF_VERSION" | tr -d '\r' | tr -d '\n')  # Remove carriage return and newlines
    if [[ "$VCF_VERSION" != "1" && "$VCF_VERSION" != "2" && "$VCF_VERSION" != "3" ]]; then
        echo "Error: Invalid VCF version ($VCF_VERSION) for pool $POOL. Skipping."
        exit 1
    fi
    
    VCF_FILE="$OUTPUT_DIR/subset_vcf_${VCF_VERSION}.vcf.gz"

    if [[ ! -f "$VCF_FILE" ]]; then
        echo "Error: Expected subset VCF does not exist: $VCF_FILE"
        exit 1
    fi

    echo "Using precomputed subset VCF: $VCF_FILE"
fi

echo "Final checkpoint before running Souporcell:"
echo "Pool: $POOL"
echo "VCF File: ${VCF_FILE:-None}"

# Create output directory for the pool
POOL_OUTPUT="$OUTPUT_DIR/$POOL"
mkdir -p "$POOL_OUTPUT"

echo "Running Souporcell pipeline for $POOL"

# Run Souporcell with Singularity, including VCF if applicable
singularity exec \
    --bind "$CELLRANGER_DIR1":"$CELLRANGER_DIR1" \
    --bind "$CELLRANGER_DIR2":"$CELLRANGER_DIR2" \
    --bind "$CELLBENDER_DIR":"$CELLBENDER_DIR" \
    --bind "$OUTPUT_DIR":"$OUTPUT_DIR" \
    --bind "$REF_FASTA":"$REF_FASTA" \
    "$SINGULARITY_IMAGE" \
    souporcell_pipeline.py \
    --bam "$BAM_FILE" \
    --barcodes "$BARCODES_FILE" \
    --fasta "$REF_FASTA/genome.fa" \
    --out_dir "$POOL_OUTPUT" \
    --threads 8 \
    --clusters 4 \
    $( [[ -n "$VCF_FILE" ]] && echo "--known_genotypes $VCF_FILE" )

echo "Souporcell pipeline completed for $POOL"
