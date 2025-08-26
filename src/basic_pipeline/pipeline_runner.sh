#!/usr/bin/env bash
set -euo pipefail

# ----------------------
# Directories (mounted by user)
# ----------------------
INPUT_DIR="/inputs"
REF_DIR="/ref"

# ----------------------
# Parse pipeline options
# ----------------------
RUN_CANCER=true
RUN_MYCO=true
RUN_PACNET=true
RUN_EKARYO=true
RUN_OUTLIERS=true

for arg in "$@"; do
    case $arg in
        --no-cancer) RUN_CANCER=false ;;
        --no-myco) RUN_MYCO=false ;;
        --no-pacnet) RUN_PACNET=false ;;
        --no-ekaryo) RUN_EKARYO=false ;;
        --no-outliers) RUN_OUTLIERS=false ;;
    esac
done

# ----------------------
# Check tools and references
# ----------------------
bash scripts/bash/check_refs_and_tools.sh "$REF_DIR"

# ----------------------
# Loop over VCF inputs
# ----------------------
VCF_FILES=("$INPUT_DIR"/*.vcf)
if [ ${#VCF_FILES[@]} -eq 0 ]; then
    echo "[ERROR] No VCF files found in $INPUT_DIR"
    exit 1
fi

# ----------------------
# Run modules
# ----------------------
for vcf in "${VCF_FILES[@]}"; do
    echo "[INFO] Processing: $vcf"

    if [ "$RUN_CANCER" = true ]; then
        echo "[INFO] Running cancer mutation calling..."
        bash modules/cancer_mutation_calling/step1_call_mutations.sh "$vcf" "$REF_DIR"
        Rscript modules/cancer_mutation_calling/step2_process_mutations.R "$vcf" "$REF_DIR"
    fi

    if [ "$RUN_MYCO" = true ]; then
        echo "[INFO] Running mycoplasma detection..."
        bash modules/mycoplasma_detection/detect_mycoplasma.sh "$vcf" "$REF_DIR"
    fi

    if [ "$RUN_PACNET" = true ]; then
        echo "[INFO] Running PacNet..."
        Rscript modules/pacnet/run_pacnet.R "$vcf" "$REF_DIR"
    fi

    if [ "$RUN_EKARYO" = true ]; then
        echo "[INFO] Running eKaryo..."
        Rscript modules/eKaryo/run_eKaryo.R "$vcf" "$REF_DIR"
    fi

    if [ "$RUN_OUTLIERS" = true ]; then
        echo "[INFO] Running outlier detection..."
        Rscript modules/outliers/detect_outliers.R "$vcf" "$REF_DIR"
    fi
done

echo "[INFO] Pipeline completed."
