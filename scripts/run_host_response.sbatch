#!/bin/bash
#SBATCH --job-name=host_response
#SBATCH --output=/scratch/sahlab/kathie/IBD_plant_modular/results/host_response_%j.out
#SBATCH --error=/scratch/sahlab/kathie/IBD_plant_modular/results/host_response_%j.err
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --partition=general
#SBATCH --account=sahlab

# Input arguments
DATA_PATH="$1"
METADATA_FILE="$2"
OUTPUT_DIR="$3"

# Source the shared environment setup (handles all conda activation)
source /scratch/sahlab/kathie/IBD_plant_modular/config/htcf_env_setup.sh host_response "$DATA_PATH" "$METADATA_FILE" "$OUTPUT_DIR"

# Change to core analysis directory
cd /scratch/sahlab/kathie/IBD_plant_modular/core_analysis/host_response

# Run the host response signatures analysis
echo "Running host response signatures analysis..."
echo "Data path: $1"
echo "Metadata: $2"
echo "Output: $3"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "Expected outputs: 13 files"
echo "=========================================="

python host_response_signatures.py \
    --data-path "$1" \
    --metadata "$2" \
    --output-dir "$3" \
    --job-id "$SLURM_JOB_ID"

echo "=========================================="
echo "Host response analysis completed: $(date)"

# Validate 13 output files were created
OUTPUT_COUNT=$(find "$3" -name "*.json" -o -name "*.tsv" -o -name "*.png" -o -name "*.gml" -o -name "*.log" -o -name "*.txt" | wc -l)
echo "Output files created: $OUTPUT_COUNT (expected: 13)"

if [ "$OUTPUT_COUNT" -eq 13 ]; then
    echo "✅ All 13 expected output files created"
else
    echo "⚠️  Expected 13 files, found $OUTPUT_COUNT"
fi