#!/bin/bash
#SBATCH --job-name=temporal_analysis
#SBATCH --output=/scratch/sahlab/kathie/IBD_plant_modular/results/temporal_analysis_%j.out
#SBATCH --error=/scratch/sahlab/kathie/IBD_plant_modular/results/temporal_analysis_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --account=sahlab

# Input arguments
DATA_PATH="$1"
METADATA_FILE="$2"
OUTPUT_DIR="$3"

# Source the shared environment setup (handles all conda activation)
source /scratch/sahlab/kathie/IBD_plant_modular/config/htcf_env_setup.sh temporal "$DATA_PATH" "$METADATA_FILE" "$OUTPUT_DIR"

# Change to core analysis directory
cd /scratch/sahlab/kathie/IBD_plant_modular/core_analysis/temporal

# Run the proven working temporal analysis
echo "Running temporal ecological fingerprinting analysis..."
echo "Data path: $1"
echo "Metadata: $2" 
echo "Output: $3"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "=========================================="

python temporal_ecological_fingerprinting.py \
    --data-path "$1" \
    --metadata "$2" \
    --output-dir "$3" \
    --job-id "$SLURM_JOB_ID"

echo "=========================================="
echo "Temporal analysis completed: $(date)"