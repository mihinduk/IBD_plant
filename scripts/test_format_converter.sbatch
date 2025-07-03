#!/bin/bash
#SBATCH --job-name=test_format_converter
#SBATCH --output=/scratch/sahlab/kathie/IBD_plant_modular/test_format_converter_%j.out
#SBATCH --error=/scratch/sahlab/kathie/IBD_plant_modular/test_format_converter_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --account=sahlab

# Test data path
TEST_DATA_PATH="$1"

# Source the shared environment setup
source /scratch/sahlab/kathie/IBD_plant_modular/config/htcf_env_setup.sh format_converter "$TEST_DATA_PATH"

# TEST_DATA_PATH is now available again as first argument ($1)

# Change to the modular directory
cd /scratch/sahlab/kathie/IBD_plant_modular

# Run the format converter test
echo "🧪 Testing Format Converter on HTCF"
echo "Test data: $TEST_DATA_PATH"
echo "Job ID: $SLURM_JOB_ID"
echo "================================"

python3 scripts/test_format_converter.py "$TEST_DATA_PATH"

echo "================================"
echo "Test completed!"