#!/bin/bash
# Shared environment setup based on proven working configuration
# From: htcf_hierarchical_temporal_v3.sbatch (job 25178500 - successful run)
# Usage: source htcf_env_setup.sh MODULE_NAME

MODULE_NAME="${1:-default}"

echo "Setting up HTCF environment for module: $MODULE_NAME"

# Save any existing arguments (for scripts that have them)
# Note: We shift off the MODULE_NAME first
shift
SAVED_ARGS=("$@")

# Clear positional parameters before conda activation (CRITICAL!)
set --

# Activate conda (same as working scripts)
echo "Activating conda environment from /ref/sahlab/software/anaconda3..."
source /ref/sahlab/software/anaconda3/bin/activate

# Restore arguments after conda activation
set -- "${SAVED_ARGS[@]}"

# Set job-specific Python environment to avoid conflicts
export PYTHONUSERBASE="/scratch/sahlab/kathie/python_${MODULE_NAME}_${SLURM_JOB_ID}"

# Verify Python location and version
echo "Python location: $(which python)"
echo "Python version: $(python --version)"

# Install base requirements (what ALL modules need)
echo "Ensuring base packages are installed..."
python -m pip install --user --quiet pandas numpy scipy > /dev/null 2>&1

# Module-specific package installations
case "$MODULE_NAME" in
    "temporal")
        echo "Installing temporal analysis packages..."
        python -m pip install --user --quiet scikit-learn matplotlib seaborn networkx statsmodels dtaidistance > /dev/null 2>&1
        ;;
    "host_response")
        echo "Installing host response analysis packages..."
        python -m pip install --user --quiet scikit-learn matplotlib seaborn networkx > /dev/null 2>&1
        ;;
    "format_converter")
        echo "Installing format converter packages..."
        # Base packages are sufficient
        ;;
    "pipeline")
        echo "Installing pipeline packages..."
        python -m pip install --user --quiet psutil openpyxl > /dev/null 2>&1
        ;;
    *)
        echo "Using base packages only for module: $MODULE_NAME"
        ;;
esac

echo "Environment ready for module: $MODULE_NAME"
echo "PYTHONUSERBASE: $PYTHONUSERBASE"
echo "----------------------------------------"