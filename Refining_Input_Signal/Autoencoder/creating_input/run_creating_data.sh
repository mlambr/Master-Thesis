#!/bin/bash

# Bash script to run the creating_data.py script
# Usage: bash run_creating_data.sh

# Set job-specific parameters
#SBATCH --job-name=creating_data     # Job name
#SBATCH --output=job_output.log      # Standard output and error log
#SBATCH --ntasks=1                   # Number of tasks (single process)
#SBATCH --cpus-per-task=4            # Number of CPU cores per task
#SBATCH --mem=250G                    # Memory per node (adjust as needed)

# Define the location of your Python script
PYTHON_SCRIPT="creating_data.py"

# Run the Python script
echo "Running Python script..."
python $PYTHON_SCRIPT

echo "Job completed."
