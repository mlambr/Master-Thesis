#!/bin/bash
#SBATCH --job-name=process_chromosomes   # Job name
#SBATCH --output=process_chromosomes.out # Standard output log
#SBATCH --error=process_chromosomes.err  # Standard error log
#SBATCH --ntasks=1                       # Number of tasks (1 for GNU Parallel)
#SBATCH --cpus-per-task=6                # Number of CPU cores per task
#SBATCH --mem=100G                        # Memory per node

# Check if a directory is provided as an argument, otherwise use the current directory
if [ -z "$1" ]; then
  input_dir="."
else
  input_dir="$1"
fi

# Find all .tsv.gz files in the specified directory (without searching subdirectories)
chromosome_files=$(find "$input_dir" -maxdepth 1 -type f -name "*.tsv.gz")

# Load GNU Parallel module if required (uncomment the line below if your system uses modules)
# module load parallel

# Run the Python script in parallel for each file
echo "$chromosome_files" | parallel -j 24 python3.8 process_chromosome.py
