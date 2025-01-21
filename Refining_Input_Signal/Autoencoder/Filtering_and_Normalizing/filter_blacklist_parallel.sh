#!/bin/bash
#SBATCH --job-name=filter_blacklist_parallel
#SBATCH --output=filter_blacklist_parallel_%j.log  # Log output file
#SBATCH --ntasks=1                           # Number of tasks (we'll handle parallelization within the script)
#SBATCH --cpus-per-task=8                    # Number of CPUs per task (adjust based on your system)
#SBATCH --mem=32G                            # Memory per node (adjust based on your file sizes and needs)

# Check if input directory is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <input_directory>"
    exit 1
fi

# Define input and output directories
input_dir="$1"
output_dir="${input_dir}/without_blacklisted"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Define the blacklist file
blacklist_file="/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/blacklisted/blacklisted_regions_hg38.bed"

# Process all tsv.gz files in the input directory
# Use GNU parallel to process files in parallel
find "$input_dir" -maxdepth 1 -name "*.tsv.gz" | parallel -j 8 --bar '
    file="{}"
    filename=$(basename "$file")
    echo "$filename"
    output_file="${output_dir}/${filename}"
    echo "$output_file"

    # Filter the file based on the blacklist
    zcat "$file" | awk -v blacklist_file="$blacklist_file" "
    BEGIN {
        # Read the blacklist file and store the blacklisted regions in an array
        while (getline < blacklist_file) {
            # Store blacklisted regions as start:end format
            blacklisted[\$1] = blacklisted[\$1] \",\" \$2 \"-\" \$3
        }
    }
    # For each line in the tsv.gz file, check if it is in the blacklist
    {
        chrom = \$1
        pos = \$2
        blacklisted_region = 0

        # Check if the position falls within any of the blacklisted ranges
        split(blacklisted[chrom], ranges, \",\")
        for (i in ranges) {
            split(ranges[i], range, \"-\")
            if (pos >= range[1] && pos <= range[2]) {
                blacklisted_region = 1
                break
            }
        }

        # If the position is not in any blacklisted region, print the line
        if (!blacklisted_region) {
            print \$0
        }
    }
    " > "$output_file"
'

echo "Filtering complete. Processed files are in $output_dir"
