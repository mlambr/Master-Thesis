#!/bin/bash
#SBATCH --job-name=downsample_bam    # Job name
#SBATCH --output=downsample_bam.out # Standard output log
#SBATCH --error=downsample_bam.err  # Standard error log
#SBATCH --ntasks=1                  # Number of tasks
#SBATCH --cpus-per-task=4           # Number of CPU cores
#SBATCH --mem=16G                   # Memory per node
#SBATCH --time=90:00:00             # Time limit (hh:mm:ss)

# Load necessary modules
#module load picard

# Define the BAM file and its original coverage
declare -A samples
samples=(
    ["../BH01_hg38_sorted.bam"]="73.1518"  # Original coverage
)

# Target coverage
target_coverage=5.0

# Function to downsample BAM files
downsample() {
    local input="$1"
    local output="$2"
    local p="$3"

    echo "Downsampling $input to proportion $p"
    java -jar $EBROOTPICARD/picard.jar DownsampleSam \
        I="$input" \
        O="$output" \
        P="$p"
    echo "Created downsampled file: $output"
}

# Loop through the samples array and process the specified BAM file
for sample in "${!samples[@]}"; do
    # Read original coverage from the array
    original_coverage="${samples[$sample]}"
    echo "Processing $sample with original coverage $original_coverage"

    # Calculate the proportion to downsample to 5x
    proportion=$(awk -v target="$target_coverage" -v original="$original_coverage" 'BEGIN {print target / original}')
    echo "Calculated proportion for $sample: $proportion"

    # Downsample to target coverage
    downsample "$sample" "${sample%.bam}_5x.bam" "$proportion"
done

echo "Downsampling completed."
