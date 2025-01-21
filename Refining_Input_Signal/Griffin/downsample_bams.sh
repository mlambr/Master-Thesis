#!/bin/bash
# Load necessary modules
#enable_modules
#module load picard  # Load the Picard module, if needed

# Define the BAM files and their corresponding original coverages and P values
declare -A samples
samples=(
    ["SeCT-26_t1.sortByCoord.bam"]="26.9723 0.37 0.11 0.033 0.01 0.0038"
    ["SeCT-26_t2.sortByCoord.bam"]="26.0296 0.37 0.11 0.033 0.01 0.0038"
    ["SeCT-26_t3.sortByCoord.bam"]="28.2068 0.37 0.11 0.033 0.01 0.0038"
    ["SeCT-58.sortByCoord.bam"]="14.3002 0.56 0.16 0.045 0.014 0.0038"
    ["SeCT-82.sortByCoord.bam"]="25.8346 0.37 0.11 0.033 0.01 0.0038"
)

# Function to downsample BAM files
downsample() {
    local input="$1"
    local output="$2"
    local p="$3"
    
    echo "Downsampling $input to $p coverage"
    java -jar $EBROOTPICARD/picard.jar DownsampleSam \
        I="$input" \
        O="$output" \
        P="$p"
    echo "Created downsampled file: $output"
}

# Loop through each sample to downsample to specific coverages
for sample in "${!samples[@]}"; do
    # Read original coverage and corresponding P values from the array
    read coverage p10 p3 p1 p03 p01 <<< "${samples[$sample]}"
    echo "Processing $sample with original coverage $coverage"

    # Run downsample commands in parallel based on coverage
    if (( $(echo "$coverage >= 10" | bc -l) )); then
        downsample "$sample" "${sample%.bam}_10x.bam" "$p10" &  # 10x coverage
    fi

    if (( $(echo "$coverage >= 3" | bc -l) )); then
        downsample "$sample" "${sample%.bam}_3x.bam" "$p3" &  # 3x coverage
    fi

    if (( $(echo "$coverage >= 1" | bc -l) )); then
        downsample "$sample" "${sample%.bam}_1x.bam" "$p1" &  # 1x coverage
    fi

    if (( $(echo "$coverage >= 0.3" | bc -l) )); then
        downsample "$sample" "${sample%.bam}_0.3x.bam" "$p03" &  # 0.3x coverage
    fi

    # Always downsample to 0.1x
    downsample "$sample" "${sample%.bam}_0.1x.bam" "$p01" &  # 0.1x coverage
done

# Wait for all background jobs to finish
wait
echo "All downsampling tasks completed."
