#!/bin/bash
#SBATCH --job-name=gc_correction
#SBATCH --output=gc_correction_%j.out
#SBATCH --error=gc_correction_%j.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G

# Check if the input BAM file, GC bias file, and reference genome are provided
if [ $# -lt 3 ]; then
    echo "Usage: $0 <input.bam> <GC_bias.tsv> <reference_genome.fa> [min_quality]"
    exit 1
fi

# Input BAM file, GC bias file, and reference genome
input_bam="$1"
gc_bias_file="$2"
reference_genome_file="$3"

# Minimum mapping quality (default: 30)
min_quality=${4:-30}

# Output directory
output_dir="${input_bam%.bam}_by_chromosome_GC_corrected_final_final"
mkdir -p "$output_dir"

# Debugging: Confirm output directory creation
if [ ! -d "$output_dir" ]; then
    echo "Error: Failed to create output directory $output_dir"
    exit 1
else
    echo "Output directory: $output_dir created successfully."
fi

# Create temporary file for GC bias lookup
temp_gc_bias_file=$(mktemp)
awk '{print $1 "-" $2, $7}' "$gc_bias_file" > "$temp_gc_bias_file"

# Function to process a single chromosome
process_chromosome() {
    chrom=$1
    min_quality=$2
    echo "Processing chromosome: $chrom with min_quality: $min_quality"

    # Fetch the reference sequence for the chromosome
    region_seq_file=$(mktemp)
    samtools faidx "$reference_genome_file" "$chrom" > "$region_seq_file"

    # Debugging: Confirm reference sequence fetch
    if [ ! -s "$region_seq_file" ]; then
        echo "Error: Failed to fetch reference sequence for $chrom"
        rm "$region_seq_file"
        exit 1
    fi

    intermediate_file=$(mktemp)

    # Process each fragment and calculate GC correction
    samtools view -q "$min_quality" -f 0x2 -F 0x410 "$input_bam" "$chrom" | \
    awk -v chrom="$chrom" -v gc_bias_file="$temp_gc_bias_file" -v ref_seq_file="$region_seq_file" '
    BEGIN {
        # Load GC bias values into an array
        while ((getline line < gc_bias_file) > 0) {
            split(line, fields, " ");
            gc_bias[fields[1]] = fields[2];
        }
        close(gc_bias_file);

        # Load the reference sequence into a single string
        while ((getline line < ref_seq_file) > 0) {
            if (line ~ /^>/) continue; # Skip FASTA headers
            ref_seq = ref_seq line;
        }
        close(ref_seq_file);
    }
    {
        fragment_length = $9;
        if (fragment_length >= 120 && fragment_length <= 200) {
            # Calculate fragment positions
            fragment_start = $4 - 1; # 0-based index
            fragment_end = fragment_start + fragment_length - 1;

            # Skip invalid indices
            if (fragment_start < 0 || fragment_end >= length(ref_seq)) next;

            # Extract reference sequence
            fragment_seq = substr(ref_seq, fragment_start + 1, fragment_length);

            # Count GC bases in the reference sequence
            gc_count = gsub(/[GCgc]/, "", fragment_seq);
            # Count ambiguous bases and add random 0 or 1 for each
            ambiguous_count = gsub(/[NRYKMBHDVnrykmbhdv]/, "", fragment_seq);
            gc_count += int(rand() * 2 * ambiguous_count); # Randomly add 0 or 1 for each ambiguous base

            # Create the key with fragment length and GC count
            key = int(fragment_length) "-" int(gc_count);

            # Check GC bias and apply only if above threshold
            if (key in gc_bias && gc_bias[key] >= 0.05) {
                gc_correction = 1.0 / gc_bias[key];

                # Calculate midpoint
                midpoint = int(($4 + $4 + fragment_length) / 2);

                # Print to intermediate file (midpoint, uncorrected count, gc-corrected count)
                print midpoint, 1, gc_correction;
            }
        }
    }' | \
    sort -k1,1n | \
    awk '{
        # Accumulate counts for each position
        if ($1 == last_pos) {
            uncorrected_count += $2;
            corrected_count += $3;
        } else {
            if (NR > 1) {
                print last_pos, uncorrected_count, corrected_count;
            }
            last_pos = $1;
            uncorrected_count = $2;
            corrected_count = $3;
        }
    }
    END {
        if (NR > 0) {
            print last_pos, uncorrected_count, corrected_count;
        }
    }' > "$intermediate_file"

    # Fill missing positions with zeros
    awk -v chrom="$chrom" 'BEGIN {OFS="\t"; pos = 1} {
        if (NR == 1) {
            while (pos < $1) {
                print chrom, pos, 0, 0;
                pos++;
            }
        }
        while (pos < $1) {
            print chrom, pos, 0, 0;
            pos++;
        }
        print chrom, $1, $2, $3;
        pos = $1 + 1;
    } END {
        while (pos <= end) {
            print chrom, pos, 0, 0;
            pos++;
        }
    }' "$intermediate_file" | gzip > "${output_dir}/${chrom}_counts.tsv.gz"

    # Clean up temporary files
    rm "$intermediate_file"
    rm "$region_seq_file"

    echo "Finished processing chromosome: $chrom"
}

export -f process_chromosome

# Export variables to make them available in the parallel environment
export input_bam temp_gc_bias_file reference_genome_file output_dir

# Run processing in parallel for each chromosome
samtools idxstats "$input_bam" | awk '$3 > 0 {print $1}' | \
parallel -j 24 process_chromosome {} "$min_quality"

# Cleanup temporary GC bias file
rm "$temp_gc_bias_file"

echo "Processing completed. Output written to $output_dir"
