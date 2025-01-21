import os
import gzip
import numpy as np
from Bio import SeqIO


def creating_data():
    # Output file and input file paths
    output_file_path = '/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/Data/chromosome8.txt'
    fullysampled_folder_path = '/cluster/work/medinfmk/cfDNA-Snyder/Snyder/maria/sorted_by_chromosome_GC_corrected_final_final/done/without_blacklisted'
    undersampled_folder_paths = [
        '/cluster/work/medinfmk/cfDNA-Snyder/Snyder/maria/BH01_hg38_5x_by_chromosome_GC_corrected_final_final/done/without_blacklisted'
    ]
    nucleotide_file_path = "/cluster/work/medinfmk/cfDNA-Snyder/Snyder/maria/hg38_maria.fa"
    #chromosomes = [1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 'X', 'Y']
    chromosomes = [8]
    chunks_length = 25_000

    # Preload the nucleotide data into memory for faster access
    nucleotide_data = {}
    for seq_record in SeqIO.parse(nucleotide_file_path, "fasta"):
        nucleotide_data[seq_record.name] = seq_record.seq

    with open(output_file_path, 'a') as output_file:
        for chromosome in chromosomes:
            chromosome = f'chr{chromosome}'
            print(f'Processing chromosome: {chromosome}')
            
            # Load undersampled data
            undersampled_data = {}
            for undersampled_folder_path in undersampled_folder_paths:
                undersampling_coeff = '5x'
                undersampled_file_path = f"{undersampled_folder_path}/{chromosome}_counts.tsv.gz"
                print(f'Loading undersampled data: {undersampled_file_path}')
                with gzip.open(undersampled_file_path, 'r') as undersampled_file:
                    undersampled_data[f'undersampled_{undersampling_coeff}'] = np.genfromtxt(undersampled_file, delimiter='\t')[:, 4]

            # Load fullysampled data
            fullysampled_file_path = f"{fullysampled_folder_path}/{chromosome}_counts.tsv.gz"
            with gzip.open(fullysampled_file_path, 'r') as fullysampled_file:
                fullysampled = np.genfromtxt(fullysampled_file, delimiter='\t')
                length = len(fullysampled)

                # Process data in chunks
                for i in range(0, length, chunks_length):
                    print(f'Chunk: {i*chunks_length} - {(i+1)*chunks_length - 1}')
                    line = ''
                    start = int(fullysampled[i, 1])
                    end = int(fullysampled[min(i + chunks_length, length) - 1, 1])
                    score_fully = fullysampled[i:i + chunks_length, 4]
                    #score_fully =  f'{",".join(map(str, score_fully))}'
                    # Start constructing the output line
                    line += f'{chromosome}&{start}&{end}&{score_fully.tolist()}&'

                    # Add undersampled data
                    for key, val in undersampled_data.items():
                        undersampled_score = val[i:i + chunks_length]
                        #undersampled_score = f'{",".join(map(str, undersampled_score))}'
                        line += f'{undersampled_score.tolist()}&'

                    # Add nucleotide data
                    sequence = nucleotide_data[chromosome][i:i + chunks_length]
                    # Convert sequence to a NumPy array of characters
                    sequence_array = np.array(list(sequence))
                    nucleotide_mapping = {'A': 0, 'a': 0, 'C': 1, 'c': 1, 'G': 2, 'g': 2, 'T': 3, 't': 3, 'N': 4, 'n': 4}

                    # Use a vectorized approach to map nucleotides to numeric values
                    mapping_func = np.vectorize(nucleotide_mapping.get)
                    sequence_mapped = mapping_func(sequence_array)
                    # sequence_mapped = [nucleotide_mapping[nuc] for nuc in sequence]
                    line += f'{sequence_mapped.tolist()}&\n'

                    # Write the line to the output file
                    output_file.write(line)
                    output_file.flush()


creating_data()
