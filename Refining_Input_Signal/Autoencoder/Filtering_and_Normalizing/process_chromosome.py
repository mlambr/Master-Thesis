import sys
import pandas as pd
from scipy.signal import medfilt
import gzip

def running_median(data, window_size):
    """Apply running median with a specified window size."""
    return medfilt(data, kernel_size=window_size)

def process_chromosome(file_path):
    """Process a single file: read, normalize the data, and save it."""
    # Read the tsv.gz file
    with gzip.open(file_path, 'rt') as f:
        df = pd.read_csv(f, sep='\t', header=None)
        
        # Ensure the file has at least 4 columns
        if df.shape[1] >= 4:
            # Extract the fourth column (index 3)
            data = df.iloc[:, 3].values
            
            # Normalize using the running median
            normalized_data = data - running_median(data, window_size=375)
            
            # Add the normalized data as a new column
            df['normalized'] = normalized_data
            
            # Write the updated DataFrame back to the original file
            with gzip.open(file_path, 'wt') as out_f:
                df.to_csv(out_f, sep='\t', header=False, index=False)

            print(f"Normalized data added to: {file_path}")
        else:
            print(f"Skipping file {file_path}, not enough columns.")

if __name__ == "__main__":
    # The first argument passed to the script will be the file path
    file_path = sys.argv[1]
    process_chromosome(file_path)
