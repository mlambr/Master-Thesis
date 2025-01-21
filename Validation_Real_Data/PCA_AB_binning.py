import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process nucleosome data with compartments and interpeak distances.")
    
    parser.add_argument("--chromosome_info", required=True, help="Path to chromosome lengths file.")
    parser.add_argument("--compartments", required=True, help="Path to compartments file.")
    parser.add_argument("--blacklist", required=True, help="Path to blacklist BED file.")
    parser.add_argument("--sample_files", required=True, nargs='+', help="Paths to sample files.")
    parser.add_argument("--output", required=True, help="Path to save output.")
    
    return parser.parse_args()

args = parse_arguments()

chromosome_information = pd.read_csv(args.chromosome_info, sep='\t', header=None)
# chromosome_information = pd.read_csv('/mnt/DATA3/SeCT/chromosome_lengths.txt', sep='\t', header = None)
chromosome_information = chromosome_information[chromosome_information[0]!='chrY']
chromosome_information = chromosome_information[chromosome_information[0]!='chrX']
chromosomes = np.array(chromosome_information[0])
chromosome_lengths = np.array(chromosome_information[1])

compartments_hg38 = pd.read_csv(args.compartments, delimiter='\t', header=None, usecols=range(4))
# compartments_hg38 = pd.read_csv("/mnt/DATA1/nucleosome/wps/GM12878_track_hg38_noheader.bed",delimiter='\t', header=None,usecols=range(4))
new_column_names2 = ['chromosome', 'start',"end", 'compartment']
compartments_hg38.columns = new_column_names2
compartments_hg38 = compartments_hg38.sort_values(by=['chromosome', 'start'])

bins = range(len(compartments_hg38))
compartments_hg38['bin'] = bins

# Create a dictionary from the lists
chromosome_lengths_dict = dict(zip(chromosomes, chromosome_lengths))

# Filter out chromosomes X and Y
compartments_hg38 = compartments_hg38[~compartments_hg38['chromosome'].isin(['chrX', 'chrY'])]

# Function to calculate global positions
def calculate_global_position(df, chromosome_lengths_dict):
    global_start = []
    global_end = []

    cumulative_length = 0
    for chromosome in df['chromosome'].unique():
        chromosome_df = df[df['chromosome'] == chromosome]
        cumulative_length = sum(chromosome_lengths_dict[chr] for chr in chromosome_lengths_dict if chr < chromosome)
        
        global_start.extend(chromosome_df['start'] + cumulative_length)
        global_end.extend(chromosome_df['end'] + cumulative_length)
    
    df['global_start'] = global_start
    df['global_end'] = global_end
    return df

# Adjust start and end positions to global positions
compartments_hg38 = calculate_global_position(compartments_hg38, chromosome_lengths_dict)


# Combine the global_start and global_end positions
global_positions = np.concatenate([compartments_hg38['global_start'].values, compartments_hg38['global_end'].values])

# Get unique positions and sort them to get the bin edges
bin_edges = np.unique(global_positions)
bin_edges = np.sort(bin_edges)

def compartment_annotation(data, compartment):
    data = data.dropna() 
    compartment = compartment.dropna()

    # Perform a merge_asof operation
    merged_df = pd.merge_asof(data.sort_values(["position","chromosome"]), compartment.sort_values(["start","chromosome"]), by='chromosome', left_on='position', right_on='start', direction='backward')

    # Filter rows where position is within the start-end range
    merged_df = merged_df[(merged_df['position'] >= merged_df['start']) & (merged_df['position'] <= merged_df['end'])]

    # Drop unnecessary columns and rearrange columns if needed
    #merged_df = merged_df[['chromosome', 'position', 'compartment']]

    return merged_df


def read_bed_file(bed_file):
    """
    Reads a BED file containing blacklisted regions.

    Args:
    bed_file (str): Path to the BED file.

    Returns:
    pd.DataFrame: DataFrame with columns ['chromosome', 'start', 'end'].
    """
    # Read the BED file with no headers and provide column names
    blacklisted_regions = pd.read_csv(bed_file, sep='\t', header=None, names=['chromosome', 'start', 'end'])
    return blacklisted_regions

blacklisted_regions = read_bed_file(args.blacklist)
# blacklisted_regions = read_bed_file('/mnt/DATA1/resources/reference_genomes/blacklist/ENCFF356LFX-hg38.bed')

def z_scale(input_dfs):
    normalized_dfs = []
    for df in input_dfs:
        mean_interpeak = df['interpeak_distance'].mean()
        std_interpeak = df['interpeak_distance'].std()
        df['z_scaled_interpeak_distance'] = (df['interpeak_distance'] - mean_interpeak) / std_interpeak
        normalized_dfs.append(df)
    return normalized_dfs

def get_mean_interpeak_distances(sample_dfs):
    mean_interpeak_distances = []
    for df in sample_dfs:
        mean_distances = df.groupby('bin')['z_scaled_interpeak_distance'].mean()
        mean_interpeak_distances.append(mean_distances)
    return mean_interpeak_distances

def get_combined_df(unique_sample_names, mean_interpeak_distances):
    combined_data = []  # Create an empty list to collect DataFrames
    
    # Loop through each Series (df) in mean_interpeak_distances
    for i, df in enumerate(mean_interpeak_distances):
        sample_name = unique_sample_names[i]  # Get the sample name
        
        # Convert the Series to a DataFrame using its index as the 'bin' column
        mean_interpeak_distances[i] = pd.DataFrame({
            'sample_name': [sample_name] * len(df),  # Repeat the sample name for all rows
            'mean_interpeak_distance': df.values,  # Extract the values from the Series
            'bin': df.index  # Use the index as the 'bin' numbers
        })
        
        # Append the resulting DataFrame to the list
        combined_data.append(mean_interpeak_distances[i])
    
    # Concatenate all the DataFrames into one combined DataFrame
    combined_df = pd.concat(combined_data, ignore_index=True)
    
    return combined_df


def filter_interpeak_distance(dfs, threshold):
    """
    Filters rows in each DataFrame where interpeak_distance is above a certain threshold.

    Args:
    dfs (list of pd.DataFrame): List of pandas DataFrames.
    threshold (float): Threshold value for interpeak_distance.

    Returns:
    list of pd.DataFrame: List of DataFrames with rows filtered based on interpeak_distance.
    """
    filtered_dfs = []
    for df in dfs:
        # Remove rows where interpeak_distance is above the threshold
        filtered_df = df[df['interpeak_distance'] <= threshold].copy()
        filtered_dfs.append(filtered_df)
    
    return filtered_dfs

sample_files_SeCT = args.sample_files
# sample_files_SeCT = ['/mnt/DATA3/nucleosome/nucleosome_centers/other_samples/SeCT-20.sortByCoord_by_chromosome/output/SeCT20_simulation_peaks_distance.bedgraph',
#                       '/mnt/DATA3/nucleosome/nucleosome_centers/other_samples/SeCT-22.sortByCoord_by_chromosome/output/SeCT22_simulation_peaks_distance.bedgraph',
#                       '/mnt/DATA3/nucleosome/nucleosome_centers/other_samples/SeCT-26_t1.sortByCoord_by_chromosome/output/SeCT26_t1_simulation_peaks_distance.bedgraph',
#                       '/mnt/DATA3/nucleosome/nucleosome_centers/other_samples/SeCT-26_t2.sortByCoord_by_chromosome/output/SeCT26_t2_simulation_peaks_distance.bedgraph',
#                       '/mnt/DATA3/nucleosome/nucleosome_centers/other_samples/SeCT-26_t3.sortByCoord_by_chromosome/output/SeCT26_t3_simulation_peaks_distance.bedgraph',
#                       '/mnt/DATA3/nucleosome/nucleosome_centers/other_samples/SeCT-47.sortByCoord_by_chromosome/output/SeCT47_simulation_peaks_distance.bedgraph',
#                       '/mnt/DATA3/nucleosome/nucleosome_centers/other_samples/SeCT-58.sortByCoord_by_chromosome/output/SeCT58_simulation_peaks_distance.bedgraph',
#                       '/mnt/DATA3/nucleosome/nucleosome_centers/other_samples/SeCT-61.sortByCoord_by_chromosome/output/SeCT61_simulation_peaks_distance.bedgraph',
#                     '/mnt/DATA3/nucleosome/nucleosome_centers/other_samples/SeCT-82.sortByCoord_by_chromosome/output/SeCT82_simulation_peaks_distance.bedgraph']  # List of sample file paths


sample_dfs_SeCT = []
for file in sample_files_SeCT:
    df = pd.read_csv(file, sep='\t', header=None, names=['chromosome','start','end', 'position', 'interpeak_distance'])
    sample_dfs_SeCT.append(df)
sample_dfs_SeCT = [df[~df['chromosome'].isin(['chrX', 'chrY'])] for df in sample_dfs_SeCT]

def convert_to_global_position(df, cumulative_lengths_mapping):
    # Add a new column with the global position
    df['global_position'] = df.apply(lambda row: row['position'] + cumulative_lengths_mapping[row['chromosome']], axis=1)
    return df

# Example chromosome lengths (replace with actual data)
chromosome_names = chromosomes  # Replace with actual chromosome names

# Create a cumulative length mapping
cumulative_lengths = np.cumsum(chromosome_lengths)
chromosome_to_cumulative = dict(zip(chromosome_names, np.concatenate(([0], cumulative_lengths[:-1]))))

for df in sample_dfs_SeCT:
    df = convert_to_global_position(df, chromosome_to_cumulative)

# Function to assign bins and filter rows
def assign_bins(df, bins_df):
    # Merge the df with bins_df based on global_position falling between global_start and global_end
    df_with_bins = pd.merge_asof(
        df.sort_values('global_position'),  # Sort the dataframe by global_position
        bins_df[['bin', 'global_start', 'global_end']].sort_values('global_start'),  # Sort bins_df by global_start
        left_on='global_position',  # The key for the DataFrame
        right_on='global_start',  # The lower bound
        direction='backward'  # We search for the nearest lower value of global_start
    )
    
    # Filter rows where global_position is within the bin range
    df_with_bins = df_with_bins[
        (df_with_bins['global_position'] >= df_with_bins['global_start']) &
        (df_with_bins['global_position'] <= df_with_bins['global_end'])
    ]
    
    # Drop unnecessary columns
    df_with_bins = df_with_bins.drop(columns=['global_start', 'global_end'])
    
    return df_with_bins

# Apply the function to each df in the list of dfs
sample_dfs_SeCT_with_bins = [assign_bins(df, compartments_hg38) for df in sample_dfs_SeCT]


# Example chromosome lengths (replace with actual data)
chromosome_names = chromosomes  # Replace with actual chromosome names

# Create a cumulative length mapping
cumulative_lengths = np.cumsum(chromosome_lengths)
chromosome_to_cumulative = dict(zip(chromosome_names, np.concatenate(([0], cumulative_lengths[:-1]))))

blacklisted_regions = blacklisted_regions[~blacklisted_regions['chromosome'].isin(['chrX', 'chrY'])]
blacklisted_regions['global_start'] = blacklisted_regions.apply(lambda row: row['start'] + chromosome_to_cumulative[row['chromosome']], axis=1)
blacklisted_regions['global_end'] = blacklisted_regions.apply(lambda row: row['end'] + chromosome_to_cumulative[row['chromosome']], axis=1)

# Function to filter out rows based on blacklisted regions without loops
def filter_blacklisted_regions_fast(dfs, blacklisted_regions):
    filtered_dfs = []

    # Convert blacklisted_regions into numpy arrays for fast operations
    blacklist_starts = blacklisted_regions['global_start'].values
    blacklist_ends = blacklisted_regions['global_end'].values

    # For each DataFrame in dfs, filter out the rows within blacklisted regions
    for df in dfs:
        positions = df['global_position'].values
        
        # Create a mask that checks if each position falls outside all blacklisted regions
        # We use broadcasting to compare positions against all ranges in blacklist simultaneously
        mask = np.ones(len(positions), dtype=bool)
        for start, end in zip(blacklist_starts, blacklist_ends):
            mask &= ~((positions >= start) & (positions <= end))
        
        # Apply the mask to filter the DataFrame
        filtered_dfs.append(df[mask])

    return filtered_dfs

# Example usage:

# Assuming `blacklisted_regions` is a DataFrame with 'global_start' and 'global_end' columns
# Assuming `dfs` is the list of DataFrames with a 'global_position' column

# Apply the filter to the list of DataFrames
without_blacklisted_sample_dfs_SeCT = filter_blacklisted_regions_fast(sample_dfs_SeCT_with_bins, blacklisted_regions)

filtered_sample_dfs_SeCT = filter_interpeak_distance(without_blacklisted_sample_dfs_SeCT, 350)
normalized_sample_dfs_SeCT = z_scale(filtered_sample_dfs_SeCT)
mean_interpeak_distances_SeCT = get_mean_interpeak_distances(normalized_sample_dfs_SeCT)
unique_sample_names_SeCT = ['SeCT-20','SeCT-22','SeCT-26_t1','SeCT-26_t2', 'SeCT-26_t3','SeCT-47','SeCT-58','SeCT-61','SeCT-82']
combined_df_SeCT = get_combined_df(unique_sample_names_SeCT, mean_interpeak_distances_SeCT)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

name = ""

# Define the sample categories
high_tumor_fraction = ['SeCT-26_t1', 'SeCT-26_t2', 'SeCT-26_t3']
benign = ['SeCT-58', 'SeCT-82']
low_tumor_fraction = ['SeCT-20', 'SeCT-22', 'SeCT-47', 'SeCT-61']

# Create a mapping of sample names to categories
sample_category = {}
for sample in high_tumor_fraction:
    sample_category[sample] = 'High Tumor Fraction'
for sample in benign:
    sample_category[sample] = 'Benign'
for sample in low_tumor_fraction:
    sample_category[sample] = 'Low Tumor Fraction'

# Pivot the DataFrame to have a row for each sample and a column for each bin
pivoted_df = combined_df_SeCT.pivot(index='sample_name', columns='bin', values='mean_interpeak_distance').fillna(0)

# Standardize the pivoted DataFrame
scaler = StandardScaler()
scaled_df = scaler.fit_transform(pivoted_df)

# Perform PCA
pca = PCA()
pca_result = pca.fit_transform(scaled_df)

# Create a DataFrame containing the PCA results
pca_df = pd.DataFrame(data=pca_result, columns=[f'PC{i+1}' for i in range(pca_result.shape[1])], index=pivoted_df.index)

# Add category information to pca_df
pca_df['category'] = pca_df.index.map(sample_category)

# Define colors for each category
colors = {
    'High Tumor Fraction': 'red',
    'Low Tumor Fraction': 'orange',
    'Benign': 'blue'
}

# Plot the PCA results
plt.figure(figsize=(4, 4), dpi=200)

# Plot each category with a different color
for category, color in colors.items():
    category_df = pca_df[pca_df['category'] == category]
    plt.scatter(category_df['PC1'], category_df['PC2'], c=color, label=category, alpha=0.5)

# title = "PCA of Mean Interpeak Distances " + name
# plt.title(title, pad=20, fontsize=15)
plt.xlabel('Principal Component 1', labelpad=15, fontsize=12)
plt.ylabel('Principal Component 2', fontsize=12)

# Add sample names
for i, sample in enumerate(pca_df.index):
    plt.text(pca_df.loc[sample, 'PC1'], pca_df.loc[sample, 'PC2'], sample, fontsize=10)

plt.grid(True)
# plt.legend(title='')
plt.legend(title='', loc='upper left', bbox_to_anchor=(1, 1))  # Position legend outside the plot
# plt.savefig('/mnt/DATA3/nucleosome/pca_plot.eps', format='eps', bbox_inches='tight')
plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.2)
ax = plt.gca()  # Get current axis

ax.spines['top'].set_visible(False)  # Hide the top spine
ax.spines['right'].set_visible(False)  # Hide the right spine
ax.spines['left'].set_linewidth(0.5)  # Customize the left spine
ax.spines['bottom'].set_linewidth(0.5)  # Customize the bottom spine

args.output
plt.savefig(args.output, format='svg',bbox_inches='tight')
# plt.savefig('/mnt/DATA3/nucleosome/images/PCA_AB_binning_fragment_centers.svg', format='svg',bbox_inches='tight')
plt.show()

