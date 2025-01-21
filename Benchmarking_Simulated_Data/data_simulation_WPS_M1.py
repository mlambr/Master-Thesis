#!/usr/bin/env python
# coding: utf-8

# In[38]:


import numpy as np
import random
from scipy.signal import savgol_filter
from whittaker_eilers import WhittakerSmoother
import pandas as pd
import gzip


# In[39]:


def save_data_as_tsv_gz(array, filename):
    # Create DataFrame
    df = pd.DataFrame({
        'Column1': [1] * len(array),
        'Column2': range(1, len(array) + 1),
        'Column3': [0] * len(array),
        'Column4': [0] * len(array),
        'Column5': array
    })

    # Save DataFrame as gzipped TSV
    with gzip.open(filename, 'wt', encoding='utf-8') as f:
        df.to_csv(f, sep='\t', index=False, header=False)


# In[40]:


def get_nucleosome_positions_phase_shift(length, num_regular_peaks, num_phase_shifted_peaks, min_interpeak_distance, max_interpeak_distance):    
    regular_peak_positions = np.linspace(0, length, num_regular_peaks, dtype=int)
    
    phase_shift = 70
    shifted_peak_positions = regular_peak_positions[0:num_phase_shifted_peaks]+phase_shift

    nucleosome_positions = np.concatenate((regular_peak_positions, shifted_peak_positions))
    nucleosome_positions = np.sort(nucleosome_positions)
    mask_regular = np.in1d(nucleosome_positions, regular_peak_positions)
    mask_shifted = ~mask_regular

    return nucleosome_positions, mask_regular, mask_shifted


# In[41]:


def get_nucleosome_map(nucleosome_position, mask_regular, mask_shifted,locus_length=20_000):
    nucleosome_position +=100
    array = np.full(locus_length+201, 0)

    # Extend the mask to include positions surrounding the peaks (-83 to +83)
    start_indices = nucleosome_position[mask_regular] - 50
    end_indices = nucleosome_position[mask_regular]+ 51
    mask_surrounding = np.zeros_like(array, dtype=bool)
    for start, end in zip(start_indices, end_indices):
            mask_surrounding[start:end] = True
    values_regular = np.arange(50, 80)
    repeated_values_regular = np.tile(values_regular, len(nucleosome_position[mask_regular]) // len(values_regular) + 1)[:len(nucleosome_position[mask_regular])]
    np.random.shuffle(repeated_values_regular)
    repeated_values_regular = np.repeat(repeated_values_regular, 101)
    array[mask_surrounding] = repeated_values_regular[:len(array[mask_surrounding])]

    values_shifted = np.arange(20, 50)
    array2 = np.full(locus_length+201, 0)

    # Extend the mask to include positions surrounding the peaks (-83 to +83)
    start_indices = nucleosome_position[mask_shifted] - 50
    end_indices = nucleosome_position[mask_shifted]+ 51
    mask_surrounding = np.zeros_like(array, dtype=bool)
    for start, end in zip(start_indices, end_indices):
            mask_surrounding[start:end] = True
    repeated_values_shifted = np.tile(values_shifted, len(nucleosome_position[mask_shifted]) // len(values_shifted) + 1)[:len(nucleosome_position[mask_shifted])]
    np.random.shuffle(repeated_values_shifted)
    repeated_values_shifted = np.repeat(repeated_values_shifted, 101)
    print("len", len( array2[mask_surrounding]))
    print("repeated_values_shifted", len(repeated_values_shifted))
    array2[mask_surrounding] = repeated_values_shifted[:len(array2[mask_surrounding])]
    
    final_array = array + array2
    return final_array


# In[42]:


def get_data_phase_shift(number,length, num_regular_peaks, num_phase_shifted_peaks, min_interpeak_distance, max_interpeak_distance):
    nucleosome_positions, regular, shifted = get_nucleosome_positions_phase_shift(length, num_regular_peaks, num_phase_shifted_peaks, min_interpeak_distance, max_interpeak_distance)
    nucleosome_map = get_nucleosome_map(nucleosome_positions,regular, shifted)
    noise = np.random.normal(loc=0, scale=100, size=len(nucleosome_map))
    nucleosome_map_noisy = nucleosome_map + noise 

    smoothed_array = savgol_filter(nucleosome_map_noisy, 100, 1)
    noise = np.random.normal(loc=0, scale=2, size=len(smoothed_array))
    noisy_data = smoothed_array + noise

    link = "/mnt/DATA1/nucleosome/wps/BH01/Monte_carlo/test_all/test_many_times/smoothed_array_phase_shift" + str(number) + ".tsv.gz"
    link2 = "/mnt/DATA1/nucleosome/wps/BH01/Monte_carlo/test_all/test_many_times/nucleosome_positions_phase_shift" + str(number) + ".tsv"
    np.savetxt(link2, nucleosome_positions, delimiter='\t', fmt='%d')
    save_data_as_tsv_gz(noisy_data,link )
    


# In[43]:


def get_nucleosome_map_only_regular(nucleosome_position, mask_regular, mask_shifted,locus_length=20_000):
    # Create an array of length 30,001 filled with -50
    nucleosome_position +=100
    array = np.full(locus_length+201, 0)
    start_indices = nucleosome_position[mask_regular] - 50
    end_indices = nucleosome_position[mask_regular]+ 51
    mask_surrounding = np.zeros_like(array, dtype=bool)
    for start, end in zip(start_indices, end_indices):
            mask_surrounding[start:end] = True
    values_regular = np.arange(50, 80)
    repeated_values_regular = np.tile(values_regular, len(nucleosome_position[mask_regular]) // len(values_regular) + 1)[:len(nucleosome_position[mask_regular])]
    np.random.shuffle(repeated_values_regular)
    repeated_values_regular = np.repeat(repeated_values_regular, 101)
    array[mask_surrounding] = repeated_values_regular[:len(array[mask_surrounding])]
    return array


# In[44]:


def get_data_phase_only_regular(number,length, num_regular_peaks, min_interpeak_distance, max_interpeak_distance):
    nucleosome_positions, regular, shifted = get_nucleosome_positions_phase_shift(length, num_regular_peaks, 0, min_interpeak_distance, max_interpeak_distance)
    nucleosome_map = get_nucleosome_map_only_regular(nucleosome_positions,regular, shifted)
    noise = np.random.normal(loc=0, scale=100, size=len(nucleosome_map))

    nucleosome_map_noisy = nucleosome_map + noise 
    smoothed_array = savgol_filter(nucleosome_map_noisy, 100, 1)
    noise = np.random.normal(loc=0, scale=2, size=len(smoothed_array))
    noisy_data = smoothed_array + noise
    
    link = "/mnt/DATA1/nucleosome/wps/BH01/Monte_carlo/test_all/test_many_times/smoothed_array_only_regular" + str(number) + ".tsv.gz"
    link2 = "/mnt/DATA1/nucleosome/wps/BH01/Monte_carlo/test_all/test_many_times/nucleosome_positions_only_regular" + str(number) + ".tsv"
    np.savetxt(link2, nucleosome_positions, delimiter='\t', fmt='%d')
    save_data_as_tsv_gz(noisy_data,link )


# In[45]:


def generate_random_nucleosome_positions(length,num_positions):
    nucleosome_positions = []
    for _ in range(num_positions):
        position = random.randint(0, length)
        nucleosome_positions.append(position)
    return np.array(nucleosome_positions)


# In[46]:


def get_nucleosome_map_only_random(nucleosome_position):
    nucleosome_position +=100
    array = np.full(20201, 0)
    values_regular = np.arange(50, 80)
    repeated_values_regular = np.tile(values_regular, len(nucleosome_position) // len(values_regular) + 1)[:len(nucleosome_position)]
    start_indices = nucleosome_position - 50
    end_indices = nucleosome_position+ 51
    mask_surrounding = np.zeros_like(array)
    for start, end in zip(start_indices, end_indices):
        np.random.shuffle(repeated_values_regular)
        value = repeated_values_regular[5]
        mask_surrounding[start:end] += value
        
    return mask_surrounding


# In[47]:


def get_data_phase_only_random(number,length, num_peaks):
    nucleosome_positions = generate_random_nucleosome_positions(length,num_peaks)
    nucleosome_map = get_nucleosome_map_only_random(nucleosome_positions)
    noise = np.random.normal(loc=0, scale=100, size=len(nucleosome_map))
    nucleosome_map_noisy = nucleosome_map + noise 

    smoothed_array = savgol_filter(nucleosome_map_noisy, 100, 1)
    noise = np.random.normal(loc=0, scale=2, size=len(smoothed_array))
    noisy_data = smoothed_array + noise
    
    link = "/mnt/DATA1/nucleosome/wps/BH01/Monte_carlo/test_all/test_many_times/smoothed_array_only_random" + str(number) + ".tsv.gz"
    link2 = "/mnt/DATA1/nucleosome/wps/BH01/Monte_carlo/test_all/test_many_times/nucleosome_positions_only_random" + str(number) + ".tsv"
    np.savetxt(link2, nucleosome_positions, delimiter='\t', fmt='%d')
    save_data_as_tsv_gz(noisy_data,link )


# In[48]:


def get_nucleosome_positions_with_random(length, num_regular_peaks, num_phase_shifted_peaks, min_interpeak_distance, max_interpeak_distance):
    regular_peak_positions = np.linspace(0, length, num_regular_peaks, dtype=int)
    random_indices = np.random.choice(num_regular_peaks, num_phase_shifted_peaks, replace=False)
    chosen_regular_peak_positions = regular_peak_positions[random_indices]
    num_extra_shifts = num_phase_shifted_peaks * 2  # Generate twice as many to ensure enough non-zero shifts
    shifts = np.random.randint(-150, 150, num_extra_shifts)
    # Filter out 0s from the shifts
    shifts = shifts[shifts != 0]
    # Take only the required number of non-zero shifts
    shifts = shifts[:num_phase_shifted_peaks]    
    phase_shifted_peaks = chosen_regular_peak_positions + shifts
    # Ensure all phase-shifted peaks are within the range of 0 and 30000
    phase_shifted_peak_positions = phase_shifted_peaks[(phase_shifted_peaks >= 0) & (phase_shifted_peaks <= 20000)]
    # Ensure we have exactly num_phase_shifted_peaks phase-shifted peaks
    phase_shifted_peak_positions = phase_shifted_peak_positions[:num_phase_shifted_peaks]
    # Concatenate regular and phase-shifted peaks
    nucleosome_positions = np.concatenate((regular_peak_positions, phase_shifted_peak_positions))
    # If the total number of nucleosomes is less than 150, randomly add more
    num_additional_nucleosomes = 150 - len(nucleosome_positions)
    if num_additional_nucleosomes > 0:
        additional_positions = np.random.randint(0, length, num_additional_nucleosomes)
        nucleosome_positions = np.concatenate((nucleosome_positions, additional_positions))
    nucleosome_positions = np.sort(nucleosome_positions)
    mask_regular = np.in1d(nucleosome_positions, regular_peak_positions)
    mask_shifted = ~mask_regular
    return nucleosome_positions, mask_regular, mask_shifted


# In[49]:


def get_data_phase_with_random(number,length, num_regular_peaks, num_phase_shifted_peaks, min_interpeak_distance, max_interpeak_distance):
    nucleosome_positions, regular, shifted = get_nucleosome_positions_with_random(length, num_regular_peaks, num_phase_shifted_peaks, min_interpeak_distance, max_interpeak_distance)
    nucleosome_map = get_nucleosome_map(nucleosome_positions,regular, shifted)
    noise = np.random.normal(loc=0, scale=100, size=len(nucleosome_map))

    nucleosome_map_noisy = nucleosome_map + noise
    
    smoothed_array = savgol_filter(nucleosome_map_noisy, 100, 1)
    noise = np.random.normal(loc=0, scale=2, size=len(smoothed_array))
    noisy_data = smoothed_array + noise
    
    link = "/mnt/DATA1/nucleosome/wps/BH01/Monte_carlo/test_all/test_many_times/smoothed_array_with_random" + str(number) + ".tsv.gz"
    link2 = "/mnt/DATA1/nucleosome/wps/BH01/Monte_carlo/test_all/test_many_times/nucleosome_positions_with_random" + str(number) + ".tsv"
    np.savetxt(link2, nucleosome_positions, delimiter='\t', fmt='%d')
    save_data_as_tsv_gz(noisy_data,link )
    


# In[50]:


def get_nucleosome_positions_with_remove(length, num_regular_peaks, num_phase_shifted_peaks, min_interpeak_distance, max_interpeak_distance,num_remove):
    regular_peak_positions = np.linspace(0, length, num_regular_peaks, dtype=int)
    random_indices = np.random.choice(num_regular_peaks, num_phase_shifted_peaks, replace=False)
    chosen_regular_peak_positions = regular_peak_positions[random_indices]
    #randomly remove number og elements
    indices_to_remove = random.sample(range(len(regular_peak_positions)), num_remove)
    regular_peak_positions = np.delete(regular_peak_positions, indices_to_remove)
   
    num_extra_shifts = num_phase_shifted_peaks * 2  # Generate twice as many to ensure enough non-zero shifts
    shifts = np.random.randint(-150, 150, num_extra_shifts)
    # Filter out 0s from the shifts
    shifts = shifts[shifts != 0]
    # Take only the required number of non-zero shifts
    shifts = shifts[:num_phase_shifted_peaks]    
    phase_shifted_peaks = chosen_regular_peak_positions + shifts
    # Ensure all phase-shifted peaks are within the range of 0 and 30000
    phase_shifted_peak_positions = phase_shifted_peaks[(phase_shifted_peaks >= 0) & (phase_shifted_peaks <= 20000)]
    # Ensure we have exactly num_phase_shifted_peaks phase-shifted peaks
    phase_shifted_peak_positions = phase_shifted_peak_positions[:num_phase_shifted_peaks]
    # Concatenate regular and phase-shifted peaks
    nucleosome_positions = np.concatenate((regular_peak_positions, phase_shifted_peak_positions))
    # If the total number of nucleosomes is less than 150, randomly add more
    num_additional_nucleosomes = num_regular_peaks+num_phase_shifted_peaks-num_remove - len(nucleosome_positions)
    if num_additional_nucleosomes > 0:
        additional_positions = np.random.randint(0, length, num_additional_nucleosomes)
        nucleosome_positions = np.concatenate((nucleosome_positions, additional_positions))
    nucleosome_positions = np.sort(nucleosome_positions)
    mask_regular = np.in1d(nucleosome_positions, regular_peak_positions)
    mask_shifted = ~mask_regular
    return nucleosome_positions, mask_regular, mask_shifted


# In[51]:


def get_data_phase_with_remove(number,length, num_regular_peaks, num_phase_shifted_peaks, min_interpeak_distance, max_interpeak_distance,num_remove):
    nucleosome_positions, regular, shifted = get_nucleosome_positions_with_remove(length, num_regular_peaks, num_phase_shifted_peaks, min_interpeak_distance, max_interpeak_distance,num_remove)
    nucleosome_map = get_nucleosome_map(nucleosome_positions,regular, shifted)
    noise = np.random.normal(loc=0, scale=100, size=len(nucleosome_map))
 
    nucleosome_map_noisy = nucleosome_map + noise 
    
    smoothed_array = savgol_filter(nucleosome_map_noisy, 100, 1)
    noise = np.random.normal(loc=0, scale=2, size=len(smoothed_array))
    noisy_data = smoothed_array + noise

    link = "/mnt/DATA1/nucleosome/wps/BH01/Monte_carlo/test_all/test_many_times/smoothed_array_with_remove" + str(number)+str(num_phase_shifted_peaks) + ".tsv.gz"
    link2 = "/mnt/DATA1/nucleosome/wps/BH01/Monte_carlo/test_all/test_many_times/nucleosome_positions_with_remove" + str(number)+str(num_phase_shifted_peaks) + ".tsv"
    np.savetxt(link2, nucleosome_positions, delimiter='\t', fmt='%d')
    save_data_as_tsv_gz(noisy_data,link )


# In[52]:


def get_all_data(times):
    length = 20_000
    num_regular_peaks = 100
    min_interpeak_distance = 180
    max_interpeak_distance = 220
    for i in range(times):
        get_data_phase_with_remove(i,length, 60, 50, min_interpeak_distance, max_interpeak_distance,10)
        get_data_phase_with_remove(i,length, num_regular_peaks, 50, min_interpeak_distance, max_interpeak_distance,10)
        get_data_phase_only_random(i,length, 100)
        get_data_phase_with_random(i,length, 66, 34, min_interpeak_distance, max_interpeak_distance)
        get_data_phase_only_regular(i,length, num_regular_peaks, min_interpeak_distance, max_interpeak_distance)
        get_data_phase_shift(i,length, 66, 34, min_interpeak_distance, max_interpeak_distance)
        

        


# In[54]:


get_all_data(1)

