import os
import torch
import gzip
import pandas as pd
from torch.utils.data import DataLoader
from torch.utils.data import Dataset
import torch.nn.functional as F
# from Model_corrected import LitAutoEncoder, custom_collate_fn

import torch.nn as nn
import pytorch_lightning as pl

class LitAutoEncoder(pl.LightningModule):
    def __init__(self, input_size=25_000, bn=160*5, k_size1=16, k_size2=4, k_size3=5, dil_1=2, dil_2=9):
        super().__init__()

        L_out1 = input_size - dil_1 * (k_size1 - 1)
        L_out2 = L_out1 - dil_2 * (k_size2 - 1)

        self.encoder = nn.Sequential(
            nn.Conv1d(2, 8, kernel_size=k_size1, dilation=dil_1),
            nn.ReLU(),
            nn.BatchNorm1d(8),
            nn.Conv1d(8, 16, kernel_size=k_size2, dilation=dil_2),
            nn.ReLU(),
            nn.Flatten(),
            nn.Dropout(0.5),
            nn.Linear(16 * L_out2, bn)
        )

        self.decoder = nn.Sequential(
            nn.Linear(bn, 16 * L_out2),
            nn.ReLU(),
            nn.Unflatten(1, (16, L_out2)),
            nn.ConvTranspose1d(16, 8, kernel_size=k_size2, dilation=dil_2),
            nn.BatchNorm1d(8),
            nn.ReLU(),
            nn.ConvTranspose1d(8, 1, kernel_size=k_size1, dilation=dil_1),
        )

    def forward(self, x1, x2):
        x = torch.stack((x1, x2), dim=1)
        x = self.encoder(x)
        x = self.decoder(x)
        return x.squeeze(1)
    
# Custom collate function to handle None values and padding
def custom_collate_fn(batch):
    batch = [b for b in batch if b is not None]
    # Pad each batch to the maximum sequence length found in the batch
    max_length = max([len(item['fullysampled']) for item in batch])  # Find max length in batch
    for item in batch:
        item['fullysampled'] = F.pad(item['fullysampled'], (0, max_length - len(item['fullysampled'])), value=0)
        item['undersampled'] = F.pad(item['undersampled'], (0, max_length - len(item['undersampled'])), value=0)
        item['nucleotide'] = F.pad(item['nucleotide'], (0, max_length - len(item['nucleotide'])), value=0)

    return torch.utils.data.dataloader.default_collate(batch)
# Dataset Class
class DataLoading(Dataset):
    def __init__(self, txt_file_path, max_length=25000):
        with open(txt_file_path, 'r') as file:
            self.data_txt = file.read().splitlines()
        self.max_length = max_length  # Set the maximum length for padding

    def __getitem__(self, idx):
        try:
            data = self.data_txt[idx]
            parts = data.split('&')
            chromosome = parts[0].strip()
            start = int(parts[1].strip())
            end = int(parts[2].strip())
            fullysampled = torch.tensor(eval(parts[3].strip()), dtype=torch.float32)
            undersampled = torch.tensor(eval(parts[4].strip()), dtype=torch.float32)
            nucleotide = torch.tensor(eval(parts[5].strip()), dtype=torch.float32)
        except Exception as e:
            print(f"Error processing line {idx}: {e}")
            return None  # Return None for invalid data

        # Pad sequences to the fixed max length
        fullysampled = F.pad(fullysampled, (0, self.max_length - len(fullysampled)), value=0)
        undersampled = F.pad(undersampled, (0, self.max_length - len(undersampled)), value=0)
        nucleotide = F.pad(nucleotide, (0, self.max_length - len(nucleotide)), value=0)

        return {
            'chromosome': chromosome,
            'start': start,
            'end': end,
            'fullysampled': fullysampled,
            'undersampled': undersampled,
            'nucleotide': nucleotide
        }

    def __len__(self):
        return len(self.data_txt)

# Path to the checkpoint
#model_checkpoint_path = '/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/checkpoint/model-epoch=39-val_loss=2.70.ckpt'  
#model_checkpoint_path = '/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/checkpoint/model-epoch=49-val_loss=0.00.ckpt'
model_checkpoint_path = '/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/checkpoint/model-epoch=29-val_loss=0.00.ckpt'

# Paths to test data and output files
test_data_path = '/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/Data/testing_bigger.txt'  # Update with your test data file path
predictions_output_path = "/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/bigger_predictions.tsv.gz"
ground_truth_output_path = "/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/bigger_ground_truth.tsv.gz"

# Load the trained model
print('Model about to be loaded')
model = LitAutoEncoder.load_from_checkpoint(model_checkpoint_path)
print('Model was loaded')
model.eval()  # Set model to evaluation mode
model.to("cuda" if torch.cuda.is_available() else "cpu")

# Dataset Class for lazy loading of batches
class DataLoading(Dataset):
    def __init__(self, txt_file_path, max_length=25000):
        self.txt_file_path = txt_file_path
        self.max_length = max_length  # Set the maximum length for padding

    def __getitem__(self, idx):
        with open(self.txt_file_path, 'r') as file:
            data = file.readlines()[idx]
        parts = data.split('&')
        chromosome = parts[0].strip()
        start = int(parts[1].strip())
        end = int(parts[2].strip())
        fullysampled = torch.tensor(eval(parts[3].strip()), dtype=torch.float32)
        undersampled = torch.tensor(eval(parts[4].strip()), dtype=torch.float32)
        nucleotide = torch.tensor(eval(parts[5].strip()), dtype=torch.float32)

        # Pad sequences to the fixed max length
        fullysampled = F.pad(fullysampled, (0, self.max_length - len(fullysampled)), value=0)
        undersampled = F.pad(undersampled, (0, self.max_length - len(undersampled)), value=0)
        nucleotide = F.pad(nucleotide, (0, self.max_length - len(nucleotide)), value=0)

        return {
            'chromosome': chromosome,
            'start': start,
            'end': end,
            'fullysampled': fullysampled,
            'undersampled': undersampled,
            'nucleotide': nucleotide
        }

    def __len__(self):
        with open(self.txt_file_path, 'r') as file:
            return len(file.readlines())

# Helper function to save data in the required format
#def save_tsv_gz(output_path, data, mode='w'):
 #   with gzip.open(output_path, mode+'t') as f:  # 't' mode for text
  #      for row in data:
   #         f.write(("\t".join(map(str, row)) + "\n").encode('utf-8'))  # Encode string to bytes

def save_tsv_gz(output_path, data, mode='w'):
    with gzip.open(output_path, mode + 't') as f:  # 't' mode for text
        for row in data:
            f.write("\t".join(map(str, row)) + "\n")  # Write directly as a string

# Prepare to load test data
test_dataset = DataLoading(test_data_path)
test_loader = DataLoader(test_dataset, batch_size=1, shuffle=False, collate_fn=custom_collate_fn)

print('data loaded')

# Prepare to store data for the current batch
device = "cuda" if torch.cuda.is_available() else "cpu"

# Iterate over the test data and process each batch
for batch_idx, batch in enumerate(test_loader):
    # Extract data from the batch
    chrom = batch['chromosome'][0]
    start = batch['start'][0]
    end = batch['end'][0]
    positions = range(start, end)
    yu_data = batch['undersampled'].to(device)
    sequence = batch['nucleotide'].to(device)
    yf_data = batch['fullysampled'][0].cpu().numpy()

    # Get predictions from the model
    with torch.no_grad():
        y_pred = model(yu_data, sequence).squeeze(0).cpu().numpy()

    # Prepare the data for the current batch
    predictions_data = []
    ground_truth_data = []
    for pos, gt, pred in zip(positions, yf_data, y_pred):
        predictions_data.append([chrom, pos, "", pred])
        ground_truth_data.append([chrom, pos, "", gt])

    # Save predictions and ground truth to separate .tsv.gz files in append mode
    save_tsv_gz(predictions_output_path, predictions_data, mode='a')
    save_tsv_gz(ground_truth_output_path, ground_truth_data, mode='a')

    # Print progress every 100 batches
    if batch_idx % 100 == 0:
        print(f"Processed {batch_idx} batches...")

print(f"Predictions saved to {predictions_output_path}")
print(f"Ground truth saved to {ground_truth_output_path}")
