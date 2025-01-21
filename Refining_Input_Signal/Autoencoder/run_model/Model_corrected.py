import os
import time
import torch
import numpy as np
import matplotlib.pyplot as plt
from torch.utils.data import Dataset, DataLoader
import torch.nn as nn
import torch.nn.functional as F
import lightning.pytorch as pl
from lightning.pytorch import Trainer
from lightning.pytorch.loggers import TensorBoardLogger
from lightning.pytorch.callbacks import ModelCheckpoint

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

# Lightning Module
class LitAutoEncoder(pl.LightningModule):
    def __init__(self, input_size=25_000, bn=160*5, k_size1=16, k_size2=4, k_size3=5, dil_1=2, dil_2=9, learning_rate=1e-3):
        super().__init__()
        self.loss = nn.MSELoss()
        self.learning_rate = learning_rate
        self.save_hyperparameters()

        L_out1 = input_size - dil_1 * (k_size1 - 1)
        L_out2 = L_out1 - dil_2 * (k_size2 - 1)
        L_out3 = (L_out2 - (k_size3 - 1) - 1) // k_size3 + 1

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

    def training_step(self, batch, batch_idx):
        yf_data, yu_data, sequence = batch['fullysampled'], batch['undersampled'], batch['nucleotide']
        y_predicted = self(yu_data, sequence)
        loss = self.loss(y_predicted, yf_data)
        self.log("train_loss", loss)

        if self.global_step % 10 == 0:
            self._log_plot(y_predicted, yf_data, "Training")

        return loss

    def validation_step(self, batch, batch_idx):
        yf_data, yu_data, sequence = batch['fullysampled'], batch['undersampled'], batch['nucleotide']
        y_predicted = self(yu_data, sequence)
        loss = self.loss(y_predicted, yf_data)
        self.log("val_loss", loss)

        if self.global_step % 10 == 0:
            self._log_plot(y_predicted, yf_data, "Validation")

    def _log_plot(self, y_predicted, yf_data, stage):
        y_pred_1000 = y_predicted[0].detach().cpu().numpy()[:1000]
        yf_data_1000 = yf_data[0].detach().cpu().numpy()[:1000]
        plt.plot(y_pred_1000, label='Predicted WPS')
        plt.plot(yf_data_1000, label='Actual WPS')
        plt.legend()
        plt.title(f'{stage} | Epoch {self.current_epoch} Step {self.global_step}')
        fig = plt.gcf()
        self.logger.experiment.add_figure(f'{stage}_WPS_Plot', fig, self.global_step)
        plt.close(fig)

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=self.learning_rate)

# Data Loading
train_dataset = DataLoading('/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/Data/training.txt')
val_dataset = DataLoading('/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/Data/validation.txt')

train_loader = DataLoader(train_dataset, batch_size=8, shuffle=True, num_workers=1, collate_fn=custom_collate_fn)
val_loader = DataLoader(val_dataset, batch_size=1, shuffle=False, collate_fn=custom_collate_fn)

# Callbacks and Trainer
checkpoint_callback = ModelCheckpoint(
    dirpath='/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/checkpoint',
    filename='model-{epoch:02d}-{val_loss:.2f}',
    save_last=True,
    every_n_epochs=10
)

logger = TensorBoardLogger(
    "/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/logs",
    name="Model1_7"
)

trainer = Trainer(max_epochs=50, logger=logger, callbacks=[checkpoint_callback])
trainer.fit(LitAutoEncoder(), train_loader, val_loader)
