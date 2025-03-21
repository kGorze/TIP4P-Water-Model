#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os

# Define paths to data files
data_298K = "/home/konrad_guest/Documents/research/cursor/tip4p_298K/analysis/rmsd_data.dat"
data_273K = "/home/konrad_guest/Documents/research/cursor/tip4p_273K/analysis/rmsd_data.dat"
data_compressibility = "/home/konrad_guest/Documents/research/cursor/iteration_5_TIP4P/analysis/rmsd_data.dat"

# Load data
data_298K_array = np.loadtxt(data_298K, comments='#')
data_273K_array = np.loadtxt(data_273K, comments='#')
data_compressibility_array = np.loadtxt(data_compressibility, comments='#')

# Convert time from ps to ns
time_298K = data_298K_array[:, 0] / 1000  # ps to ns
rmsd_298K = data_298K_array[:, 1]

time_273K = data_273K_array[:, 0] / 1000  # ps to ns
rmsd_273K = data_273K_array[:, 1]

time_compressibility = data_compressibility_array[:, 0] / 1000  # ps to ns
rmsd_compressibility = data_compressibility_array[:, 1]

# Create the plot
plt.figure(figsize=(12, 8))

# Plot data
plt.plot(time_298K, rmsd_298K, 'b-', linewidth=2, label='TIP4P 298K')
plt.plot(time_273K, rmsd_273K, 'g-', linewidth=2, label='TIP4P 273K')
plt.plot(time_compressibility, rmsd_compressibility, 'r-', linewidth=2, label='TIP4P Compressibility')

# Calculate mean RMSD values for each dataset
mean_298K = np.mean(rmsd_298K)
mean_273K = np.mean(rmsd_273K)
mean_compressibility = np.mean(rmsd_compressibility)

# Add horizontal lines for mean values
plt.axhline(y=mean_298K, color='b', linestyle='--', alpha=0.7)
plt.axhline(y=mean_273K, color='g', linestyle='--', alpha=0.7)
plt.axhline(y=mean_compressibility, color='r', linestyle='--', alpha=0.7)

# Add annotations for mean values
plt.annotate(f'Mean: {mean_298K:.2f} Å', xy=(0.02, 0.92), xycoords='axes fraction', 
             color='b', fontsize=10, bbox=dict(facecolor='white', alpha=0.7))
plt.annotate(f'Mean: {mean_273K:.2f} Å', xy=(0.02, 0.88), xycoords='axes fraction', 
             color='g', fontsize=10, bbox=dict(facecolor='white', alpha=0.7))
plt.annotate(f'Mean: {mean_compressibility:.2f} Å', xy=(0.02, 0.84), xycoords='axes fraction', 
             color='r', fontsize=10, bbox=dict(facecolor='white', alpha=0.7))

# Add title and labels
plt.title('TIP4P Water RMSD Comparison', fontsize=16)
plt.xlabel('Simulation Time (ns)', fontsize=14)
plt.ylabel('RMSD (Å)', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(loc='best', fontsize=12)

# Add standard deviations
std_298K = np.std(rmsd_298K)
std_273K = np.std(rmsd_273K)
std_compressibility = np.std(rmsd_compressibility)

# Add text box with statistics
stats_text = (
    f'TIP4P 298K: Mean = {mean_298K:.2f} Å, Std Dev = {std_298K:.2f} Å\n'
    f'TIP4P 273K: Mean = {mean_273K:.2f} Å, Std Dev = {std_273K:.2f} Å\n'
    f'TIP4P Compressibility: Mean = {mean_compressibility:.2f} Å, Std Dev = {std_compressibility:.2f} Å'
)
plt.figtext(0.5, 0.01, stats_text, ha='center', fontsize=12, 
            bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))

# Adjust layout and save
plt.tight_layout(rect=[0, 0.05, 1, 0.95])
plt.savefig('combined_rmsd_comparison.png', dpi=300)
plt.close()

print("Combined RMSD plot has been saved as 'combined_rmsd_comparison.png'") 