#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.patches as patches

# Define the simulation directories and their labels
data_paths = [
    "tip4p_273K/analysis/rmsd_data.dat",
    "tip4p_298K/analysis/rmsd_data.dat",
    "iteration_5_TIP4P/analysis/rmsd_data.dat"
]

labels = [
    "TIP4P @ 273 K",
    "TIP4P @ 298 K",
    "TIP4P Compressibility test"
]

# Define colors for each simulation (matching the color scheme in rmsd_comparison.py)
colors = [
    "#1f77b4",  # Blue for tip4p_273K
    "#ff7f0e",  # Orange for tip4p_298K
    "#2ca02c"   # Green for iteration_5_TIP4P
]

# Create the comparison directory if it doesn't exist
comparison_dir = "comparison"
os.makedirs(comparison_dir, exist_ok=True)

# Define simulation phases (approximate times in ns)
# These values should be adjusted based on your actual simulation protocol
phases = {
    "EM": (0, 0.05),    # Energy minimization (first 50 ps)
    "NVT": (0.05, 3.0),  # NVT equilibration (50 ps to 3 ns)
    "MD": (3.0, 8.0)     # Production MD (3 ns to 8 ns)
}

# Create a figure with a larger size for better visibility
plt.figure(figsize=(14, 10))

# Create a light red background for the entire plot
ax = plt.gca()
ax.add_patch(patches.Rectangle((0, 0), 8, 80, color='red', alpha=0.1, zorder=-100))

# Add vertical lines for phase boundaries
for phase, (start, end) in phases.items():
    if start > 0:  # Don't draw a line at the beginning
        plt.axvline(x=start, color='gray', linestyle='-', linewidth=1.5, alpha=0.7)
    
    # Add phase label in a box
    if phase != "EM":  # EM phase is usually too small to label clearly
        plt.text((start + end) / 2, 45, phase, 
                 bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'),
                 ha='center', va='center', fontsize=12, fontweight='bold')

# Dictionary to store statistics for each simulation
stats = {}

# Loop through each simulation's data file
for path, label, color in zip(data_paths, labels, colors):
    try:
        # Load the RMSD data (skipping header if present)
        data = np.loadtxt(path, comments='#')
        
        # Extract time and RMSD columns
        time_ps = data[:, 0]
        rmsd_A = data[:, 1]
        
        # Convert time from ps to ns for better readability
        time_ns = time_ps / 1000.0
        
        # Plot on the same figure with consistent colors
        plt.plot(time_ns, rmsd_A, label=label, color=color, linewidth=2)
        
        # Calculate statistics
        mean_rmsd = np.mean(rmsd_A)
        std_rmsd = np.std(rmsd_A)
        min_rmsd = np.min(rmsd_A)
        max_rmsd = np.max(rmsd_A)
        
        # Store statistics
        stats[label] = {
            "mean": mean_rmsd,
            "std": std_rmsd,
            "min": min_rmsd,
            "max": max_rmsd
        }
        
        # Add a horizontal line for the mean RMSD value
        plt.axhline(y=mean_rmsd, color=color, linestyle='--', alpha=0.7)
        
        print(f"Loaded data from {path}")
        print(f"  Time range: {time_ns.min():.2f} - {time_ns.max():.2f} ns")
        print(f"  RMSD range: {rmsd_A.min():.2f} - {rmsd_A.max():.2f} Å")
        print(f"  Mean RMSD: {mean_rmsd:.2f} Å")
        print(f"  Std Dev: {std_rmsd:.2f} Å")
        print()
        
    except Exception as e:
        print(f"Error loading data from {path}: {e}")

# Add labels and title
plt.xlabel("Simulation Time (ns)", fontsize=14, fontweight='bold')
plt.ylabel("RMSD (Å)", fontsize=14, fontweight='bold')
plt.title("TIP4P Water RMSD vs. Time", fontsize=16, fontweight='bold')

# Set axis limits
plt.xlim(0, 8)  # Adjust based on your simulation length
plt.ylim(0, 50)  # Adjust based on your RMSD values

# Add a grid for better readability
plt.grid(True, alpha=0.3)

# Add a legend with a semi-transparent background
plt.legend(loc="upper right", framealpha=0.7, fontsize=12)

# Add statistics boxes for each simulation
for i, (label, stat) in enumerate(stats.items()):
    # Create statistics text
    stat_text = f"Mean RMSD: {stat['mean']:.2f} Å\n"
    stat_text += f"Std Dev: {stat['std']:.2f} Å\n"
    stat_text += f"Min RMSD: {stat['min']:.2f} Å\n"
    stat_text += f"Max RMSD: {stat['max']:.2f} Å"
    
    # Add statistics box
    plt.text(7.8, 45 - i*10, stat_text, 
             bbox=dict(facecolor='white', edgecolor=colors[i], boxstyle='round,pad=0.5'),
             ha='right', va='top', fontsize=10, color=colors[i])

# Add dotted horizontal line at 30 Å for reference
plt.axhline(y=30, color='brown', linestyle=':', linewidth=1.5, alpha=0.7)

# Adjust layout
plt.tight_layout()

# Save the figure with high resolution
output_file = os.path.join(comparison_dir, "combined_rmsd_phases.png")
plt.savefig(output_file, dpi=300)
print(f"Saved combined RMSD plot with phase labels to {output_file}")

# Create a second figure focusing on the equilibrated region
plt.figure(figsize=(14, 10))

# Create a light red background for the entire plot
ax = plt.gca()
ax.add_patch(patches.Rectangle((0, 30), 8, 15, color='red', alpha=0.1, zorder=-100))

# Loop through each simulation's data file again for the zoomed plot
for path, label, color in zip(data_paths, labels, colors):
    try:
        # Load the RMSD data
        data = np.loadtxt(path, comments='#')
        
        # Extract time and RMSD columns
        time_ps = data[:, 0]
        rmsd_A = data[:, 1]
        
        # Convert time from ps to ns
        time_ns = time_ps / 1000.0
        
        # Find the equilibration point (assuming it's around 250 ps based on previous analysis)
        equil_idx = np.where(time_ps >= 250)[0][0] if len(np.where(time_ps >= 250)[0]) > 0 else 0
        
        # Plot only the equilibrated portion
        plt.plot(time_ns[equil_idx:], rmsd_A[equil_idx:], label=label, color=color, linewidth=2)
        
        # Add a horizontal line for the mean of the equilibrated RMSD
        equil_mean = np.mean(rmsd_A[equil_idx:])
        plt.axhline(y=equil_mean, color=color, linestyle='--', alpha=0.7)
        
        # Add phase boundaries
        for phase, (start, end) in phases.items():
            if start > 0 and start >= time_ns[equil_idx]:
                plt.axvline(x=start, color='gray', linestyle='-', linewidth=1.5, alpha=0.7)
            
            # Add phase label in a box
            if phase != "EM" and end >= time_ns[equil_idx]:
                phase_start = max(start, time_ns[equil_idx])
                phase_mid = (phase_start + end) / 2
                plt.text(phase_mid, 45, phase, 
                         bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'),
                         ha='center', va='center', fontsize=12, fontweight='bold')
        
    except Exception as e:
        print(f"Error loading data for zoomed plot from {path}: {e}")

# Add labels and title for the zoomed plot
plt.xlabel("Simulation Time (ns)", fontsize=14, fontweight='bold')
plt.ylabel("RMSD (Å)", fontsize=14, fontweight='bold')
plt.title("Equilibrated TIP4P Water RMSD Comparison", fontsize=16, fontweight='bold')

# Set axis limits for better focus on the equilibrated region
plt.xlim(0, 8)
plt.ylim(30, 45)

# Add a grid for better readability
plt.grid(True, alpha=0.3)

# Add a legend with a semi-transparent background
plt.legend(loc="upper right", framealpha=0.7, fontsize=12)

# Add statistics boxes for equilibrated data
for i, (label, color) in enumerate(zip(labels, colors)):
    try:
        data = np.loadtxt(data_paths[i], comments='#')
        time_ps = data[:, 0]
        rmsd_A = data[:, 1]
        time_ns = time_ps / 1000.0
        
        equil_idx = np.where(time_ps >= 250)[0][0] if len(np.where(time_ps >= 250)[0]) > 0 else 0
        equil_data = rmsd_A[equil_idx:]
        
        equil_mean = np.mean(equil_data)
        equil_std = np.std(equil_data)
        equil_min = np.min(equil_data)
        equil_max = np.max(equil_data)
        
        # Create statistics text
        stat_text = f"Mean RMSD: {equil_mean:.2f} Å\n"
        stat_text += f"Std Dev: {equil_std:.2f} Å\n"
        stat_text += f"Min RMSD: {equil_min:.2f} Å\n"
        stat_text += f"Max RMSD: {equil_max:.2f} Å"
        
        # Add statistics box
        plt.text(7.8, 44 - i*3, stat_text, 
                 bbox=dict(facecolor='white', edgecolor=color, boxstyle='round,pad=0.5'),
                 ha='right', va='top', fontsize=9, color=color)
    except Exception as e:
        print(f"Error calculating equilibrated statistics for {label}: {e}")

# Adjust layout
plt.tight_layout()

# Save the zoomed figure
zoomed_output_file = os.path.join(comparison_dir, "equilibrated_rmsd_phases.png")
plt.savefig(zoomed_output_file, dpi=300)
print(f"Saved equilibrated RMSD plot with phase labels to {zoomed_output_file}")

# Show the plots if running interactively
# plt.show() 