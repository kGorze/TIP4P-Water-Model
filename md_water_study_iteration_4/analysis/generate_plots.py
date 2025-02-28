#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# Directory setup
analysis_dir = sys.argv[1]
plots_dir = sys.argv[2]

# Function to read xvg files
def read_xvg(filename):
    x = []
    y = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("#") or line.startswith("@"):
                continue
            values = line.strip().split()
            if len(values) >= 2:
                x.append(float(values[0]))
                y.append(float(values[1]))
    return np.array(x), np.array(y)

# Create plots directory if it does not exist
os.makedirs(plots_dir, exist_ok=True)

# Plot RDF O-O
print("Creating O-O RDF plot...")
try:
    x, y = read_xvg(os.path.join(analysis_dir, "rdf_OO.xvg"))
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, "b-", linewidth=2)
    plt.xlabel("r (nm)", fontsize=14)
    plt.ylabel("g(r)", fontsize=14)
    plt.title("Oxygen-Oxygen Radial Distribution Function", fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(plots_dir, "rdf_oo_plot.png"), dpi=300, bbox_inches="tight")
    plt.close()
except Exception as e:
    print(f"Error plotting O-O RDF: {e}")

# Plot Temperature
print("Creating temperature plot...")
try:
    x, y = read_xvg(os.path.join(analysis_dir, "temperature.xvg"))
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, "r-", linewidth=2)
    plt.xlabel("Time (ps)", fontsize=14)
    plt.ylabel("Temperature (K)", fontsize=14)
    plt.title("Temperature vs Time", fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(plots_dir, "temperature_plot.png"), dpi=300, bbox_inches="tight")
    plt.close()
except Exception as e:
    print(f"Error plotting temperature: {e}")

# Plot RMSD
print("Creating RMSD plot...")
try:
    x, y = read_xvg(os.path.join(analysis_dir, "rmsd.xvg"))
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, "b-", linewidth=2)
    plt.xlabel("Time (ns)", fontsize=14)
    plt.ylabel("RMSD (nm)", fontsize=14)
    plt.title("Root Mean Square Deviation vs Time", fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(plots_dir, "rmsd_plot.png"), dpi=300, bbox_inches="tight")
    plt.close()
except Exception as e:
    print(f"Error plotting RMSD: {e}")

print("Basic plots created successfully!")
