#!/usr/bin/env python3
"""
This script reads the histogram data from 'distance_histogram.dat',
where each line contains:
    lower_bound   upper_bound   count
It computes the midpoint for each bin and plots the data as a bar chart.
"""

import matplotlib
matplotlib.use('TkAgg')   # Switch backend to TkAgg to avoid GTK4 issues

import numpy as np
import matplotlib.pyplot as plt

# Load the data, skipping any comment lines (starting with '#')
data = np.loadtxt("distance_histogram.dat", comments='#')

if data.size == 0:
    raise ValueError("No data found in 'distance_histogram.dat'.")

# Extract columns: lower bound, upper bound, and count
lower_bounds = data[:, 0]
upper_bounds = data[:, 1]
counts = data[:, 2]

# Compute midpoints for each bin
midpoints = (lower_bounds + upper_bounds) / 2.0

# Determine the bin width (assuming uniform bin width)
bin_width = upper_bounds[0] - lower_bounds[0]

# Create a bar chart
plt.figure(figsize=(8, 6))
plt.bar(midpoints, counts, width=bin_width, align='center', edgecolor='black')
plt.xlabel("Distance (Å)")
plt.ylabel("Count")
plt.title("Distance Histogram")
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()

