#!/usr/bin/env python3
"""
Hydrogen Bond Lifetime Comparison

This script compares hydrogen bond lifetimes between different water models:
- TIP4P at 273K
- TIP4P at 298K
- Iteration 5 TIP4P model

The script generates comparative plots of the hydrogen bond lifetime distributions
and cumulative distributions, similar to the individual hblife_detailed_plot.png files.

Output is saved to the comparison directory.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import expon
import pandas as pd

# Create the hydrogen bond comparison directory if it doesn't exist
comparison_dir = "comparison/hydrogen_bonds"
os.makedirs(comparison_dir, exist_ok=True)

# Define the directories containing the different models
directories = [
    "tip4p_273K",
    "tip4p_298K",
    "iteration_5_TIP4P"
]

# Define a consistent color scheme and labels for each system
COLOR_MAP = {
    "tip4p_273K": "#1f77b4",      # Blue
    "tip4p_298K": "#ff7f0e",      # Orange
    "iteration_5_TIP4P": "#2ca02c" # Green
}

LABELS = {
    "tip4p_273K": "TIP4P (273K)",
    "tip4p_298K": "TIP4P (298K)",
    "iteration_5_TIP4P": "Iteration 5 TIP4P"
}

def load_hbond_lifetime_data(directory):
    """
    Load hydrogen bond lifetime data and statistics from a directory.
    
    Parameters:
    -----------
    directory : str
        Path to the directory containing the analysis data
        
    Returns:
    --------
    tuple
        (lifetimes_array, lifetime_stats)
    """
    # Initialize variables
    lifetimes = None
    stats = {}
    
    # Path to the hydrogen bond lifetime data file
    lifetimes_file = os.path.join(directory, "analysis", "hbond_lifetimes.dat")
    
    # Path to the hydrogen bond lifetime statistics file
    stats_file = os.path.join(directory, "analysis", "hbond_lifetime_stats.txt")
    
    # Load hydrogen bond lifetimes
    if os.path.exists(lifetimes_file):
        try:
            lifetimes = np.loadtxt(lifetimes_file, comments='#')
            print(f"Loaded {len(lifetimes)} hydrogen bond lifetimes from {lifetimes_file}")
        except Exception as e:
            print(f"Error loading hydrogen bond lifetimes from {lifetimes_file}: {e}")
    else:
        print(f"Warning: Hydrogen bond lifetime data file not found at {lifetimes_file}")
    
    # Load hydrogen bond lifetime statistics
    if os.path.exists(stats_file):
        try:
            with open(stats_file, 'r') as f:
                for line in f:
                    if "Mean lifetime:" in line:
                        stats['mean'] = float(line.split(':')[1].split()[0])
                    elif "Median lifetime:" in line:
                        stats['median'] = float(line.split(':')[1].split()[0])
                    elif "Standard deviation:" in line:
                        stats['std'] = float(line.split(':')[1].split()[0])
                    elif "Minimum lifetime:" in line:
                        stats['min'] = float(line.split(':')[1].split()[0])
                    elif "Maximum lifetime:" in line:
                        stats['max'] = float(line.split(':')[1].split()[0])
                    elif "25th percentile:" in line:
                        stats['p25'] = float(line.split(':')[1].split()[0])
                    elif "75th percentile:" in line:
                        stats['p75'] = float(line.split(':')[1].split()[0])
            print(f"Loaded hydrogen bond lifetime statistics from {stats_file}")
        except Exception as e:
            print(f"Error loading hydrogen bond lifetime statistics from {stats_file}: {e}")
    else:
        print(f"Warning: Hydrogen bond lifetime statistics file not found at {stats_file}")
    
    return lifetimes, stats

def generate_comparison_plots(data_dict, stats_dict):
    """
    Generate comparison plots for hydrogen bond lifetimes.
    
    Parameters:
    -----------
    data_dict : dict
        Dictionary containing hydrogen bond lifetime arrays for each system
    stats_dict : dict
        Dictionary containing hydrogen bond lifetime statistics for each system
    """
    # Create a figure for the comparison plots
    fig = plt.figure(figsize=(14, 12))
    
    # Create a 2x1 grid for the two plots
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.3)
    
    # Find the global maximum lifetime across all systems for consistent x-axis
    global_max_lifetime = 0
    for system in directories:
        if system in stats_dict:
            global_max_lifetime = max(global_max_lifetime, stats_dict[system]['max'])
    
    # Round up to the nearest whole number for a cleaner plot
    global_max_lifetime = np.ceil(global_max_lifetime)
    print(f"Using global maximum lifetime of {global_max_lifetime} ps for x-axis scaling")
    
    # Plot 1: Histograms of hydrogen bond lifetimes
    ax1 = fig.add_subplot(gs[0])
    
    # Plot histograms for each system
    for system in directories:
        if system in data_dict and data_dict[system] is not None:
            # Get the lifetime data
            lifetimes = data_dict[system]
            
            # Calculate histogram bins based on global maximum
            bins = np.linspace(0, global_max_lifetime, 20)
            
            # Plot histogram
            ax1.hist(lifetimes, bins=bins, alpha=0.6, 
                    color=COLOR_MAP[system], 
                    label=f"{LABELS[system]}", 
                    density=True, histtype='bar', edgecolor='black', linewidth=1)
            
            # Plot exponential fit
            x = np.linspace(0, global_max_lifetime, 1000)
            loc, scale = expon.fit(lifetimes)
            ax1.plot(x, expon.pdf(x, loc=loc, scale=scale), 
                    color=COLOR_MAP[system], linestyle='-', linewidth=2,
                    label=f"{LABELS[system]} Exponential fit (τ = {scale:.2f} ps)")
            
            # Add vertical lines for mean and median
            mean = stats_dict[system]['mean']
            median = stats_dict[system]['median']
            
            ax1.axvline(mean, color=COLOR_MAP[system], linestyle='--', alpha=0.8, 
                        label=f"{LABELS[system]} Mean: {mean:.2f} ps")
    
    # Customize the histogram plot
    ax1.set_xlabel('Hydrogen bond lifetime (ps)', fontweight='bold')
    ax1.set_ylabel('Probability density', fontweight='bold')
    ax1.set_title('Distribution of Hydrogen Bond Lifetimes', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, global_max_lifetime)
    ax1.legend(loc='upper right', fontsize=10)
    
    # Plot 2: Cumulative distributions
    ax2 = fig.add_subplot(gs[1])
    
    # Plot cumulative distributions for each system
    for system in directories:
        if system in data_dict and data_dict[system] is not None:
            # Get the lifetime data
            lifetimes = data_dict[system]
            
            # Calculate bins for the cumulative histogram based on global maximum
            bins = np.linspace(0, global_max_lifetime, 30)
            
            # Plot cumulative histogram
            ax2.hist(lifetimes, bins=bins, alpha=0.7, color=COLOR_MAP[system], 
                    density=True, cumulative=True, histtype='step', linewidth=2,
                    label=f"{LABELS[system]} Cumulative distribution")
            
            # Plot exponential CDF
            x = np.linspace(0, global_max_lifetime, 1000)
            loc, scale = expon.fit(lifetimes)
            ax2.plot(x, expon.cdf(x, loc=loc, scale=scale), 
                     color=COLOR_MAP[system], linestyle='--', linewidth=1.5,
                     label=f"{LABELS[system]} Exponential CDF")
            
            # Add markers for percentiles
            percentiles = [25, 50, 75, 90]
            for p in percentiles:
                if p == 25:
                    p_value = stats_dict[system]['p25'] if 'p25' in stats_dict[system] else np.percentile(lifetimes, p)
                elif p == 50:
                    p_value = stats_dict[system]['median']
                elif p == 75:
                    p_value = stats_dict[system]['p75'] if 'p75' in stats_dict[system] else np.percentile(lifetimes, p)
                else:
                    p_value = np.percentile(lifetimes, p)
                
                ax2.plot(p_value, p/100, 'o', color=COLOR_MAP[system], markersize=8)
                
                # Add percentile annotations only for the 50% point to avoid cluttering
                if p == 50:
                    ax2.annotate(f"{p}%: {p_value:.2f} ps", 
                              xy=(p_value, p/100),
                              xytext=(p_value+0.4, p/100+0.05),
                              color=COLOR_MAP[system],
                              arrowprops=dict(arrowstyle="->", color=COLOR_MAP[system]))
    
    # Customize the cumulative distribution plot
    ax2.set_xlabel('Hydrogen bond lifetime (ps)', fontweight='bold')
    ax2.set_ylabel('Cumulative probability', fontweight='bold')
    ax2.set_title('Cumulative Distribution of Hydrogen Bond Lifetimes', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, global_max_lifetime)
    ax2.legend(loc='lower right', fontsize=10)
    
    # Create a summary table of statistics
    stats_df = pd.DataFrame()
    for system in directories:
        if system in stats_dict:
            stats_df[LABELS[system]] = pd.Series({
                'Mean (ps)': stats_dict[system]['mean'],
                'Median (ps)': stats_dict[system]['median'],
                'Std Dev (ps)': stats_dict[system]['std'],
                '25% (ps)': stats_dict[system]['p25'] if 'p25' in stats_dict[system] else np.nan,
                '75% (ps)': stats_dict[system]['p75'] if 'p75' in stats_dict[system] else np.nan,
                'Min (ps)': stats_dict[system]['min'],
                'Max (ps)': stats_dict[system]['max']
            })
    
    # Add a row for percent difference from TIP4P 298K (room temperature) as reference
    if 'TIP4P (298K)' in stats_df.columns:
        reference_mean = stats_df.loc['Mean (ps)', 'TIP4P (298K)']
        reference_median = stats_df.loc['Median (ps)', 'TIP4P (298K)']
        
        percent_diff_means = {}
        percent_diff_medians = {}
        
        for column in stats_df.columns:
            if column != 'TIP4P (298K)':
                mean_val = stats_df.loc['Mean (ps)', column]
                median_val = stats_df.loc['Median (ps)', column]
                
                percent_diff_mean = ((mean_val - reference_mean) / reference_mean) * 100
                percent_diff_median = ((median_val - reference_median) / reference_median) * 100
                
                percent_diff_means[column] = percent_diff_mean
                percent_diff_medians[column] = percent_diff_median
        
        # Add reference value for TIP4P 298K
        percent_diff_means['TIP4P (298K)'] = 0.0
        percent_diff_medians['TIP4P (298K)'] = 0.0
        
        # Add to DataFrame
        stats_df.loc['% Diff from 298K (Mean)', :] = pd.Series(percent_diff_means)
        stats_df.loc['% Diff from 298K (Median)', :] = pd.Series(percent_diff_medians)
    
    # Add the table to the figure
    ax_table = fig.add_axes([0.15, 0.01, 0.7, 0.1], frame_on=False)
    ax_table.axis('off')
    
    # Format the table data for display
    table_data = np.round(stats_df.values, 2)
    
    # Create the table
    table = ax_table.table(
        cellText=table_data,
        rowLabels=stats_df.index,
        colLabels=stats_df.columns,
        cellLoc='center',
        loc='center'
    )
    
    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)
    
    # Color code the cells in the percent difference rows
    if 'TIP4P (298K)' in stats_df.columns:
        for (row, col), cell in table.get_celld().items():
            if row > 0 and stats_df.index[row-1] in ['% Diff from 298K (Mean)', '% Diff from 298K (Median)']:
                value = table_data[row-1, col-1]
                if abs(value) < 5:
                    cell.set_facecolor('#e6ffe6')  # Light green for small differences
                elif abs(value) < 15:
                    cell.set_facecolor('#ffffcc')  # Light yellow for moderate differences
                else:
                    cell.set_facecolor('#ffe6e6')  # Light red for large differences
    
    # Add title to the table
    ax_table.set_title('Hydrogen Bond Lifetime Statistics', fontsize=12, fontweight='bold', pad=20)
    
    # Add overall figure title
    fig.suptitle('Comparison of Hydrogen Bond Lifetimes between Water Models', 
                fontsize=16, fontweight='bold', y=0.98)
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.12, 1, 0.96])
    
    # Save the figure
    plt.savefig(os.path.join(comparison_dir, 'hbond_lifetime_comparison.png'), dpi=300, bbox_inches='tight')
    print(f"Saved hydrogen bond lifetime comparison plot to {os.path.join(comparison_dir, 'hbond_lifetime_comparison.png')}")
    
    # Also save the statistics table as a CSV file
    stats_df.to_csv(os.path.join(comparison_dir, 'hbond_lifetime_stats_comparison.csv'))
    print(f"Saved hydrogen bond lifetime statistics to {os.path.join(comparison_dir, 'hbond_lifetime_stats_comparison.csv')}")

def main():
    print("Loading hydrogen bond lifetime data from all directories...")
    
    # Load data from each directory
    data_dict = {}
    stats_dict = {}
    
    for directory in directories:
        lifetimes, stats = load_hbond_lifetime_data(directory)
        if lifetimes is not None:
            data_dict[directory] = lifetimes
            stats_dict[directory] = stats
    
    if not data_dict:
        print("No valid hydrogen bond lifetime data found. Please run hbond_analysis.py in each directory first.")
        return
    
    print("Generating hydrogen bond lifetime comparison plots...")
    generate_comparison_plots(data_dict, stats_dict)
    
    print("Hydrogen bond lifetime comparison complete!")

if __name__ == "__main__":
    main() 