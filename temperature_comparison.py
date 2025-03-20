#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from scipy import stats as scipy_stats
from matplotlib.ticker import MultipleLocator

# Define directories to compare
directories = [
    "tip4p_273K",
    "tip4p_298K",
    "iteration_5_TIP4P"
]

# Create comparison directory if it doesn't exist
comparison_dir = "comparison"
os.makedirs(comparison_dir, exist_ok=True)

# Define a consistent color scheme for each system
COLOR_MAP = {
    "tip4p_273K": "#1f77b4",      # Blue
    "tip4p_298K": "#ff7f0e",      # Orange
    "iteration_5_TIP4P": "#2ca02c", # Green - Keep directory name as-is in the color map
    "tip4p_273K_compressability": "#2ca02c" # Green - Also add display name to the color map
}

# Target temperature values (K)
TARGET_TEMPERATURES = {
    "tip4p_273K": 273.15,      # 0°C
    "tip4p_298K": 298.15,      # 25°C
    "iteration_5_TIP4P": 273.15,  # Original directory name
    "tip4p_273K_compressability": 273.15  # Display name
}

def load_temperature_data(directory):
    """Load temperature data from a directory's analysis folder."""
    data_path = os.path.join(directory, "analysis", "temperature_stats.dat")
    
    # First try to load the statistics file to get mean and std
    temp_stats = {}
    if os.path.exists(data_path):
        try:
            with open(data_path, 'r') as f:
                for line in f:
                    if "Mean temperature:" in line:
                        temp_stats['mean'] = float(line.split(':')[1].split()[0])
                    elif "Std deviation:" in line:
                        temp_stats['std'] = float(line.split(':')[1].split()[0])
        except Exception as e:
            print(f"Error reading temperature stats from {data_path}: {e}")
    
    # Now try to load the raw data if available
    # Check if there's a temperature data file from thermodynamics_analysis.py
    # This could be in different formats, so we'll try a few possibilities
    possible_data_files = [
        os.path.join(directory, "analysis", "temperature_data.dat"),
        os.path.join(directory, "analysis", "temperature_data.txt")
    ]
    
    data_file = None
    for file_path in possible_data_files:
        if os.path.exists(file_path):
            data_file = file_path
            break
    
    if data_file:
        try:
            # Try to load the data file
            data = pd.read_csv(data_file, sep='\t', comment='#')
            
            # Check if the data has the expected format
            if 'Time (ps)' in data.columns and 'Temperature' in data.columns:
                # Convert time from ps to ns
                data['Time (ns)'] = data['Time (ps)'] / 1000.0
                
                # Add directory name as a label, but rename iteration_5_TIP4P for display
                system_name = os.path.basename(directory)
                if system_name == "iteration_5_TIP4P":
                    data['System'] = "tip4p_273K_compressability"
                else:
                    data['System'] = system_name
                
                return data, temp_stats
        except Exception as e:
            print(f"Error loading temperature data from {data_file}: {e}")
    
    # If we couldn't load the data file, try to extract from the EDR file
    # This is a fallback and requires the thermodynamics_analysis.py script
    print(f"No temperature data file found for {directory}, trying to extract from EDR file...")
    
    # Try to import the parse_edr_file function from the thermodynamics_analysis.py script
    try:
        sys.path.append(os.path.join(directory, "scripts"))
        from utils import parse_edr_file
        
        # Parse temperature data from EDR file
        edr_file = os.path.join(directory, "md.edr")
        if os.path.exists(edr_file):
            temp_data = parse_edr_file(edr_file, properties=['Temperature'])
            if 'Temperature' in temp_data and temp_data['Temperature'] is not None:
                data = temp_data['Temperature']
                
                # Convert time from ps to ns
                data['Time (ns)'] = data['Time (ps)'] / 1000.0
                
                # Rename iteration_5_TIP4P for display
                system_name = os.path.basename(directory)
                if system_name == "iteration_5_TIP4P":
                    data['System'] = "tip4p_273K_compressability"
                else:
                    data['System'] = system_name
                
                return data, temp_stats
    except Exception as e:
        print(f"Error extracting temperature data from EDR file: {e}")
    
    # If we still don't have data, return None
    print(f"Warning: Could not load temperature data for {directory}")
    return None, temp_stats

def add_significance_bars(ax, data, systems, y_position):
    """Add significance bars between systems based on t-test."""
    pairs = []
    for i in range(len(systems)):
        for j in range(i+1, len(systems)):
            pairs.append((systems[i], systems[j]))
    
    for pair in pairs:
        sys1_data = data[data['System'] == pair[0]]['Temperature'].values
        sys2_data = data[data['System'] == pair[1]]['Temperature'].values
        
        # Perform t-test
        t_stat, p_value = scipy_stats.ttest_ind(sys1_data, sys2_data, equal_var=False)
        
        # Get x positions
        x1 = systems.index(pair[0])
        x2 = systems.index(pair[1])
        
        # Add significance bar
        if p_value < 0.001:
            sig_symbol = '***'
        elif p_value < 0.01:
            sig_symbol = '**'
        elif p_value < 0.05:
            sig_symbol = '*'
        else:
            sig_symbol = 'ns'
            
        ax.plot([x1, x2], [y_position, y_position], 'k-', linewidth=1)
        ax.text((x1 + x2) / 2, y_position + 0.05, sig_symbol, ha='center')

def calculate_running_average(x, y, window=50):
    """Calculate running average with specified window size."""
    avg_y = np.zeros_like(y)
    n = len(y)
    
    for i in range(n):
        start = max(0, i - window // 2)
        end = min(n, i + window // 2 + 1)
        avg_y[i] = np.mean(y[start:end])
    
    return x, avg_y

def main():
    print("Loading temperature data from all directories...")
    
    # Load data from each directory
    all_data = []
    all_stats = {}
    
    for directory in directories:
        data, stats = load_temperature_data(directory)
        if data is not None:
            all_data.append(data)
            all_stats[directory] = stats
    
    if not all_data:
        print("No valid temperature data found. Please run thermodynamics_analysis.py in each directory first.")
        return
    
    # Combine all data into a single DataFrame
    combined_data = pd.concat(all_data, ignore_index=True)
    
    print("Generating temperature comparison plots...")
    
    # Set a nice style for plots
    sns.set_style("whitegrid")
    plt.rcParams.update({
        'font.size': 12,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'figure.figsize': (12, 8),
        'figure.dpi': 100,
    })
    
    # Get unique systems from the combined data
    systems = combined_data['System'].unique().tolist()
    
    # Ensure the systems are in the same order as the directories
    ordered_systems = []
    for dir_name in directories:
        if dir_name == "iteration_5_TIP4P":
            ordered_systems.append("tip4p_273K_compressability")
        else:
            ordered_systems.append(dir_name)
    
    # Make sure all systems are in combined_data
    systems = [sys for sys in ordered_systems if sys in combined_data['System'].unique()]
    
    # Plot 1: Temperature vs Time for all systems
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Add vertical grid lines at major time intervals (converted from 1000 ps to 1 ns)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.grid(True, which='major', axis='both', linestyle='--', alpha=0.7)
    
    # Plot raw data
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        ax.plot(system_data['Time (ns)'], system_data['Temperature'], 
                label=system, color=COLOR_MAP[system], alpha=0.5, linewidth=1)
    
    # Add smoothed trend lines
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        # Sort by time to ensure proper running average
        system_data = system_data.sort_values('Time (ns)')
        # Calculate running average
        smooth_time, smooth_temp = calculate_running_average(
            system_data['Time (ns)'].values, system_data['Temperature'].values, window=100)
        ax.plot(smooth_time, smooth_temp, color=COLOR_MAP[system], 
                linewidth=2.5, alpha=0.8, linestyle='-')
    
    # Add horizontal lines for target temperatures
    for system in systems:
        if system in TARGET_TEMPERATURES:
            ax.axhline(y=TARGET_TEMPERATURES[system], color=COLOR_MAP[system], 
                      linestyle='--', alpha=0.5, linewidth=1)
    
    # Create a clean, bold legend
    legend = ax.legend(title="System", loc='upper right', framealpha=0.9, fontsize=12)
    plt.setp(legend.get_title(), fontweight='bold')
    plt.setp(legend.get_texts(), fontweight='bold')
    
    ax.set_xlabel('Time (ns)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Temperature (K)', fontweight='bold', fontsize=12)
    ax.set_title('Temperature Comparison', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, 'temperature_comparison.png'), dpi=300)
    plt.close()
    
    # Plot 2: Temperature distribution (KDE plot) for all systems
    fig, ax = plt.subplots(figsize=(12, 8))
    
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        sns.kdeplot(system_data['Temperature'], label=system, fill=True, 
                   alpha=0.3, color=COLOR_MAP[system], ax=ax)
        
        # Add vertical line for target temperature (more subtle)
        if system in TARGET_TEMPERATURES:
            ax.axvline(x=TARGET_TEMPERATURES[system], color=COLOR_MAP[system], 
                      linestyle='--', alpha=0.5, linewidth=1)
    
    # Create a clean, bold legend
    legend = ax.legend(title="System", loc='upper right', framealpha=0.9, fontsize=12)
    plt.setp(legend.get_title(), fontweight='bold')
    plt.setp(legend.get_texts(), fontweight='bold')
    
    ax.set_xlabel('Temperature (K)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Probability Density', fontweight='bold', fontsize=12)
    ax.set_title('Temperature Distribution', fontsize=14, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, 'temperature_distribution_comparison.png'), dpi=300)
    plt.close()
    
    # Plot 3: Boxplot of temperature values with statistical significance
    fig, ax = plt.subplots(figsize=(10, 8))
    
    try:
        # Try with hue parameter (newer seaborn versions)
        sns.boxplot(x='System', y='Temperature', data=combined_data, ax=ax,
                   palette=COLOR_MAP, hue='System', legend=False, notch=True, showfliers=False)
    except:
        # Fallback for older seaborn versions
        sns.boxplot(x='System', y='Temperature', data=combined_data, ax=ax,
                   palette=COLOR_MAP, notch=True, showfliers=False)
    
    # Add individual data points with jitter and transparency
    sns.stripplot(x='System', y='Temperature', data=combined_data, ax=ax,
                 size=4, color='black', alpha=0.3, jitter=True)
    
    # Add target temperature lines (more subtle)
    for i, system in enumerate(systems):
        if system in TARGET_TEMPERATURES:
            ax.hlines(y=TARGET_TEMPERATURES[system], xmin=i-0.4, xmax=i+0.4, 
                     colors=COLOR_MAP[system], linestyles='dashed', linewidth=1, alpha=0.5)
    
    ax.set_xlabel('System', fontweight='bold', fontsize=12)
    ax.set_ylabel('Temperature (K)', fontweight='bold', fontsize=12)
    ax.set_title('Temperature Distribution by System', fontsize=14, fontweight='bold')
    ax.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, 'temperature_boxplot.png'), dpi=300)
    plt.close()
    
    # Plot 4: Temperature fluctuations (standard deviation comparison)
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Calculate standard deviations for each system
    std_data = []
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        std = system_data['Temperature'].std()
        mean = system_data['Temperature'].mean()
        target = TARGET_TEMPERATURES.get(system, mean)
        rel_std = (std / target) * 100  # Relative std as percentage
        
        std_data.append({
            'System': system,
            'Standard Deviation (K)': std,
            'Relative Std (%)': rel_std
        })
    
    std_df = pd.DataFrame(std_data)
    
    # Create a bar plot for standard deviations
    bar_positions = np.arange(len(systems))
    bar_width = 0.35
    
    # Plot absolute standard deviations
    bars1 = ax.bar(bar_positions - bar_width/2, std_df['Standard Deviation (K)'], 
                  bar_width, label='Absolute Std (K)', alpha=0.7)
    
    # Plot relative standard deviations
    ax2 = ax.twinx()
    bars2 = ax2.bar(bar_positions + bar_width/2, std_df['Relative Std (%)'], 
                   bar_width, label='Relative Std (%)', alpha=0.7, color='orange')
    
    # Add labels and title
    ax.set_xlabel('System', fontweight='bold', fontsize=12)
    ax.set_ylabel('Standard Deviation (K)', fontweight='bold', fontsize=12)
    ax2.set_ylabel('Relative Standard Deviation (%)', fontweight='bold', fontsize=12)
    ax.set_title('Temperature Fluctuations', fontsize=14, fontweight='bold')
    
    # Set x-tick labels
    ax.set_xticks(bar_positions)
    ax.set_xticklabels(systems, fontweight='bold')
    
    # Add a legend with bold text
    legend = ax.legend(handles=[bars1, bars2], loc='upper left')
    plt.setp(legend.get_texts(), fontweight='bold')
    
    # Add value labels on bars
    for i, bar in enumerate(bars1):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.05,
               f'{height:.2f}', ha='center', va='bottom', fontweight='bold', fontsize=9)
    
    for i, bar in enumerate(bars2):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                f'{height:.2f}%', ha='center', va='bottom', fontweight='bold', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, 'temperature_fluctuations.png'), dpi=300)
    plt.close()
    
    # Calculate and save statistics
    temp_stats = combined_data.groupby('System')['Temperature'].agg(['mean', 'std', 'min', 'max', 'count']).reset_index()
    temp_stats.columns = ['System', 'Mean Temperature (K)', 'Std Dev (K)', 
                         'Min Temperature (K)', 'Max Temperature (K)', 'Sample Count']
    
    # Add target comparison
    temp_stats['Target Temperature (K)'] = temp_stats['System'].map(
        lambda x: TARGET_TEMPERATURES.get(x, 'N/A'))
    temp_stats['Deviation from Target (K)'] = temp_stats.apply(
        lambda row: (row['Mean Temperature (K)'] - float(row['Target Temperature (K)'])) 
        if row['Target Temperature (K)'] != 'N/A' else 'N/A', axis=1)
    temp_stats['Relative Deviation (%)'] = temp_stats.apply(
        lambda row: ((row['Mean Temperature (K)'] - float(row['Target Temperature (K)'])) / 
                    float(row['Target Temperature (K)']) * 100) 
        if row['Target Temperature (K)'] != 'N/A' else 'N/A', axis=1)
    
    stats_file = os.path.join(comparison_dir, 'temperature_statistics.txt')
    
    with open(stats_file, 'w') as f:
        f.write("Temperature Statistics Comparison\n")
        f.write("================================\n\n")
        
        # Write statistical significance test results
        f.write("Statistical Significance (p-values from t-test):\n")
        f.write("-----------------------------------------------\n")
        for i, sys1 in enumerate(systems):
            for j, sys2 in enumerate(systems):
                if i < j:
                    sys1_data = combined_data[combined_data['System'] == sys1]['Temperature'].values
                    sys2_data = combined_data[combined_data['System'] == sys2]['Temperature'].values
                    t_stat, p_value = scipy_stats.ttest_ind(sys1_data, sys2_data, equal_var=False)
                    sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
                    f.write(f"{sys1} vs {sys2}: p={p_value:.6f} {sig}\n")
        f.write("\n")
        
        # Write detailed statistics for each system
        for _, row in temp_stats.iterrows():
            f.write(f"System: {row['System']}\n")
            f.write(f"  Mean Temperature: {row['Mean Temperature (K)']:.4f} K\n")
            f.write(f"  Std Deviation: {row['Std Dev (K)']:.4f} K\n")
            f.write(f"  Min Temperature: {row['Min Temperature (K)']:.4f} K\n")
            f.write(f"  Max Temperature: {row['Max Temperature (K)']:.4f} K\n")
            f.write(f"  Sample Count: {row['Sample Count']}\n")
            f.write(f"  Target Temperature: {row['Target Temperature (K)']} K\n")
            if row['Deviation from Target (K)'] != 'N/A':
                f.write(f"  Deviation from Target: {row['Deviation from Target (K)']:.4f} K\n")
                f.write(f"  Relative Deviation: {row['Relative Deviation (%)']:.4f}%\n")
            else:
                f.write(f"  Deviation from Target: N/A\n")
                f.write(f"  Relative Deviation: N/A\n")
            f.write("\n")
    
    print(f"Saved temperature comparison plots and statistics to {comparison_dir}/")
    print(f"Statistics saved to {stats_file}")

if __name__ == "__main__":
    import sys
    main() 