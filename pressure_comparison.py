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

# Target pressure values (now in atm)
TARGET_PRESSURE = {
    "tip4p_273K": 0.986923,      # 1 bar in atm
    "tip4p_298K": 0.986923,      # 1 bar in atm
    "iteration_5_TIP4P": 0.986923,  # 1 bar in atm (original name)
    "tip4p_273K_compressability": 0.986923  # 1 bar in atm (display name)
}

# Conversion factor from bar to atm
BAR_TO_ATM = 0.986923

def load_pressure_data(directory):
    """Load pressure data from a directory's analysis folder."""
    data_path = os.path.join(directory, "analysis", "pressure_stats.dat")
    
    # First try to load the statistics file to get mean and std
    pressure_stats = {}
    if os.path.exists(data_path):
        try:
            with open(data_path, 'r') as f:
                for line in f:
                    if "Mean pressure:" in line:
                        pressure_stats['mean'] = float(line.split(':')[1].split()[0])
                    elif "Standard deviation:" in line:
                        pressure_stats['std'] = float(line.split(':')[1].split()[0])
        except Exception as e:
            print(f"Error reading pressure stats from {data_path}: {e}")
    
    # Now try to load the raw data if available
    # Check if there's a pressure data file from thermodynamics_analysis.py
    data_file = os.path.join(directory, "analysis", "pressure_data.dat")
    
    if os.path.exists(data_file):
        try:
            # Try to load the data file
            data = pd.read_csv(data_file, sep='\t', comment='#')
            
            # Check if the data has the expected format
            if 'Time (ps)' in data.columns and 'Pressure' in data.columns:
                # Convert time from ps to ns
                data['Time (ns)'] = data['Time (ps)'] / 1000.0
                
                # Convert pressure from bar to atm
                data['Pressure (atm)'] = data['Pressure'] * BAR_TO_ATM
                
                # Add directory name as a label, but rename iteration_5_TIP4P for display
                system_name = os.path.basename(directory)
                if system_name == "iteration_5_TIP4P":
                    data['System'] = "tip4p_273K_compressability"
                else:
                    data['System'] = system_name
                
                return data, pressure_stats
        except Exception as e:
            print(f"Error loading pressure data from {data_file}: {e}")
    
    # If we couldn't load the data file, try to extract from the EDR file
    # This is a fallback and requires the thermodynamics_analysis.py script
    print(f"No pressure data file found for {directory}, trying to extract from EDR file...")
    
    # Try to import the parse_edr_file function from the thermodynamics_analysis.py script
    try:
        sys.path.append(os.path.join(directory, "scripts"))
        from utils import parse_edr_file
        
        # Parse pressure data from EDR file
        edr_file = os.path.join(directory, "md.edr")
        if os.path.exists(edr_file):
            pressure_data = parse_edr_file(edr_file, properties=['Pressure'])
            if 'Pressure' in pressure_data and pressure_data['Pressure'] is not None:
                data = pressure_data['Pressure']
                
                # Convert time from ps to ns
                data['Time (ns)'] = data['Time (ps)'] / 1000.0
                
                # Convert pressure from bar to atm
                data['Pressure (atm)'] = data['Pressure'] * BAR_TO_ATM
                
                # Rename iteration_5_TIP4P for display
                system_name = os.path.basename(directory)
                if system_name == "iteration_5_TIP4P":
                    data['System'] = "tip4p_273K_compressability"
                else:
                    data['System'] = system_name
                
                return data, pressure_stats
    except Exception as e:
        print(f"Error extracting pressure data from EDR file: {e}")
    
    # If we still don't have data, return None
    print(f"Warning: Could not load pressure data for {directory}")
    return None, pressure_stats

def add_significance_bars(ax, data, systems, y_position):
    """Add significance bars between systems based on t-test."""
    pairs = []
    for i in range(len(systems)):
        for j in range(i+1, len(systems)):
            pairs.append((systems[i], systems[j]))
    
    for pair in pairs:
        sys1_data = data[data['System'] == pair[0]]['Pressure'].values
        sys2_data = data[data['System'] == pair[1]]['Pressure'].values
        
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
        ax.text((x1 + x2) / 2, y_position + 5, sig_symbol, ha='center')

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
    print("Loading pressure data from all directories...")
    
    # Load data from each directory
    all_data = []
    all_stats = {}
    
    for directory in directories:
        data, stats = load_pressure_data(directory)
        if data is not None:
            all_data.append(data)
            all_stats[directory] = stats
    
    if not all_data:
        print("No valid pressure data found. Please run thermodynamics_analysis.py in each directory first.")
        return
    
    # Combine all data into a single DataFrame
    combined_data = pd.concat(all_data, ignore_index=True)
    
    print("Generating pressure comparison plots...")
    
    # Set a nice style for plots
    sns.set_style("whitegrid")
    plt.rcParams.update({
        'font.size': 12,
        'font.family': 'Times New Roman',
        'font.weight': 'bold',
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
    
    # Plot 1: Pressure vs Time for all systems
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Add vertical grid lines at major time intervals (converted from 1000 ps to 1 ns)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.grid(True, which='major', axis='both', linestyle='--', alpha=0.7)
    
    # Plot raw data with low alpha for clarity (pressure data can be noisy)
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        ax.plot(system_data['Time (ns)'], system_data['Pressure (atm)'], 
                label=system, color=COLOR_MAP[system], alpha=0.3, linewidth=0.5)
    
    # Add smoothed trend lines with larger window due to pressure fluctuations
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        # Sort by time to ensure proper running average
        system_data = system_data.sort_values('Time (ns)')
        # Calculate running average with larger window for pressure
        smooth_time, smooth_pressure = calculate_running_average(
            system_data['Time (ns)'].values, system_data['Pressure (atm)'].values, window=200)
        ax.plot(smooth_time, smooth_pressure, color=COLOR_MAP[system], 
                linewidth=2.5, alpha=0.8, linestyle='-')
    
    # Add horizontal line for target pressure (converted to atm)
    ax.axhline(y=TARGET_PRESSURE["tip4p_273K"], color='black', linestyle='--', alpha=0.5, linewidth=1)
    
    # Create a clean, bold legend
    legend = ax.legend(title="System", loc='upper right', framealpha=0.9, fontsize=12)
    plt.setp(legend.get_title(), fontweight='bold')
    plt.setp(legend.get_texts(), fontweight='bold')
    
    ax.set_xlabel('Time (ns)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Pressure (atm)', fontweight='bold', fontsize=12)
    ax.set_title('Pressure Comparison', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, 'pressure_comparison.png'), dpi=300)
    plt.close()
    
    # Plot 2: Pressure distribution (KDE plot) for all systems
    fig, ax = plt.subplots(figsize=(12, 8))
    
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        sns.kdeplot(system_data['Pressure (atm)'], label=system, fill=True, 
                   alpha=0.3, color=COLOR_MAP[system], ax=ax)
        
        # Add vertical line for mean pressure
        mean_pressure = system_data['Pressure (atm)'].mean()
        ax.axvline(x=mean_pressure, color=COLOR_MAP[system], 
                  linestyle='-', alpha=0.7, linewidth=1.5)
    
    # Add vertical line for target pressure
    ax.axvline(x=TARGET_PRESSURE["tip4p_273K"], color='black', linestyle='--', alpha=0.5, linewidth=1)
    
    ax.set_xlabel('Pressure (atm)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Probability Density', fontweight='bold', fontsize=12)
    ax.set_title('Pressure Distribution', fontsize=14, fontweight='bold')
    
    # Create a clean, bold legend
    legend = ax.legend(title="System", loc='upper right', framealpha=0.9, fontsize=12)
    plt.setp(legend.get_title(), fontweight='bold')
    plt.setp(legend.get_texts(), fontweight='bold')
    
    ax.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, 'pressure_distribution_comparison.png'), dpi=300)
    plt.close()
    
    # Plot 3: Boxplot of pressure values with statistical significance
    fig, ax = plt.subplots(figsize=(10, 8))
    
    try:
        # Try with hue parameter (newer seaborn versions)
        sns.boxplot(x='System', y='Pressure (atm)', data=combined_data, ax=ax,
                   palette=COLOR_MAP, hue='System', legend=False, notch=True, showfliers=False)
    except:
        # Fallback for older seaborn versions
        sns.boxplot(x='System', y='Pressure (atm)', data=combined_data, ax=ax,
                   palette=COLOR_MAP, notch=True, showfliers=False)
    
    # Add individual data points with jitter and transparency
    # For pressure, we'll sample points to avoid overcrowding
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        # Sample every 10th point to avoid overcrowding
        sampled_data = system_data.iloc[::10]
        sns.stripplot(x='System', y='Pressure (atm)', data=sampled_data, ax=ax,
                     size=3, color='black', alpha=0.2, jitter=True)
    
    # Add target pressure line
    ax.axhline(y=TARGET_PRESSURE["tip4p_273K"], color='black', linestyle='--', alpha=0.5, linewidth=1)
    
    ax.set_xlabel('System', fontweight='bold', fontsize=12)
    ax.set_ylabel('Pressure (atm)', fontweight='bold', fontsize=12)
    ax.set_title('Pressure Distribution by System', fontsize=14, fontweight='bold')
    ax.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, 'pressure_boxplot.png'), dpi=300)
    plt.close()
    
    # Plot 4: Pressure convergence (running average)
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Calculate running average for each system
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        # Sort by time to ensure proper running average
        system_data = system_data.sort_values('Time (ns)')
        
        # Calculate cumulative average (running average from start)
        time = system_data['Time (ns)'].values
        pressure = system_data['Pressure (atm)'].values
        cumulative_avg = np.cumsum(pressure) / np.arange(1, len(pressure) + 1)
        
        ax.plot(time, cumulative_avg, color=COLOR_MAP[system], 
                linewidth=2, alpha=0.8, label=f"{system}")
    
    # Add horizontal line for target pressure
    ax.axhline(y=TARGET_PRESSURE["tip4p_273K"], color='black', linestyle='--', alpha=0.5, linewidth=1)
    
    # Add vertical grid lines
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.grid(True, which='major', axis='both', linestyle='--', alpha=0.7)
    
    ax.set_xlabel('Time (ns)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Cumulative Average Pressure (atm)', fontweight='bold', fontsize=12)
    ax.set_title('Pressure Convergence', fontsize=14, fontweight='bold')
    
    # Create a clean, bold legend
    legend = ax.legend(title="System", loc='best', framealpha=0.9, fontsize=12)
    plt.setp(legend.get_title(), fontweight='bold')
    plt.setp(legend.get_texts(), fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, 'pressure_convergence.png'), dpi=300)
    plt.close()
    
    # Calculate and save statistics
    pressure_stats = combined_data.groupby('System')['Pressure (atm)'].agg(['mean', 'std', 'min', 'max', 'count']).reset_index()
    pressure_stats.columns = ['System', 'Mean Pressure (atm)', 'Std Dev (atm)', 
                             'Min Pressure (atm)', 'Max Pressure (atm)', 'Sample Count']
    
    # Add target comparison
    pressure_stats['Target Pressure (atm)'] = pressure_stats['System'].map(
        lambda x: TARGET_PRESSURE.get(x, BAR_TO_ATM))
    pressure_stats['Deviation from Target (atm)'] = pressure_stats.apply(
        lambda row: (row['Mean Pressure (atm)'] - row['Target Pressure (atm)']), axis=1)
    pressure_stats['Relative Deviation (%)'] = pressure_stats.apply(
        lambda row: ((row['Mean Pressure (atm)'] - row['Target Pressure (atm)']) / 
                    row['Target Pressure (atm)'] * 100), axis=1)
    
    stats_file = os.path.join(comparison_dir, 'pressure_statistics.txt')
    
    with open(stats_file, 'w') as f:
        f.write("Pressure Statistics Comparison\n")
        f.write("=============================\n\n")
        
        # Write statistical significance test results
        f.write("Statistical Significance (p-values from t-test):\n")
        f.write("-----------------------------------------------\n")
        for i, sys1 in enumerate(systems):
            for j, sys2 in enumerate(systems):
                if i < j:
                    sys1_data = combined_data[combined_data['System'] == sys1]['Pressure (atm)'].values
                    sys2_data = combined_data[combined_data['System'] == sys2]['Pressure (atm)'].values
                    t_stat, p_value = scipy_stats.ttest_ind(sys1_data, sys2_data, equal_var=False)
                    sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
                    f.write(f"{sys1} vs {sys2}: p={p_value:.6f} {sig}\n")
        f.write("\n")
        
        # Write detailed statistics for each system
        for _, row in pressure_stats.iterrows():
            f.write(f"System: {row['System']}\n")
            f.write(f"  Mean Pressure: {row['Mean Pressure (atm)']:.4f} atm\n")
            f.write(f"  Std Deviation: {row['Std Dev (atm)']:.4f} atm\n")
            f.write(f"  Min Pressure: {row['Min Pressure (atm)']:.4f} atm\n")
            f.write(f"  Max Pressure: {row['Max Pressure (atm)']:.4f} atm\n")
            f.write(f"  Sample Count: {row['Sample Count']}\n")
            f.write(f"  Target Pressure: {row['Target Pressure (atm)']} atm\n")
            f.write(f"  Deviation from Target: {row['Deviation from Target (atm)']:.4f} atm\n")
            f.write(f"  Relative Deviation: {row['Relative Deviation (%)']:.4f}%\n")
            
            # Add note about pressure fluctuations
            rel_std = (row['Std Dev (atm)'] / abs(row['Mean Pressure (atm)'])) * 100
            f.write(f"  Relative Fluctuation: {rel_std:.2f}%\n")
            f.write("\n")
        
        # Add note about pressure fluctuations in NPT simulations
        f.write("\nNote about Pressure Fluctuations in NPT Simulations:\n")
        f.write("------------------------------------------------\n")
        f.write("Large pressure fluctuations are normal in NPT simulations, especially for small systems.\n")
        f.write("The fluctuations scale as 1/√N, where N is the number of particles.\n")
        f.write("Only the average pressure should converge to the target pressure (0.987 atm).\n")
        f.write("Standard deviations of hundreds of atm are typical for water systems.\n")
    
    print(f"Saved pressure comparison plots and statistics to {comparison_dir}/")
    print(f"Statistics saved to {stats_file}")

if __name__ == "__main__":
    import sys
    main() 