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
    "iteration_5_TIP4P": "#2ca02c" # Green
}

# Experimental reference values (g/cm³)
EXPERIMENTAL_DENSITIES = {
    "tip4p_273K": 0.9998,      # Water at 273K (0°C)
    "tip4p_298K": 0.9970,      # Water at 298K (25°C)
    "iteration_5_TIP4P": 0.9970  # Assuming same reference as 298K
}

def load_density_data(directory):
    """Load density data from a directory's analysis folder."""
    data_path = os.path.join(directory, "analysis", "density_vs_time.txt")
    
    if not os.path.exists(data_path):
        print(f"Warning: Data file not found at {data_path}")
        return None
    
    try:
        # Load data with pandas for easier handling
        data = pd.read_csv(data_path, sep='\t', comment='#', header=None)
        
        # Check if the data has the expected format (at least 2 columns)
        if data.shape[1] < 2:
            print(f"Warning: Data in {data_path} doesn't have enough columns")
            return None
            
        # Assign column names based on expected format
        if data.shape[1] >= 5:  # Full format with box dimensions
            data.columns = ['Time', 'Density', 'Box_X', 'Box_Y', 'Box_Z']
        else:  # Minimal format with just time and density
            data.columns = ['Time', 'Density']
            
        # Add directory name as a label
        data['System'] = os.path.basename(directory)
        
        return data
    except Exception as e:
        print(f"Error loading data from {data_path}: {e}")
        return None

def add_significance_bars(ax, data, systems, y_position):
    """Add significance bars between systems based on t-test."""
    pairs = []
    for i in range(len(systems)):
        for j in range(i+1, len(systems)):
            pairs.append((systems[i], systems[j]))
    
    for pair in pairs:
        sys1_data = data[data['System'] == pair[0]]['Density'].values
        sys2_data = data[data['System'] == pair[1]]['Density'].values
        
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
        ax.text((x1 + x2) / 2, y_position + 0.0005, sig_symbol, ha='center')

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
    print("Loading density data from all directories...")
    
    # Load data from each directory
    all_data = []
    for directory in directories:
        data = load_density_data(directory)
        if data is not None:
            all_data.append(data)
    
    if not all_data:
        print("No valid data found. Please run calculate_density.py in each directory first.")
        return
    
    # Combine all data into a single DataFrame
    combined_data = pd.concat(all_data, ignore_index=True)
    
    print("Generating comparison plots...")
    
    # Set a nice style for plots
    sns.set_style("whitegrid")
    plt.rcParams.update({
        'font.size': 12,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'figure.figsize': (12, 8),
        'figure.dpi': 100,
    })
    
    # Get unique systems in the order they appear in directories
    systems = [sys for sys in directories if sys in combined_data['System'].unique()]
    
    # Plot 1: Density vs Time for all systems
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Add vertical grid lines at major time intervals (every 1000 ps)
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.grid(True, which='major', axis='both', linestyle='--', alpha=0.7)
    
    # Plot raw data
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        ax.plot(system_data['Time'], system_data['Density'], 
                label=system, color=COLOR_MAP[system], alpha=0.8, linewidth=1)
    
    # Add smoothed trend lines
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        # Sort by time to ensure proper running average
        system_data = system_data.sort_values('Time')
        # Calculate running average
        smooth_time, smooth_density = calculate_running_average(
            system_data['Time'].values, system_data['Density'].values, window=100)
        ax.plot(smooth_time, smooth_density, color=COLOR_MAP[system], 
                linewidth=2.5, alpha=0.7, linestyle='-')
    
    # Add horizontal lines for experimental values
    for system in systems:
        if system in EXPERIMENTAL_DENSITIES:
            ax.axhline(y=EXPERIMENTAL_DENSITIES[system], color=COLOR_MAP[system], 
                      linestyle='--', alpha=0.7, linewidth=1.5)
    
    # Add annotations for experimental values with consistent decimal places
    legend_elements = []
    for system in systems:
        if system in EXPERIMENTAL_DENSITIES:
            # Format with consistent 4 decimal places
            exp_value = f"{EXPERIMENTAL_DENSITIES[system]:.4f}"
            legend_elements.append(
                Patch(facecolor='none', edgecolor=COLOR_MAP[system], linestyle='--',
                      label=f"{system} (Exp: {exp_value} g/cm³)")
            )
    
    # Create a second legend for experimental values
    if legend_elements:
        second_legend = ax.legend(handles=legend_elements, loc='lower right', 
                                 title="Experimental References", framealpha=0.7)
        ax.add_artist(second_legend)
    
    # Add primary legend for simulation data
    ax.legend(title="Simulation Data", loc='upper right', framealpha=0.7)
    
    ax.set_xlabel('Time (ps)', fontweight='bold')
    ax.set_ylabel('Density (g/cm³)', fontweight='bold')
    ax.set_title('Water Density Comparison Over Time', fontsize=14, fontweight='bold')
    
    # Add an inset to zoom in on a region of interest
    # Find a good region to zoom in - last 20% of the simulation
    time_max = combined_data['Time'].max()
    zoom_start = time_max * 0.8
    
    # Create inset axes with proper parameters
    axins = ax.inset_axes([0.15, 0.15, 0.3, 0.3])
    
    # Add subtle background shading to the inset
    axins.patch.set_facecolor('lightgray')
    axins.patch.set_alpha(0.1)
    
    # Add a border to the inset
    for spine in axins.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)
    
    # Plot data in the inset
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        zoom_data = system_data[system_data['Time'] >= zoom_start]
        axins.plot(zoom_data['Time'], zoom_data['Density'], 
                  color=COLOR_MAP[system], linewidth=1.5)
    
    axins.set_title('Last 20% of Simulation', fontsize=10, fontweight='bold')
    axins.grid(True, linestyle='--', alpha=0.5)
    
    # Add annotation to highlight the inset with a more prominent arrow
    ax.annotate('Zoomed region', xy=(zoom_start, combined_data['Density'].min()),
               xytext=(zoom_start * 0.7, combined_data['Density'].min()),
               arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=10),
               fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, 'density_comparison.png'), dpi=300)
    plt.close()
    
    # Plot 2: Density distribution (KDE plot) for all systems
    fig, ax = plt.subplots(figsize=(12, 8))
    
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        sns.kdeplot(system_data['Density'], label=system, fill=True, 
                   alpha=0.3, color=COLOR_MAP[system], ax=ax)
        
        # Add vertical line for experimental value
        if system in EXPERIMENTAL_DENSITIES:
            ax.axvline(x=EXPERIMENTAL_DENSITIES[system], color=COLOR_MAP[system], 
                      linestyle='--', alpha=0.7, linewidth=1.5)
            
            # Add annotation for experimental value
            ax.annotate(f"Exp: {EXPERIMENTAL_DENSITIES[system]:.4f}", 
                       xy=(EXPERIMENTAL_DENSITIES[system], 0), 
                       xytext=(EXPERIMENTAL_DENSITIES[system], 1),
                       rotation=90, ha='right', va='bottom',
                       color=COLOR_MAP[system], fontsize=10)
    
    ax.set_xlabel('Density (g/cm³)', fontweight='bold')
    ax.set_ylabel('Probability Density', fontweight='bold')
    ax.set_title('Water Density Distribution Comparison', fontsize=14, fontweight='bold')
    ax.legend(title="System", loc='upper right', framealpha=0.7)
    ax.grid(True, linestyle='--', alpha=0.7)
    
    # Add text with statistical information
    stats_text = ""
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        mean = system_data['Density'].mean()
        std = system_data['Density'].std()
        stats_text += f"{system}: {mean:.4f} ± {std:.4f} g/cm³\n"
    
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
           verticalalignment='top', horizontalalignment='left',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, 'density_distribution_comparison.png'), dpi=300)
    plt.close()
    
    # Plot 3: Box volume comparison (if available)
    if 'Box_X' in combined_data.columns:
        # Create figure with shared y-axis
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8), 
                                      gridspec_kw={'width_ratios': [3, 1]},
                                      sharey=True)  # Share y-axis for better comparison
        
        # Calculate volumes for all systems
        volume_data = []
        mean_volumes = {}
        std_volumes = {}
        
        for system in systems:
            system_data = combined_data[combined_data['System'] == system]
            # Calculate volume in nm³
            system_volume = system_data['Box_X'] * system_data['Box_Y'] * system_data['Box_Z']
            
            # Store mean and std for each system
            mean_volumes[system] = system_volume.mean()
            std_volumes[system] = system_volume.std()
            
            # Store data for violin plot
            temp_df = pd.DataFrame({
                'System': system,
                'Volume': system_volume
            })
            volume_data.append(temp_df)
        
        # Combine all volume data
        volume_df = pd.concat(volume_data, ignore_index=True)
        
        # Add vertical grid lines at major time intervals (every 1000 ps)
        ax1.xaxis.set_major_locator(MultipleLocator(1000))
        ax1.grid(True, which='major', axis='both', linestyle='--', alpha=0.7)
        
        # Time series plot
        for system in systems:
            system_data = combined_data[combined_data['System'] == system]
            # Calculate volume in nm³
            volume = system_data['Box_X'] * system_data['Box_Y'] * system_data['Box_Z']
            ax1.plot(system_data['Time'], volume, label=system, 
                    color=COLOR_MAP[system], alpha=0.8, linewidth=1.5)
            
            # Add horizontal line for mean volume
            ax1.axhline(y=mean_volumes[system], color=COLOR_MAP[system], 
                       linestyle='-', alpha=0.3, linewidth=2)
        
        # Add legend for time series
        ax1.legend(title="System", loc='upper right', framealpha=0.7)
        
        # Set labels for time series plot
        ax1.set_xlabel('Time (ps)', fontweight='bold')
        ax1.set_ylabel('Box Volume (nm³)', fontweight='bold')
        ax1.set_title('Simulation Box Volume Over Time', fontsize=14, fontweight='bold')
        
        # Violin plot for volume distribution
        try:
            # Try with hue parameter (newer seaborn versions)
            sns.violinplot(x='System', y='Volume', data=volume_df, ax=ax2, 
                          palette=COLOR_MAP, hue='System', legend=False, inner='quartile')
        except:
            # Fallback for older seaborn versions
            sns.violinplot(x='System', y='Volume', data=volume_df, ax=ax2, 
                          palette=COLOR_MAP, inner='quartile')
        
        # Set labels for distribution plot
        ax2.set_title('Volume Distribution', fontsize=14, fontweight='bold')
        ax2.set_xlabel('')
        ax2.set_ylabel('')  # No need for y-label since it's shared
        ax2.grid(True, axis='y', linestyle='--', alpha=0.7)
        
        # Add horizontal lines for mean volumes on violin plot
        for i, system in enumerate(systems):
            ax2.axhline(y=mean_volumes[system], xmin=i/len(systems), xmax=(i+1)/len(systems),
                       color=COLOR_MAP[system], linestyle='-', alpha=0.3, linewidth=2)
        
        # Add annotations with mean volumes and standard deviations
        for i, system in enumerate(systems):
            mean_vol = mean_volumes[system]
            std_vol = std_volumes[system]
            
            # Position the text to avoid truncation
            y_offset = 15 if i == 0 else 10  # Larger offset for first system to avoid truncation
            
            # Create annotation with mean and std
            ax2.annotate(f"Mean: {mean_vol:.2f} nm³\nStd: {std_vol:.2f} nm³", 
                        xy=(i, mean_vol),
                        xytext=(0, y_offset), textcoords='offset points',
                        ha='center', va='bottom',
                        bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.9))
        
        # Adjust layout to prevent text truncation
        plt.tight_layout(pad=2.0)
        plt.savefig(os.path.join(comparison_dir, 'volume_comparison.png'), dpi=300)
        plt.close()
    
    # Plot 4: Boxplot of density values with statistical significance
    fig, ax = plt.subplots(figsize=(10, 8))
    
    try:
        # Try with hue parameter (newer seaborn versions)
        sns.boxplot(x='System', y='Density', data=combined_data, ax=ax,
                   palette=COLOR_MAP, hue='System', legend=False, notch=True, showfliers=False)
    except:
        # Fallback for older seaborn versions
        sns.boxplot(x='System', y='Density', data=combined_data, ax=ax,
                   palette=COLOR_MAP, notch=True, showfliers=False)
    
    # Add individual data points with jitter and transparency
    sns.stripplot(x='System', y='Density', data=combined_data, ax=ax,
                 size=4, color='black', alpha=0.3, jitter=True)
    
    # Add significance bars
    max_density = combined_data['Density'].max()
    add_significance_bars(ax, combined_data, systems, max_density * 1.01)
    
    # Add experimental reference lines
    for i, system in enumerate(systems):
        if system in EXPERIMENTAL_DENSITIES:
            ax.hlines(y=EXPERIMENTAL_DENSITIES[system], xmin=i-0.4, xmax=i+0.4, 
                     colors=COLOR_MAP[system], linestyles='dashed', linewidth=2, alpha=0.7)
            ax.text(i, EXPERIMENTAL_DENSITIES[system] - 0.001, 
                   f"Exp: {EXPERIMENTAL_DENSITIES[system]:.4f}", 
                   ha='center', va='top', fontsize=9, 
                   bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.2'))
    
    ax.set_xlabel('System', fontweight='bold')
    ax.set_ylabel('Density (g/cm³)', fontweight='bold')
    ax.set_title('Density Distribution by System', fontsize=14, fontweight='bold')
    ax.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add legend for significance
    legend_elements = [
        Patch(facecolor='none', edgecolor='black', label='* p < 0.05'),
        Patch(facecolor='none', edgecolor='black', label='** p < 0.01'),
        Patch(facecolor='none', edgecolor='black', label='*** p < 0.001'),
        Patch(facecolor='none', edgecolor='black', label='ns: not significant')
    ]
    ax.legend(handles=legend_elements, loc='upper right', 
             title="Statistical Significance", framealpha=0.7)
    
    plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, 'density_boxplot.png'), dpi=300)
    plt.close()
    
    # Calculate and save statistics
    density_stats = combined_data.groupby('System')['Density'].agg(['mean', 'std', 'min', 'max', 'count']).reset_index()
    density_stats.columns = ['System', 'Mean Density (g/cm³)', 'Std Dev (g/cm³)', 
                    'Min Density (g/cm³)', 'Max Density (g/cm³)', 'Sample Count']
    
    # Add experimental comparison
    density_stats['Experimental Value (g/cm³)'] = density_stats['System'].map(
        lambda x: EXPERIMENTAL_DENSITIES.get(x, 'N/A'))
    density_stats['Deviation from Exp (%)'] = density_stats.apply(
        lambda row: ((row['Mean Density (g/cm³)'] - float(row['Experimental Value (g/cm³)'])) / 
                    float(row['Experimental Value (g/cm³)']) * 100) 
        if row['Experimental Value (g/cm³)'] != 'N/A' else 'N/A', axis=1)
    
    stats_file = os.path.join(comparison_dir, 'density_statistics.txt')
    
    with open(stats_file, 'w') as f:
        f.write("Density Statistics Comparison\n")
        f.write("============================\n\n")
        
        # Write statistical significance test results
        f.write("Statistical Significance (p-values from t-test):\n")
        f.write("-----------------------------------------------\n")
        for i, sys1 in enumerate(systems):
            for j, sys2 in enumerate(systems):
                if i < j:
                    sys1_data = combined_data[combined_data['System'] == sys1]['Density'].values
                    sys2_data = combined_data[combined_data['System'] == sys2]['Density'].values
                    t_stat, p_value = scipy_stats.ttest_ind(sys1_data, sys2_data, equal_var=False)
                    sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
                    f.write(f"{sys1} vs {sys2}: p={p_value:.6f} {sig}\n")
        f.write("\n")
        
        # Write detailed statistics for each system
        for _, row in density_stats.iterrows():
            f.write(f"System: {row['System']}\n")
            f.write(f"  Mean Density: {row['Mean Density (g/cm³)']:.4f} g/cm³\n")
            f.write(f"  Std Deviation: {row['Std Dev (g/cm³)']:.4f} g/cm³\n")
            f.write(f"  Min Density: {row['Min Density (g/cm³)']:.4f} g/cm³\n")
            f.write(f"  Max Density: {row['Max Density (g/cm³)']:.4f} g/cm³\n")
            f.write(f"  Sample Count: {row['Sample Count']}\n")
            f.write(f"  Experimental Value: {row['Experimental Value (g/cm³)']} g/cm³\n")
            if row['Deviation from Exp (%)'] != 'N/A':
                f.write(f"  Deviation from Experimental: {row['Deviation from Exp (%)']:.4f}%\n")
            else:
                f.write(f"  Deviation from Experimental: N/A\n")
            f.write("\n")
    
    print(f"Saved enhanced comparison plots and statistics to {comparison_dir}/")
    print(f"Statistics saved to {stats_file}")

if __name__ == "__main__":
    main() 