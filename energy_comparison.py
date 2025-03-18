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

# Energy components to analyze
ENERGY_COMPONENTS = ['Potential', 'Kinetic', 'Total-Energy']

def load_energy_data(directory):
    """Load energy data from a directory's analysis folder."""
    # First try to load the statistics file to get mean and std
    stats_path = os.path.join(directory, "analysis", "energy_stats.dat")
    energy_stats = {}
    
    if os.path.exists(stats_path):
        try:
            with open(stats_path, 'r') as f:
                for line in f:
                    if "Potential energy:" in line:
                        parts = line.split(':')[1].strip().split('±')
                        energy_stats['potential_mean'] = float(parts[0].strip())
                        energy_stats['potential_std'] = float(parts[1].strip().split()[0])
                    elif "Kinetic energy:" in line:
                        parts = line.split(':')[1].strip().split('±')
                        energy_stats['kinetic_mean'] = float(parts[0].strip())
                        energy_stats['kinetic_std'] = float(parts[1].strip().split()[0])
                    elif "Total energy:" in line:
                        parts = line.split(':')[1].strip().split('±')
                        energy_stats['total_mean'] = float(parts[0].strip())
                        energy_stats['total_std'] = float(parts[1].strip().split()[0])
                    elif "Energy drift:" in line:
                        energy_stats['drift'] = float(line.split(':')[1].split()[0])
        except Exception as e:
            print(f"Error reading energy stats from {stats_path}: {e}")
    
    # Now try to load the raw data if available
    data_path = os.path.join(directory, "analysis", "energy_data.dat")
    
    if os.path.exists(data_path):
        try:
            # Try to load the data file
            data = pd.read_csv(data_path, sep='\t', comment='#')
            
            # Check if the data has the expected format
            expected_columns = ['Time (ps)', 'Potential', 'Kinetic', 'Total']
            if all(col in data.columns for col in expected_columns):
                # Add directory name as a label
                data['System'] = os.path.basename(directory)
                return data, energy_stats
        except Exception as e:
            print(f"Error loading energy data from {data_path}: {e}")
    
    # If we couldn't load the data file, try to extract from the EDR file
    # This is a fallback and requires the thermodynamics_analysis.py script
    print(f"No energy data file found for {directory}, trying to extract from EDR file...")
    
    # Try to import the parse_edr_file function from the thermodynamics_analysis.py script
    try:
        sys.path.append(os.path.join(directory, "scripts"))
        from utils import parse_edr_file
        
        # Parse energy data from EDR file
        edr_file = os.path.join(directory, "md.edr")
        if os.path.exists(edr_file):
            energy_data = parse_edr_file(edr_file, properties=ENERGY_COMPONENTS)
            
            if all(comp in energy_data for comp in ENERGY_COMPONENTS):
                # Combine data into a single DataFrame
                data = pd.DataFrame({
                    'Time (ps)': energy_data['Potential']['Time (ps)'],
                    'Potential': energy_data['Potential']['Potential'],
                    'Kinetic': energy_data['Kinetic']['Kinetic'],
                    'Total': energy_data['Total-Energy']['Total-Energy']
                })
                
                # Add directory name as a label
                data['System'] = os.path.basename(directory)
                return data, energy_stats
    except Exception as e:
        print(f"Error extracting energy data from EDR file: {e}")
    
    # If we still don't have data, return None
    print(f"Warning: Could not load energy data for {directory}")
    return None, energy_stats

def calculate_running_average(x, y, window=50):
    """Calculate running average with specified window size."""
    avg_y = np.zeros_like(y)
    n = len(y)
    
    for i in range(n):
        start = max(0, i - window // 2)
        end = min(n, i + window // 2 + 1)
        avg_y[i] = np.mean(y[start:end])
    
    return x, avg_y

def plot_energy_component_comparison(combined_data, systems, component, output_file):
    """Plot comparison of a specific energy component across systems."""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Add vertical grid lines at major time intervals (every 1000 ps)
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.grid(True, which='major', axis='both', linestyle='--', alpha=0.7)
    
    # Plot raw data with low alpha for clarity
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        ax.plot(system_data['Time (ps)'], system_data[component], 
                label=system, color=COLOR_MAP[system], alpha=0.3, linewidth=0.5)
    
    # Add smoothed trend lines
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        # Sort by time to ensure proper running average
        system_data = system_data.sort_values('Time (ps)')
        # Calculate running average
        smooth_time, smooth_energy = calculate_running_average(
            system_data['Time (ps)'].values, system_data[component].values, window=100)
        ax.plot(smooth_time, smooth_energy, color=COLOR_MAP[system], 
                linewidth=2.5, alpha=0.8, linestyle='-')
    
    # Add legend
    ax.legend(title="System", loc='best', framealpha=0.7)
    
    ax.set_xlabel('Time (ps)', fontweight='bold')
    ax.set_ylabel(f'{component} Energy (kJ/mol)', fontweight='bold')
    ax.set_title(f'{component} Energy Comparison Over Time', fontsize=14, fontweight='bold')
    
    # Add an inset to zoom in on a region of interest
    # Find a good region to zoom in - last 20% of the simulation
    time_max = combined_data['Time (ps)'].max()
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
        zoom_data = system_data[system_data['Time (ps)'] >= zoom_start]
        axins.plot(zoom_data['Time (ps)'], zoom_data[component], 
                  color=COLOR_MAP[system], linewidth=1.5)
    
    axins.set_title('Last 20% of Simulation', fontsize=10, fontweight='bold')
    axins.grid(True, linestyle='--', alpha=0.5)
    
    # Add annotation to highlight the inset with a more prominent arrow
    y_min = combined_data[component].min()
    ax.annotate('Zoomed region', xy=(zoom_start, y_min),
               xytext=(zoom_start * 0.7, y_min),
               arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=10),
               fontsize=11, fontweight='bold')
    
    # Add text with statistical information
    stats_text = ""
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        mean = system_data[component].mean()
        std = system_data[component].std()
        stats_text += f"{system}: {mean:.2f} ± {std:.2f} kJ/mol\n"
    
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
           verticalalignment='top', horizontalalignment='left',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def plot_energy_distribution_comparison(combined_data, systems, component, output_file):
    """Plot distribution comparison of a specific energy component across systems."""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        sns.kdeplot(system_data[component], label=system, fill=True, 
                   alpha=0.3, color=COLOR_MAP[system], ax=ax)
        
        # Add vertical line for mean energy
        mean_energy = system_data[component].mean()
        ax.axvline(x=mean_energy, color=COLOR_MAP[system], 
                  linestyle='-', alpha=0.7, linewidth=1.5)
        
        # Add annotation for mean energy
        ax.annotate(f"Mean: {mean_energy:.2f} kJ/mol", 
                   xy=(mean_energy, 0), 
                   xytext=(mean_energy, 0.0001),
                   rotation=90, ha='right', va='bottom',
                   color=COLOR_MAP[system], fontsize=10)
    
    ax.set_xlabel(f'{component} Energy (kJ/mol)', fontweight='bold')
    ax.set_ylabel('Probability Density', fontweight='bold')
    ax.set_title(f'{component} Energy Distribution Comparison', fontsize=14, fontweight='bold')
    ax.legend(title="System", loc='best', framealpha=0.7)
    ax.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def plot_energy_conservation(combined_data, systems, output_file):
    """Plot energy conservation (drift) comparison across systems."""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Add vertical grid lines at major time intervals (every 1000 ps)
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.grid(True, which='major', axis='both', linestyle='--', alpha=0.7)
    
    # Calculate and plot energy drift for each system
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        # Sort by time
        system_data = system_data.sort_values('Time (ps)')
        
        # Calculate drift as deviation from initial total energy
        if len(system_data) > 0:
            time = system_data['Time (ps)'].values
            total_energy = system_data['Total'].values
            initial_energy = total_energy[0]
            drift = total_energy - initial_energy
            
            # Convert time to ns for drift calculation
            time_ns = time / 1000
            
            # Linear regression to find the slope (drift rate)
            slope, intercept, r_value, p_value, std_err = scipy_stats.linregress(time_ns, total_energy)
            
            # Plot drift
            ax.plot(time, drift, color=COLOR_MAP[system], alpha=0.3, linewidth=0.5)
            
            # Plot drift trend line
            drift_line = (intercept - initial_energy) + slope * time_ns
            ax.plot(time, drift_line, color=COLOR_MAP[system], 
                    linewidth=2, alpha=0.8, linestyle='-',
                    label=f"{system}: {slope:.4f} kJ/mol/ns")
    
    ax.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    
    ax.set_xlabel('Time (ps)', fontweight='bold')
    ax.set_ylabel('Energy Drift (kJ/mol)', fontweight='bold')
    ax.set_title('Energy Conservation Comparison', fontsize=14, fontweight='bold')
    ax.legend(title="Drift Rate", loc='best', framealpha=0.7)
    
    # Add annotation about energy conservation
    ax.text(0.02, 0.02, 
           "Note: Good energy conservation is indicated by minimal drift.\n"
           "The slope of the trend line represents the drift rate in kJ/mol/ns.",
           transform=ax.transAxes, 
           verticalalignment='bottom', horizontalalignment='left',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def plot_energy_fluctuations(combined_data, systems, output_file):
    """Plot energy fluctuations comparison across systems."""
    # Create figure with 3 subplots (one for each energy component)
    fig, axs = plt.subplots(3, 1, figsize=(12, 15), sharex=True)
    
    components = ['Potential', 'Kinetic', 'Total']
    titles = ['Potential Energy Fluctuations', 'Kinetic Energy Fluctuations', 'Total Energy Fluctuations']
    
    for i, (component, title) in enumerate(zip(components, titles)):
        ax = axs[i]
        
        # Add vertical grid lines
        ax.xaxis.set_major_locator(MultipleLocator(1000))
        ax.grid(True, which='major', axis='both', linestyle='--', alpha=0.7)
        
        # Calculate and plot fluctuations for each system
        for system in systems:
            system_data = combined_data[combined_data['System'] == system]
            # Sort by time
            system_data = system_data.sort_values('Time (ps)')
            
            if len(system_data) > 0:
                time = system_data['Time (ps)'].values
                energy = system_data[component].values
                
                # For total energy, detrend to show fluctuations around the drift line
                if component == 'Total':
                    # Calculate drift line
                    time_ns = time / 1000
                    slope, intercept, _, _, _ = scipy_stats.linregress(time_ns, energy)
                    drift_line = intercept + slope * time_ns
                    # Detrend
                    energy = energy - drift_line
                
                # Calculate running average
                _, smooth_energy = calculate_running_average(time, energy, window=100)
                
                # Calculate fluctuations (deviations from running average)
                fluctuations = energy - smooth_energy
                
                # Plot fluctuations
                ax.plot(time, fluctuations, color=COLOR_MAP[system], alpha=0.5, linewidth=0.5)
                
                # Calculate and plot standard deviation bands
                std = np.std(fluctuations)
                ax.axhline(y=std, color=COLOR_MAP[system], linestyle='--', alpha=0.5)
                ax.axhline(y=-std, color=COLOR_MAP[system], linestyle='--', alpha=0.5)
                
                # Fill between std bands
                ax.fill_between(time, -std, std, color=COLOR_MAP[system], alpha=0.1)
                
                # Add label with std info
                if i == 0:  # Only add to legend in first subplot
                    ax.plot([], [], color=COLOR_MAP[system], linewidth=2, 
                           label=f"{system}: σ = {std:.2f} kJ/mol")
        
        ax.set_ylabel('Fluctuation (kJ/mol)', fontweight='bold')
        ax.set_title(title, fontsize=12, fontweight='bold')
        
        if i == 0:  # Only add legend to first subplot
            ax.legend(title="Standard Deviation", loc='upper right', framealpha=0.7)
    
    # Set common x label
    axs[2].set_xlabel('Time (ps)', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def plot_energy_components_stacked(combined_data, systems, output_file):
    """Plot stacked energy components for each system."""
    # Create a figure with subplots for each system
    fig, axs = plt.subplots(len(systems), 1, figsize=(12, 4*len(systems)), sharex=True)
    
    # If only one system, axs will not be an array
    if len(systems) == 1:
        axs = [axs]
    
    for i, system in enumerate(systems):
        ax = axs[i]
        system_data = combined_data[combined_data['System'] == system]
        
        # Sort by time
        system_data = system_data.sort_values('Time (ps)')
        
        if len(system_data) > 0:
            time = system_data['Time (ps)'].values
            potential = system_data['Potential'].values
            kinetic = system_data['Kinetic'].values
            total = system_data['Total'].values
            
            # Plot components
            ax.plot(time, potential, 'r-', alpha=0.7, label='Potential')
            ax.plot(time, kinetic, 'b-', alpha=0.7, label='Kinetic')
            ax.plot(time, total, 'g-', alpha=0.7, label='Total')
            
            # Add horizontal lines for means
            ax.axhline(y=np.mean(potential), color='r', linestyle='--', alpha=0.5)
            ax.axhline(y=np.mean(kinetic), color='b', linestyle='--', alpha=0.5)
            ax.axhline(y=np.mean(total), color='g', linestyle='--', alpha=0.5)
            
            # Add grid
            ax.grid(True, linestyle='--', alpha=0.7)
            
            # Add title and legend
            ax.set_title(f'{system} Energy Components', fontsize=12, fontweight='bold')
            ax.legend(loc='upper right')
            
            # Add text with statistics
            stats_text = (f"Potential: {np.mean(potential):.2f} ± {np.std(potential):.2f} kJ/mol\n"
                         f"Kinetic: {np.mean(kinetic):.2f} ± {np.std(kinetic):.2f} kJ/mol\n"
                         f"Total: {np.mean(total):.2f} ± {np.std(total):.2f} kJ/mol")
            
            ax.text(0.02, 0.95, stats_text, transform=ax.transAxes, 
                   verticalalignment='top', horizontalalignment='left',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            # Set y-label
            ax.set_ylabel('Energy (kJ/mol)', fontweight='bold')
    
    # Set common x-label
    axs[-1].set_xlabel('Time (ps)', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def main():
    print("Loading energy data from all directories...")
    
    # Load data from each directory
    all_data = []
    all_stats = {}
    
    for directory in directories:
        data, stats = load_energy_data(directory)
        if data is not None:
            all_data.append(data)
            all_stats[directory] = stats
    
    if not all_data:
        print("No valid energy data found. Please run thermodynamics_analysis.py in each directory first.")
        return
    
    # Combine all data into a single DataFrame
    combined_data = pd.concat(all_data, ignore_index=True)
    
    print("Generating energy comparison plots...")
    
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
    
    # Plot 1: Potential energy comparison
    plot_energy_component_comparison(
        combined_data, systems, 'Potential', 
        os.path.join(comparison_dir, 'potential_energy_comparison.png')
    )
    
    # Plot 2: Kinetic energy comparison
    plot_energy_component_comparison(
        combined_data, systems, 'Kinetic', 
        os.path.join(comparison_dir, 'kinetic_energy_comparison.png')
    )
    
    # Plot 3: Total energy comparison
    plot_energy_component_comparison(
        combined_data, systems, 'Total', 
        os.path.join(comparison_dir, 'total_energy_comparison.png')
    )
    
    # Plot 4: Potential energy distribution
    plot_energy_distribution_comparison(
        combined_data, systems, 'Potential', 
        os.path.join(comparison_dir, 'potential_energy_distribution.png')
    )
    
    # Plot 5: Kinetic energy distribution
    plot_energy_distribution_comparison(
        combined_data, systems, 'Kinetic', 
        os.path.join(comparison_dir, 'kinetic_energy_distribution.png')
    )
    
    # Plot 6: Total energy distribution
    plot_energy_distribution_comparison(
        combined_data, systems, 'Total', 
        os.path.join(comparison_dir, 'total_energy_distribution.png')
    )
    
    # Plot 7: Energy conservation comparison
    plot_energy_conservation(
        combined_data, systems, 
        os.path.join(comparison_dir, 'energy_conservation_comparison.png')
    )
    
    # Plot 8: Energy fluctuations comparison
    plot_energy_fluctuations(
        combined_data, systems, 
        os.path.join(comparison_dir, 'energy_fluctuations_comparison.png')
    )
    
    # Plot 9: Stacked energy components for each system
    plot_energy_components_stacked(
        combined_data, systems, 
        os.path.join(comparison_dir, 'energy_components_stacked.png')
    )
    
    # Calculate and save statistics
    energy_stats = {}
    
    for component in ['Potential', 'Kinetic', 'Total']:
        # Calculate statistics for each component
        component_stats = combined_data.groupby('System')[component].agg(
            ['mean', 'std', 'min', 'max', 'count']).reset_index()
        component_stats.columns = ['System', f'Mean {component} (kJ/mol)', f'Std Dev {component} (kJ/mol)', 
                                  f'Min {component} (kJ/mol)', f'Max {component} (kJ/mol)', 'Sample Count']
        energy_stats[component] = component_stats
    
    # Calculate energy drift for each system
    drift_stats = []
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        # Sort by time
        system_data = system_data.sort_values('Time (ps)')
        
        if len(system_data) > 0:
            time = system_data['Time (ps)'].values
            total_energy = system_data['Total'].values
            
            # Convert time to ns for drift calculation
            time_ns = time / 1000
            
            # Linear regression to find the slope (drift rate)
            slope, intercept, r_value, p_value, std_err = scipy_stats.linregress(time_ns, total_energy)
            
            # Calculate relative drift
            mean_energy = np.mean(total_energy)
            rel_drift = (slope / abs(mean_energy)) * 100  # percent per ns
            
            drift_stats.append({
                'System': system,
                'Drift Rate (kJ/mol/ns)': slope,
                'R-squared': r_value**2,
                'Relative Drift (%/ns)': rel_drift
            })
    
    drift_df = pd.DataFrame(drift_stats)
    
    # Save statistics to file
    stats_file = os.path.join(comparison_dir, 'energy_statistics.txt')
    
    with open(stats_file, 'w') as f:
        f.write("Energy Statistics Comparison\n")
        f.write("===========================\n\n")
        
        # Write energy drift statistics
        f.write("Energy Drift Statistics:\n")
        f.write("----------------------\n")
        for _, row in drift_df.iterrows():
            f.write(f"System: {row['System']}\n")
            f.write(f"  Drift Rate: {row['Drift Rate (kJ/mol/ns)']:.6f} kJ/mol/ns\n")
            f.write(f"  R-squared: {row['R-squared']:.6f}\n")
            f.write(f"  Relative Drift: {row['Relative Drift (%/ns)']:.6f}%/ns\n")
            
            # Add qualitative assessment of energy conservation
            rel_drift = abs(row['Relative Drift (%/ns)'])
            if rel_drift < 0.001:
                quality = "Excellent"
            elif rel_drift < 0.01:
                quality = "Very good"
            elif rel_drift < 0.1:
                quality = "Good"
            elif rel_drift < 1.0:
                quality = "Acceptable"
            else:
                quality = "Poor"
            
            f.write(f"  Energy Conservation Quality: {quality}\n\n")
        
        # Write detailed statistics for each energy component
        for component in ['Potential', 'Kinetic', 'Total']:
            f.write(f"\n{component} Energy Statistics:\n")
            f.write("-" * (len(component) + 19) + "\n")
            
            for _, row in energy_stats[component].iterrows():
                f.write(f"System: {row['System']}\n")
                f.write(f"  Mean: {row[f'Mean {component} (kJ/mol)']:.4f} kJ/mol\n")
                f.write(f"  Std Deviation: {row[f'Std Dev {component} (kJ/mol)']:.4f} kJ/mol\n")
                f.write(f"  Min: {row[f'Min {component} (kJ/mol)']:.4f} kJ/mol\n")
                f.write(f"  Max: {row[f'Max {component} (kJ/mol)']:.4f} kJ/mol\n")
                f.write(f"  Sample Count: {row['Sample Count']}\n")
                
                # Calculate relative fluctuation
                rel_fluct = (row[f'Std Dev {component} (kJ/mol)'] / 
                            abs(row[f'Mean {component} (kJ/mol)'])) * 100
                f.write(f"  Relative Fluctuation: {rel_fluct:.2f}%\n\n")
    
    print(f"Saved energy comparison plots and statistics to {comparison_dir}/")
    print(f"Statistics saved to {stats_file}")

if __name__ == "__main__":
    import sys
    main() 