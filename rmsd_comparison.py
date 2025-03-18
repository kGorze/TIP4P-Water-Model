#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from scipy import stats as scipy_stats
from matplotlib.ticker import MultipleLocator
import sys
import MDAnalysis as mda
from MDAnalysis.analysis import rms

# Define directories to compare
directories = [
    "tip4p_273K",
    "tip4p_298K",
    "iteration_5_TIP4P"
]

# Create comparison directory if it doesn't exist
comparison_dir = "comparison"
rmsd_dir = os.path.join(comparison_dir, "rmsd")
os.makedirs(rmsd_dir, exist_ok=True)

# Define a consistent color scheme for each system
COLOR_MAP = {
    "tip4p_273K": "#1f77b4",      # Blue
    "tip4p_298K": "#ff7f0e",      # Orange
    "iteration_5_TIP4P": "#2ca02c" # Green
}

def calculate_rmsd_direct(universe, reference_frame=0, selection='name OW'):
    """
    Calculate RMSD directly without saving to a file.
    This is a simplified version of the calculate_rmsd function from rmsd_analysis.py.
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory
    reference_frame : int
        Frame to use as reference
    selection : str
        Atom selection string
        
    Returns:
    --------
    tuple
        (time_array, rmsd_array)
    """
    print("Calculating RMSD...")
    
    # Select atoms
    atoms = universe.select_atoms(selection)
    print(f"Selected {len(atoms)} atoms for RMSD calculation")
    
    # Set reference coordinates
    universe.trajectory[reference_frame]
    reference_coords = atoms.positions.copy()
    
    # Calculate RMSD for each frame
    time_array = []
    rmsd_array = []
    
    for ts in universe.trajectory:
        # Print progress every 100 frames
        if ts.frame % 100 == 0:
            print(f"Processing frame {ts.frame}/{len(universe.trajectory)}")
            
        time_array.append(ts.time)
        rmsd_array.append(rms.rmsd(atoms.positions, reference_coords, superposition=True))
    
    return np.array(time_array), np.array(rmsd_array)

def load_rmsd_data(directory):
    """Load RMSD data from a directory's analysis folder or calculate it if not available."""
    # First try to load pre-calculated RMSD data if it exists
    rmsd_file = os.path.join(directory, "analysis", "rmsd_data.dat")
    
    if os.path.exists(rmsd_file):
        try:
            # Try to load the data file
            data = pd.read_csv(rmsd_file, sep='\t', comment='#', header=None)
            
            # Check if the data has the expected format (at least 2 columns)
            if data.shape[1] >= 2:
                # Assign column names based on expected format
                data.columns = ['Time', 'RMSD']
                
                # Add directory name as a label
                data['System'] = os.path.basename(directory)
                print(f"Successfully loaded RMSD data for {directory}")
                return data
            else:
                print(f"Warning: RMSD data in {rmsd_file} doesn't have enough columns")
        except Exception as e:
            print(f"Error loading RMSD data from {rmsd_file}: {e}")
    
    # If we couldn't load pre-calculated data, calculate it now
    print(f"No pre-calculated RMSD data found for {directory}, calculating now...")
    
    # Load the trajectory
    tpr_file = os.path.join(directory, "md.tpr")
    trr_file = os.path.join(directory, "md.trr")
    
    if not os.path.exists(tpr_file) or not os.path.exists(trr_file):
        print(f"Trajectory files not found for {directory}")
        return None
    
    try:
        # Load the universe
        u = mda.Universe(tpr_file, trr_file)
        
        # Calculate RMSD directly
        time_array, rmsd_array = calculate_rmsd_direct(u)
        
        # Create a DataFrame
        data = pd.DataFrame({
            'Time': time_array,
            'RMSD': rmsd_array
        })
        
        # Add directory name as a label
        data['System'] = os.path.basename(directory)
        
        return data
    except Exception as e:
        print(f"Error calculating RMSD for {directory}: {e}")
    
    # If we still don't have data, return None
    print(f"Warning: Could not load or calculate RMSD data for {directory}")
    return None

def detect_equilibration(time_array, rmsd_array, window_size=100):
    """
    Detect equilibration point in RMSD data.
    
    Parameters:
    -----------
    time_array : numpy.ndarray
        Array of time values
    rmsd_array : numpy.ndarray
        Array of RMSD values
    window_size : int
        Size of the sliding window for variance calculation
    
    Returns:
    --------
    tuple
        (equilibration_time, equilibration_index)
    """
    # Ensure arrays are numpy arrays
    time_array = np.array(time_array)
    rmsd_array = np.array(rmsd_array)
    
    # Calculate variance in sliding windows
    variances = []
    for i in range(len(rmsd_array) - window_size):
        window = rmsd_array[i:i+window_size]
        variances.append(np.var(window))
    
    # Find the point where variance stabilizes
    # We'll use a simple approach: find where the variance drops below
    # a threshold relative to the maximum variance
    if len(variances) > 0:
        max_var = max(variances)
        threshold = max_var * 0.1  # 10% of maximum variance
        
        for i, var in enumerate(variances):
            if var < threshold:
                # Found equilibration point
                equil_index = i
                equil_time = time_array[i] if i < len(time_array) else time_array[-1]
                return equil_time, equil_index
    
    # If no clear equilibration point is found, return the last 20% of the trajectory
    equil_index = int(len(time_array) * 0.8)
    equil_time = time_array[equil_index] if equil_index < len(time_array) else time_array[-1]
    return equil_time, equil_index

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
    print("Loading RMSD data from all directories...")
    
    # Load data from each directory
    all_data = []
    equilibration_points = {}
    
    for directory in directories:
        data = load_rmsd_data(directory)
        if data is not None:
            # Detect equilibration point
            equil_time, equil_index = detect_equilibration(
                data['Time'].values, data['RMSD'].values
            )
            equilibration_points[directory] = (equil_time, equil_index)
            
            all_data.append(data)
    
    if not all_data:
        print("No valid RMSD data found. Please run rmsd_analysis.py in each directory first.")
        return
    
    # Combine all data into a single DataFrame
    combined_data = pd.concat(all_data, ignore_index=True)
    
    print("Generating RMSD comparison plots...")
    
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
    
    # Plot 1: RMSD vs Time for all systems
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Add vertical grid lines at major time intervals (every 1000 ps)
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.grid(True, which='major', axis='both', linestyle='--', alpha=0.7)
    
    # Plot raw data
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        ax.plot(system_data['Time'], system_data['RMSD'], 
                label=system, color=COLOR_MAP[system], alpha=0.5, linewidth=1)
    
    # Add smoothed trend lines
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        # Sort by time to ensure proper running average
        system_data = system_data.sort_values('Time')
        # Calculate running average
        smooth_time, smooth_rmsd = calculate_running_average(
            system_data['Time'].values, system_data['RMSD'].values, window=100)
        ax.plot(smooth_time, smooth_rmsd, color=COLOR_MAP[system], 
                linewidth=2.5, alpha=0.8, linestyle='-')
    
    # Add vertical lines for equilibration points
    for system in systems:
        if system in equilibration_points:
            equil_time, _ = equilibration_points[system]
            ax.axvline(x=equil_time, color=COLOR_MAP[system], 
                      linestyle='--', alpha=0.7, linewidth=1.5)
            
            # Add annotation for equilibration point
            ax.annotate(f"Equilibration: {equil_time:.0f} ps", 
                       xy=(equil_time, ax.get_ylim()[0]), 
                       xytext=(equil_time, ax.get_ylim()[0] * 1.1),
                       rotation=90, ha='right', va='bottom',
                       color=COLOR_MAP[system], fontsize=10,
                       arrowprops=dict(arrowstyle="->", color=COLOR_MAP[system]))
    
    # Add legend
    ax.legend(title="System", loc='upper right', framealpha=0.7)
    
    ax.set_xlabel('Time (ps)', fontweight='bold')
    ax.set_ylabel('RMSD (Å)', fontweight='bold')
    ax.set_title('Water Oxygen RMSD Comparison Over Time', fontsize=14, fontweight='bold')
    
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
        axins.plot(zoom_data['Time'], zoom_data['RMSD'], 
                  color=COLOR_MAP[system], linewidth=1.5)
    
    axins.set_title('Last 20% of Simulation', fontsize=10, fontweight='bold')
    axins.grid(True, linestyle='--', alpha=0.5)
    
    # Add annotation to highlight the inset with a more prominent arrow
    ax.annotate('Zoomed region', xy=(zoom_start, combined_data['RMSD'].min()),
               xytext=(zoom_start * 0.7, combined_data['RMSD'].min()),
               arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=10),
               fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(rmsd_dir, 'rmsd_comparison.png'), dpi=300)
    plt.close()
    
    # Plot 2: RMSD distribution (KDE plot) for all systems
    fig, ax = plt.subplots(figsize=(12, 8))
    
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        
        # Get equilibrated data (after equilibration point)
        if system in equilibration_points:
            _, equil_index = equilibration_points[system]
            equil_data = system_data.iloc[equil_index:]
        else:
            equil_data = system_data
        
        # Plot KDE for equilibrated data
        sns.kdeplot(equil_data['RMSD'], label=f"{system} (Equilibrated)", fill=True, 
                   alpha=0.3, color=COLOR_MAP[system], ax=ax)
        
        # Add vertical line for mean RMSD
        mean_rmsd = equil_data['RMSD'].mean()
        ax.axvline(x=mean_rmsd, color=COLOR_MAP[system], 
                  linestyle='-', alpha=0.7, linewidth=1.5)
        
        # Add annotation for mean RMSD
        ax.annotate(f"Mean: {mean_rmsd:.4f} Å", 
                   xy=(mean_rmsd, 0), 
                   xytext=(mean_rmsd, 0.1),
                   rotation=90, ha='right', va='bottom',
                   color=COLOR_MAP[system], fontsize=10)
    
    ax.set_xlabel('RMSD (Å)', fontweight='bold')
    ax.set_ylabel('Probability Density', fontweight='bold')
    ax.set_title('Water Oxygen RMSD Distribution Comparison (Equilibrated)', fontsize=14, fontweight='bold')
    ax.legend(title="System", loc='upper right', framealpha=0.7)
    ax.grid(True, linestyle='--', alpha=0.7)
    
    # Add text with statistical information
    stats_text = ""
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        
        # Get equilibrated data
        if system in equilibration_points:
            _, equil_index = equilibration_points[system]
            equil_data = system_data.iloc[equil_index:]
        else:
            equil_data = system_data
            
        mean = equil_data['RMSD'].mean()
        std = equil_data['RMSD'].std()
        stats_text += f"{system}: {mean:.4f} ± {std:.4f} Å\n"
    
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
           verticalalignment='top', horizontalalignment='left',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(os.path.join(rmsd_dir, 'rmsd_distribution_comparison.png'), dpi=300)
    plt.close()
    
    # Plot 3: Boxplot of RMSD values with statistical significance
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create a DataFrame with equilibrated data for boxplot
    boxplot_data = []
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        
        # Get equilibrated data
        if system in equilibration_points:
            _, equil_index = equilibration_points[system]
            equil_data = system_data.iloc[equil_index:]
        else:
            equil_data = system_data
            
        # Add to boxplot data
        boxplot_data.append(pd.DataFrame({
            'System': system,
            'RMSD': equil_data['RMSD'].values
        }))
    
    boxplot_df = pd.concat(boxplot_data, ignore_index=True)
    
    try:
        # Try with hue parameter (newer seaborn versions)
        sns.boxplot(x='System', y='RMSD', data=boxplot_df, ax=ax,
                   palette=COLOR_MAP, hue='System', legend=False, notch=True, showfliers=False)
    except:
        # Fallback for older seaborn versions
        sns.boxplot(x='System', y='RMSD', data=boxplot_df, ax=ax,
                   palette=COLOR_MAP, notch=True, showfliers=False)
    
    # Add individual data points with jitter and transparency
    # Sample points to avoid overcrowding
    for system in systems:
        system_data = boxplot_df[boxplot_df['System'] == system]
        # Sample every 10th point to avoid overcrowding
        sampled_data = system_data.iloc[::10]
        sns.stripplot(x='System', y='RMSD', data=sampled_data, ax=ax,
                     size=3, color='black', alpha=0.2, jitter=True)
    
    # Add significance bars
    # Function to add significance bars
    def add_significance_bars(ax, data, systems, y_position):
        """Add significance bars between systems based on t-test."""
        pairs = []
        for i in range(len(systems)):
            for j in range(i+1, len(systems)):
                pairs.append((systems[i], systems[j]))
        
        for pair in pairs:
            sys1_data = data[data['System'] == pair[0]]['RMSD'].values
            sys2_data = data[data['System'] == pair[1]]['RMSD'].values
            
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
            ax.text((x1 + x2) / 2, y_position + 0.005, sig_symbol, ha='center')
    
    # Add significance bars
    max_rmsd = boxplot_df['RMSD'].max()
    add_significance_bars(ax, boxplot_df, systems, max_rmsd * 1.05)
    
    ax.set_xlabel('System', fontweight='bold')
    ax.set_ylabel('RMSD (Å)', fontweight='bold')
    ax.set_title('RMSD Distribution by System (Equilibrated)', fontsize=14, fontweight='bold')
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
    plt.savefig(os.path.join(rmsd_dir, 'rmsd_boxplot.png'), dpi=300)
    plt.close()
    
    # Plot 4: RMSD fluctuations (standard deviation comparison)
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Calculate standard deviations for each system
    std_data = []
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        
        # Get equilibrated data
        if system in equilibration_points:
            _, equil_index = equilibration_points[system]
            equil_data = system_data.iloc[equil_index:]
        else:
            equil_data = system_data
            
        std = equil_data['RMSD'].std()
        mean = equil_data['RMSD'].mean()
        rel_std = (std / mean) * 100  # Relative std as percentage
        
        std_data.append({
            'System': system,
            'Standard Deviation (Å)': std,
            'Relative Std (%)': rel_std
        })
    
    std_df = pd.DataFrame(std_data)
    
    # Create a bar plot for standard deviations
    bar_positions = np.arange(len(systems))
    bar_width = 0.35
    
    # Plot absolute standard deviations
    bars1 = ax.bar(bar_positions - bar_width/2, std_df['Standard Deviation (Å)'], 
                  bar_width, label='Absolute Std (Å)', alpha=0.7)
    
    # Plot relative standard deviations
    ax2 = ax.twinx()
    bars2 = ax2.bar(bar_positions + bar_width/2, std_df['Relative Std (%)'], 
                   bar_width, label='Relative Std (%)', alpha=0.7, color='orange')
    
    # Add labels and title
    ax.set_xlabel('System', fontweight='bold')
    ax.set_ylabel('Standard Deviation (Å)', fontweight='bold')
    ax2.set_ylabel('Relative Standard Deviation (%)', fontweight='bold')
    ax.set_title('RMSD Fluctuations Comparison', fontsize=14, fontweight='bold')
    
    # Set x-tick labels
    ax.set_xticks(bar_positions)
    ax.set_xticklabels(systems)
    
    # Add a legend
    ax.legend(handles=[bars1, bars2], loc='upper left')
    
    # Add value labels on bars
    for i, bar in enumerate(bars1):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.001,
               f'{height:.4f} Å', ha='center', va='bottom', fontsize=9)
    
    for i, bar in enumerate(bars2):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{height:.2f}%', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(os.path.join(rmsd_dir, 'rmsd_fluctuations.png'), dpi=300)
    plt.close()
    
    # Calculate and save statistics
    # Create a DataFrame with equilibrated data for statistics
    stats_data = []
    for system in systems:
        system_data = combined_data[combined_data['System'] == system]
        
        # Get equilibrated data
        if system in equilibration_points:
            equil_time, equil_index = equilibration_points[system]
            equil_data = system_data.iloc[equil_index:]
        else:
            equil_time = "N/A"
            equil_data = system_data
            
        # Calculate statistics
        mean = equil_data['RMSD'].mean()
        std = equil_data['RMSD'].std()
        min_val = equil_data['RMSD'].min()
        max_val = equil_data['RMSD'].max()
        count = len(equil_data)
        
        stats_data.append({
            'System': system,
            'Equilibration Time (ps)': equil_time if equil_time != "N/A" else np.nan,
            'Mean RMSD (Å)': mean,
            'Std Dev (Å)': std,
            'Min RMSD (Å)': min_val,
            'Max RMSD (Å)': max_val,
            'Sample Count': count,
            'Relative Std (%)': (std / mean) * 100
        })
    
    rmsd_stats = pd.DataFrame(stats_data)
    
    # Save statistics to file
    stats_file = os.path.join(rmsd_dir, 'rmsd_statistics.txt')
    
    with open(stats_file, 'w') as f:
        f.write("RMSD Statistics Comparison\n")
        f.write("=========================\n\n")
        
        # Write statistical significance test results
        f.write("Statistical Significance (p-values from t-test):\n")
        f.write("-----------------------------------------------\n")
        for i, sys1 in enumerate(systems):
            for j, sys2 in enumerate(systems):
                if i < j:
                    sys1_data = boxplot_df[boxplot_df['System'] == sys1]['RMSD'].values
                    sys2_data = boxplot_df[boxplot_df['System'] == sys2]['RMSD'].values
                    t_stat, p_value = scipy_stats.ttest_ind(sys1_data, sys2_data, equal_var=False)
                    sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
                    f.write(f"{sys1} vs {sys2}: p={p_value:.6f} {sig}\n")
        f.write("\n")
        
        # Write detailed statistics for each system
        for _, row in rmsd_stats.iterrows():
            f.write(f"System: {row['System']}\n")
            if not np.isnan(row['Equilibration Time (ps)']):
                f.write(f"  Equilibration Time: {row['Equilibration Time (ps)']:.1f} ps\n")
            else:
                f.write(f"  Equilibration Time: N/A\n")
            f.write(f"  Mean RMSD: {row['Mean RMSD (Å)']:.6f} Å\n")
            f.write(f"  Std Deviation: {row['Std Dev (Å)']:.6f} Å\n")
            f.write(f"  Min RMSD: {row['Min RMSD (Å)']:.6f} Å\n")
            f.write(f"  Max RMSD: {row['Max RMSD (Å)']:.6f} Å\n")
            f.write(f"  Sample Count: {row['Sample Count']}\n")
            f.write(f"  Relative Fluctuation: {row['Relative Std (%)']:.2f}%\n")
            f.write("\n")
        
        # Add note about RMSD interpretation
        f.write("\nNote about RMSD Interpretation:\n")
        f.write("-----------------------------\n")
        f.write("RMSD measures structural deviation from a reference structure (first frame).\n")
        f.write("Lower RMSD values indicate structures closer to the reference.\n")
        f.write("Equilibration is detected when RMSD fluctuations stabilize.\n")
        f.write("For water systems, RMSD primarily reflects diffusion and reorganization of water molecules.\n")
        f.write("Differences between systems may indicate varying degrees of structural stability or mobility.\n")
    
    print(f"Saved RMSD comparison plots and statistics to {rmsd_dir}/")
    print(f"Statistics saved to {stats_file}")

if __name__ == "__main__":
    main() 