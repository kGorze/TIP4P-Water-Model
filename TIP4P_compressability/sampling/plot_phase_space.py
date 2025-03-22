#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import sys
import pandas as pd
from scipy.signal import savgol_filter

def run_gmx_energy(edr_file, tpr_file):
    """
    Extract pressure and temperature data from GROMACS energy file
    """
    if not os.path.exists(edr_file):
        raise FileNotFoundError(f"EDR file {edr_file} not found!")
    if not os.path.exists(tpr_file):
        raise FileNotFoundError(f"TPR file {tpr_file} not found!")

    # Create input for gmx energy
    with open('temp_input.txt', 'w') as f:
        f.write("Temperature\nPressure\n0\n")  # Changed order to Temperature first

    # Run gmx energy command
    cmd = f"gmx energy -f {edr_file} -s {tpr_file} -o energy.xvg < temp_input.txt"
    try:
        process = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running gmx energy: {e}")
        print(f"Error output: {e.stderr}")
        sys.exit(1)
    finally:
        if os.path.exists('temp_input.txt'):
            os.remove('temp_input.txt')

    # Read the generated xvg file
    if not os.path.exists('energy.xvg'):
        raise FileNotFoundError("energy.xvg was not created!")
    
    try:
        # Skip comment lines (those starting with # or @)
        with open('energy.xvg', 'r') as f:
            lines = f.readlines()
        data_lines = [l for l in lines if not (l.startswith('#') or l.startswith('@'))]
        
        # Convert to DataFrame with correct column order
        data = pd.DataFrame([line.split() for line in data_lines], 
                          columns=['Time', 'Temperature', 'Pressure'],  # Changed order
                          dtype=float)
        
        # Convert pressure from bar to atm and time to ns
        data['Pressure'] = data['Pressure'] * 0.986923  # bar to atm
        data['Time'] = data['Time'] / 1000  # ps to ns
        
        # Print some basic statistics to verify the data
        print("\nRaw Data Statistics:")
        print(f"Temperature: mean = {data['Temperature'].mean():.2f} K, std = {data['Temperature'].std():.2f} K")
        print(f"Min: {data['Temperature'].min():.2f} K, Max: {data['Temperature'].max():.2f} K")
        print(f"Pressure: mean = {data['Pressure'].mean():.2f} atm, std = {data['Pressure'].std():.2f} atm")
        print(f"Min: {data['Pressure'].min():.2f} atm, Max: {data['Pressure'].max():.2f} atm")
        
        return data
    except Exception as e:
        print(f"Error processing energy.xvg: {e}")
        sys.exit(1)
    finally:
        if os.path.exists('energy.xvg'):
            os.remove('energy.xvg')

def preprocess_data(data, downsample_factor=5, filter_outliers=True, smooth=True):
    """Preprocess the data for clearer visualization"""
    # Create a copy to avoid modifying original data
    processed_data = data.copy()
    
    # Filter outliers based on pressure
    if filter_outliers:
        p_q1, p_q3 = np.percentile(processed_data['Pressure'], [1, 99])
        p_iqr = p_q3 - p_q1
        p_lower = p_q1 - 1.5 * p_iqr
        p_upper = p_q3 + 1.5 * p_iqr
        
        # Filter out extreme outliers
        mask = (processed_data['Pressure'] >= p_lower) & (processed_data['Pressure'] <= p_upper)
        filtered_data = processed_data[mask]
        
        # Only use filtered data if we haven't lost too many points
        if len(filtered_data) > 0.9 * len(processed_data):
            processed_data = filtered_data
            print(f"Filtered outliers: removed {len(data) - len(processed_data)} points")
    
    # Downsample to reduce clutter
    processed_data = processed_data.iloc[::downsample_factor].reset_index(drop=True)
    print(f"Downsampled data from {len(data)} to {len(processed_data)} points")
    
    # Apply Savitzky-Golay filter for smoothing if needed
    if smooth and len(processed_data) > 10:
        # Make sure we have enough data points for the filter
        window_size = min(11, len(processed_data) // 2 * 2 - 1)  # must be odd
        if window_size >= 3:
            try:
                processed_data['Pressure_smooth'] = savgol_filter(processed_data['Pressure'], window_size, 3)
                processed_data['Temperature_smooth'] = savgol_filter(processed_data['Temperature'], window_size, 3)
                print(f"Applied smoothing with window size {window_size}")
            except Exception as e:
                print(f"Could not apply smoothing: {e}")
                processed_data['Pressure_smooth'] = processed_data['Pressure']
                processed_data['Temperature_smooth'] = processed_data['Temperature']
        else:
            processed_data['Pressure_smooth'] = processed_data['Pressure']
            processed_data['Temperature_smooth'] = processed_data['Temperature']
    else:
        processed_data['Pressure_smooth'] = processed_data['Pressure']
        processed_data['Temperature_smooth'] = processed_data['Temperature']
    
    return processed_data

def create_phase_space_plot(data, output_suffix="", connect_points=False, downsample_factor=5, smooth=True):
    """Create phase space trajectory plot with time-based coloring"""
    plt.style.use('seaborn-v0_8-darkgrid')
    fig = plt.figure(figsize=(15, 10))
    
    # Preprocess data for visualization
    if connect_points:
        plot_data = preprocess_data(data, downsample_factor=downsample_factor, smooth=smooth)
    else:
        plot_data = data
    
    # Calculate temperature range centered around mean with reasonable std
    temp_mean = plot_data['Temperature'].mean()
    temp_std = plot_data['Temperature'].std()
    t_min = temp_mean - 3 * temp_std
    t_max = temp_mean + 3 * temp_std
    
    # Calculate pressure range using percentiles to avoid outliers
    p_min, p_max = plot_data['Pressure'].quantile([0.01, 0.99])
    p_range = p_max - p_min
    p_min -= 0.05 * p_range
    p_max += 0.05 * p_range
    
    # Create grid for subplots
    gs = plt.GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 3])
    
    # Temperature distribution (top)
    ax_temp_hist = fig.add_subplot(gs[0, 0])
    ax_temp_hist.hist(plot_data['Temperature'], bins=50, color='royalblue', alpha=0.7)
    ax_temp_hist.set_xlabel('Temperature (K)')
    ax_temp_hist.set_ylabel('Frequency')
    ax_temp_hist.set_xlim(t_min, t_max)
    
    # Main phase space plot
    ax_main = fig.add_subplot(gs[1, 0])
    
    if connect_points:
        # For connected trajectory, use the smoothed data
        t_col = 'Temperature_smooth'
        p_col = 'Pressure_smooth'
        
        # Set up colormap for time progression
        cmap = plt.cm.viridis
        norm = plt.Normalize(plot_data['Time'].min(), plot_data['Time'].max())
        
        # Plot segments with color gradient for time
        points = np.array([plot_data[t_col], plot_data[p_col]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        # Create a LineCollection with varying linewidth
        from matplotlib.collections import LineCollection
        lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=1.5, alpha=0.7)
        lc.set_array(plot_data['Time'][:-1])
        line = ax_main.add_collection(lc)
        
        # Add scatter at regular intervals for time indication
        interval = max(1, len(plot_data) // 50)  # At most 50 markers
        markers = plot_data.iloc[::interval]
        scatter = ax_main.scatter(markers[t_col], markers[p_col],
                                c=markers['Time'], cmap=cmap, norm=norm,
                                s=30, alpha=0.8, marker='o', edgecolors='k', linewidths=0.5)
        
        # Add special markers for start and end
        ax_main.scatter(plot_data[t_col].iloc[0], plot_data[p_col].iloc[0], 
                       c='blue', s=100, marker='*', label='Start', edgecolors='k', zorder=10)
        ax_main.scatter(plot_data[t_col].iloc[-1], plot_data[p_col].iloc[-1], 
                       c='red', s=100, marker='*', label='End', edgecolors='k', zorder=10)
        
        # Add arrows to show direction
        arrow_idx = np.linspace(0, len(plot_data)-2, 10, dtype=int)
        for i in arrow_idx:
            ax_main.annotate('', 
                            xy=(plot_data[t_col].iloc[i+1], plot_data[p_col].iloc[i+1]),
                            xytext=(plot_data[t_col].iloc[i], plot_data[p_col].iloc[i]),
                            arrowprops=dict(arrowstyle='->', color=cmap(norm(plot_data['Time'].iloc[i])), 
                                            lw=1.5, alpha=0.8),
                            zorder=5)
    else:
        # Original scatter plot with all data points
        scatter = ax_main.scatter(data['Temperature'], data['Pressure'],
                                c=data['Time'], cmap='viridis',
                                s=2, alpha=0.6)
    
    ax_main.set_xlabel('Temperature (K)')
    ax_main.set_ylabel('Pressure (atm)')
    ax_main.set_xlim(t_min, t_max)
    ax_main.set_ylim(p_min, p_max)
    
    # Add target lines
    ax_main.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, label='Target P (1.0 atm)')
    ax_main.axvline(x=273.15, color='red', linestyle='--', alpha=0.5, label='Target T (273.15 K)')
    
    # Pressure distribution (right)
    ax_press_hist = fig.add_subplot(gs[1, 1])
    ax_press_hist.hist(plot_data['Pressure'], bins=50, orientation='horizontal',
                      color='forestgreen', alpha=0.7)
    ax_press_hist.set_ylabel('Pressure (atm)')
    ax_press_hist.set_xlabel('Frequency')
    ax_press_hist.set_ylim(p_min, p_max)
    
    # Print some statistics
    print(f"\nPlotting ranges:")
    print(f"Temperature range: {t_min:.2f} K to {t_max:.2f} K")
    print(f"Pressure range: {p_min:.2f} atm to {p_max:.2f} atm")
    
    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax_main)
    cbar.set_label('Time (ns)')
    
    # Add legend and title
    ax_main.legend(loc='upper left')
    
    title = 'Phase Space Trajectory and Distributions'
    if connect_points:
        title += f' (Connected Points, n={len(plot_data)})'
    else:
        title += f' (Scatter, n={len(data)})'
    plt.suptitle(title, fontsize=14)
    
    # Add info about preprocessing
    if connect_points:
        info_text = f"Downsampled by factor {downsample_factor}"
        if smooth:
            info_text += ", with smoothing applied"
        ax_main.text(0.02, 0.02, info_text, transform=ax_main.transAxes, 
                    fontsize=8, alpha=0.7, bbox=dict(facecolor='white', alpha=0.5))
    
    # Adjust layout and save
    plt.tight_layout()
    output_file = f'phase_space_trajectory{output_suffix}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved as: {output_file}")
    plt.close()

def main():
    # Set paths
    edr_file = "../md.edr"  # Path to your EDR file
    tpr_file = "../md.tpr"  # Path to your TPR file
    
    try:
        # Read the energy data
        print("Reading energy data from GROMACS files...")
        data = run_gmx_energy(edr_file, tpr_file)
        
        # Create both versions of the plot
        print("\nCreating scatter plot version...")
        create_phase_space_plot(data, "_scatter", connect_points=False)
        
        print("\nCreating connected trajectory version...")
        # Test different downsampling factors based on data size
        if len(data) > 5000:
            downsample_factor = 20
        elif len(data) > 1000:
            downsample_factor = 10
        else:
            downsample_factor = 5
            
        create_phase_space_plot(data, "_trajectory", connect_points=True, 
                              downsample_factor=downsample_factor, smooth=True)
        
        # Create a more detailed version with less downsampling
        print("\nCreating detailed trajectory version...")
        create_phase_space_plot(data, "_trajectory_detailed", connect_points=True, 
                              downsample_factor=max(1, downsample_factor//2), smooth=True)
        
        print("Done!")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 