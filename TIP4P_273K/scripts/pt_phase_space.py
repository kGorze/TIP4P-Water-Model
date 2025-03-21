#!/usr/bin/env python3
"""
Pressure-Temperature Phase Space Visualization for TIP4P Water

This script creates visualizations of the pressure-temperature phase space
trajectory during molecular dynamics simulations.

Usage:
    python pt_phase_space.py

Output:
    - Phase space trajectory plot
    - 2D histogram/topology map of visited states
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import pandas as pd
from utils import parse_edr_file, save_plot

def create_pt_phase_space_plot(edr_file='../md.edr', output_dir='../analysis'):
    """
    Create a visualization of pressure-temperature phase space
    
    Parameters:
    -----------
    edr_file : str
        Path to the GROMACS energy file (.edr)
    output_dir : str
        Directory to save output files
    """
    print("Analyzing pressure-temperature phase space...")
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Extract temperature and pressure data from the energy file
    data = parse_edr_file(edr_file, properties=['Temperature', 'Pressure'])
    
    # Check if data extraction was successful
    if 'Temperature' not in data or 'Pressure' not in data:
        print("Error: Could not extract temperature or pressure data from EDR file")
        return None
    
    # Extract data
    temp_data = data['Temperature']
    press_data = data['Pressure']
    
    # Merge the dataframes on 'Time (ps)'
    merged_df = pd.merge(temp_data, press_data, on='Time (ps)', how='inner')
    
    # Convert pressure from bar to atm (1 bar = 0.986923 atm)
    conversion_factor = 0.986923
    merged_df['Pressure'] = merged_df['Pressure'] * conversion_factor
    
    # Calculate statistics
    temp_mean = merged_df['Temperature'].mean()
    temp_std = merged_df['Temperature'].std()
    press_mean = merged_df['Pressure'].mean()
    press_std = merged_df['Pressure'].std()
    
    # Save combined data
    merged_df.to_csv(f'{output_dir}/pt_phase_space_data.dat', sep='\t', index=False, 
                    float_format='%.6f', header=True)
    
    # Save statistics
    with open(f'{output_dir}/pt_phase_space_stats.dat', 'w') as f:
        f.write(f"# Pressure-Temperature Phase Space Statistics\n")
        f.write(f"# Temperature mean: {temp_mean:.4f} K\n")
        f.write(f"# Temperature std: {temp_std:.4f} K\n")
        f.write(f"# Pressure mean: {press_mean:.4f} atm\n")
        f.write(f"# Pressure std: {press_std:.4f} atm\n")
        f.write(f"# Total points: {len(merged_df)}\n")
    
    # Create figure with subplots
    fig = plt.figure(figsize=(15, 12))
    gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 3])
    
    # Main scatter plot (bottom left)
    ax_main = plt.subplot(gs[1, 0])
    scatter = ax_main.scatter(
        merged_df['Temperature'], 
        merged_df['Pressure'],
        c=merged_df['Time (ps)'],
        cmap='viridis',
        alpha=0.7,
        s=5,
        label='PT Trajectory'
    )
    
    # Add a colorbar for time
    cb = plt.colorbar(scatter, ax=ax_main)
    cb.set_label('Time (ps)')
    
    # Add markers for the start and end points
    ax_main.scatter(
        merged_df['Temperature'].iloc[0], 
        merged_df['Pressure'].iloc[0],
        c='red', 
        s=100, 
        marker='o', 
        label='Start',
        edgecolors='black'
    )
    ax_main.scatter(
        merged_df['Temperature'].iloc[-1], 
        merged_df['Pressure'].iloc[-1],
        c='blue', 
        s=100, 
        marker='o', 
        label='End',
        edgecolors='black'
    )
    
    # Add mean point
    ax_main.scatter(
        temp_mean, 
        press_mean,
        c='yellow', 
        s=100, 
        marker='*', 
        label='Mean',
        edgecolors='black'
    )
    
    # Add target values (convert 1 bar to atm)
    target_pressure_atm = 1.0 * conversion_factor
    ax_main.axhline(y=target_pressure_atm, color='k', linestyle='--', alpha=0.5, label=f'Target Pressure ({target_pressure_atm:.2f} atm)')
    ax_main.axvline(x=273.15, color='k', linestyle='--', alpha=0.5, label='Target Temp (273.15 K)')
    
    # Set labels and title
    ax_main.set_xlabel('Temperature (K)')
    ax_main.set_ylabel('Pressure (atm)')
    ax_main.set_title('Pressure-Temperature Phase Space Trajectory')
    ax_main.legend(loc='upper left')
    ax_main.grid(True, alpha=0.3)
    
    # Set axis limits to zoom on the relevant region
    temp_padding = max(3.0, temp_std * 3)
    press_padding = max(250.0 * conversion_factor, press_std * 1.5)
    ax_main.set_xlim(temp_mean - temp_padding, temp_mean + temp_padding)
    ax_main.set_ylim(press_mean - press_padding, press_mean + press_padding)
    
    # Temperature histogram (top left)
    ax_temp = plt.subplot(gs[0, 0], sharex=ax_main)
    ax_temp.hist(merged_df['Temperature'], bins=50, color='blue', alpha=0.7)
    ax_temp.axvline(x=temp_mean, color='r', linestyle='--', label=f'Mean: {temp_mean:.2f} K')
    ax_temp.axvline(x=273.15, color='k', linestyle='--', alpha=0.5, label='Target: 273.15 K')
    ax_temp.set_ylabel('Frequency')
    ax_temp.set_title('Temperature Distribution')
    ax_temp.grid(True, alpha=0.3)
    plt.setp(ax_temp.get_xticklabels(), visible=False)
    
    # Pressure histogram (bottom right)
    ax_press = plt.subplot(gs[1, 1], sharey=ax_main)
    ax_press.hist(merged_df['Pressure'], bins=50, color='green', alpha=0.7, orientation='horizontal')
    ax_press.axhline(y=press_mean, color='r', linestyle='--', label=f'Mean: {press_mean:.2f} atm')
    ax_press.axhline(y=target_pressure_atm, color='k', linestyle='--', alpha=0.5, label=f'Target: {target_pressure_atm:.2f} atm')
    ax_press.set_xlabel('Frequency')
    ax_press.set_title('Pressure Distribution')
    ax_press.grid(True, alpha=0.3)
    plt.setp(ax_press.get_yticklabels(), visible=False)
    
    # Empty axes for top right
    plt.subplot(gs[0, 1]).set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/pt_phase_space_trajectory.png', dpi=300, bbox_inches='tight')
    print(f"Saved phase space trajectory plot to {output_dir}/pt_phase_space_trajectory.png")
    
    # Create a 2D histogram (topology map) of the phase space
    plt.figure(figsize=(12, 10))
    
    # Calculate histogram
    h, xedges, yedges, im = plt.hist2d(
        merged_df['Temperature'], 
        merged_df['Pressure'],
        bins=50,
        cmap='viridis',
        norm=LogNorm(),
        alpha=0.8
    )
    
    # Add colorbar
    cbar = plt.colorbar(im)
    cbar.set_label('Count (log scale)')
    
    # Add markers
    plt.scatter(
        temp_mean, 
        press_mean,
        c='yellow', 
        s=100, 
        marker='*', 
        label='Mean',
        edgecolors='black'
    )
    plt.scatter(
        merged_df['Temperature'].iloc[0], 
        merged_df['Pressure'].iloc[0],
        c='red', 
        s=80, 
        marker='o', 
        label='Start',
        edgecolors='black'
    )
    plt.scatter(
        merged_df['Temperature'].iloc[-1], 
        merged_df['Pressure'].iloc[-1],
        c='blue', 
        s=80, 
        marker='o', 
        label='End',
        edgecolors='black'
    )
    
    # Add target values
    plt.axhline(y=target_pressure_atm, color='w', linestyle='--', alpha=0.7, label=f'Target Pressure ({target_pressure_atm:.2f} atm)')
    plt.axvline(x=273.15, color='w', linestyle='--', alpha=0.7, label='Target Temp (273.15 K)')
    
    # Set labels and title
    plt.xlabel('Temperature (K)')
    plt.ylabel('Pressure (atm)')
    plt.title('Pressure-Temperature Phase Space Density Map')
    plt.legend(loc='upper right')
    plt.grid(False)
    
    # Set axis limits to zoom on the relevant region
    plt.xlim(temp_mean - temp_padding, temp_mean + temp_padding)
    plt.ylim(press_mean - press_padding, press_mean + press_padding)
    
    # Add a text box with statistics
    stats_text = '\n'.join([
        f'Mean Temperature: {temp_mean:.2f} K (σ={temp_std:.2f} K)',
        f'Mean Pressure: {press_mean:.2f} atm (σ={press_std:.2f} atm)',
        f'Target: 273.15 K, {target_pressure_atm:.2f} atm',
        f'Most visited region: {xedges[np.unravel_index(np.argmax(h), h.shape)[0]]:.2f} K, ' +
        f'{yedges[np.unravel_index(np.argmax(h), h.shape)[1]]:.2f} atm'
    ])
    plt.text(0.02, 0.02, stats_text, transform=plt.gca().transAxes, fontsize=10,
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/pt_phase_space_density.png', dpi=300, bbox_inches='tight')
    print(f"Saved phase space density map to {output_dir}/pt_phase_space_density.png")
    
    # Create a 3D line plot showing the PT trajectory over time
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the trajectory
    ax.plot(
        merged_df['Time (ps)'], 
        merged_df['Temperature'], 
        merged_df['Pressure'], 
        color='blue', 
        alpha=0.7,
        linewidth=1
    )
    
    # Mark the start and end points
    ax.scatter(
        merged_df['Time (ps)'].iloc[0], 
        merged_df['Temperature'].iloc[0], 
        merged_df['Pressure'].iloc[0],
        color='red',
        s=100,
        label='Start'
    )
    ax.scatter(
        merged_df['Time (ps)'].iloc[-1], 
        merged_df['Temperature'].iloc[-1], 
        merged_df['Pressure'].iloc[-1],
        color='blue',
        s=100,
        label='End'
    )
    
    # Add reference planes for target values
    time_max = merged_df['Time (ps)'].max()
    xgrid, zgrid = np.meshgrid(
        np.linspace(0, time_max, 10),
        np.linspace(press_mean - press_padding, press_mean + press_padding, 10)
    )
    ygrid = np.ones_like(xgrid) * 273.15
    ax.plot_surface(xgrid, ygrid, zgrid, color='gray', alpha=0.2)
    
    xgrid, ygrid = np.meshgrid(
        np.linspace(0, time_max, 10),
        np.linspace(temp_mean - temp_padding, temp_mean + temp_padding, 10)
    )
    zgrid = np.ones_like(xgrid) * target_pressure_atm
    ax.plot_surface(xgrid, ygrid, zgrid, color='gray', alpha=0.2)
    
    # Set labels and title
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Temperature (K)')
    ax.set_zlabel('Pressure (atm)')
    ax.set_title('Pressure-Temperature Trajectory Over Time')
    ax.legend()
    
    # Set view angle
    ax.view_init(elev=30, azim=45)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/pt_phase_space_3d.png', dpi=300, bbox_inches='tight')
    print(f"Saved 3D phase space trajectory plot to {output_dir}/pt_phase_space_3d.png")
    
    return merged_df

def main():
    """Main function to run the analysis"""
    # Get the absolute path of the script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Set up paths
    edr_file = os.path.join(script_dir, '..', 'md.edr')
    output_dir = os.path.join(script_dir, '..', 'analysis')
    
    # Run the analysis
    create_pt_phase_space_plot(edr_file, output_dir)
    
    print("Pressure-Temperature phase space analysis complete!")

if __name__ == "__main__":
    main() 