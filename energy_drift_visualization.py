#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from scipy import stats as scipy_stats
from matplotlib.ticker import MultipleLocator
from matplotlib.gridspec import GridSpec
import matplotlib.transforms as mtransforms

# Ensure output directory exists
output_dir = "visualization_outputs"
os.makedirs(output_dir, exist_ok=True)

# Set better visual style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

# Define simulation parameters for illustration
np.random.seed(42)  # For reproducibility
simulation_time = 10  # ns
time_step = 0.002  # ps
steps_per_ns = int(1 / time_step)
time_points = np.linspace(0, simulation_time, simulation_time * steps_per_ns)

def generate_synthetic_energy_data():
    """
    Generate synthetic energy data to illustrate energy drift concepts.
    Returns a dictionary with different simulation energy profiles.
    """
    # Create time array in ps
    time_ps = time_points * 1000  # Convert ns to ps
    
    # Ideal theoretical data (no drift, perfect conservation)
    ideal_energy = np.ones_like(time_points) * -12500
    
    # Typical drift rates mentioned (21-25 kJ/mol/ns)
    typical_drift_rate_low = 21  # kJ/mol/ns
    typical_drift_rate_high = 25  # kJ/mol/ns
    
    # Generate energy profiles with different drift rates
    # Positive drift
    positive_drift = ideal_energy + typical_drift_rate_high * time_points + np.random.normal(0, 3, size=len(time_points))
    
    # Negative drift
    negative_drift = ideal_energy - typical_drift_rate_low * time_points + np.random.normal(0, 3, size=len(time_points))
    
    # Small drift (good conservation) - removed to simplify visualization
    # small_drift = ideal_energy + 5 * time_points + np.random.normal(0, 3, size=len(time_points))
    
    # Large drift (poor conservation) - removed to simplify visualization
    # large_drift = ideal_energy + 50 * time_points + np.random.normal(0, 8, size=len(time_points))
    
    return {
        'time_ps': time_ps,
        'ideal': ideal_energy,
        'positive_drift': positive_drift,
        'negative_drift': negative_drift
    }

def plot_simplified_energy_drift(data, output_file):
    """
    Create a simplified, cleaner visualization of energy conservation and drift.
    Focus on key concepts rather than showing too many different scenarios.
    """
    fig = plt.figure(figsize=(15, 10))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[3, 2])
    
    # Main plot: Absolute energy values (Top panel)
    ax1 = fig.add_subplot(gs[0, :])
    
    # Plot energy trajectories - only show key examples to reduce clutter
    ax1.plot(data['time_ps']/1000, data['ideal'], 'k--', alpha=0.8, linewidth=2, label='Ideal (No Drift)')
    ax1.plot(data['time_ps']/1000, data['positive_drift'], 'r-', alpha=0.8, linewidth=2.5, label='Positive Drift (~25 kJ/mol/ns)')
    ax1.plot(data['time_ps']/1000, data['negative_drift'], 'b-', alpha=0.8, linewidth=2.5, label='Negative Drift (~21 kJ/mol/ns)')
    
    ax1.set_xlabel('Simulation Time (ns)', fontweight='bold')
    ax1.set_ylabel('Total Energy (kJ/mol)', fontweight='bold')
    ax1.set_title('Energy Conservation in Molecular Dynamics Simulations', fontsize=18, fontweight='bold')
    ax1.legend(loc='best', frameon=True, fancybox=True, shadow=True, fontsize=12)
    
    # Add annotation about numerical integration
    ax1.text(0.02, 0.04, 
           "Total energy should be conserved in isolated systems, but numerical integration methods\n"
           "introduce small errors in each time step that accumulate over the simulation.",
           transform=ax1.transAxes, 
           verticalalignment='bottom', horizontalalignment='left',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.9), fontsize=12)
    
    # Left bottom plot: Illustration of error accumulation
    ax3 = fig.add_subplot(gs[1, 0])
    
    # Time points for step-wise illustration
    step_times = np.linspace(0, 0.1, 30)  # Fewer points for clearer visualization
    
    # Generate step-wise error accumulation
    step_energy = -12500 + np.zeros_like(step_times)
    step_errors = np.random.normal(0.05, 0.02, size=len(step_times))  # Small random errors
    
    for i in range(1, len(step_times)):
        step_energy[i] = step_energy[i-1] + step_errors[i]
    
    # Plot step-wise error accumulation
    ax3.plot(step_times, step_energy, 'k-', alpha=0.8, linewidth=2)
    ax3.scatter(step_times, step_energy, c='r', s=40, alpha=0.7)
    
    # Add arrows to highlight errors - only add a few strategic ones
    key_indices = [5, 10, 15, 20, 25]
    for i in key_indices:
        ax3.annotate('', 
                     xy=(step_times[i], step_energy[i]), 
                     xytext=(step_times[i], step_energy[i-1]),
                     arrowprops=dict(arrowstyle='->', color='blue', alpha=0.8, linewidth=1.5))
    
    ax3.set_xlabel('Time (ns)', fontweight='bold')
    ax3.set_ylabel('Total Energy (kJ/mol)', fontweight='bold')
    ax3.set_title('Error Accumulation in Integration Steps', fontsize=16, fontweight='bold')
    
    # Add text explaining error accumulation
    ax3.text(0.05, 0.05, 
            "Each integration step introduces\n"
            "a small numerical error.\n"
            "These errors accumulate over time,\n"
            "resulting in energy drift.",
            transform=ax3.transAxes, 
            verticalalignment='bottom', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9), fontsize=12)
    
    # Right bottom plot: Acceptable drift ranges
    ax4 = fig.add_subplot(gs[1, 1])
    
    # Create bar chart of different drift rates
    drift_categories = ['Excellent', 'Good', 'Poor']  # Simplified categories
    drift_values = [5, 25, 100]  # Simplified values
    drift_colors = ['#2ca02c', '#1f77b4', '#d62728']  # Green, Blue, Red
    
    # Highlight TIP4P water model range
    tip4p_min = 21
    tip4p_max = 25
    
    # Create bar chart
    bars = ax4.bar(drift_categories, drift_values, color=drift_colors, alpha=0.7, width=0.6)
    
    # Add a horizontal band for TIP4P range
    ax4.axhline(y=tip4p_min, color='#ff7f0e', linestyle='--', alpha=0.8, linewidth=2)
    ax4.axhline(y=tip4p_max, color='#ff7f0e', linestyle='--', alpha=0.8, linewidth=2)
    ax4.axvspan(drift_categories.index('Good') - 0.4, drift_categories.index('Good') + 0.4, 
               alpha=0.2, color='#ff7f0e', label='TIP4P Range (21-25 kJ/mol/ns)')
    
    ax4.set_xlabel('Energy Conservation Quality', fontweight='bold')
    ax4.set_ylabel('Drift Rate (kJ/mol/ns)', fontweight='bold')
    ax4.set_title('Acceptable Ranges for Energy Drift', fontsize=16, fontweight='bold')
    
    # Add annotation about TIP4P
    ax4.text(0.5, 0.95, 
            "TIP4P water models typically show\n"
            "drift rates of 21-25 kJ/mol/ns,\n"
            "which falls within the good range\n"
            "for MD simulations.",
            transform=ax4.transAxes, 
            verticalalignment='top', horizontalalignment='center',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9), fontsize=12)
    
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.3, wspace=0.25)
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()

def main():
    """Main function to run the visualizations."""
    # Generate synthetic data
    data = generate_synthetic_energy_data()
    
    # Create simplified visualization
    plot_simplified_energy_drift(data, os.path.join(output_dir, "energy_drift_simplified.png"))
    
    print("Simplified energy drift visualization created successfully!")
    print(f"Output file saved to {output_dir}/energy_drift_simplified.png")

if __name__ == "__main__":
    main() 