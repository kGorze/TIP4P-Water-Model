#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import stats as scipy_stats
from matplotlib.gridspec import GridSpec

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
simulation_time = 20  # ns (longer to clearly show plateau)
time_step = 0.002  # ps
steps_per_ns = int(1 / time_step)
time_points = np.linspace(0, simulation_time, simulation_time * steps_per_ns)

def generate_synthetic_data():
    """Generate synthetic data for RMSD and energy drift illustration."""
    # Create time array in ns
    time_ns = time_points
    
    # Base energy
    base_energy = -12500
    
    # Create RMSD data - faster equilibration for 298K, slower for 273K
    rmsd_273K = 0.55 * (1 - np.exp(-time_points/2.5)) + 0.1
    rmsd_273K += np.random.normal(0, 0.015, size=len(time_points))
    
    rmsd_298K = 0.58 * (1 - np.exp(-time_points/1.5)) + 0.1
    rmsd_298K += np.random.normal(0, 0.015, size=len(time_points))
    
    # Create energy drift data - consistent with TIP4P water model rates
    energy_273K = base_energy - 21 * time_points + np.random.normal(0, 3, size=len(time_points))
    energy_298K = base_energy + 25 * time_points + np.random.normal(0, 3, size=len(time_points))
    
    return {
        'time_ns': time_ns,
        'rmsd_273K': rmsd_273K,
        'rmsd_298K': rmsd_298K,
        'energy_273K': energy_273K,
        'energy_298K': energy_298K
    }

def plot_simplified_rmsd_energy(data, output_file):
    """Create a simplified visualization showing RMSD plateaus despite energy drift."""
    # Create figure with just two subplots
    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 1])
    
    # Define colors
    colors = {
        '273K': '#1f77b4',  # Blue
        '298K': '#ff7f0e'   # Orange
    }
    
    # Top plot: RMSD over time with equilibration phase highlighted
    ax1 = fig.add_subplot(gs[0])
    
    # Create shaded region for equilibration phase
    equilibration_end = 5  # ns
    ax1.axvspan(0, equilibration_end, alpha=0.15, color='gray', label='Equilibration Phase')
    ax1.axvspan(equilibration_end, simulation_time, alpha=0.05, color='green', label='Production Phase')
    
    # Add vertical line at equilibration end
    ax1.axvline(x=equilibration_end, color='k', linestyle='--', alpha=0.5)
    
    # Plot RMSD data
    ax1.plot(data['time_ns'], data['rmsd_273K'], color=colors['273K'], 
             linewidth=2.5, label='TIP4P 273K')
    ax1.plot(data['time_ns'], data['rmsd_298K'], color=colors['298K'], 
             linewidth=2.5, label='TIP4P 298K')
    
    # Add plateau indicators - horizontal lines at stabilized values
    plateau_time = 10  # ns
    plateau_idx = np.argmin(np.abs(data['time_ns'] - plateau_time))
    
    plateau_273K = data['rmsd_273K'][plateau_idx:].mean()
    plateau_298K = data['rmsd_298K'][plateau_idx:].mean()
    
    ax1.axhline(y=plateau_273K, color=colors['273K'], linestyle='--', alpha=0.7)
    ax1.axhline(y=plateau_298K, color=colors['298K'], linestyle='--', alpha=0.7)
    
    # Add annotation about higher temperature causing faster equilibration
    ax1.annotate('Faster equilibration at 298K\ndue to higher molecular mobility', 
                xy=(2.5, data['rmsd_298K'][int(2.5*steps_per_ns)]),
                xytext=(4, data['rmsd_298K'][int(2.5*steps_per_ns)] + 0.15),
                arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
                fontsize=12)
    
    # Add annotation about plateaus
    ax1.text(equilibration_end + 0.5, data['rmsd_298K'][plateau_idx] + 0.03, 
            'RMSD stabilizes at plateau\ndespite continued energy drift',
            fontsize=12, ha='left', va='center',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.7))
    
    ax1.set_xlabel('Simulation Time (ns)', fontweight='bold')
    ax1.set_ylabel('RMSD (nm)', fontweight='bold')
    ax1.set_title('Root Mean Square Deviation (RMSD) Stabilization', fontsize=16, fontweight='bold')
    ax1.legend(loc='lower right', frameon=True)
    
    # Bottom plot: Energy drift - making it simple but clear
    ax2 = fig.add_subplot(gs[1])
    
    # Same shaded regions for consistency
    ax2.axvspan(0, equilibration_end, alpha=0.15, color='gray')
    ax2.axvspan(equilibration_end, simulation_time, alpha=0.05, color='green')
    ax2.axvline(x=equilibration_end, color='k', linestyle='--', alpha=0.5)
    
    # Plot energy data
    ax2.plot(data['time_ns'], data['energy_273K'], color=colors['273K'], 
             linewidth=2.5, label='TIP4P 273K')
    ax2.plot(data['time_ns'], data['energy_298K'], color=colors['298K'], 
             linewidth=2.5, label='TIP4P 298K')
    
    # Calculate and annotate drift rates
    slope_273K, _, _, _, _ = scipy_stats.linregress(data['time_ns'], data['energy_273K'])
    slope_298K, _, _, _, _ = scipy_stats.linregress(data['time_ns'], data['energy_298K'])
    
    # Add a single clear annotation about energy drift
    ax2.text(0.97, 0.05, 
            f"Energy drift rates:\n"
            f"273K: {abs(slope_273K):.1f} kJ/mol/ns\n"
            f"298K: {abs(slope_298K):.1f} kJ/mol/ns\n\n"
            f"Despite continuous energy drift,\n"
            f"structural integrity is maintained\n"
            f"as shown by stable RMSD plateaus.",
            transform=ax2.transAxes, 
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9), fontsize=12)
    
    ax2.set_xlabel('Simulation Time (ns)', fontweight='bold')
    ax2.set_ylabel('Total Energy (kJ/mol)', fontweight='bold')
    ax2.set_title('Continuous Energy Drift in MD Simulations', fontsize=16, fontweight='bold')
    ax2.legend(loc='upper right')
    
    # Finalize layout
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.3)
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()

def main():
    """Main function to run the visualizations."""
    try:
        print("Starting simplified RMSD stability visualization...")
        
        # Generate synthetic data
        data = generate_synthetic_data()
        
        # Create visualization
        output_path = os.path.join(output_dir, "rmsd_stability_simple.png")
        plot_simplified_rmsd_energy(data, output_path)
        
        print(f"Simplified RMSD stability visualization created successfully!")
        print(f"Output file saved to {output_path}")
    except Exception as e:
        print(f"Error occurred: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 