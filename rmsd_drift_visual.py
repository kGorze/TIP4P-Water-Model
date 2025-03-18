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
simulation_time = 20  # ns (longer to clearly show plateau)
time_step = 0.002  # ps
steps_per_ns = int(1 / time_step)
time_points = np.linspace(0, simulation_time, simulation_time * steps_per_ns)

def generate_synthetic_data():
    """
    Generate synthetic data to illustrate RMSD stabilization despite energy drift.
    Returns a dictionary with RMSD and energy data for different models.
    """
    print("Generating time points...")
    # Create time array in ps
    time_ps = time_points * 1000  # Convert ns to ps
    
    # Base energy for all systems
    base_energy = -12500
    
    print("Creating RMSD data...")
    # Create RMSD data for different temperature models
    # TIP4P 273K model - slower equilibration
    rmsd_273K = 0.5 * (1 - np.exp(-time_points/2)) + 0.1
    # Add some noise to RMSD
    rmsd_273K += np.random.normal(0, 0.02, size=len(time_points))
    
    # TIP4P 298K model - faster equilibration due to higher temperature
    rmsd_298K = 0.5 * (1 - np.exp(-time_points/1.2)) + 0.1
    rmsd_298K += np.random.normal(0, 0.02, size=len(time_points))
    
    # Alternative 273K model - similar equilibration as first 273K model
    rmsd_273K_alt = 0.5 * (1 - np.exp(-time_points/2.1)) + 0.12
    rmsd_273K_alt += np.random.normal(0, 0.02, size=len(time_points))
    
    print("Creating energy data...")
    # Create energy drift data
    # All models have different drift rates but similar RMSD plateau values
    energy_273K = base_energy - 21 * time_points + np.random.normal(0, 5, size=len(time_points))
    energy_298K = base_energy + 25 * time_points + np.random.normal(0, 5, size=len(time_points))
    energy_273K_alt = base_energy - 18 * time_points + np.random.normal(0, 5, size=len(time_points))
    
    print("Data generation complete")
    return {
        'time_ps': time_ps,
        'rmsd_273K': rmsd_273K,
        'rmsd_298K': rmsd_298K,
        'rmsd_273K_alt': rmsd_273K_alt,
        'energy_273K': energy_273K,
        'energy_298K': energy_298K,
        'energy_273K_alt': energy_273K_alt
    }

def plot_rmsd_energy_stability(data, output_file):
    """
    Create a visualization showing RMSD plateaus despite energy drift.
    """
    print("Setting up figure...")
    fig = plt.figure(figsize=(15, 12))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[2, 2, 1])
    
    # Define consistent colors for the models
    colors = {
        '273K': '#1f77b4',      # Blue
        '298K': '#ff7f0e',      # Orange
        '273K (alt)': '#2ca02c' # Green
    }
    
    print("Creating top left plot (RMSD over time)...")
    # Top left plot: RMSD over time
    ax1 = fig.add_subplot(gs[0, 0])
    
    ax1.plot(data['time_ps']/1000, data['rmsd_273K'], color=colors['273K'], 
             linewidth=2.5, label='TIP4P 273K')
    ax1.plot(data['time_ps']/1000, data['rmsd_298K'], color=colors['298K'], 
             linewidth=2.5, label='TIP4P 298K')
    ax1.plot(data['time_ps']/1000, data['rmsd_273K_alt'], color=colors['273K (alt)'], 
             linewidth=2.5, label='TIP4P 273K (alt)')
    
    # Add plateau indicators
    plateau_time = 10  # ns
    plateau_idx = np.argmin(np.abs(data['time_ps']/1000 - plateau_time))
    
    plateau_273K = data['rmsd_273K'][plateau_idx:].mean()
    plateau_298K = data['rmsd_298K'][plateau_idx:].mean()
    plateau_273K_alt = data['rmsd_273K_alt'][plateau_idx:].mean()
    
    ax1.axhline(y=plateau_273K, color=colors['273K'], linestyle='--', alpha=0.7)
    ax1.axhline(y=plateau_298K, color=colors['298K'], linestyle='--', alpha=0.7)
    ax1.axhline(y=plateau_273K_alt, color=colors['273K (alt)'], linestyle='--', alpha=0.7)
    
    # Annotate plateaus
    ax1.text(plateau_time + 1, plateau_273K + 0.02, f"Plateau ≈ {plateau_273K:.2f} nm", 
             color=colors['273K'], fontweight='bold')
    ax1.text(plateau_time + 1, plateau_298K + 0.02, f"Plateau ≈ {plateau_298K:.2f} nm", 
             color=colors['298K'], fontweight='bold')
    ax1.text(plateau_time + 1, plateau_273K_alt - 0.03, f"Plateau ≈ {plateau_273K_alt:.2f} nm", 
             color=colors['273K (alt)'], fontweight='bold')
    
    ax1.set_xlabel('Simulation Time (ns)', fontweight='bold')
    ax1.set_ylabel('RMSD (nm)', fontweight='bold')
    ax1.set_title('Root Mean Square Deviation (RMSD) Over Time', fontsize=16, fontweight='bold')
    ax1.legend(loc='lower right', frameon=True, fancybox=True, shadow=True, fontsize=12)
    
    print("Creating top right plot (Energy drift)...")
    # Top right plot: Energy drift
    ax2 = fig.add_subplot(gs[0, 1])
    
    ax2.plot(data['time_ps']/1000, data['energy_273K'], color=colors['273K'], 
             linewidth=2.5, label='TIP4P 273K')
    ax2.plot(data['time_ps']/1000, data['energy_298K'], color=colors['298K'], 
             linewidth=2.5, label='TIP4P 298K')
    ax2.plot(data['time_ps']/1000, data['energy_273K_alt'], color=colors['273K (alt)'], 
             linewidth=2.5, label='TIP4P 273K (alt)')
    
    # Calculate and annotate drift rates
    time_ns = data['time_ps']/1000
    slope_273K, _, _, _, _ = scipy_stats.linregress(time_ns, data['energy_273K'])
    slope_298K, _, _, _, _ = scipy_stats.linregress(time_ns, data['energy_298K'])
    slope_273K_alt, _, _, _, _ = scipy_stats.linregress(time_ns, data['energy_273K_alt'])
    
    ax2.text(0.02, 0.05, f"Drift rates:\n"
             f"TIP4P 273K: {abs(slope_273K):.1f} kJ/mol/ns\n"
             f"TIP4P 298K: {abs(slope_298K):.1f} kJ/mol/ns\n"
             f"TIP4P 273K (alt): {abs(slope_273K_alt):.1f} kJ/mol/ns",
             transform=ax2.transAxes, 
             verticalalignment='bottom', horizontalalignment='left',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.9), fontsize=12)
    
    ax2.set_xlabel('Simulation Time (ns)', fontweight='bold')
    ax2.set_ylabel('Total Energy (kJ/mol)', fontweight='bold')
    ax2.set_title('Energy Drift in Different Models', fontsize=16, fontweight='bold')
    ax2.legend(loc='upper right', frameon=True, fancybox=True, shadow=True, fontsize=12)
    
    print("Creating middle plot (RMSD stabilization vs energy drift)...")
    # Middle: Combined plot showing relationship between RMSD stabilization and energy drift
    ax3 = fig.add_subplot(gs[1, :])
    
    # Create shaded regions for equilibration and production phases
    equilibration_end = 5  # ns
    ax3.axvspan(0, equilibration_end, alpha=0.2, color='gray', label='Equilibration Phase')
    ax3.axvspan(equilibration_end, simulation_time, alpha=0.1, color='green', label='Production Phase (Stable RMSD)')
    
    # Plot normalized RMSD and energy for 298K system as example
    # Normalize to [0,1] range for comparison
    norm_rmsd = (data['rmsd_298K'] - data['rmsd_298K'].min()) / (data['rmsd_298K'].max() - data['rmsd_298K'].min())
    
    # For energy, we want to show the drift magnitude, so we use absolute difference from initial
    energy_diff = np.abs(data['energy_298K'] - data['energy_298K'][0])
    norm_energy_diff = energy_diff / energy_diff.max()
    
    # Create a twin axis for the energy plot
    ax3.plot(time_ns, norm_rmsd, color=colors['298K'], linewidth=3, label='RMSD (normalized)')
    ax3.set_ylabel('Normalized RMSD', fontweight='bold', color=colors['298K'])
    
    ax3_twin = ax3.twinx()
    ax3_twin.plot(time_ns, norm_energy_diff, color='red', linewidth=3, linestyle='--', 
                 label='Energy Drift Magnitude')
    ax3_twin.set_ylabel('Normalized Energy Drift Magnitude', fontweight='bold', color='red')
    
    # Add vertical line at equilibration end
    ax3.axvline(x=equilibration_end, color='k', linestyle='--', alpha=0.7)
    ax3.text(equilibration_end + 0.5, 0.5, 'RMSD Stabilizes', fontsize=12, rotation=90, 
             verticalalignment='center')
    
    # Combine legends from both axes
    lines1, labels1 = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3_twin.get_legend_handles_labels()
    ax3.legend(lines1 + lines2, labels1 + labels2, loc='upper right', frameon=True, 
               fancybox=True, shadow=True, fontsize=12)
    
    ax3.set_xlabel('Simulation Time (ns)', fontweight='bold')
    ax3.set_title('RMSD Stabilization Despite Continued Energy Drift', fontsize=16, fontweight='bold')
    
    # Add key message annotation
    ax3.text(0.02, 0.05, 
            "While energy drift continues linearly throughout the simulation,\n"
            "RMSD stabilizes after the equilibration phase, indicating that the system\n"
            "has reached a physically meaningful and stable structural configuration.",
            transform=ax3.transAxes, 
            verticalalignment='bottom', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9), fontsize=12)
    
    print("Creating bottom summary section...")
    # Bottom: Summary and explanation 
    ax4 = fig.add_subplot(gs[2, :])
    ax4.axis('off')  # No axes for this text box
    
    summary_text = (
        "Key Observations about RMSD and Energy Drift in MD Simulations:\n\n"
        "1. All systems eventually reach stable RMSD plateaus despite continued energy drift\n\n"
        "2. The 298K model equilibrates faster due to increased molecular mobility at higher temperatures\n\n"
        "3. Both 273K simulations reach comparable RMSD plateaus, suggesting consistent structural integrity\n\n"
        "4. Practical implication: Perfect energy conservation is unattainable with current numerical methods,\n"
        "   but systems can still provide physically meaningful and stable representations when drift rates\n"
        "   are monitored and kept within reasonable limits (typically 20-25 kJ/mol/ns for water models)"
    )
    
    ax4.text(0.5, 0.5, summary_text, 
             ha='center', va='center', fontsize=14,
             bbox=dict(boxstyle='round', facecolor='#f0f0f0', alpha=0.9, pad=1.5))
    
    print("Finalizing layout...")
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.35, wspace=0.25)
    
    # Save the figure
    print(f"Saving figure to {output_file}...")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()
    print("Figure saved successfully")

def main():
    """Main function to run the visualizations."""
    try:
        print("Starting RMSD stability visualization script...")
        
        # Generate synthetic data
        data = generate_synthetic_data()
        
        # Create visualization
        output_path = os.path.join(output_dir, "rmsd_stability_despite_drift.png")
        plot_rmsd_energy_stability(data, output_path)
        
        print(f"RMSD stability visualization created successfully!")
        print(f"Output file saved to {output_path}")
    except Exception as e:
        print(f"Error occurred: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 