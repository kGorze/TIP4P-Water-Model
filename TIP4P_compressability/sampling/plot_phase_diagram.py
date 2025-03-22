#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import sys
import pandas as pd
from scipy.signal import savgol_filter
from scipy.constants import Boltzmann
from scipy.ndimage import gaussian_filter
from scipy.spatial import ConvexHull
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import matplotlib.colors as mcolors

# Stała Boltzmanna w kcal/(mol·K)
k_B_kcal = 0.0019872041

# Referencyjne punkty fazowe wody (MPa, °C)
WATER_PHASE_DATA = {
    'triple_points': [
        # [ciśnienie (MPa), temperatura (°C), opis]
        [0.000611657, 0.010, "Gas-Liquid-Ice Ih (Triple Point)"],
        [209.9, -21.985, "Liquid-Ice Ih-Ice III"],
        [350.1, -16.986, "Liquid-Ice III-Ice V"],
        [632.4, 0.16, "Liquid-Ice V-Ice VI"],
        [2216, 81.85, "Liquid-Ice VI-Ice VII"],
    ],
    'phase_lines': [
        # [ciśnienie start (MPa), temp start (°C), ciśnienie end (MPa), temp end (°C), opis]
        [0.000611657, 0.010, 209.9, -21.985, "Ice Ih-Liquid"],
        [0.000611657, 0.010, 0.000611657, 100, "Liquid-Gas"],
        [209.9, -21.985, 350.1, -16.986, "Ice III-Liquid"],
        [350.1, -16.986, 632.4, 0.16, "Ice V-Liquid"],
        [632.4, 0.16, 2216, 81.85, "Ice VI-Liquid"],
    ]
}

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
        
        # Convert pressure from bar to MPa and temperature from K to °C for consistency with reference data
        data['Pressure'] = data['Pressure'] * 0.1  # bar to MPa
        data['Temperature'] = data['Temperature'] - 273.15  # K to °C
        data['Time'] = data['Time'] / 1000  # ps to ns
        
        # Print some basic statistics to verify the data
        print("\nRaw Data Statistics:")
        print(f"Temperature: mean = {data['Temperature'].mean():.2f} °C, std = {data['Temperature'].std():.2f} °C")
        print(f"Min: {data['Temperature'].min():.2f} °C, Max: {data['Temperature'].max():.2f} °C")
        print(f"Pressure: mean = {data['Pressure'].mean():.2f} MPa, std = {data['Pressure'].std():.2f} MPa")
        print(f"Min: {data['Pressure'].min():.2f} MPa, Max: {data['Pressure'].max():.2f} MPa")
        print(f"Total number of data points: {len(data)}")
        
        return data
    except Exception as e:
        print(f"Error processing energy.xvg: {e}")
        sys.exit(1)
    finally:
        if os.path.exists('energy.xvg'):
            os.remove('energy.xvg')

def estimate_phase_boundaries(data, bins=80, sigma=1.0):
    """
    Estimate phase boundaries from simulation data using free energy analysis
    
    Returns dict with estimated phase regions and boundaries
    """
    # Calculate temperature and pressure ranges for histogram
    t_min, t_max = data['Temperature'].min(), data['Temperature'].max()
    p_min, p_max = data['Pressure'].min(), data['Pressure'].max()
    
    # Add padding to ranges
    t_padding = 0.05 * (t_max - t_min)
    p_padding = 0.05 * (p_max - p_min)
    
    t_min -= t_padding
    t_max += t_padding
    p_min -= p_padding
    p_max += p_padding
    
    # Create 2D histogram
    H, xedges, yedges = np.histogram2d(
        data['Temperature'], 
        data['Pressure'], 
        bins=bins, 
        range=[[t_min, t_max], [p_min, p_max]],
        density=True
    )
    
    # Apply Gaussian smoothing to reduce noise - increased sigma for better boundaries
    H_smooth = gaussian_filter(H, sigma=sigma)
    
    # Calculate free energy landscape
    # Use average temperature (in K) for kT
    avg_temp_K = data['Temperature'].mean() + 273.15
    kT = k_B_kcal * avg_temp_K
    
    # Calculate free energy: F = -kT*ln(P)
    epsilon = 1e-10  # Small constant to avoid log(0)
    F = -kT * np.log(H_smooth.T + epsilon)
    
    # Normalize free energy to have minimum at 0
    F -= np.min(F)
    
    # Get temperature and pressure coordinates
    temp_coords = 0.5 * (xedges[1:] + xedges[:-1])
    pres_coords = 0.5 * (yedges[1:] + yedges[:-1])
    
    # Find regions of low free energy (potential phases)
    # Use a threshold to identify separate basins
    # Adjusted threshold for better phase separation
    threshold = 0.15 * np.max(F)  
    low_energy_regions = (F < threshold)
    
    # For now, just return the computed data
    # In a more sophisticated analysis, we would identify distinct phases and their boundaries
    return {
        'temperature': temp_coords,
        'pressure': pres_coords,
        'free_energy': F,
        'low_energy_regions': low_energy_regions,
        'histogram': H_smooth,
        'threshold': threshold
    }

def plot_phase_diagram(data, phase_data, ref_data=None, output_file='phase_diagram.png'):
    """
    Plot phase diagram from simulation data and compare with reference data if provided
    """
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Extract coordinates and free energy data
    temp = phase_data['temperature']
    pres = phase_data['pressure']
    F = phase_data['free_energy']
    
    # Get the sigma value used for smoothing
    sigma = 1.0  # Default value
    if 'sigma' in globals():
        sigma = globals()['sigma']
    
    # Set up an improved custom colormap for better contrast
    # Use a diverging colormap with stronger contrast for better visualization of phase boundaries
    # Deep blue for stable phases, bright red for boundaries
    colors = [(0, 0, 0.8), (0.5, 0.5, 1), (1, 1, 1), (1, 0.5, 0.5), (0.8, 0, 0)]
    cmap = LinearSegmentedColormap.from_list('enhanced_phase_cmap', colors, N=256)
    
    # Normalize colormap with enhanced contrast
    vmin, vmax = np.percentile(F, [5, 95])  # Focus on the middle 90% of values
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    
    # Plot the free energy landscape with enhanced color contrast
    contourf = ax.contourf(temp, pres, F, levels=50, cmap=cmap, norm=norm)
    cbar = plt.colorbar(contourf, label='Free Energy (kcal/mol)')
    
    # Add contour lines with labels for better boundary visualization
    contour = ax.contour(temp, pres, F, levels=5, colors='k', alpha=0.6, linewidths=0.7)
    plt.clabel(contour, inline=True, fontsize=8, fmt='%.1f')
    
    # Set axis limits based on the data range
    t_min, t_max = temp.min(), temp.max()
    p_min, p_max = pres.min(), pres.max()
    
    # Implement split scale for pressure axis to better visualize the full range
    # Keep linear scale but focus on data range where points exist
    # Calculate where 95% of the data points fall
    p_range_min, p_range_max = np.percentile(data['Pressure'], [2.5, 97.5])
    
    # Set reasonable limits with padding
    t_padding = 0.1 * (t_max - t_min)
    p_padding = 0.1 * (p_range_max - p_range_min)
    
    ax.set_xlim(t_min - t_padding, t_max + t_padding)
    ax.set_ylim(p_range_min - p_padding, p_range_max + p_padding)
    
    # Overlay the data points with better transparency
    if len(data) > 5000:
        # Downsample for clarity
        sample_data = data.iloc[::len(data)//1000]
    else:
        sample_data = data
    ax.scatter(sample_data['Temperature'], sample_data['Pressure'], 
              color='grey', alpha=0.05, s=2, label='Simulation points')
    
    # If reference data is provided, add it to the plot
    if ref_data:
        # Add triple points
        triple_points = ref_data['triple_points']
        for point in triple_points:
            p, t, desc = point
            if t_min - t_padding <= t <= t_max + t_padding and p_range_min - p_padding <= p <= p_range_max + p_padding:
                ax.scatter(t, p, color='red', s=60, marker='o', edgecolors='black', zorder=5)
                ax.annotate(desc, (t, p), xytext=(10, 10), textcoords='offset points', 
                           fontsize=8, bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
        
        # Add phase lines
        phase_lines = ref_data['phase_lines']
        for line in phase_lines:
            p1, t1, p2, t2, desc = line
            # Only plot lines that fall within our visualization window
            if ((t_min - t_padding <= t1 <= t_max + t_padding or t_min - t_padding <= t2 <= t_max + t_padding) and
                (p_range_min - p_padding <= p1 <= p_range_max + p_padding or p_range_min - p_padding <= p2 <= p_range_max + p_padding)):
                ax.plot([t1, t2], [p1, p2], 'r--', linewidth=1.8, label=desc)
    
    # Add zoomed inset for the triple point region
    # Find region around the triple point (0.01°C, 0.000611657 MPa)
    triple_point_temp = 0.01
    triple_point_pres = 0.000611657
    
    # Create inset axes for zoomed region
    # Focus on region near triple point within actual data range
    inset_t_min = max(t_min, -5)
    inset_t_max = min(t_max, 5)
    inset_p_min = max(p_range_min, -10)
    inset_p_max = min(p_range_max, 10)
    
    # Only add inset if the triple point region is within our data
    if (inset_t_min < triple_point_temp < inset_t_max or 
        data['Temperature'].min() < triple_point_temp < data['Temperature'].max()):
        
        # Create inset
        axins = inset_axes(ax, width="40%", height="40%", loc='lower right',
                           bbox_to_anchor=(0.95, 0.05, 1, 1), bbox_transform=ax.transAxes)
        
        # Plot data in inset
        axins.contourf(temp, pres, F, levels=50, cmap=cmap, norm=norm)
        
        # Set inset limits
        axins.set_xlim(inset_t_min, inset_t_max)
        axins.set_ylim(inset_p_min, inset_p_max)
        
        # Add triple point and phase lines
        if ref_data:
            # Add triple point if it falls within our inset range
            if inset_t_min <= triple_point_temp <= inset_t_max and inset_p_min <= triple_point_pres <= inset_p_max:
                axins.scatter(triple_point_temp, triple_point_pres, 
                             color='red', s=100, marker='o', edgecolors='black', zorder=5)
                axins.annotate("Triple Point", (triple_point_temp, triple_point_pres), 
                              xytext=(20, 20), textcoords='offset points',
                              fontsize=8, bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7),
                              arrowprops=dict(arrowstyle='->', color='black'))
            
            # Add relevant phase lines
            for line in phase_lines:
                p1, t1, p2, t2, desc = line
                if ((inset_t_min <= t1 <= inset_t_max or inset_t_min <= t2 <= inset_t_max) and
                    (inset_p_min <= p1 <= inset_p_max or inset_p_min <= p2 <= inset_p_max)):
                    axins.plot([t1, t2], [p1, p2], 'r--', linewidth=1.5)
        
        # Add connecting lines to show the zoomed region
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
        axins.set_xlabel('T (°C)', fontsize=8)
        axins.set_ylabel('P (MPa)', fontsize=8)
        axins.tick_params(axis='both', which='major', labelsize=6)
    
    # Add phase labels for clearer interpretation
    # Identify potential phases based on free energy minima
    F_smooth = gaussian_filter(F, sigma=2.0)  # Extra smoothing to find major minima
    local_min = (F_smooth < np.roll(F_smooth, 1, axis=0)) & (F_smooth < np.roll(F_smooth, -1, axis=0)) & \
                (F_smooth < np.roll(F_smooth, 1, axis=1)) & (F_smooth < np.roll(F_smooth, -1, axis=1))
    
    # Find the indices of the local minima
    min_indices = np.where(local_min)
    
    # Label the potential phases if we have multiple local minima
    if len(min_indices[0]) > 0:
        # Get values at minima
        min_values = F_smooth[min_indices]
        # Sort by depth of minima (lowest free energy first)
        sorted_idx = np.argsort(min_values)
        
        # Add labels for the top 3 minima if we have that many
        phase_labels = ["Potential Phase 1", "Potential Phase 2", "Potential Phase 3"]
        for i in range(min(len(sorted_idx), 3)):
            idx = sorted_idx[i]
            x_idx, y_idx = min_indices[1][idx], min_indices[0][idx]
            t_val, p_val = temp[x_idx], pres[y_idx]
            
            # Only add label if within visible range
            if ax.get_xlim()[0] <= t_val <= ax.get_xlim()[1] and ax.get_ylim()[0] <= p_val <= ax.get_ylim()[1]:
                ax.annotate(phase_labels[i], (t_val, p_val), 
                          xytext=(0, 10), textcoords='offset points',
                          fontsize=10, ha='center',
                          bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))
    
    # Add labels and title
    ax.set_xlabel('Temperature (°C)', fontsize=12)
    ax.set_ylabel('Pressure (MPa)', fontsize=12)
    ax.set_title('Phase Diagram from Simulation Data', fontsize=14)
    
    # Handle duplicate legend entries (if any phase lines have the same description)
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper left', fontsize=9)
    
    # Add info text with more details
    info_text = (f"Data points: {len(data)}\n"
                f"Temperature range: {data['Temperature'].min():.1f} to {data['Temperature'].max():.1f} °C\n"
                f"Pressure range: {data['Pressure'].min():.1f} to {data['Pressure'].max():.1f} MPa\n"
                f"Free energy threshold: {phase_data['threshold']:.2f} kcal/mol")
    
    plt.figtext(0.02, 0.02, info_text, fontsize=9, 
               bbox=dict(facecolor='white', alpha=0.8))
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nEnhanced phase diagram saved as: {output_file}")
    plt.close()

def plot_detailed_comparison(data, phase_data, ref_data, output_file='phase_comparison.png'):
    """
    Create a detailed comparison of simulation phase diagram with reference data
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9))
    
    # Extract coordinates and free energy data
    temp = phase_data['temperature']
    pres = phase_data['pressure']
    F = phase_data['free_energy']
    
    # Set up colormaps
    # For simulation data: Low energy (stable phases) in blue, high energy (boundaries) in red
    colors = [(0, 0, 1), (1, 1, 1), (1, 0, 0)]
    cmap = LinearSegmentedColormap.from_list('phase_cmap', colors, N=256)
    
    # Plot simulation-based phase diagram (left panel)
    cf1 = ax1.contourf(temp, pres, F, levels=50, cmap=cmap)
    ax1.contour(temp, pres, F, levels=10, colors='k', alpha=0.3, linewidths=0.5)
    
    # Overlay sparse sample of data points
    if len(data) > 5000:
        sample_data = data.iloc[::len(data)//500]
    else:
        sample_data = data
    ax1.scatter(sample_data['Temperature'], sample_data['Pressure'], 
                color='grey', alpha=0.1, s=1)
    
    # Add reference to left panel for comparison
    if ref_data:
        # Add triple points
        triple_points = ref_data['triple_points']
        for point in triple_points:
            p, t, desc = point
            if t > -50 and t < 150 and p < 1000:  # Filter to relevant range
                ax1.scatter(t, p, color='red', s=50, marker='o', edgecolors='black', zorder=5)
        
        # Add phase lines (limited to visible range)
        phase_lines = ref_data['phase_lines']
        for line in phase_lines:
            p1, t1, p2, t2, desc = line
            if (t1 > -50 and t1 < 150 and p1 < 1000) or (t2 > -50 and t2 < 150 and p2 < 1000):
                ax1.plot([t1, t2], [p1, p2], 'r--', linewidth=1.5)
    
    # Add colorbar
    cbar1 = plt.colorbar(cf1, ax=ax1)
    cbar1.set_label('Free Energy (kcal/mol)')
    
    # Set axis labels and title
    ax1.set_xlabel('Temperature (°C)')
    ax1.set_ylabel('Pressure (MPa)')
    ax1.set_title('Phase Diagram from Simulation')
    
    # Plot reference phase diagram (right panel)
    # Create a cleaner visualization of the reference data
    
    # Set reasonable limits for visualization
    t_min_ref, t_max_ref = -30, 100
    p_min_ref, p_max_ref = 0, 700
    
    ax2.set_xlim(t_min_ref, t_max_ref)
    ax2.set_ylim(p_min_ref, p_max_ref)
    
    # Set log scale for pressure to better visualize the full range
    ax2.set_yscale('log')
    
    # Shade different phase regions with distinct colors
    # This is simplified - a proper diagram would require exact phase boundaries
    
    # Define some colors for different phases
    phase_colors = {
        'Ice Ih': 'lightblue',
        'Ice III': 'skyblue',
        'Ice V': 'deepskyblue',
        'Ice VI': 'royalblue',
        'Liquid': 'blue',
        'Gas': 'lightgrey'
    }
    
    # Helper function to create a polygon from points defining phase boundary
    def create_phase_polygon(points, closed=True):
        if closed and points[0] != points[-1]:
            points.append(points[0])  # Close the polygon
        return Polygon(points, alpha=0.3, edgecolor='black', linewidth=0.5)
    
    # Add phase lines
    for line in ref_data['phase_lines']:
        p1, t1, p2, t2, desc = line
        ax2.plot([t1, t2], [p1, p2], 'k-', linewidth=2, label=desc)
    
    # Add triple points with annotations
    for point in ref_data['triple_points']:
        p, t, desc = point
        if t_min_ref <= t <= t_max_ref and p_min_ref <= p <= p_max_ref:
            ax2.scatter(t, p, color='red', s=80, marker='o', edgecolors='black', zorder=5)
            ax2.annotate(desc, (t, p), xytext=(10, 5), textcoords='offset points',
                        fontsize=9, bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
    
    # Add labels and annotations for major phases
    ax2.text(50, 0.1, "Gas", fontsize=14, ha='center')
    ax2.text(50, 10, "Liquid", fontsize=14, ha='center')
    ax2.text(-20, 100, "Ice Ih", fontsize=14, ha='center')
    ax2.text(-15, 300, "Ice III", fontsize=12, ha='center')
    ax2.text(0, 500, "Ice V/VI", fontsize=12, ha='center')
    
    # Add critical point (for water: T ≈ 374°C, P ≈ 22.1 MPa)
    ax2.scatter(374, 22.1, color='black', s=100, marker='*', edgecolors='white', zorder=5)
    ax2.annotate("Critical Point\n(374°C, 22.1 MPa)", (374, 22.1), xytext=(-40, -40),
                textcoords='offset points', fontsize=10, 
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
    
    # Set axis labels and title
    ax2.set_xlabel('Temperature (°C)')
    ax2.set_ylabel('Pressure (MPa)')
    ax2.set_title('Reference Phase Diagram of Water')
    
    # Handle duplicate legend entries
    handles, labels = ax2.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax2.legend(by_label.values(), by_label.keys(), loc='upper right', fontsize=9)
    
    # Add descriptive text
    fig.text(0.5, 0.01, 
             "Comparison of simulation-derived phase behavior with reference water phase diagram.\n"
             "Note: Pressure scales differ between plots. The simulation data covers a limited region of the phase space.",
             ha='center', fontsize=10)
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nDetailed comparison saved as: {output_file}")
    plt.close()

def main():
    # Set paths
    edr_file = "../md.edr"  # Path to your EDR file
    tpr_file = "../md.tpr"  # Path to your TPR file
    
    try:
        # Read the energy data
        print("Reading energy data from GROMACS files...")
        data = run_gmx_energy(edr_file, tpr_file)
        
        # Estimate phase boundaries from simulation data
        print("\nAnalyzing phase behavior...")
        phase_data = estimate_phase_boundaries(data, bins=80, sigma=1.0)
        
        # Create simple phase diagram
        print("\nCreating phase diagram...")
        plot_phase_diagram(data, phase_data, ref_data=WATER_PHASE_DATA)
        
        # Create detailed comparison with reference data
        print("\nCreating detailed comparison with reference phase diagram...")
        plot_detailed_comparison(data, phase_data, WATER_PHASE_DATA)
        
        print("\nDone!")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 