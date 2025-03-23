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
from scipy import stats

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

def estimate_phase_boundaries(data, bins=80, sigma=1.5):
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
    
    # Apply Gaussian smoothing to reduce noise - increased sigma for cleaner boundaries
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
    # Improve threshold selection - use percentile instead of fixed fraction
    threshold = np.percentile(F, 15)  # Use lowest 15% of free energy values
    low_energy_regions = (F < threshold)
    
    # For now, just return the computed data
    # In a more sophisticated analysis, we would identify distinct phases and their boundaries
    return {
        'temperature': temp_coords,
        'pressure': pres_coords,
        'free_energy': F,
        'low_energy_regions': low_energy_regions,
        'histogram': H_smooth,
        'threshold': threshold,
        'sigma': sigma  # Store sigma value for later reference
    }

def plot_phase_diagram(data, phase_data, ref_data=None, output_file='phase_diagram.png'):
    """
    Plot publication-quality phase diagram from simulation data and compare with reference data if provided
    """
    # Use a clean, professional style
    plt.style.use('seaborn-v0_8-white')
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Extract coordinates and free energy data
    temp = phase_data['temperature']
    pres = phase_data['pressure']
    F = phase_data['free_energy']
    
    # Get the sigma value used for smoothing
    sigma = phase_data.get('sigma', 1.5)  # Use value from phase_data or default to 1.5
    
    # Create a custom scientific colormap with better contrast between phases
    # Blue for stable phases, white for intermediate, red for boundaries
    from matplotlib.colors import ListedColormap
    colors = [(0, 0, 0.7), (0.1, 0.1, 0.9), (0.4, 0.4, 1), (0.8, 0.8, 1), 
              (1, 1, 1), 
              (1, 0.8, 0.8), (1, 0.4, 0.4), (0.9, 0.1, 0.1), (0.7, 0, 0)]
    cmap = LinearSegmentedColormap.from_list('enhanced_phase_cmap', colors, N=256)
    
    # Normalize colormap with enhanced contrast - adaptive percentile approach
    vmin, vmax = np.percentile(F, [2, 98])  # Focus on the middle 96% of values for better contrast
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    
    # Plot the free energy landscape with enhanced color contrast and more levels for smoother boundaries
    contourf = ax.contourf(temp, pres, F, levels=100, cmap=cmap, norm=norm, antialiased=True)
    
    # Create a more elegant, scientifically-styled colorbar
    cbar = plt.colorbar(contourf, ax=ax, pad=0.02)
    cbar.set_label('Free Energy (kcal/mol)', fontsize=12, fontweight='bold')
    cbar.ax.tick_params(labelsize=10)
    
    # Add contour lines with labels for better boundary visualization - fewer levels for clarity
    contour = ax.contour(temp, pres, F, levels=4, colors='k', alpha=0.5, linewidths=0.7)
    plt.clabel(contour, inline=True, fontsize=8, fmt='%.1f')
    
    # Calculate optimal axis limits - focus on the meaningful data range
    t_min, t_max = temp.min(), temp.max()
    p_min, p_max = pres.min(), pres.max()
    
    # Use percentile approach to avoid outliers skewing the limits
    p_range_min, p_range_max = np.percentile(data['Pressure'], [1, 99])
    
    # Set reasonable limits with padding
    t_padding = 0.1 * (t_max - t_min)
    p_padding = 0.05 * (p_range_max - p_range_min)  # Reduced padding for better focus
    
    ax.set_xlim(t_min - t_padding, t_max + t_padding)
    ax.set_ylim(p_range_min - p_padding, p_range_max + p_padding)
    
    # Overlay the data points with better transparency - just enough to show the sampling
    if len(data) > 5000:
        # Downsample for clarity
        sample_data = data.iloc[::len(data)//2000]  # More points for better representation
    else:
        sample_data = data
    
    # Use very small and transparent points to show data without overwhelming the plot
    ax.scatter(sample_data['Temperature'], sample_data['Pressure'], 
              color='grey', alpha=0.03, s=1, label='Simulation points', rasterized=True)
    
    # If reference data is provided, add it to the plot with enhanced styling
    if ref_data:
        # Add triple points
        triple_points = ref_data['triple_points']
        for point in triple_points:
            p, t, desc = point
            if t_min - t_padding <= t <= t_max + t_padding and p_range_min - p_padding <= p <= p_range_max + p_padding:
                ax.scatter(t, p, color='red', s=70, marker='o', edgecolors='black', linewidth=1.5, zorder=5)
                ax.annotate(desc.split(" (")[0], (t, p), xytext=(10, 10), textcoords='offset points', 
                           fontsize=9, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9, ec="grey"))
        
        # Add phase lines with improved styling
        phase_lines = ref_data['phase_lines']
        for line in phase_lines:
            p1, t1, p2, t2, desc = line
            # Only plot lines that fall within our visualization window
            if ((t_min - t_padding <= t1 <= t_max + t_padding or t_min - t_padding <= t2 <= t_max + t_padding) and
                (p_range_min - p_padding <= p1 <= p_range_max + p_padding or p_range_min - p_padding <= p2 <= p_range_max + p_padding)):
                ax.plot([t1, t2], [p1, p2], 'r--', linewidth=2, label=desc)
    
    # Add zoomed inset for the triple point region
    # Find region around the triple point (0.01°C, 0.000611657 MPa)
    triple_point_temp = 0.01
    triple_point_pres = 0.000611657
    
    # Create inset axes for zoomed region if relevant to our data
    # Focus on region near triple point within actual data range
    inset_t_min = max(t_min, -3)  # More focused range
    inset_t_max = min(t_max, 3)
    inset_p_min = max(p_range_min, -5)
    inset_p_max = min(p_range_max, 5)
    
    # Only add inset if the triple point region is within our data and adds useful information
    in_range = (inset_t_min < triple_point_temp < inset_t_max or 
                data['Temperature'].min() < triple_point_temp < data['Temperature'].max())
    near_pressure = np.min(np.abs(data['Pressure'] - triple_point_pres)) < 5
    
    if in_range and near_pressure:
        # Create inset with better style
        axins = inset_axes(ax, width="35%", height="35%", loc='lower right',
                           bbox_to_anchor=(0.97, 0.03, 1, 1), bbox_transform=ax.transAxes)
        
        # Plot data in inset with same styling
        axins.contourf(temp, pres, F, levels=100, cmap=cmap, norm=norm)
        
        # Set inset limits
        axins.set_xlim(inset_t_min, inset_t_max)
        axins.set_ylim(inset_p_min, inset_p_max)
        
        # Add triple point and phase lines
        if ref_data:
            # Add triple point if it falls within our inset range
            if inset_t_min <= triple_point_temp <= inset_t_max and inset_p_min <= triple_point_pres <= inset_p_max:
                axins.scatter(triple_point_temp, triple_point_pres, 
                             color='red', s=120, marker='o', edgecolors='black', zorder=5)
                axins.annotate("Triple Point", (triple_point_temp, triple_point_pres), 
                              xytext=(20, 20), textcoords='offset points',
                              fontsize=9, fontweight='bold', 
                              bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9, ec="grey"),
                              arrowprops=dict(arrowstyle='->', color='black', linewidth=1.5))
            
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
        for i in range(min(len(sorted_idx), 3)):
            idx = sorted_idx[i]
            x_idx, y_idx = min_indices[1][idx], min_indices[0][idx]
            t_val, p_val = temp[x_idx], pres[y_idx]
            
            # Use physically meaningful labels
            phase_label = assign_phase_label(t_val, p_val)
            
            # Only add label if within visible range and avoid overlapping labels
            if ax.get_xlim()[0] <= t_val <= ax.get_xlim()[1] and ax.get_ylim()[0] <= p_val <= ax.get_ylim()[1]:
                # Enhanced label styling
                ax.annotate(phase_label, (t_val, p_val), 
                          xytext=(0, 10), textcoords='offset points',
                          fontsize=11, fontweight='bold', ha='center',
                          bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9, ec="grey"))
    
    # Add labels and title with enhanced styling
    ax.set_xlabel('Temperature (°C)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Pressure (MPa)', fontsize=14, fontweight='bold')
    ax.set_title('Phase Diagram from Simulation Data', fontsize=16, fontweight='bold')
    
    # Handle duplicate legend entries (if any phase lines have the same description)
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    legend = ax.legend(by_label.values(), by_label.keys(), loc='upper left', fontsize=10, framealpha=0.9)
    legend.get_frame().set_edgecolor('grey')
    
    # Add info text with more details - styled as a scientific footnote
    info_text = (f"Data points: {len(data)}\n"
                f"Temperature range: {data['Temperature'].min():.1f} to {data['Temperature'].max():.1f} °C\n"
                f"Pressure range: {data['Pressure'].min():.1f} to {data['Pressure'].max():.1f} MPa\n"
                f"Free energy threshold: {phase_data['threshold']:.2f} kcal/mol\n"
                f"Smoothing sigma: {sigma:.1f}")
    
    plt.figtext(0.02, 0.02, info_text, fontsize=9, 
               bbox=dict(facecolor='white', alpha=0.9, ec="grey", boxstyle="round,pad=0.5"))
    
    # Add grid for better readability with scientific styling
    ax.grid(True, linestyle='--', alpha=0.3)
    
    # Improve tick labels
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    # Adjust layout and save with high resolution
    plt.tight_layout()
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    print(f"\nPublication-quality phase diagram saved as: {output_file}")
    plt.close()

def plot_detailed_comparison(data, phase_data, ref_data, output_file='phase_comparison.png'):
    """
    Create a publication-quality detailed comparison of simulation phase diagram with reference data
    """
    # Use a clean, professional style
    plt.style.use('seaborn-v0_8-white')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Extract coordinates and free energy data
    temp = phase_data['temperature']
    pres = phase_data['pressure']
    F = phase_data['free_energy']
    
    # Get the sigma value used for smoothing
    sigma = phase_data.get('sigma', 1.5)
    
    # Set reasonable limits for visualization
    t_min_ref, t_max_ref = -30, 100
    p_min_ref, p_max_ref = 0, 700
    
    # Calculate reasonable limits for simulation data
    t_sim_min, t_sim_max = data['Temperature'].min(), data['Temperature'].max()
    p_sim_min, p_sim_max = data['Pressure'].min(), data['Pressure'].max()
    
    # Add padding to simulation ranges
    t_padding = 0.1 * (t_sim_max - t_sim_min)
    p_padding = 0.1 * (p_sim_max - p_sim_min)
    
    t_sim_min -= t_padding
    t_sim_max += t_padding
    p_sim_min -= p_padding if p_sim_min > 0 else 0  # Only add negative padding if p_min is positive
    p_sim_max += p_padding
    
    # Set up enhanced colormaps for simulation data
    # Scientific colormap: Blue for stable phases, white for intermediate, red for boundaries
    colors = [(0, 0, 0.7), (0.1, 0.1, 0.9), (0.4, 0.4, 1), (0.8, 0.8, 1), 
              (1, 1, 1), 
              (1, 0.8, 0.8), (1, 0.4, 0.4), (0.9, 0.1, 0.1), (0.7, 0, 0)]
    cmap = LinearSegmentedColormap.from_list('enhanced_phase_cmap', colors, N=256)
    
    # Normalize colormap with enhanced contrast
    vmin, vmax = np.percentile(F, [2, 98])
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    
    # Plot simulation-based phase diagram (left panel)
    cf1 = ax1.contourf(temp, pres, F, levels=100, cmap=cmap, norm=norm, antialiased=True)
    ax1.contour(temp, pres, F, levels=5, colors='k', alpha=0.3, linewidths=0.5)
    
    # Overlay sparse sample of data points - more subtle to highlight the phases
    if len(data) > 5000:
        sample_data = data.iloc[::len(data)//2000]  # More points but still manageable
    else:
        sample_data = data
    ax1.scatter(sample_data['Temperature'], sample_data['Pressure'], 
                color='grey', alpha=0.03, s=1, rasterized=True)
    
    # Add reference to left panel for comparison with enhanced styling
    if ref_data:
        # Add triple points with better styling
        triple_points = ref_data['triple_points']
        for point in triple_points:
            p, t, desc = point
            if t_sim_min <= t <= t_sim_max and p_sim_min <= p <= p_sim_max:
                ax1.scatter(t, p, color='red', s=70, marker='o', edgecolors='black', linewidth=1.5, zorder=5)
                # Shorter annotation to avoid clutter
                short_desc = desc.split(" (")[0]
                ax1.annotate(short_desc, (t, p), xytext=(10, 5), textcoords='offset points',
                           fontsize=9, fontweight='bold', 
                           bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9, ec="grey"))
        
        # Add phase lines (limited to visible range) with improved styling
        phase_lines = ref_data['phase_lines']
        for line in phase_lines:
            p1, t1, p2, t2, desc = line
            if ((t_sim_min <= t1 <= t_sim_max or t_sim_min <= t2 <= t_sim_max) and
                (p_sim_min <= p1 <= p_sim_max or p_sim_min <= p2 <= p_sim_max)):
                ax1.plot([t1, t2], [p1, p2], 'r--', linewidth=2)
    
    # Add colorbar with improved styling
    cbar1 = plt.colorbar(cf1, ax=ax1, pad=0.01)
    cbar1.set_label('Free Energy (kcal/mol)', fontsize=12, fontweight='bold')
    cbar1.ax.tick_params(labelsize=10)
    
    # Add phase labels for simulation data with improved styling
    # Identify potential phases based on free energy minima
    F_smooth = gaussian_filter(F, sigma=2.0)
    local_min = (F_smooth < np.roll(F_smooth, 1, axis=0)) & (F_smooth < np.roll(F_smooth, -1, axis=0)) & \
                (F_smooth < np.roll(F_smooth, 1, axis=1)) & (F_smooth < np.roll(F_smooth, -1, axis=1))
    
    # Find and label the minima
    min_indices = np.where(local_min)
    sim_phase_labels = []
    
    if len(min_indices[0]) > 0:
        min_values = F_smooth[min_indices]
        sorted_idx = np.argsort(min_values)
        
        for i in range(min(len(sorted_idx), 3)):
            idx = sorted_idx[i]
            x_idx, y_idx = min_indices[1][idx], min_indices[0][idx]
            t_val, p_val = temp[x_idx], pres[y_idx]
            
            # Use physically meaningful labels
            phase_label = assign_phase_label(t_val, p_val)
            sim_phase_labels.append((t_val, p_val, phase_label))
            
            # Only add label if within visible range
            if ax1.get_xlim()[0] <= t_val <= ax1.get_xlim()[1] and ax1.get_ylim()[0] <= p_val <= ax1.get_ylim()[1]:
                ax1.annotate(phase_label, (t_val, p_val), 
                           xytext=(0, 10), textcoords='offset points',
                           fontsize=11, fontweight='bold', ha='center',
                           bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9, ec="grey"))
    
    # Set axis labels and title with enhanced styling
    ax1.set_xlabel('Temperature (°C)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Pressure (MPa)', fontsize=12, fontweight='bold')
    ax1.set_title('Phase Diagram from Simulation', fontsize=14, fontweight='bold')
    
    # Add grid for better readability
    ax1.grid(True, linestyle='--', alpha=0.3)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    
    # Plot reference phase diagram (right panel) with enhanced styling
    ax2.set_xlim(t_min_ref, t_max_ref)
    ax2.set_ylim(p_min_ref, p_max_ref)
    
    # Use log scale for pressure to better visualize the full range, with minor ticks
    ax2.set_yscale('log')
    ax2.yaxis.set_minor_formatter(plt.FormatStrFormatter("%.1f"))
    ax2.yaxis.set_minor_locator(plt.LogLocator(subs=[0.2, 0.4, 0.6, 0.8]))
    
    # Create a better visualization of the different phase regions
    # Define colors for different phases - use color blind friendly palette
    phase_colors = {
        'Ice Ih': '#bad0ff',  # Light blue
        'Ice III': '#8fbaff', # Medium blue
        'Ice V': '#5b9aff',   # Darker blue 
        'Ice VI': '#366fda',  # Deep blue
        'Liquid': '#0048c0',  # Very deep blue
        'Gas': '#f0f0f0'      # Light grey
    }
    
    # Helper function to create a polygon from points defining phase boundary
    def create_phase_polygon(points, closed=True):
        if closed and points[0] != points[-1]:
            points.append(points[0])  # Close the polygon
        return Polygon(points, alpha=0.4, edgecolor='black', linewidth=0.7)
    
    # Add simplified phase regions using polygon shapes where possible
    # These are approximations for visualization purposes
    
    # Gas phase region (approximate)
    gas_points = [
        [0, 0.0006], 
        [100, 0.0006], 
        [100, 0.1], 
        [0, 0.1]
    ]
    gas_poly = create_phase_polygon(gas_points)
    ax2.add_patch(gas_poly)
    
    # Liquid region (approximate)
    liquid_points = [
        [0, 0.0006],
        [0, 200],
        [-21.985, 209.9],
        [-21.985, 0.0006]
    ]
    liquid_poly = create_phase_polygon(liquid_points)
    ax2.add_patch(liquid_poly)
    
    # Ice Ih region (approximate)
    ice_Ih_points = [
        [-21.985, 0.0006],
        [-21.985, 209.9],
        [-50, 209.9],
        [-50, 0.0006]
    ]
    ice_Ih_poly = create_phase_polygon(ice_Ih_points)
    ax2.add_patch(ice_Ih_poly)
    
    # Add phase lines with improved styling - bolder lines
    for line in ref_data['phase_lines']:
        p1, t1, p2, t2, desc = line
        ax2.plot([t1, t2], [p1, p2], 'k-', linewidth=2, label=desc)
    
    # Add triple points with enhanced annotations
    for point in ref_data['triple_points']:
        p, t, desc = point
        if t_min_ref <= t <= t_max_ref and p_min_ref <= p <= p_max_ref:
            ax2.scatter(t, p, color='red', s=80, marker='o', edgecolors='black', linewidth=1.5, zorder=5)
            # Shorter, cleaner labels
            short_desc = desc.split(" (")[0]
            ax2.annotate(short_desc, (t, p), xytext=(10, 5), textcoords='offset points',
                        fontsize=9, fontweight='bold',
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9, ec="grey"))
    
    # Add labels and annotations for major phases with improved styling
    ax2.text(50, 0.05, "Gas", fontsize=14, ha='center', fontweight='bold')
    ax2.text(50, 10, "Liquid", fontsize=14, ha='center', fontweight='bold')
    ax2.text(-20, 100, "Ice Ih", fontsize=14, ha='center', fontweight='bold')
    ax2.text(-15, 300, "Ice III", fontsize=12, ha='center', fontweight='bold')
    ax2.text(0, 500, "Ice V/VI", fontsize=12, ha='center', fontweight='bold')
    
    # Add critical point (for water: T ≈ 374°C, P ≈ 22.1 MPa) with improved styling
    ax2.scatter(374, 22.1, color='black', s=120, marker='*', edgecolors='white', linewidth=1.5, zorder=5)
    ax2.annotate("Critical Point\n(374°C, 22.1 MPa)", (374, 22.1), xytext=(-40, -40),
                textcoords='offset points', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9, ec="grey"),
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2", linewidth=1.5))
    
    # Set axis labels and title with improved styling
    ax2.set_xlabel('Temperature (°C)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Pressure (MPa)', fontsize=12, fontweight='bold')
    ax2.set_title('Reference Phase Diagram of Water', fontsize=14, fontweight='bold')
    
    # Add grid for better readability
    ax2.grid(True, linestyle='--', alpha=0.3)
    ax2.tick_params(axis='both', which='major', labelsize=10)
    
    # Handle duplicate legend entries with improved styling
    handles, labels = ax2.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    legend = ax2.legend(by_label.values(), by_label.keys(), loc='upper right', fontsize=9, framealpha=0.9)
    legend.get_frame().set_edgecolor('grey')
    
    # Generate phase mapping text based on identified phases in simulation
    mapping_text = "Simulation phase regions map to reference as follows:\n"
    
    for i, (t, p, label) in enumerate(sim_phase_labels[:3]):  # Limit to top 3 phases
        ref_phase = "Unknown"
        if label == "Liquid Water":
            ref_phase = "Liquid"
        elif label == "Ice Ih-like":
            ref_phase = "Ice Ih"
        elif label == "Ice III/V/VI-like":
            ref_phase = "Ice III/V/VI"
        elif label == "Stretched Liquid":
            ref_phase = "Metastable liquid (negative P)"
        elif label == "Liquid-Gas Transition":
            ref_phase = "Liquid-Gas coexistence"
            
        mapping_text += f"- {label} (T={t:.1f}°C, P={p:.1f}MPa) → {ref_phase}\n"
    
    # Add descriptive text with improved styling
    fig.text(0.5, 0.01, mapping_text +
             "\nNote: Pressure scales differ between plots. The simulation data covers a limited region of phase space.\n"
             "Negative pressures in simulation may correspond to metastable stretched states not shown in reference diagram.",
             ha='center', fontsize=10, bbox=dict(facecolor='white', alpha=0.9, ec="grey", boxstyle="round,pad=0.5"))
    
    # Add connecting arrows to show corresponding phases between the two diagrams
    if sim_phase_labels:
        for t, p, label in sim_phase_labels:
            if label == "Liquid Water":
                # Draw arrow from simulation liquid to reference liquid
                fig.patches.extend([plt.Arrow(0.3, 0.6, 0.2, 0, width=0.05, facecolor='grey', alpha=0.5)])
            elif label == "Ice Ih-like":
                # Draw arrow from simulation ice to reference ice
                fig.patches.extend([plt.Arrow(0.3, 0.4, 0.2, -0.1, width=0.05, facecolor='grey', alpha=0.5)])
    
    # Adjust layout and save with high resolution
    plt.tight_layout()
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    print(f"\nPublication-quality comparison saved as: {output_file}")
    plt.close()

def filter_data(data, p_threshold=-20, z_threshold=3.0):
    """
    Filter out data points with unrealistic negative pressures and statistical outliers
    
    Parameters:
    -----------
    data : pandas DataFrame
        The raw simulation data
    p_threshold : float
        Threshold for filtering extreme negative pressures (in MPa)
    z_threshold : float
        Z-score threshold for identifying statistical outliers
        
    Returns:
    --------
    data_filtered : pandas DataFrame
        The filtered simulation data
    """
    # Step 1: Filter out extreme negative pressures
    reasonable_data = data[data['Pressure'] > p_threshold].copy()
    
    # Step 2: Apply statistical methods to identify and remove outliers
    z_scores = stats.zscore(reasonable_data[['Temperature', 'Pressure']])
    abs_z_scores = np.abs(z_scores)
    filtered_entries = (abs_z_scores < z_threshold).all(axis=1)
    data_filtered = reasonable_data[filtered_entries]
    
    # Print statistics about filtered data
    print("\nFiltering Results:")
    print(f"Original data points: {len(data)}")
    print(f"After negative pressure filter (P > {p_threshold} MPa): {len(reasonable_data)}")
    print(f"After outlier filter (|z| < {z_threshold}): {len(data_filtered)}")
    print(f"Temperature range after filtering: {data_filtered['Temperature'].min():.1f} to {data_filtered['Temperature'].max():.1f} °C")
    print(f"Pressure range after filtering: {data_filtered['Pressure'].min():.1f} to {data_filtered['Pressure'].max():.1f} MPa")
    
    return data_filtered

def assign_phase_label(t_val, p_val):
    """
    Assign physically meaningful phase labels based on temperature and pressure values
    
    Parameters:
    -----------
    t_val : float
        Temperature value in °C
    p_val : float
        Pressure value in MPa
        
    Returns:
    --------
    label : str
        Physically meaningful phase label
    """
    # Major phase regions with more precise criteria
    
    # Triple point of water
    triple_point_dist = np.sqrt((t_val - 0.01)**2 + (p_val - 0.000611657)**2)
    if triple_point_dist < 0.5:  # Near triple point
        return "Gas-Liquid-Ice Ih\n(Triple Point)"
        
    # Negative pressure regimes
    if p_val < -10:
        if t_val < 0:
            return "Stretched Supercooled\nLiquid"
        else:
            return "Stretched Liquid"
    elif p_val < 0:
        if t_val < 0:
            return "Metastable\nSupercooled Liquid"
        else:
            return "Metastable Liquid"
            
    # Low pressure regimes
    if p_val < 0.1:
        if t_val > 0.1 and t_val < 100:
            return "Liquid-Gas\nTransition"
        elif t_val <= 0:
            return "Ice Ih"
    
    # Normal and high pressure regimes        
    if p_val < 200:
        if t_val < 0:
            return "Ice Ih"
        else:
            return "Liquid Water"
    elif p_val < 350:
        if t_val < -10:
            return "Ice III"
        elif t_val < 0:
            return "Ice Ih/III\nTransition"
        else:
            return "Liquid (High P)"
    elif p_val < 632:
        if t_val < -5:
            return "Ice V"
        else:
            return "Liquid (Very High P)"
    else:
        if t_val < 0:
            return "Ice VI"
        else:
            return "Liquid (Extreme P)"

def main():
    # Set paths
    edr_file = "../md.edr"  # Path to your EDR file
    tpr_file = "../md.tpr"  # Path to your TPR file
    
    try:
        # Read the energy data
        print("Reading energy data from GROMACS files...")
        data = run_gmx_energy(edr_file, tpr_file)
        
        # Filter data to remove unrealistic points and outliers
        print("\nFiltering data to remove unrealistic points and outliers...")
        data_filtered = filter_data(data, p_threshold=-20, z_threshold=3.0)
        
        # Estimate phase boundaries from simulation data
        print("\nAnalyzing phase behavior...")
        # Increased bins and sigma for better boundary detection
        phase_data = estimate_phase_boundaries(data_filtered, bins=100, sigma=1.5)
        
        # Create simple phase diagram
        print("\nCreating phase diagram...")
        plot_phase_diagram(data_filtered, phase_data, ref_data=WATER_PHASE_DATA)
        
        # Create detailed comparison with reference data
        print("\nCreating detailed comparison with reference phase diagram...")
        plot_detailed_comparison(data_filtered, phase_data, WATER_PHASE_DATA)
        
        # Also create a version with the original unfiltered data for comparison
        print("\nCreating phase diagram with unfiltered data for comparison...")
        phase_data_unfiltered = estimate_phase_boundaries(data, bins=80, sigma=1.0)
        plot_phase_diagram(data, phase_data_unfiltered, ref_data=WATER_PHASE_DATA, 
                         output_file='phase_diagram_unfiltered.png')
        
        print("\nDone!")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 