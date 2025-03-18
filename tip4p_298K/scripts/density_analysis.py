#!/usr/bin/env python3
"""
Density Analysis for TIP4P Water

This script calculates various density profiles and maps from GROMACS trajectory files.
It generates: 
- Density histograms
- 1D density profiles 
- 2D density maps
- Radial density maps

Usage:
    python density_analysis.py [tpr_file] [trajectory_file]

Default:
    Uses md.tpr and md.xtc in the current directory
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from utils import load_universe, save_plot
import MDAnalysis as mda
# Only import density if needed, with error handling
try:
    from MDAnalysis.analysis import density
except ImportError:
    print("Warning: density module not available in MDAnalysis")
    density = None

def get_water_selection(universe):
    """
    Get a selection of water atoms from the universe.
    Tries different common water residue names and atom names.
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        The universe containing the water molecules
        
    Returns:
    --------
    MDAnalysis.AtomGroup, int
        Selection of water atoms and atoms per water molecule (3 for TIP3P, 4 for TIP4P with virtual site)
    """
    # Try different common water residue names
    try:
        water = universe.select_atoms("resname SOL")
        if len(water) == 0:
            water = universe.select_atoms("resname TIP4 or resname TIP4P or resname WAT")
        if len(water) == 0:
            water = universe.select_atoms("name OW or name O")  # Try selecting oxygens
    except:
        print("Warning: Could not select water molecules by residue name. Trying atom names.")
        try:
            water = universe.select_atoms("name OW or name O")  # Try selecting oxygens
        except:
            print("Error: Could not select water atoms. Using all atoms as fallback.")
            water = universe.atoms  # Fallback to all atoms
    
    # Try to determine atoms per water molecule
    # Check if this is TIP4P (4 atoms) or TIP3P (3 atoms)
    atoms_per_molecule = 3  # Default to 3 (TIP3P)
    
    # Check for virtual sites (M or EP sites in TIP4P)
    virtual_sites = universe.select_atoms("name M or name EP or name LP")
    if len(virtual_sites) > 0:
        print(f"Detected virtual sites ({len(virtual_sites)} atoms). Assuming TIP4P model with 4 atoms per molecule.")
        atoms_per_molecule = 4
    else:
        print("No virtual sites detected. Assuming 3 atoms per water molecule.")
    
    return water, atoms_per_molecule

def calculate_density_histogram(universe, n_frames=None, phase_frames=None, phase_names=None, 
                               atoms_per_molecule=None, water_molar_mass=18.01528, exp_density=0.9998, temperature=273.0):
    """
    Calculate the density histogram of the system over time
    """
    # Get water selection and determine atoms per molecule
    water, detected_atoms_per_molecule = get_water_selection(universe)
    
    # Use provided atoms_per_molecule if specified, otherwise use detected value
    if atoms_per_molecule is None:
        atoms_per_molecule = detected_atoms_per_molecule
    
    # Calculate number of water molecules
    n_water_molecules = len(water) / atoms_per_molecule
    
    print(f"Found {len(water)} water atoms")
    print(f"Using {atoms_per_molecule} atoms per water molecule")
    print(f"Number of water molecules: {n_water_molecules}")
    print(f"Using water molar mass of {water_molar_mass} g/mol")
    print(f"Using experimental density of {exp_density} g/cm³ at {temperature}K")
    
    # Calculate system volume and densities for each frame
    densities = []
    volumes = []
    times = []
    
    # Iterate through trajectory
    if n_frames is None:
        n_frames = len(universe.trajectory)
    
    # Constants for density calculation
    avogadro = 6.02214076e23  # Avogadro's number
    nm3_to_cm3 = 1e-21        # Conversion from nm³ to cm³
    amu_to_g_per_cm3 = 1.660539067  # Conversion factor for amu/nm³ to g/cm³
    
    print("\nStarting density calculations...")
    print(f"Processing {n_frames} frames...")
    
    # Debug: Print initial box dimensions
    ts = universe.trajectory[0]
    print(f"\nInitial box dimensions:")
    print(f"Box vectors (nm):")  # Note: GROMACS uses nm
    print(f"a: {ts.dimensions[0]:.3f} nm")
    print(f"b: {ts.dimensions[1]:.3f} nm")
    print(f"c: {ts.dimensions[2]:.3f} nm")
    print(f"α: {ts.dimensions[3]:.2f}°")
    print(f"β: {ts.dimensions[4]:.2f}°")
    print(f"γ: {ts.dimensions[5]:.2f}°")
    
    # Store initial volume for NPT check
    initial_volume = ts.volume  # Volume in nm³
    volume_changes = []
    
    for frame_idx, ts in enumerate(universe.trajectory[:n_frames]):
        # Calculate box volume directly from dimensions for each frame
        box_length_a = ts.dimensions[0]  # nm
        box_length_b = ts.dimensions[1]  # nm
        box_length_c = ts.dimensions[2]  # nm
        volume = box_length_a * box_length_b * box_length_c  # nm³
        
        # Calculate density using direct conversion
        # Method 1: Using the conversion factor
        density = (water.masses.sum() / volume) * amu_to_g_per_cm3
        
        # Alternative method (commented out):
        # Method 2: Using Avogadro's number
        # mass_g = water.masses.sum() / avogadro  # Convert amu to g
        # volume_cm3 = volume * nm3_to_cm3  # Convert nm³ to cm³
        # density = mass_g / volume_cm3  # g/cm³
        
        # Store volume in m³ for output
        volumes.append(volume * 1e-27)  # Convert nm³ to m³
        densities.append(density)
        
        # Store time in ps
        times.append(ts.time)
        
        # Calculate volume change for NPT check
        volume_changes.append(abs(volume - initial_volume) / initial_volume)
        
        # Print progress and values every 50 frames
        if frame_idx % 50 == 0:
            print(f"Frame {frame_idx}/{n_frames}:")
            print(f"  Box dimensions: {box_length_a:.3f} × {box_length_b:.3f} × {box_length_c:.3f} nm³")
            print(f"  Volume: {volume:.2f} nm³")
            print(f"  Mass: {water.masses.sum():.6e} amu")
            print(f"  Density: {density:.6f} g/cm³")
    
    # Calculate statistics
    mean_density = np.mean(densities)
    std_density = np.std(densities)
    min_density = np.min(densities)
    max_density = np.max(densities)
    
    # Check if this is likely an NPT simulation
    max_volume_change = max(volume_changes) * 100  # Convert to percentage
    
    print("\nDensity calculation complete!")
    print(f"Mean density: {mean_density:.6f} g/cm³")
    print(f"Standard deviation: {std_density:.6f} g/cm³")
    print(f"Range: {min_density:.6f} - {max_density:.6f} g/cm³")
    print(f"Maximum volume change: {max_volume_change:.2f}%")
    
    if max_volume_change < 0.01:
        print("\nWARNING: Box volume shows no significant changes.")
        print("This might be an NVT simulation or there might be an issue with the trajectory.")
        print("For NPT simulations, expect some volume fluctuations.")
    
    # Ensure analysis directory exists
    analysis_dir = '../analysis'
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    
    # Save data with header containing calculation parameters
    header = (
        f"# Density calculation parameters:\n"
        f"# Water atoms: {len(water)}\n"
        f"# Atoms per molecule: {atoms_per_molecule}\n"
        f"# Water molecules: {n_water_molecules}\n"
        f"# Water molar mass: {water_molar_mass} g/mol\n"
        f"# Frames analyzed: {len(densities)}\n"
        f"# Mean density: {mean_density:.6f} g/cm³\n"
        f"# Standard deviation: {std_density:.6f} g/cm³\n"
        f"#\n"
        f"# Frame\tTime (ps)\tDensity (g/cm³)\tVolume (m³)"
    )
    
    np.savetxt(os.path.join(analysis_dir, 'density_vs_time.dat'),
              np.column_stack((np.arange(len(densities)), times, densities, volumes)),
              header=header,
              comments='')
    
    # Calculate percent difference from experimental
    percent_diff = ((mean_density - exp_density) / exp_density) * 100
    
    # Create figure for density vs time plot
    plt.figure(figsize=(12, 7))
    
    # Plot the actual density data with higher transparency
    plt.plot(times, densities, 'b-', linewidth=0.5, label='Instantaneous', alpha=0.3)
    
    # Calculate and plot running average (over 50 frames) for smoother visualization
    window = min(50, len(densities) // 10)  # Use 10% of total frames or 50, whichever is smaller
    if window > 1:
        running_avg = np.convolve(densities, np.ones(window)/window, mode='valid')
        running_avg_times = times[window-1:]
        plt.plot(running_avg_times, running_avg, 'b-', linewidth=2.0, 
                label=f'{window}-frame Average', alpha=1.0)
    
    # Add mean line
    plt.axhline(mean_density, color='red', linestyle='dashed', 
                linewidth=2, label=f'Mean: {mean_density:.4f} g/cm³')
    
    # Add standard deviation bands
    plt.axhline(mean_density + std_density, color='red', linestyle='dotted', linewidth=1)
    plt.axhline(mean_density - std_density, color='red', linestyle='dotted', linewidth=1)
    plt.fill_between(times, mean_density - std_density, mean_density + std_density, 
                    color='red', alpha=0.1, label=f'±1σ: {std_density:.4f} g/cm³')
    
    # Add experimental density reference
    plt.axhline(exp_density, color='green', linestyle='dotted', 
                linewidth=2, label=f'Exp. at {temperature}K: {exp_density:.4f} g/cm³')
    
    # Add TIP4P literature value reference
    tip4p_lit_density = 1.002  # g/cm^3 at 298K
    plt.axhline(tip4p_lit_density, color='purple', linestyle='-.', 
                linewidth=1.5, label=f'TIP4P Lit: ~{tip4p_lit_density:.3f} g/cm³')
    
    # Add text box with statistics
    stats_text = (
        f"Mean Density: {mean_density:.4f} g/cm³\n"
        f"Standard Dev: {std_density:.4f} g/cm³\n"
        f"Diff from Exp: {percent_diff:.1f}%"
    )
    
    plt.text(0.97, 0.97, stats_text, transform=plt.gca().transAxes,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(facecolor='white', edgecolor='gray', alpha=0.9, 
                     boxstyle='round,pad=0.4', linewidth=0.5),
            fontsize=10)
    
    # Set axis limits with padding
    density_range = max_density - min_density
    plt.ylim(min_density - density_range*0.1, max_density + density_range*0.1)
    plt.xlim(min(times), max(times))
    
    # Improve grid and ticks
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.minorticks_on()
    
    # Set labels with proper units
    plt.xlabel('Time (ps)', fontsize=12, fontweight='bold')
    plt.ylabel('Density (g/cm³)', fontsize=12, fontweight='bold')
    plt.title('TIP4P Water Density vs. Time', fontsize=14, fontweight='bold')
    
    # Optimize legend position and style
    plt.legend(loc='upper left', framealpha=0.9, edgecolor='gray', 
              fontsize=9, ncol=1, frameon=True, title="Legend")
    
    # Add simulation conditions note
    plt.figtext(0.5, 0.01, 
               f"TIP4P water model simulation at ~{temperature}K and 1 bar | {len(densities)} frames analyzed",
               ha='center', fontsize=10, style='italic')
    
    # Improve overall style
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    
    # Save the plot
    save_plot(plt, 'density_vs_time.png')
    
    return mean_density

def calculate_1d_density_profile(universe, axis='z', bins=100, n_frames=None):
    """
    Calculate 1D density profile along a specified axis
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory
    axis : str
        Axis along which to calculate the profile ('x', 'y', or 'z')
    bins : int
        Number of bins for histogram
    n_frames : int or None
        Number of frames to analyze (None = all frames)
        
    Returns:
    --------
    tuple
        (bin_centers, hist) - the x and y values for the density profile
    """
    water_oxygen = universe.select_atoms("name OW")
    
    axis_num = {'x': 0, 'y': 1, 'z': 2}[axis]
    
    if n_frames is None:
        n_frames = len(universe.trajectory)
    
    # Initialize arrays
    all_positions = np.zeros((n_frames * len(water_oxygen), 3))
    
    # Collect positions from each frame
    for i, ts in enumerate(universe.trajectory[:n_frames]):
        start_idx = i * len(water_oxygen)
        end_idx = start_idx + len(water_oxygen)
        all_positions[start_idx:end_idx] = water_oxygen.positions
    
    # Get the values along the selected axis
    axis_values = all_positions[:, axis_num]
    
    # Calculate histogram
    hist, bin_edges = np.histogram(axis_values, bins=bins, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Ensure analysis directory exists
    analysis_dir = '../analysis'
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    
    # Save data
    np.savetxt(os.path.join(analysis_dir, f'density_profile_{axis}.dat'), 
               np.column_stack((bin_centers, hist)),
               header=f'{axis} (Å)\tDensity (normalized)', 
               comments='# ')
    
    # Plot clean density profile
    plt.figure(figsize=(10, 6))
    plt.plot(bin_centers, hist, 'b-', linewidth=1.5)
    
    # Add simple annotation
    plt.text(0.02, 0.95, 
             f"Oscillations show water layering structure",
             transform=plt.gca().transAxes, 
             verticalalignment='top', horizontalalignment='left',
             bbox=dict(facecolor='white', alpha=0.7))
    
    plt.xlabel(f'{axis} coordinate (Å)')
    plt.ylabel('Density (normalized)')
    plt.title(f'TIP4P Water Density Profile along {axis}-axis')
    plt.grid(True, alpha=0.3)
    save_plot(plt, f'density_profile_{axis}_plot.png')
    
    return bin_centers, hist

def create_combined_density_profiles(profiles_data, output_dir='../analysis'):
    """
    Create a combined plot of density profiles along x, y, and z axes
    
    Parameters:
    -----------
    profiles_data : dict
        Dictionary containing the profile data for each axis
        Format: {'x': (bin_centers_x, hist_x), 'y': (bin_centers_y, hist_y), 'z': (bin_centers_z, hist_z)}
    output_dir : str
        Directory to save the output plot
    """
    # Create figure with 3 subplots (vertically stacked)
    fig, axes = plt.subplots(3, 1, figsize=(12, 12), sharex=True)
    
    # Plot each profile in its own subplot
    for i, axis in enumerate(['x', 'y', 'z']):
        if axis in profiles_data:
            bin_centers, hist = profiles_data[axis]
            axes[i].plot(bin_centers, hist, 'b-', linewidth=1.5)
            
            # Add annotation
            axes[i].text(0.02, 0.95, 
                        f"Oscillations show water layering structure",
                        transform=axes[i].transAxes, 
                        verticalalignment='top', horizontalalignment='left',
                        bbox=dict(facecolor='white', alpha=0.7))
            
            axes[i].set_ylabel('Density (normalized)')
            axes[i].set_title(f'Density Profile along {axis}-axis')
            axes[i].grid(True, alpha=0.3)
            
            # Set y-axis limits to be consistent across all plots
            y_min = 0
            y_max = max(hist) * 1.05  # Add 5% margin
            axes[i].set_ylim(y_min, y_max)
    
    # Set common x-axis label
    axes[2].set_xlabel('Coordinate (Å)')
    
    # Add overall title
    plt.suptitle('TIP4P Water Density Profiles', fontsize=16)
    
    # Adjust spacing between subplots
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)  # Make room for the suptitle
    
    # Save the combined plot
    plt.savefig(os.path.join(output_dir, 'combined_density_profiles.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Combined density profiles plot saved to", os.path.join(output_dir, 'combined_density_profiles.png'))

def calculate_2d_density_map(universe, plane='xy', nbins=100, n_frames=None):
    """
    Calculate 2D density map on a specified plane
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory
    plane : str
        Plane for the density map ('xy', 'yz', or 'xz')
    nbins : int
        Number of bins for each dimension (increased for better resolution)
    n_frames : int or None
        Number of frames to analyze (None = all frames)
        
    Returns:
    --------
    tuple
        (X, Y, hist) - the meshgrid X, Y and the 2D histogram values
    """
    water_oxygen = universe.select_atoms("name OW")
    
    axis_map = {'xy': (0, 1), 'yz': (1, 2), 'xz': (0, 2)}
    axis1, axis2 = axis_map[plane]
    axis_labels = ['x', 'y', 'z']
    
    if n_frames is None:
        n_frames = len(universe.trajectory)
    
    # Initialize arrays for positions
    all_positions = np.zeros((n_frames * len(water_oxygen), 3))
    
    # Collect positions from each frame
    for i, ts in enumerate(universe.trajectory[:n_frames]):
        start_idx = i * len(water_oxygen)
        end_idx = start_idx + len(water_oxygen)
        all_positions[start_idx:end_idx] = water_oxygen.positions
    
    # Extract coordinates for the specified plane
    x = all_positions[:, axis1]
    y = all_positions[:, axis2]
    
    # Calculate 2D histogram
    hist, xedges, yedges = np.histogram2d(x, y, bins=nbins, density=True)
    
    # Apply Gaussian smoothing for better visualization (optional)
    from scipy.ndimage import gaussian_filter
    hist_smooth = gaussian_filter(hist, sigma=1.0)
    
    # Create a meshgrid for plotting
    X, Y = np.meshgrid((xedges[:-1] + xedges[1:]) / 2, (yedges[:-1] + yedges[1:]) / 2)
    
    # Ensure analysis directory exists
    analysis_dir = '../analysis'
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    
    # Save data
    np.savetxt(os.path.join(analysis_dir, f'density_map_{plane}.dat'), hist_smooth,
              header=f'2D Density map on {plane} plane', 
              comments='# ')
    
    # Plot improved 2D density map
    plt.figure(figsize=(10, 8))
    
    # Use a better colormap with improved contrast
    cmap = plt.cm.viridis
    
    # Plot density map with improved contrast
    plt.pcolormesh(X, Y, hist_smooth.T, cmap=cmap, shading='auto')
    cbar = plt.colorbar(label='Density (normalized)')
    
    # Adjust colorbar font size
    cbar.ax.tick_params(labelsize=10)
    
    # Add annotation about density distribution
    plt.text(0.02, 0.98, 
             f"Homogeneous density distribution\nindicates proper equilibration",
             transform=plt.gca().transAxes, 
             verticalalignment='top', horizontalalignment='left',
             bbox=dict(facecolor='white', alpha=0.7))
    
    plt.xlabel(f'{axis_labels[axis1]} coordinate (Å)')
    plt.ylabel(f'{axis_labels[axis2]} coordinate (Å)')
    plt.title(f'TIP4P Water Density Map on {plane}-plane')
    plt.axis('equal')
    save_plot(plt, f'density_map_{plane}.png')
    
    return X, Y, hist_smooth

def create_combined_density_maps(density_maps_data, output_dir='../analysis'):
    """
    Create a combined plot of density maps on xy, xz, and yz planes
    
    Parameters:
    -----------
    density_maps_data : dict
        Dictionary containing the density map data for each plane
        Format: {'xy': (X, Y, hist_xy), 'xz': (X, Y, hist_xz), 'yz': (X, Y, hist_yz)}
    output_dir : str
        Directory to save the output plot
    """
    # Create figure with 2x2 grid (one empty subplot)
    fig = plt.figure(figsize=(15, 12))
    
    # Define subplot positions
    positions = {
        'xy': 1,  # top-left
        'xz': 2,  # top-right
        'yz': 3   # bottom-left
    }
    
    # Common colormap for all plots
    cmap = plt.cm.viridis
    
    # Find global min and max for consistent color scaling
    all_values = []
    for plane, (_, _, hist) in density_maps_data.items():
        all_values.extend(hist.flatten())
    
    vmin = np.min(all_values)
    vmax = np.max(all_values)
    
    # Plot each density map
    for plane, (X, Y, hist) in density_maps_data.items():
        axis_map = {'xy': (0, 1), 'yz': (1, 2), 'xz': (0, 2)}
        axis1, axis2 = axis_map[plane]
        axis_labels = ['x', 'y', 'z']
        
        # Create subplot
        ax = fig.add_subplot(2, 2, positions[plane])
        
        # Plot density map
        im = ax.pcolormesh(X, Y, hist.T, cmap=cmap, shading='auto', vmin=vmin, vmax=vmax)
        
        # Add labels and title
        ax.set_xlabel(f'{axis_labels[axis1]} coordinate (Å)')
        ax.set_ylabel(f'{axis_labels[axis2]} coordinate (Å)')
        ax.set_title(f'Density Map on {plane}-plane')
        
        # Add grid
        ax.grid(False)
        
        # Set aspect ratio to equal
        ax.set_aspect('equal')
        
        # Add annotation
        ax.text(0.02, 0.98, 
                f"Homogeneous distribution\nindicates proper equilibration",
                transform=ax.transAxes, 
                verticalalignment='top', horizontalalignment='left',
                fontsize=9,
                bbox=dict(facecolor='white', alpha=0.7))

    # Add a colorbar in the empty subplot position
    cbar_ax = fig.add_subplot(2, 2, 4)
    cbar_ax.axis('off')
    cbar = fig.colorbar(im, ax=cbar_ax, orientation='vertical', fraction=0.8)
    cbar.set_label('Density (normalized)')
    
    # Add overall title
    plt.suptitle('TIP4P Water Density Maps', fontsize=16)
    
    # Adjust spacing
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)  # Make room for the suptitle
    
    # Save the combined plot
    plt.savefig(os.path.join(output_dir, 'combined_density_maps.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Combined density maps plot saved to", os.path.join(output_dir, 'combined_density_maps.png'))

def calculate_radial_density(universe, center_selection="name OW", nbins=50, max_radius=10.0, n_frames=None):
    """
    Calculate radial density around selected atoms
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory
    center_selection : str
        Selection string for center atoms
    nbins : int
        Number of radial bins
    max_radius : float
        Maximum radius for the density calculation (in Å)
    n_frames : int or None
        Number of frames to analyze (None = all frames)
    """
    from scipy.spatial.distance import cdist
    from scipy.signal import find_peaks
    import scipy.integrate as integrate
    
    # Select atoms
    center_atoms = universe.select_atoms(center_selection)
    water_oxygen = universe.select_atoms("name OW")
    
    if n_frames is None:
        n_frames = len(universe.trajectory)
    
    # Initialize arrays
    radial_density = np.zeros(nbins)
    bin_edges = np.linspace(0, max_radius, nbins+1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_volumes = (4/3) * np.pi * (bin_edges[1:]**3 - bin_edges[:-1]**3)
    
    # Create mask to exclude very small distances (r < 1.0 Å)
    min_radius = 1.0  # Minimum distance to consider
    valid_mask = bin_centers >= min_radius
    
    # Take a subset of center atoms to speed up calculation
    if len(center_atoms) > 10:
        center_atoms = center_atoms[:10]
    
    frame_count = 0
    
    # Loop through trajectory
    for ts in universe.trajectory[:n_frames]:
        if frame_count % 10 == 0:
            print(f"Processing frame {frame_count}/{n_frames}")
        
        all_hist = np.zeros(nbins)
        
        # Loop through center atoms
        for center_pos in center_atoms.positions:
            # Calculate distances to all water oxygens
            distances = cdist([center_pos], water_oxygen.positions)[0]
            
            # Create histogram
            hist, _ = np.histogram(distances, bins=bin_edges)
            all_hist += hist
        
        # Normalize by bin volume to get density
        all_hist = all_hist / bin_volumes
        
        # Accumulate
        radial_density += all_hist
        frame_count += 1
    
    # Normalize by number of frames and center atoms
    radial_density = radial_density / (frame_count * len(center_atoms))
    
    # Ensure analysis directory exists
    analysis_dir = '../analysis'
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    
    # Save data
    np.savetxt(os.path.join(analysis_dir, 'radial_density.dat'), 
               np.column_stack((bin_centers, radial_density)),
               header='r (Å)\tDensity (normalized)', 
               comments='# ')
    
    # Find peaks in the radial density to identify solvation shells
    # Only consider data beyond min_radius
    valid_data = radial_density[valid_mask]
    valid_r = bin_centers[valid_mask]
    
    # Find peaks with prominence at least 10% of the maximum density
    peak_prominence = 0.1 * np.max(valid_data)
    peaks, peak_properties = find_peaks(valid_data, prominence=peak_prominence)
    
    # Get peak positions in Å
    peak_positions = valid_r[peaks]
    
    # Calculate coordination numbers (number of water molecules in each shell)
    # by integrating 4πr²ρ(r) from 0 to r
    coordination_numbers = []
    shell_boundaries = []
    
    if len(peak_positions) >= 1:
        # For the first shell, find the minimum after the first peak
        first_peak_idx = peaks[0]
        # Look for the minimum after the first peak
        min_after_peak = np.argmin(valid_data[first_peak_idx:first_peak_idx+10]) + first_peak_idx
        first_shell_boundary = valid_r[min_after_peak]
        shell_boundaries.append(first_shell_boundary)
        
        # Calculate coordination number for first shell
        # Integrate 4πr²ρ(r) from 0 to first_shell_boundary
        r_values = bin_centers
        integrand = 4 * np.pi * r_values**2 * radial_density
        # Use trapezoidal rule to integrate up to the shell boundary
        boundary_idx = np.searchsorted(r_values, first_shell_boundary)
        cn_first_shell = np.trapz(integrand[:boundary_idx], r_values[:boundary_idx])
        coordination_numbers.append(cn_first_shell)
    
    if len(peak_positions) >= 2:
        # For the second shell, find the minimum after the second peak
        second_peak_idx = peaks[1]
        # Convert to index in valid_data
        valid_second_peak_idx = second_peak_idx
        # Look for the minimum after the second peak
        if valid_second_peak_idx + 10 < len(valid_data):
            min_after_second_peak = np.argmin(valid_data[valid_second_peak_idx:valid_second_peak_idx+10]) + valid_second_peak_idx
            second_shell_boundary = valid_r[min_after_second_peak]
            shell_boundaries.append(second_shell_boundary)
            
            # Calculate coordination number for second shell
            boundary_idx = np.searchsorted(r_values, second_shell_boundary)
            cn_second_shell = np.trapz(integrand[:boundary_idx], r_values[:boundary_idx]) - cn_first_shell
            coordination_numbers.append(cn_second_shell)
    
    # Create an enhanced plot with more information
    plt.figure(figsize=(12, 8))
    
    # Plot the radial density with improved styling
    plt.plot(bin_centers[valid_mask], radial_density[valid_mask], 'b-', linewidth=2.5)
    
    # Fill area under the curve with light blue
    plt.fill_between(bin_centers[valid_mask], 0, radial_density[valid_mask], color='blue', alpha=0.1)
    
    # Mark the peaks and add annotations
    if len(peak_positions) > 0:
        for i, peak_pos in enumerate(peak_positions):
            peak_idx = np.where(bin_centers >= peak_pos)[0][0]
            peak_height = radial_density[peak_idx]
            
            # Mark peak with a point
            plt.plot(peak_pos, peak_height, 'ro', markersize=8)
            
            # Add annotation for each peak
            if i == 0:
                plt.annotate(f'First Shell\n({peak_pos:.2f} Å)',
                            xy=(peak_pos, peak_height),
                            xytext=(peak_pos + 0.5, peak_height + 0.01),
                            arrowprops=dict(facecolor='black', shrink=0.05, width=1.5),
                            fontsize=10,
                            bbox=dict(facecolor='white', alpha=0.8))
            elif i == 1:
                plt.annotate(f'Second Shell\n({peak_pos:.2f} Å)',
                            xy=(peak_pos, peak_height),
                            xytext=(peak_pos + 0.5, peak_height + 0.01),
                            arrowprops=dict(facecolor='black', shrink=0.05, width=1.5),
                            fontsize=10,
                            bbox=dict(facecolor='white', alpha=0.8))
    
    # Mark shell boundaries if found
    for i, boundary in enumerate(shell_boundaries):
        boundary_idx = np.where(bin_centers >= boundary)[0][0]
        boundary_height = radial_density[boundary_idx]
        
        # Draw vertical line at shell boundary
        plt.axvline(boundary, color='gray', linestyle='--', alpha=0.7)
        
        # Add annotation for shell boundary
        plt.annotate(f'Shell {i+1} Boundary\n({boundary:.2f} Å)',
                    xy=(boundary, boundary_height / 2),
                    xytext=(boundary + 0.3, boundary_height / 2),
                    arrowprops=dict(facecolor='gray', shrink=0.05, width=1, alpha=0.7),
                    fontsize=9,
                    color='gray')
    
    # Add information box with coordination numbers
    info_text = "Water Solvation Structure:\n"
    if len(coordination_numbers) >= 1:
        info_text += f"First Shell Coordination: {coordination_numbers[0]:.2f} molecules\n"
    if len(coordination_numbers) >= 2:
        info_text += f"Second Shell Coordination: {coordination_numbers[1]:.2f} molecules\n"
    
    # Add bulk density reference
    bulk_idx = np.where(bin_centers >= 8.0)[0][0]  # Use r > 8Å as bulk
    bulk_density = np.mean(radial_density[bulk_idx:])
    info_text += f"Bulk Density Ratio: {bulk_density:.2f}"
    
    # Position the text box in the upper right
    plt.text(0.97, 0.97, info_text, transform=plt.gca().transAxes,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(facecolor='white', edgecolor='gray', alpha=0.9, boxstyle='round,pad=0.5'),
            fontsize=10)
    
    # Add explanation of what the plot shows
    plt.figtext(0.5, 0.01, 
               "Radial density distribution shows how water molecules organize in shells around reference atoms.\n"
               "Peaks represent preferred distances for water molecules, indicating solvation structure.",
               ha='center', fontsize=10, style='italic')
    
    # Improve axis labels and title
    plt.xlabel('Distance from reference atom (Å)', fontsize=12, fontweight='bold')
    plt.ylabel('Density (normalized)', fontsize=12, fontweight='bold')
    plt.title('TIP4P Water Radial Density Distribution', fontsize=14, fontweight='bold')
    
    # Improve grid and ticks
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.minorticks_on()
    
    # Set axis limits to focus on the relevant range
    plt.xlim(min_radius, max_radius)
    plt.ylim(0, np.max(radial_density[valid_mask]) * 1.1)
    
    # Improve overall style
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])  # Adjust layout to make room for the note
    
    # Save the enhanced plot
    save_plot(plt, 'radial_density_map.png')

def calculate_density_histogram_multi(universes, phase_frames=None, phase_names=None, 
                                atoms_per_molecule=None, water_molar_mass=18.01528,
                                exp_density=0.9998, temperature=273.0):
    """
    Calculate density histogram for multiple trajectories
    
    Parameters:
    -----------
    universes : list of (MDAnalysis.Universe, int, int)
        List of (universe, start_frame, end_frame) tuples for each trajectory
    phase_frames : list of int
        Number of frames in each phase
    phase_names : list of str
        Names of each phase
    atoms_per_molecule : int or None
        Number of atoms per water molecule (None = auto-detect)
    water_molar_mass : float
        Molar mass of water in g/mol (default: 18.01528)
    exp_density : float
        Experimental density of water in g/cm³ (default: 0.9998 at 273K)
    temperature : float
        Simulation temperature in K (default: 273.0)
        
    Returns:
    --------
    float
        Average density across all trajectories
    """
    # Initialize arrays to store all densities and times
    all_densities = []
    all_times = []
    
    # Track cumulative time
    cumulative_time = 0.0
    
    # Define phase durations in ns (adjust these values as needed)
    # These are approximate durations for each phase
    phase_durations = [0.1, 3.0, 3.0, 5.0]  # EM, NVT, NPT, MD in ns
    
    # Process each universe
    for i, (universe, start_frame, end_frame) in enumerate(universes):
        phase_name = phase_names[i] if phase_names and i < len(phase_names) else f"Phase {i+1}"
        print(f"Processing {phase_name} trajectory...")
        
        # Get water selection
        water_selection, detected_atoms_per_molecule = get_water_selection(universe)
        n_waters = len(water_selection) // detected_atoms_per_molecule
        print(f"  Found {len(water_selection)} water atoms")
        
        # Use provided atoms_per_molecule if specified, otherwise use detected value
        if atoms_per_molecule is None:
            atoms_per_molecule = detected_atoms_per_molecule
        
        print(f"  Using {atoms_per_molecule} atoms per water molecule")
        print(f"  Using water molar mass of {water_molar_mass} g/mol")
        print(f"  Using experimental density of {exp_density} g/cm³ at {temperature}K")
        
        # Calculate density for each frame
        densities = []
        phase_times = []
        
        # Get number of frames in this phase
        if start_frame is None:
            start_frame = 0
        if end_frame is None:
            end_frame = len(universe.trajectory)
        n_frames = end_frame - start_frame
        
        # Get phase duration (use predefined or calculate based on frames)
        if i < len(phase_durations):
            phase_duration = phase_durations[i]
        else:
            phase_duration = 1.0  # Default 1 ns if not specified
        
        # Iterate through trajectory frames
        for frame_idx, ts in enumerate(universe.trajectory[start_frame:end_frame]):
            # Calculate density using MDAnalysis masses
            # Convert amu to g (1 amu = 1 g/mol)
            # For N molecules, we need to divide by Avogadro's number to get mass in g
            mass_g = water_selection.masses.sum() / 6.02214076e23  # Convert amu to g
            volume_cm3 = ts.volume * 1e-21  # Convert nm^3 to cm^3
            density = mass_g / volume_cm3  # g/cm^3
            
            densities.append(density)
            
            # Calculate relative time within this phase (in ns)
            if n_frames > 1:
                relative_time = (frame_idx / (n_frames - 1)) * phase_duration
            else:
                relative_time = 0
            
            # Add to cumulative time (in ns)
            current_time = cumulative_time + relative_time
            phase_times.append(current_time)
        
        # Update cumulative time for next phase
        cumulative_time += phase_duration
        
        # Add to overall arrays
        all_densities.extend(densities)
        all_times.extend(phase_times)
        
        # Print average density for this phase
        avg_density = np.mean(densities)
        print(f"  Added {len(densities)} frames with average density: {avg_density:.4f} g/cm³")
    
    # Calculate overall statistics
    mean_density = np.mean(all_densities)
    std_density = np.std(all_densities)
    
    # Create a new figure with a clean, simple design
    plt.figure(figsize=(12, 7))
    
    # Print time range for debugging
    print(f"Time range: {min(all_times):.2f} ns to {max(all_times):.2f} ns")
    
    # Plot a single continuous line for density
    plt.plot(all_times, all_densities, color='blue', linewidth=1.5, label='Density')
    
    # Add mean line
    plt.axhline(mean_density, color='red', linestyle='dashed', 
                linewidth=2, label=f'Mean: {mean_density:.4f} g/cm³')
    
    # Add standard deviation bands
    plt.axhline(mean_density + std_density, color='red', linestyle='dotted', linewidth=1)
    plt.axhline(mean_density - std_density, color='red', linestyle='dotted', linewidth=1)
    plt.fill_between(all_times, mean_density - std_density, mean_density + std_density, 
                    color='red', alpha=0.1, label=f'±1σ: {std_density:.4f} g/cm³')
    
    # Add experimental density reference
    plt.axhline(exp_density, color='green', linestyle='dotted', 
                linewidth=2, label=f'Exp. at {temperature}K: {exp_density:.4f} g/cm³')
    
    # Add TIP4P literature value reference
    tip4p_lit_density = 1.002  # g/cm^3 at 298K
    plt.axhline(tip4p_lit_density, color='purple', linestyle='-.', 
                linewidth=1.5, label=f'TIP4P Lit: ~{tip4p_lit_density:.3f} g/cm³')
    
    # Calculate percent difference from experimental value
    percent_diff = abs((mean_density - exp_density) / exp_density) * 100
    
    # Add text box with statistics
    stats_text = (
        f"Mean Density: {mean_density:.4f} g/cm³\n"
        f"Standard Dev: {std_density:.4f} g/cm³\n"
        f"Diff from Exp: {percent_diff:.1f}%"
    )
    
    # Position the text box in the upper right with a slight offset from the corner
    plt.text(0.97, 0.97, stats_text, transform=plt.gca().transAxes,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(facecolor='white', edgecolor='gray', alpha=0.9, 
                     boxstyle='round,pad=0.4', linewidth=0.5),
            fontsize=10)
    
    # Set y-axis limits with some padding
    # Use a reasonable range that shows the important data
    y_min = max(0, min(exp_density * 0.5, min(all_densities) * 0.9))
    y_max = max(all_densities) * 1.1
    plt.ylim(y_min, y_max)
    
    # Set x-axis limits to show full time range
    plt.xlim(0, max(all_times))
    
    # Improve grid and ticks
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.minorticks_on()
    
    # Set labels
    plt.xlabel('Simulation Time (ns)', fontsize=12, fontweight='bold')
    plt.ylabel('Density (g/cm³)', fontsize=12, fontweight='bold')
    plt.title('TIP4P Water Density vs. Time', fontsize=14, fontweight='bold')
    
    # Optimize legend position and style
    plt.legend(loc='upper left', framealpha=0.9, edgecolor='gray', 
              fontsize=9, ncol=1, frameon=True, title="Legend")
    
    # Add a note about the simulation conditions
    plt.figtext(0.5, 0.01, 
               f"TIP4P water model simulation at ~273K and 1 bar | {len(all_densities)} frames analyzed",
               ha='center', fontsize=10, style='italic')
    
    # Improve overall style
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])  # Adjust layout to make room for the note
    
    # Save the plot
    save_plot(plt, 'density_vs_time.png')
    
    return mean_density

def main():
    # Get command line arguments if provided
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze density from GROMACS trajectories')
    parser.add_argument('--tpr', nargs='+', default=['em.tpr', 'nvt.tpr', 'npt.tpr', 'md.tpr'],
                        help='TPR files for each simulation phase')
    parser.add_argument('--trj', nargs='+', default=['em.trr', 'nvt.trr', 'npt.trr', 'md.trr'],
                        help='Trajectory files for each simulation phase')
    parser.add_argument('--phases', nargs='+', default=['EM', 'NVT', 'NPT', 'MD'],
                        help='Names of simulation phases')
    parser.add_argument('--start-frames', nargs='+', type=int, 
                        help='Starting frame for each phase')
    parser.add_argument('--end-frames', nargs='+', type=int,
                        help='Ending frame for each phase')
    parser.add_argument('--single-phase', action='store_true',
                        help='Analyze only a single phase (md.tpr and md.xtc)')
    parser.add_argument('--frames', type=int, default=None,
                        help='Number of frames to analyze (None = all frames)')
    # Add water model parameters
    parser.add_argument('--atoms-per-molecule', type=int, default=None,
                        help='Number of atoms per water molecule (default: auto-detect)')
    parser.add_argument('--water-molar-mass', type=float, default=18.01528,
                        help='Molar mass of water in g/mol (default: 18.01528)')
    parser.add_argument('--exp-density', type=float, default=0.9998,
                        help='Experimental density of water in g/cm³ (default: 0.9998 at 273K)')
    parser.add_argument('--temperature', type=float, default=273.0,
                        help='Simulation temperature in K (default: 273.0)')
    # Add option to specify TIP4P model variant
    parser.add_argument('--water-model', type=str, default=None, choices=['TIP3P', 'TIP4P', 'TIP4P/2005', 'TIP4P/Ew', 'TIP4P/Ice'],
                        help='Water model type (default: auto-detect)')
    
    args = parser.parse_args()
    
    # If water model is specified, set atoms per molecule accordingly
    if args.water_model:
        if args.water_model == 'TIP3P':
            if args.atoms_per_molecule is None:
                args.atoms_per_molecule = 3
            print(f"Using {args.water_model} water model with {args.atoms_per_molecule} atoms per molecule")
        elif args.water_model.startswith('TIP4P'):
            if args.atoms_per_molecule is None:
                args.atoms_per_molecule = 4
            print(f"Using {args.water_model} water model with {args.atoms_per_molecule} atoms per molecule")
    
    # If temperature is 298K, use appropriate experimental density
    if abs(args.temperature - 298.0) < 1.0 and args.exp_density == 0.9998:
        args.exp_density = 0.9970  # g/cm³ at 298K
        print(f"Temperature is close to 298K, using experimental density of {args.exp_density} g/cm³")
    
    # If single-phase flag is set, use the traditional approach with just md files
    if args.single_phase or len(sys.argv) == 3:  # Backward compatibility with old command line format
        if len(sys.argv) == 3:
            tpr_file = sys.argv[1]
            trajectory_file = sys.argv[2]
        else:
            tpr_file = 'md.tpr'
            trajectory_file = 'md.trr'
        
        print(f"Analyzing single phase with {tpr_file} and {trajectory_file}")
        
        # Load the trajectory
        universe = load_universe(tpr_file, trajectory_file)
        
        # Calculate and plot density histogram
        print("Calculating density histogram...")
        avg_density = calculate_density_histogram(
            universe, 
            n_frames=args.frames,
            atoms_per_molecule=args.atoms_per_molecule,
            water_molar_mass=args.water_molar_mass,
            exp_density=args.exp_density,
            temperature=args.temperature
        )
        print(f"Average density: {avg_density:.4f} g/cm³")
        
        # Calculate 1D density profiles along each axis and store the data
        profiles_data = {}
        for axis in ['x', 'y', 'z']:
            print(f"Calculating 1D density profile along {axis}...")
            bin_centers, hist = calculate_1d_density_profile(universe, axis=axis, n_frames=args.frames)
            profiles_data[axis] = (bin_centers, hist)
        
        # Create combined density profiles plot
        print("Creating combined density profiles plot...")
        create_combined_density_profiles(profiles_data)
        
        # Calculate 2D density maps and store the data
        density_maps_data = {}
        for plane in ['xy', 'yz', 'xz']:
            print(f"Calculating 2D density map on {plane} plane...")
            X, Y, hist = calculate_2d_density_map(universe, plane=plane, n_frames=min(50, args.frames if args.frames else 50))
            density_maps_data[plane] = (X, Y, hist)
        
        # Create combined density maps plot
        print("Creating combined density maps plot...")
        create_combined_density_maps(density_maps_data)
        
        # Calculate radial density
        print("Calculating radial density...")
        calculate_radial_density(universe, n_frames=min(50, args.frames if args.frames else 50))
    else:
        # Multi-phase analysis
        # Ensure the number of TPR files, trajectory files, and phase names match
        if len(args.tpr) != len(args.trj) or (args.phases and len(args.tpr) != len(args.phases)):
            print("Error: Number of TPR files, trajectory files, and phase names must match")
            sys.exit(1)
        
        # Convert start and end frames to lists of integers if provided
        start_frames = args.start_frames if args.start_frames else [0] * len(args.tpr)
        end_frames = args.end_frames if args.end_frames else [None] * len(args.tpr)
        
        print(f"Analyzing {len(args.tpr)} phases: {', '.join(args.phases)}")
        
        # Load trajectories
        from utils import load_combined_universe
        universes, frame_counts = load_combined_universe(
            tpr_files=args.tpr,
            trajectory_files=args.trj,
            start_frames=start_frames,
            end_frames=end_frames
        )
        
        # Calculate and plot density histogram with phase information
        print("Calculating density histogram...")
        avg_density = calculate_density_histogram_multi(
            universes=universes,
            phase_frames=frame_counts,
            phase_names=args.phases,
            atoms_per_molecule=args.atoms_per_molecule,
            water_molar_mass=args.water_molar_mass,
            exp_density=args.exp_density,
            temperature=args.temperature
        )
        print(f"Average density: {avg_density:.4f} g/cm³")
        
        # For other analyses, use the last universe (MD phase)
        # This is a simplification - ideally, we would analyze all phases for all metrics
        # but that would require significant refactoring of all analysis functions
        print("Using MD phase for remaining analyses...")
        universe, start_frame, end_frame = universes[-1]
        
        # Calculate 1D density profiles along each axis and store the data
        profiles_data = {}
        for axis in ['x', 'y', 'z']:
            print(f"Calculating 1D density profile along {axis}...")
            bin_centers, hist = calculate_1d_density_profile(universe, axis=axis, n_frames=args.frames)
            profiles_data[axis] = (bin_centers, hist)
        
        # Create combined density profiles plot
        print("Creating combined density profiles plot...")
        create_combined_density_profiles(profiles_data)
        
        # Calculate 2D density maps and store the data
        density_maps_data = {}
        for plane in ['xy', 'yz', 'xz']:
            print(f"Calculating 2D density map on {plane} plane...")
            X, Y, hist = calculate_2d_density_map(universe, plane=plane, n_frames=min(50, args.frames if args.frames else 50))
            density_maps_data[plane] = (X, Y, hist)
        
        # Create combined density maps plot
        print("Creating combined density maps plot...")
        create_combined_density_maps(density_maps_data)
        
        # Calculate radial density
        print("Calculating radial density...")
        calculate_radial_density(universe, n_frames=min(50, args.frames if args.frames else 50))
    
    print("Density analysis complete!")

if __name__ == '__main__':
    main() 