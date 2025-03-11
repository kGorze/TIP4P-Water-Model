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
# Import matplotlib colormap without Axes3D to avoid warnings
from matplotlib import cm
from utils import load_universe, save_plot
import MDAnalysis as mda
# Only import density if needed, with error handling
try:
    from MDAnalysis.analysis import density
except ImportError:
    print("Warning: density module not available in MDAnalysis")
    density = None

def calculate_density_histogram(universe, n_frames=None):
    """
    Calculate the density histogram of the system over time
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory
    n_frames : int or None
        Number of frames to analyze (None = all frames)
    """
    # Check for resname SOL, if not try TIP4 or other common water names
    try:
        water = universe.select_atoms("resname SOL")
        if len(water) == 0:
            water = universe.select_atoms("resname TIP4 or resname TIP4P or resname WAT")
    except:
        print("Warning: Could not select water molecules. Trying a simple selection.")
        water = universe.select_atoms("name OW or name O")  # Try selecting oxygens
    
    print(f"Found {len(water)} water atoms")
    
    # Calculate system volume and densities for each frame
    densities = []
    volumes = []
    
    # Iterate through trajectory
    if n_frames is None:
        n_frames = len(universe.trajectory)
    
    for ts in universe.trajectory[:n_frames]:
        # Calculate box volume in nm^3 (MDAnalysis uses Å)
        volume = universe.dimensions[0] * universe.dimensions[1] * universe.dimensions[2] / 1000.0
        volumes.append(volume)
        
        # Calculate density in g/cm^3
        # For water, each molecule is ~18 g/mol, and 1 mol is 6.022e23 molecules
        # 1 nm^3 = 1e-21 cm^3
        n_water_molecules = len(water) / 3  # Each water has 3 atoms (1 O, 2 H)
        mass_g = n_water_molecules * 18.0 / 6.022e23
        density_g_cm3 = mass_g / (volume * 1e-21)
        densities.append(density_g_cm3)
    
    # Ensure analysis directory exists
    analysis_dir = '../analysis'
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
        print(f"Created analysis directory: {analysis_dir}")
    
    # Save data
    np.savetxt(os.path.join(analysis_dir, 'density_vs_time.dat'), 
               np.column_stack((np.arange(n_frames), densities, volumes)),
               header='Frame\tDensity (g/cm^3)\tVolume (nm^3)', 
               comments='# ')
    
    # Plot histogram
    plt.figure(figsize=(10, 6))
    plt.hist(densities, bins=30, alpha=0.7, color='blue')
    plt.axvline(np.mean(densities), color='red', linestyle='dashed', 
                linewidth=2, label=f'Mean: {np.mean(densities):.4f} g/cm³')
    plt.xlabel('Density (g/cm³)')
    plt.ylabel('Frequency')
    plt.title('TIP4P Water Density Distribution')
    plt.legend()
    plt.grid(True, alpha=0.3)
    save_plot(plt, 'density_histogram.png')
    
    # Plot density vs time
    plt.figure(figsize=(10, 6))
    plt.plot(densities, color='blue')
    plt.axhline(np.mean(densities), color='red', linestyle='dashed', 
                linewidth=2, label=f'Mean: {np.mean(densities):.4f} g/cm³')
    plt.xlabel('Frame')
    plt.ylabel('Density (g/cm³)')
    plt.title('TIP4P Water Density vs. Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    save_plot(plt, 'density_vs_time.png')
    
    return np.mean(densities)

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
    
    # Plot density profile
    plt.figure(figsize=(10, 6))
    plt.plot(bin_centers, hist, 'b-', linewidth=2)
    plt.xlabel(f'{axis} coordinate (Å)')
    plt.ylabel('Density (normalized)')
    plt.title(f'TIP4P Water Density Profile along {axis}-axis')
    plt.grid(True, alpha=0.3)
    save_plot(plt, f'density_profile_{axis}_plot.png')

def calculate_2d_density_map(universe, plane='xy', nbins=50, n_frames=None):
    """
    Calculate 2D density map on a specified plane
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory
    plane : str
        Plane for the density map ('xy', 'yz', or 'xz')
    nbins : int
        Number of bins for each dimension
    n_frames : int or None
        Number of frames to analyze (None = all frames)
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
    
    # Create a meshgrid for plotting
    X, Y = np.meshgrid((xedges[:-1] + xedges[1:]) / 2, (yedges[:-1] + yedges[1:]) / 2)
    
    # Ensure analysis directory exists
    analysis_dir = '../analysis'
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    
    # Save data
    np.savetxt(os.path.join(analysis_dir, f'density_map_{plane}.dat'), hist,
              header=f'2D Density map on {plane} plane', 
              comments='# ')
    
    # Plot 2D density map
    plt.figure(figsize=(10, 8))
    plt.pcolormesh(X, Y, hist.T, cmap=cm.viridis, shading='auto')
    plt.colorbar(label='Density (normalized)')
    plt.xlabel(f'{axis_labels[axis1]} coordinate (Å)')
    plt.ylabel(f'{axis_labels[axis2]} coordinate (Å)')
    plt.title(f'TIP4P Water Density Map on {plane}-plane')
    plt.axis('equal')
    save_plot(plt, f'density_map_{plane}.png')

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
    
    # Plot radial density
    plt.figure(figsize=(10, 6))
    plt.plot(bin_centers, radial_density, 'b-', linewidth=2)
    plt.xlabel('Distance from reference atom (Å)')
    plt.ylabel('Density (normalized)')
    plt.title('Radial Density Distribution')
    plt.grid(True, alpha=0.3)
    save_plot(plt, 'radial_density_map.png')

def main():
    # Get command line arguments if provided
    if len(sys.argv) > 2:
        tpr_file = sys.argv[1]
        trajectory_file = sys.argv[2]
    else:
        tpr_file = 'md.tpr'
        trajectory_file = 'md.xtc'
    
    # Load the trajectory
    universe = load_universe(tpr_file, trajectory_file)
    
    # Calculate and plot density histogram
    print("Calculating density histogram...")
    avg_density = calculate_density_histogram(universe, n_frames=100)
    print(f"Average density: {avg_density:.4f} g/cm³")
    
    # Calculate 1D density profiles along each axis
    for axis in ['x', 'y', 'z']:
        print(f"Calculating 1D density profile along {axis}...")
        calculate_1d_density_profile(universe, axis=axis, n_frames=100)
    
    # Calculate 2D density maps
    for plane in ['xy', 'yz', 'xz']:
        print(f"Calculating 2D density map on {plane} plane...")
        calculate_2d_density_map(universe, plane=plane, n_frames=50)
    
    # Calculate radial density
    print("Calculating radial density...")
    calculate_radial_density(universe, n_frames=50)
    
    print("Density analysis complete!")

if __name__ == '__main__':
    main() 