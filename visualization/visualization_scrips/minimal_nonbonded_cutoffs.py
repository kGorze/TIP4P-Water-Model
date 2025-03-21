#!/usr/bin/env python3

import numpy as np
import matplotlib
# Use Agg backend which is more stable
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import os

def create_output_dir():
    """Create output directory if it doesn't exist."""
    output_dir = "visualization_outputs"
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def lennard_jones_potential(r, epsilon=0.65, sigma=0.3):
    """Calculate Lennard-Jones potential."""
    return 4 * epsilon * ((sigma/r)**12 - (sigma/r)**6)

def coulomb_potential(r, q1=0.5, q2=-0.5, epsilon_r=1):
    """Calculate Coulomb potential."""
    k = 8.9875517923e9  # Coulomb constant in N·m²/C²
    # Convert to kJ/mol scale for comparison with LJ
    k_scaled = k * 6.022e-4
    return k_scaled * q1 * q2 / (epsilon_r * r)

def plot_potentials():
    """Plot potential energy curves."""
    output_dir = create_output_dir()
    
    # Create figure with minimal settings
    plt.figure(figsize=(8, 6), dpi=80)
    
    # Distance range (nm)
    r = np.linspace(0.2, 2.0, 500)  # Reduced number of points
    
    # Calculate potentials
    lj = lennard_jones_potential(r)
    coulomb = coulomb_potential(r)
    total = lj + coulomb
    
    # Define cutoff
    cutoff = 1.0  # nm
    cutoff_idx = np.abs(r - cutoff).argmin()
    
    # Plot potentials
    plt.plot(r, lj, 'b-', label='Lennard-Jones')
    plt.plot(r, coulomb, 'r-', label='Coulomb')
    plt.plot(r, total, 'k-', label='Total')
    
    # Add cutoff line
    plt.axvline(x=cutoff, color='green', linestyle='--', label=f'Cutoff ({cutoff} nm)')
    
    # Set labels and title
    plt.xlabel('Distance (nm)')
    plt.ylabel('Potential Energy')
    plt.title('Non-bonded Interaction Potentials')
    
    # Add legend
    plt.legend(loc='upper right')
    
    # Set y-axis limits to focus on the relevant part
    plt.ylim(-2, 5)
    
    # Save figure with minimal settings
    plt.savefig(os.path.join(output_dir, 'nonbonded_potentials.png'), dpi=100)
    plt.close()
    
    print(f"Potential energy curves saved to {output_dir}")

def plot_2d_cutoff():
    """Plot 2D representation of cutoff sphere."""
    output_dir = create_output_dir()
    
    # Create figure with minimal settings
    plt.figure(figsize=(6, 6), dpi=80)
    
    # Define cutoff
    cutoff = 1.0  # nm
    
    # Create a circle representing the cutoff
    cutoff_circle = Circle((0, 0), cutoff, fill=False, edgecolor='green', 
                          linestyle='--', label=f'Cutoff ({cutoff} nm)')
    
    # Create a collection of water molecules (reduced number)
    np.random.seed(42)  # For reproducibility
    n_molecules = 30
    positions = np.random.uniform(-1.5, 1.5, (n_molecules, 2))
    
    # Calculate distances from center
    distances = np.sqrt(positions[:, 0]**2 + positions[:, 1]**2)
    
    # Color molecules based on distance (inside/outside cutoff)
    colors = ['blue' if d <= cutoff else 'lightgray' for d in distances]
    sizes = [30 if d <= cutoff else 20 for d in distances]
    
    # Plot molecules
    plt.scatter(positions[:, 0], positions[:, 1], c=colors, s=sizes)
    
    # Add central molecule
    plt.scatter(0, 0, c='red', s=50, label='Central Molecule')
    
    # Add cutoff circle
    plt.gca().add_patch(cutoff_circle)
    
    # Set equal aspect ratio
    plt.gca().set_aspect('equal')
    
    # Set labels and title
    plt.xlabel('X (nm)')
    plt.ylabel('Y (nm)')
    plt.title('2D Representation of Cutoff Sphere')
    
    # Add legend
    plt.legend(loc='upper right')
    
    # Set limits
    plt.xlim(-1.5, 1.5)
    plt.ylim(-1.5, 1.5)
    
    # Save figure with minimal settings
    plt.savefig(os.path.join(output_dir, 'cutoff_2d.png'), dpi=100)
    plt.close()
    
    print(f"2D cutoff representation saved to {output_dir}")

def plot_pme_grid():
    """Plot PME grid representation."""
    output_dir = create_output_dir()
    
    # Create figure with minimal settings
    plt.figure(figsize=(6, 6), dpi=80)
    
    # Define cutoff
    cutoff = 1.0  # nm
    
    # Create a collection of water molecules (reduced number)
    np.random.seed(42)  # For reproducibility
    n_molecules = 30
    positions = np.random.uniform(-1.5, 1.5, (n_molecules, 2))
    
    # Calculate distances from center
    distances = np.sqrt(positions[:, 0]**2 + positions[:, 1]**2)
    
    # Create a grid (reduced resolution)
    grid_size = 0.4  # Larger grid size for fewer points
    x_grid = np.arange(-1.5, 1.5 + grid_size, grid_size)
    y_grid = np.arange(-1.5, 1.5 + grid_size, grid_size)
    X, Y = np.meshgrid(x_grid, y_grid)
    
    # Create a potential field (simplified for visualization)
    Z = np.zeros_like(X)
    for i in range(len(positions)):
        px, py = positions[i]
        # Only consider molecules within 2x cutoff for the field
        if distances[i] <= 2*cutoff:
            # Create a Gaussian-like potential field
            Z += np.exp(-((X - px)**2 + (Y - py)**2) / (0.4**2))
    
    # Plot the PME grid
    mesh = plt.pcolormesh(X, Y, Z, cmap='coolwarm', shading='auto')
    
    # Add cutoff circle
    cutoff_circle = Circle((0, 0), cutoff, fill=False, edgecolor='green', 
                          linestyle='--')
    plt.gca().add_patch(cutoff_circle)
    
    # Add central molecule
    plt.scatter(0, 0, c='red', s=50)
    
    # Set equal aspect ratio
    plt.gca().set_aspect('equal')
    
    # Set labels and title
    plt.xlabel('X (nm)')
    plt.ylabel('Y (nm)')
    plt.title('PME Grid for Long-Range Electrostatics')
    
    # Add colorbar
    plt.colorbar(mesh, label='Electrostatic Potential')
    
    # Set limits
    plt.xlim(-1.5, 1.5)
    plt.ylim(-1.5, 1.5)
    
    # Save figure with minimal settings
    plt.savefig(os.path.join(output_dir, 'pme_grid.png'), dpi=100)
    plt.close()
    
    print(f"PME grid representation saved to {output_dir}")

def main():
    """Run all visualizations."""
    # Create separate plots to avoid the image size issue
    plot_potentials()
    plot_2d_cutoff()
    plot_pme_grid()
    
    print("All non-bonded interaction visualizations completed successfully.")

if __name__ == "__main__":
    main() 