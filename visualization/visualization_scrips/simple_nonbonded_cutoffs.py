#!/usr/bin/env python3

import numpy as np
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
    
    # Create figure
    plt.figure(figsize=(10, 6))
    
    # Distance range (nm)
    r = np.linspace(0.2, 2.0, 1000)
    
    # Calculate potentials
    lj = lennard_jones_potential(r)
    coulomb = coulomb_potential(r)
    total = lj + coulomb
    
    # Define cutoff
    cutoff = 1.0  # nm
    cutoff_idx = np.abs(r - cutoff).argmin()
    
    # Plot potentials
    plt.plot(r, lj, 'b-', label='Lennard-Jones', linewidth=2)
    plt.plot(r, coulomb, 'r-', label='Coulomb', linewidth=2)
    plt.plot(r, total, 'k-', label='Total', linewidth=2.5)
    
    # Add cutoff line
    plt.axvline(x=cutoff, color='green', linestyle='--', linewidth=2, 
               label=f'Cutoff ({cutoff} nm)')
    
    # Shade the region beyond cutoff
    plt.fill_between(r[cutoff_idx:], total[cutoff_idx:], 0, 
                    color='lightgray', alpha=0.5, hatch='///')
    
    # Add text for PME
    plt.text(1.5, coulomb[cutoff_idx]/2, 'Long-range electrostatics\nhandled by PME', 
            fontsize=10, ha='center', va='center',
            bbox=dict(facecolor='lightyellow', alpha=0.8, boxstyle='round,pad=0.5'))
    
    # Set labels and title
    plt.xlabel('Distance (nm)', fontweight='bold')
    plt.ylabel('Potential Energy', fontweight='bold')
    plt.title('Non-bonded Interaction Potentials', fontweight='bold')
    
    # Add grid and legend
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='upper right')
    
    # Set y-axis limits to focus on the relevant part
    plt.ylim(-2, 5)
    
    # Add explanation text
    explanation = (
        "Non-Bonded Interactions in MD Simulations:\n"
        "• Short-range interactions calculated directly within cutoff radius (typically 1.0-1.2 nm)\n"
        "• Long-range electrostatics beyond cutoff handled by Particle Mesh Ewald (PME) method"
    )
    
    plt.figtext(0.5, 0.01, explanation, ha='center', fontsize=10, 
               bbox=dict(facecolor='lightyellow', alpha=0.5, boxstyle='round,pad=0.5'))
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.15, 1, 0.98])
    
    # Save figure
    plt.savefig(os.path.join(output_dir, 'nonbonded_potentials.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'nonbonded_potentials.pdf'), bbox_inches='tight')
    plt.close()
    
    print(f"Potential energy curves saved to {output_dir}")

def plot_2d_cutoff():
    """Plot 2D representation of cutoff sphere."""
    output_dir = create_output_dir()
    
    # Create figure
    plt.figure(figsize=(8, 8))
    
    # Define cutoff
    cutoff = 1.0  # nm
    
    # Create a circle representing the cutoff
    cutoff_circle = Circle((0, 0), cutoff, fill=False, edgecolor='green', 
                          linestyle='--', linewidth=2, label=f'Cutoff ({cutoff} nm)')
    
    # Create a collection of water molecules
    np.random.seed(42)  # For reproducibility
    n_molecules = 50
    positions = np.random.uniform(-1.5, 1.5, (n_molecules, 2))
    
    # Calculate distances from center
    distances = np.sqrt(positions[:, 0]**2 + positions[:, 1]**2)
    
    # Color molecules based on distance (inside/outside cutoff)
    colors = ['blue' if d <= cutoff else 'lightgray' for d in distances]
    sizes = [50 if d <= cutoff else 30 for d in distances]
    
    # Plot molecules
    plt.scatter(positions[:, 0], positions[:, 1], c=colors, s=sizes, alpha=0.7)
    
    # Add central molecule
    plt.scatter(0, 0, c='red', s=100, edgecolor='black', linewidth=1.5, label='Central Molecule')
    
    # Add cutoff circle
    plt.gca().add_patch(cutoff_circle)
    
    # Add PME grid in background
    grid_size = 0.2
    x_grid = np.arange(-1.5, 1.5 + grid_size, grid_size)
    y_grid = np.arange(-1.5, 1.5 + grid_size, grid_size)
    
    for x in x_grid:
        plt.axvline(x=x, color='lightgray', linestyle='-', linewidth=0.5, alpha=0.3)
    for y in y_grid:
        plt.axhline(y=y, color='lightgray', linestyle='-', linewidth=0.5, alpha=0.3)
    
    # Add PME annotation
    plt.text(1.2, 1.2, 'PME Grid', fontsize=10, ha='center', va='center',
            bbox=dict(facecolor='lightyellow', alpha=0.8, boxstyle='round,pad=0.5'))
    
    # Set equal aspect ratio
    plt.gca().set_aspect('equal')
    
    # Set labels and title
    plt.xlabel('X (nm)', fontweight='bold')
    plt.ylabel('Y (nm)', fontweight='bold')
    plt.title('2D Representation of Cutoff Sphere', fontweight='bold')
    
    # Add legend
    plt.legend(loc='upper right')
    
    # Set limits
    plt.xlim(-1.5, 1.5)
    plt.ylim(-1.5, 1.5)
    
    # Add explanation text
    explanation = (
        "Cutoff-based Calculation in MD Simulations:\n"
        "• Interactions within cutoff (blue) calculated directly\n"
        "• Long-range electrostatics (gray) handled by PME\n"
        "• PME grid used for efficient calculation of long-range forces"
    )
    
    plt.figtext(0.5, 0.01, explanation, ha='center', fontsize=10, 
               bbox=dict(facecolor='lightyellow', alpha=0.5, boxstyle='round,pad=0.5'))
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.15, 1, 0.98])
    
    # Save figure
    plt.savefig(os.path.join(output_dir, 'cutoff_2d.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'cutoff_2d.pdf'), bbox_inches='tight')
    plt.close()
    
    print(f"2D cutoff representation saved to {output_dir}")

def plot_pme_grid():
    """Plot PME grid representation."""
    output_dir = create_output_dir()
    
    # Create figure
    plt.figure(figsize=(8, 8))
    
    # Define cutoff
    cutoff = 1.0  # nm
    
    # Create a collection of water molecules
    np.random.seed(42)  # For reproducibility
    n_molecules = 50
    positions = np.random.uniform(-1.5, 1.5, (n_molecules, 2))
    
    # Calculate distances from center
    distances = np.sqrt(positions[:, 0]**2 + positions[:, 1]**2)
    
    # Create a grid
    grid_size = 0.3
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
            Z += np.exp(-((X - px)**2 + (Y - py)**2) / (0.3**2))
    
    # Plot the PME grid
    mesh = plt.pcolormesh(X, Y, Z, cmap='coolwarm', alpha=0.7, shading='auto')
    
    # Add grid lines
    for x in x_grid[::2]:
        plt.axvline(x=x, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
    for y in y_grid[::2]:
        plt.axhline(y=y, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
    
    # Add cutoff circle
    cutoff_circle = Circle((0, 0), cutoff, fill=False, edgecolor='green', 
                          linestyle='--', linewidth=2)
    plt.gca().add_patch(cutoff_circle)
    
    # Add central molecule
    plt.scatter(0, 0, c='red', s=100, edgecolor='black', linewidth=1.5)
    
    # Color molecules based on distance (inside/outside cutoff)
    colors = ['blue' if d <= cutoff else 'lightgray' for d in distances]
    sizes = [50 if d <= cutoff else 30 for d in distances]
    
    # Add molecules (fewer for clarity)
    plt.scatter(positions[::2, 0], positions[::2, 1], 
               c=[colors[i] for i in range(0, len(colors), 2)], 
               s=[sizes[i] for i in range(0, len(sizes), 2)], 
               alpha=0.7, edgecolor='black')
    
    # Set equal aspect ratio
    plt.gca().set_aspect('equal')
    
    # Set labels and title
    plt.xlabel('X (nm)', fontweight='bold')
    plt.ylabel('Y (nm)', fontweight='bold')
    plt.title('PME Grid for Long-Range Electrostatics', fontweight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(mesh)
    cbar.set_label('Electrostatic Potential', fontweight='bold')
    
    # Set limits
    plt.xlim(-1.5, 1.5)
    plt.ylim(-1.5, 1.5)
    
    # Add explanation text
    explanation = (
        "Particle Mesh Ewald (PME) Method:\n"
        "• Efficiently calculates long-range electrostatics\n"
        "• Uses a grid and Fast Fourier Transforms\n"
        "• Complements direct calculation within cutoff"
    )
    
    plt.figtext(0.5, 0.01, explanation, ha='center', fontsize=10, 
               bbox=dict(facecolor='lightyellow', alpha=0.5, boxstyle='round,pad=0.5'))
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.15, 1, 0.98])
    
    # Save figure
    plt.savefig(os.path.join(output_dir, 'pme_grid.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'pme_grid.pdf'), bbox_inches='tight')
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