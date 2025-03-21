#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.collections import PatchCollection
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

def visualize_nonbonded_cutoffs():
    """Create visualizations for non-bonded interaction cutoffs."""
    output_dir = create_output_dir()
    
    # Create two separate figures instead of one large figure
    # Figure 1: Potential energy curves and 2D representation
    fig1 = plt.figure(figsize=(12, 10))
    
    # 1. Potential energy curves
    ax1 = fig1.add_subplot(211)
    
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
    ax1.plot(r, lj, 'b-', label='Lennard-Jones', linewidth=2)
    ax1.plot(r, coulomb, 'r-', label='Coulomb', linewidth=2)
    ax1.plot(r, total, 'k-', label='Total', linewidth=2.5)
    
    # Add cutoff line
    ax1.axvline(x=cutoff, color='green', linestyle='--', linewidth=2, 
               label=f'Cutoff ({cutoff} nm)')
    
    # Shade the region beyond cutoff
    ax1.fill_between(r[cutoff_idx:], total[cutoff_idx:], 0, 
                    color='lightgray', alpha=0.5, hatch='///')
    
    # Add text for PME
    ax1.text(1.5, coulomb[cutoff_idx]/2, 'Long-range electrostatics\nhandled by PME', 
            fontsize=10, ha='center', va='center',
            bbox=dict(facecolor='lightyellow', alpha=0.8, boxstyle='round,pad=0.5'))
    
    # Set labels and title
    ax1.set_xlabel('Distance (nm)', fontweight='bold')
    ax1.set_ylabel('Potential Energy', fontweight='bold')
    ax1.set_title('Non-bonded Interaction Potentials', fontweight='bold')
    
    # Add grid and legend
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend(loc='upper right')
    
    # Set y-axis limits to focus on the relevant part
    ax1.set_ylim(-2, 5)
    
    # 2. 2D representation of cutoff sphere
    ax2 = fig1.add_subplot(212)
    
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
    ax2.scatter(positions[:, 0], positions[:, 1], c=colors, s=sizes, alpha=0.7)
    
    # Add central molecule
    ax2.scatter(0, 0, c='red', s=100, edgecolor='black', linewidth=1.5, label='Central Molecule')
    
    # Add cutoff circle
    ax2.add_patch(cutoff_circle)
    
    # Add PME grid in background
    grid_size = 0.2
    x_grid = np.arange(-1.5, 1.5 + grid_size, grid_size)
    y_grid = np.arange(-1.5, 1.5 + grid_size, grid_size)
    
    for x in x_grid:
        ax2.axvline(x=x, color='lightgray', linestyle='-', linewidth=0.5, alpha=0.3)
    for y in y_grid:
        ax2.axhline(y=y, color='lightgray', linestyle='-', linewidth=0.5, alpha=0.3)
    
    # Add PME annotation
    ax2.text(1.2, 1.2, 'PME Grid', fontsize=10, ha='center', va='center',
            bbox=dict(facecolor='lightyellow', alpha=0.8, boxstyle='round,pad=0.5'))
    
    # Set equal aspect ratio
    ax2.set_aspect('equal')
    
    # Set labels and title
    ax2.set_xlabel('X (nm)', fontweight='bold')
    ax2.set_ylabel('Y (nm)', fontweight='bold')
    ax2.set_title('2D Representation of Cutoff Sphere', fontweight='bold')
    
    # Add legend
    ax2.legend(loc='upper right')
    
    # Set limits
    ax2.set_xlim(-1.5, 1.5)
    ax2.set_ylim(-1.5, 1.5)
    
    # Add explanation text for figure 1
    explanation1 = (
        "Non-Bonded Interactions in MD Simulations (Part 1):\n\n"
        "• Short-range interactions (Lennard-Jones, Coulomb) calculated directly within cutoff radius (typically 1.0-1.2 nm)\n"
        "• Long-range electrostatics beyond cutoff handled by Particle Mesh Ewald (PME) method"
    )
    
    plt.figtext(0.5, 0.01, explanation1, ha='center', fontsize=12, 
               bbox=dict(facecolor='lightyellow', alpha=0.5, boxstyle='round,pad=0.5'))
    
    # Adjust layout for figure 1
    plt.tight_layout(rect=[0, 0.08, 1, 0.98])
    
    # Save figure 1
    plt.savefig(os.path.join(output_dir, 'nonbonded_cutoffs_part1.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'nonbonded_cutoffs_part1.pdf'), bbox_inches='tight')
    plt.close(fig1)
    
    # Figure 2: 3D representation and PME grid
    fig2 = plt.figure(figsize=(12, 10))
    
    # 3. 3D representation of cutoff sphere
    ax3 = fig2.add_subplot(211, projection='3d')
    
    # Create points on a sphere with fewer points to reduce complexity
    u = np.linspace(0, 2 * np.pi, 20)  # Reduced from 30
    v = np.linspace(0, np.pi, 20)      # Reduced from 30
    x = cutoff * np.outer(np.cos(u), np.sin(v))
    y = cutoff * np.outer(np.sin(u), np.sin(v))
    z = cutoff * np.outer(np.ones(np.size(u)), np.cos(v))
    
    # Plot cutoff sphere with reduced resolution
    ax3.plot_surface(x, y, z, color='green', alpha=0.2, linewidth=0.5, 
                    edgecolor='green', label='Cutoff Sphere', rstride=2, cstride=2)
    
    # Create 3D water molecules with fewer points
    np.random.seed(42)  # For reproducibility
    n_molecules_3d = 50  # Reduced from 100
    positions_3d = np.random.uniform(-1.5, 1.5, (n_molecules_3d, 3))
    
    # Calculate distances from center
    distances_3d = np.sqrt(positions_3d[:, 0]**2 + positions_3d[:, 1]**2 + positions_3d[:, 2]**2)
    
    # Color molecules based on distance (inside/outside cutoff)
    colors_3d = ['blue' if d <= cutoff else 'lightgray' for d in distances_3d]
    sizes_3d = [30 if d <= cutoff else 15 for d in distances_3d]
    
    # Plot molecules
    ax3.scatter(positions_3d[:, 0], positions_3d[:, 1], positions_3d[:, 2], 
               c=colors_3d, s=sizes_3d, alpha=0.7)
    
    # Add central molecule
    ax3.scatter(0, 0, 0, c='red', s=100, edgecolor='black', linewidth=1.5)
    
    # Set labels and title
    ax3.set_xlabel('X (nm)', fontweight='bold')
    ax3.set_ylabel('Y (nm)', fontweight='bold')
    ax3.set_zlabel('Z (nm)', fontweight='bold')
    ax3.set_title('3D Cutoff Sphere', fontweight='bold')
    
    # Set equal aspect ratio
    ax3.set_box_aspect([1, 1, 1])
    
    # Set limits
    ax3.set_xlim(-1.5, 1.5)
    ax3.set_ylim(-1.5, 1.5)
    ax3.set_zlim(-1.5, 1.5)
    
    # 4. PME grid representation
    ax4 = fig2.add_subplot(212)
    
    # Create a grid with reduced resolution
    grid_size = 0.3  # Increased from 0.25 to reduce grid points
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
    
    # Plot the PME grid - store the mesh object for the colorbar
    mesh = ax4.pcolormesh(X, Y, Z, cmap='coolwarm', alpha=0.7, shading='auto')
    
    # Add grid lines (fewer lines)
    for x in x_grid[::2]:  # Plot every other grid line
        ax4.axvline(x=x, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
    for y in y_grid[::2]:  # Plot every other grid line
        ax4.axhline(y=y, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
    
    # Add cutoff circle
    cutoff_circle2 = Circle((0, 0), cutoff, fill=False, edgecolor='green', 
                           linestyle='--', linewidth=2)
    ax4.add_patch(cutoff_circle2)
    
    # Add central molecule
    ax4.scatter(0, 0, c='red', s=100, edgecolor='black', linewidth=1.5)
    
    # Add molecules (fewer for clarity)
    ax4.scatter(positions[::2, 0], positions[::2, 1], c=[colors[i] for i in range(0, len(colors), 2)], 
               s=[sizes[i] for i in range(0, len(sizes), 2)], alpha=0.7, edgecolor='black')
    
    # Set equal aspect ratio
    ax4.set_aspect('equal')
    
    # Set labels and title
    ax4.set_xlabel('X (nm)', fontweight='bold')
    ax4.set_ylabel('Y (nm)', fontweight='bold')
    ax4.set_title('PME Grid for Long-Range Electrostatics', fontweight='bold')
    
    # Add colorbar - use the existing mesh object
    cbar = plt.colorbar(mesh, ax=ax4)
    cbar.set_label('Electrostatic Potential', fontweight='bold')
    
    # Set limits
    ax4.set_xlim(-1.5, 1.5)
    ax4.set_ylim(-1.5, 1.5)
    
    # Add explanation text for figure 2
    explanation2 = (
        "Non-Bonded Interactions in MD Simulations (Part 2):\n\n"
        "• PME uses a 3D grid and Fast Fourier Transforms to efficiently calculate long-range interactions\n"
        "• Cutoffs reduce computational cost while maintaining accuracy through PME for electrostatics"
    )
    
    plt.figtext(0.5, 0.01, explanation2, ha='center', fontsize=12, 
               bbox=dict(facecolor='lightyellow', alpha=0.5, boxstyle='round,pad=0.5'))
    
    # Adjust layout for figure 2
    plt.tight_layout(rect=[0, 0.08, 1, 0.98])
    
    # Save figure 2
    plt.savefig(os.path.join(output_dir, 'nonbonded_cutoffs_part2.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'nonbonded_cutoffs_part2.pdf'), bbox_inches='tight')
    plt.close(fig2)
    
    # Create a combined figure with all plots (simplified versions)
    fig_combined = plt.figure(figsize=(10, 8))
    
    # Add title for the combined figure
    fig_combined.suptitle('Non-Bonded Interactions in MD Simulations', fontsize=16, fontweight='bold')
    
    # 1. Potential energy curves (simplified)
    ax1_combined = fig_combined.add_subplot(221)
    ax1_combined.plot(r, lj, 'b-', label='LJ', linewidth=1.5)
    ax1_combined.plot(r, coulomb, 'r-', label='Coulomb', linewidth=1.5)
    ax1_combined.axvline(x=cutoff, color='green', linestyle='--', linewidth=1.5, label=f'Cutoff')
    ax1_combined.set_xlabel('Distance (nm)')
    ax1_combined.set_ylabel('Energy')
    ax1_combined.set_title('Interaction Potentials')
    ax1_combined.legend(loc='upper right', fontsize='small')
    ax1_combined.set_ylim(-2, 5)
    
    # 2. 2D cutoff (simplified)
    ax2_combined = fig_combined.add_subplot(222)
    ax2_combined.add_patch(Circle((0, 0), cutoff, fill=False, edgecolor='green', linestyle='--'))
    ax2_combined.scatter(positions[:, 0], positions[:, 1], c=colors, s=[s/2 for s in sizes], alpha=0.7)
    ax2_combined.scatter(0, 0, c='red', s=50, edgecolor='black')
    ax2_combined.set_aspect('equal')
    ax2_combined.set_xlabel('X (nm)')
    ax2_combined.set_ylabel('Y (nm)')
    ax2_combined.set_title('2D Cutoff Representation')
    ax2_combined.set_xlim(-1.5, 1.5)
    ax2_combined.set_ylim(-1.5, 1.5)
    
    # 3. 3D cutoff (simplified)
    ax3_combined = fig_combined.add_subplot(223, projection='3d')
    ax3_combined.plot_surface(x, y, z, color='green', alpha=0.2, linewidth=0.5, 
                             edgecolor='green', rstride=4, cstride=4)
    ax3_combined.scatter(positions_3d[::3, 0], positions_3d[::3, 1], positions_3d[::3, 2], 
                        c=[colors_3d[i] for i in range(0, len(colors_3d), 3)], 
                        s=[sizes_3d[i]/2 for i in range(0, len(sizes_3d), 3)], alpha=0.7)
    ax3_combined.scatter(0, 0, 0, c='red', s=50)
    ax3_combined.set_xlabel('X')
    ax3_combined.set_ylabel('Y')
    ax3_combined.set_zlabel('Z')
    ax3_combined.set_title('3D Cutoff Sphere')
    ax3_combined.set_box_aspect([1, 1, 1])
    
    # 4. PME grid (simplified)
    ax4_combined = fig_combined.add_subplot(224)
    mesh_combined = ax4_combined.pcolormesh(X, Y, Z, cmap='coolwarm', alpha=0.7, shading='auto')
    ax4_combined.add_patch(Circle((0, 0), cutoff, fill=False, edgecolor='green', linestyle='--'))
    ax4_combined.set_aspect('equal')
    ax4_combined.set_xlabel('X (nm)')
    ax4_combined.set_ylabel('Y (nm)')
    ax4_combined.set_title('PME Grid')
    ax4_combined.set_xlim(-1.5, 1.5)
    ax4_combined.set_ylim(-1.5, 1.5)
    
    # Add explanation text for combined figure
    explanation_combined = (
        "Key concepts: Direct calculation within cutoff, PME for long-range electrostatics"
    )
    
    plt.figtext(0.5, 0.01, explanation_combined, ha='center', fontsize=10, 
               bbox=dict(facecolor='lightyellow', alpha=0.5, boxstyle='round,pad=0.5'))
    
    # Adjust layout for combined figure
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    
    # Save combined figure
    plt.savefig(os.path.join(output_dir, 'nonbonded_cutoffs.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'nonbonded_cutoffs.pdf'), bbox_inches='tight')
    plt.close(fig_combined)
    
    print(f"Non-bonded interactions cutoff visualizations saved to {output_dir}")

if __name__ == "__main__":
    visualize_nonbonded_cutoffs() 