#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Arc, FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
import matplotlib.colors as mcolors
import os

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)

def create_output_dir():
    """Create output directory if it doesn't exist."""
    output_dir = "visualization_outputs"
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def visualize_water_constraints():
    """Create a 3D visualization of water molecule with constraints."""
    output_dir = create_output_dir()
    
    # Create figure with two subplots: 3D model and 2D schematic
    fig = plt.figure(figsize=(15, 7))
    
    # 3D Model
    ax1 = fig.add_subplot(121, projection='3d')
    
    # Define water molecule coordinates (TIP4P model)
    # Oxygen at origin
    o_pos = np.array([0, 0, 0])
    
    # Hydrogens - positioned to create the correct HOH angle (104.52°)
    angle_rad = np.deg2rad(104.52)
    bond_length = 0.9572  # OH bond length in Å
    
    h1_pos = np.array([bond_length * np.sin(angle_rad/2), 0, bond_length * np.cos(angle_rad/2)])
    h2_pos = np.array([-bond_length * np.sin(angle_rad/2), 0, bond_length * np.cos(angle_rad/2)])
    
    # Virtual site (M-site) for TIP4P
    m_pos = np.array([0, 0, 0.15])
    
    # Plot atoms
    ax1.scatter(*o_pos, color='red', s=200, label='Oxygen')
    ax1.scatter(*h1_pos, color='white', s=100, edgecolors='black', label='Hydrogen')
    ax1.scatter(*h2_pos, color='white', s=100, edgecolors='black')
    ax1.scatter(*m_pos, color='blue', s=50, alpha=0.5, label='Virtual Site (M)')
    
    # Plot bonds
    ax1.plot([o_pos[0], h1_pos[0]], [o_pos[1], h1_pos[1]], [o_pos[2], h1_pos[2]], 'k-', linewidth=2)
    ax1.plot([o_pos[0], h2_pos[0]], [o_pos[1], h2_pos[1]], [o_pos[2], h2_pos[2]], 'k-', linewidth=2)
    ax1.plot([o_pos[0], m_pos[0]], [o_pos[1], m_pos[1]], [o_pos[2], m_pos[2]], 'b--', linewidth=1, alpha=0.5)
    
    # Add bond length annotations
    midpoint1 = (o_pos + h1_pos) / 2
    midpoint2 = (o_pos + h2_pos) / 2
    
    ax1.text(midpoint1[0], midpoint1[1], midpoint1[2], f"{bond_length} Å", color='blue', fontweight='bold')
    ax1.text(midpoint2[0], midpoint2[1], midpoint2[2], f"{bond_length} Å", color='blue', fontweight='bold')
    
    # Add angle annotation
    # Create an arrow to indicate the angle
    arrow1 = Arrow3D([h1_pos[0], o_pos[0]], [h1_pos[1], o_pos[1]], [h1_pos[2], o_pos[2]],
                    mutation_scale=10, lw=1, arrowstyle='->', color='green')
    arrow2 = Arrow3D([h2_pos[0], o_pos[0]], [h2_pos[1], o_pos[1]], [h2_pos[2], o_pos[2]],
                    mutation_scale=10, lw=1, arrowstyle='->', color='green')
    
    ax1.add_artist(arrow1)
    ax1.add_artist(arrow2)
    
    # Add angle text
    angle_pos = o_pos + np.array([0, 0.5, 0])
    ax1.text(angle_pos[0], angle_pos[1], angle_pos[2], f"HOH Angle: {104.52}°", 
             color='green', fontweight='bold')
    
    # Set equal aspect ratio
    ax1.set_box_aspect([1, 1, 1])
    
    # Set labels and title
    ax1.set_xlabel('X (Å)')
    ax1.set_ylabel('Y (Å)')
    ax1.set_zlabel('Z (Å)')
    ax1.set_title('3D Water Molecule with Constraints', fontweight='bold')
    
    # Add legend
    ax1.legend(loc='upper right')
    
    # 2D Schematic
    ax2 = fig.add_subplot(122)
    
    # Draw 2D water molecule
    o_2d = np.array([0, 0])
    h1_2d = np.array([np.sin(angle_rad/2) * bond_length, np.cos(angle_rad/2) * bond_length])
    h2_2d = np.array([-np.sin(angle_rad/2) * bond_length, np.cos(angle_rad/2) * bond_length])
    
    # Plot atoms
    ax2.scatter(*o_2d, color='red', s=300, label='O')
    ax2.scatter(*h1_2d, color='white', s=150, edgecolors='black', label='H')
    ax2.scatter(*h2_2d, color='white', s=150, edgecolors='black')
    
    # Plot bonds
    ax2.plot([o_2d[0], h1_2d[0]], [o_2d[1], h1_2d[1]], 'k-', linewidth=3)
    ax2.plot([o_2d[0], h2_2d[0]], [o_2d[1], h2_2d[1]], 'k-', linewidth=3)
    
    # Add bond length annotations with arrows
    ax2.annotate(f"{bond_length} Å", 
                xy=((o_2d[0] + h1_2d[0])/2, (o_2d[1] + h1_2d[1])/2),
                xytext=((o_2d[0] + h1_2d[0])/2 + 0.3, (o_2d[1] + h1_2d[1])/2 + 0.3),
                arrowprops=dict(facecolor='blue', shrink=0.05, width=1.5),
                fontweight='bold', color='blue')
    
    ax2.annotate(f"{bond_length} Å", 
                xy=((o_2d[0] + h2_2d[0])/2, (o_2d[1] + h2_2d[1])/2),
                xytext=((o_2d[0] + h2_2d[0])/2 - 0.3, (o_2d[1] + h2_2d[1])/2 + 0.3),
                arrowprops=dict(facecolor='blue', shrink=0.05, width=1.5),
                fontweight='bold', color='blue')
    
    # Add angle annotation
    angle = Arc((0, 0), 0.5, 0.5, 
               theta1=np.rad2deg(np.arctan2(h2_2d[1], h2_2d[0])), 
               theta2=np.rad2deg(np.arctan2(h1_2d[1], h1_2d[0])), 
               color='green', linewidth=2)
    ax2.add_patch(angle)
    
    # Add angle text
    ax2.annotate(f"{104.52}°", 
                xy=(0, 0.1),
                xytext=(0, -0.5),
                arrowprops=dict(facecolor='green', shrink=0.05, width=1.5),
                fontweight='bold', color='green', ha='center')
    
    # Set equal aspect ratio
    ax2.set_aspect('equal')
    
    # Set labels and title
    ax2.set_xlabel('X (Å)')
    ax2.set_ylabel('Y (Å)')
    ax2.set_title('2D Water Molecule Schematic', fontweight='bold')
    
    # Add grid
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    # Set limits
    ax2.set_xlim(-1.5, 1.5)
    ax2.set_ylim(-1, 1.5)
    
    # Add explanation text
    explanation = (
        "Water Molecule Constraints in MD Simulations:\n"
        "• OH Bond Length: Fixed at 0.9572 Å\n"
        "• HOH Angle: Fixed at 104.52°\n\n"
        "Benefits of Constraints:\n"
        "• Eliminates high-frequency bond vibrations\n"
        "• Allows larger integration time steps (2 fs)\n"
        "• Maintains correct water geometry\n"
        "• Implemented using SHAKE/LINCS algorithms"
    )
    
    plt.figtext(0.5, 0.01, explanation, ha='center', fontsize=12, 
               bbox=dict(facecolor='lightyellow', alpha=0.5, boxstyle='round,pad=0.5'))
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.08, 1, 0.98])
    
    # Save figure
    plt.savefig(os.path.join(output_dir, 'water_molecule_constraints.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'water_molecule_constraints.pdf'), bbox_inches='tight')
    
    print(f"Water molecule constraint visualization saved to {output_dir}")
    plt.close()

if __name__ == "__main__":
    visualize_water_constraints() 