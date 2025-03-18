#!/usr/bin/env python3

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.base import AnalysisFromFunction
import sys
import os

def calculate_density(atomgroup):
    """Calculate density for a given atom group."""
    # Get total mass in amu
    total_mass = sum(atomgroup.masses)
    # Get volume in nm^3
    volume = atomgroup.dimensions[0] * atomgroup.dimensions[1] * atomgroup.dimensions[2]
    # Calculate density in g/cm^3
    # 1 amu/nm^3 = 1.660539067 g/cm^3
    density = (total_mass / volume) * 1.660539067
    return density

def main():
    # Ensure analysis directory exists
    analysis_dir = "../analysis"
    os.makedirs(analysis_dir, exist_ok=True)
    
    print("Loading trajectory...")
    # Load the system
    u = mda.Universe("../md.tpr", "../md.trr")
    
    print("Calculating density and tracking box dimensions...")
    densities = []
    times = []
    box_x = []
    box_y = []
    box_z = []
    
    # Iterate through trajectory
    for ts in u.trajectory:
        if ts.frame % 10 == 0:  # Print progress every 10 frames
            sys.stdout.write(f"\rProcessing frame {ts.frame}/{len(u.trajectory)}")
            sys.stdout.flush()
        
        density = calculate_density(u.atoms)
        densities.append(density)
        times.append(ts.time)
        box_x.append(ts.dimensions[0])
        box_y.append(ts.dimensions[1])
        box_z.append(ts.dimensions[2])
    
    print("\nGenerating plots...")
    # Convert to numpy arrays
    times = np.array(times)
    densities = np.array(densities)
    box_x = np.array(box_x)
    box_y = np.array(box_y)
    box_z = np.array(box_z)
    
    # Plot 1: Density vs Time (standalone)
    plt.figure(figsize=(10, 6))
    plt.plot(times, densities)
    plt.xlabel('Time (ps)')
    plt.ylabel('Density (g/cm³)')
    plt.title('Water Density vs Time')
    plt.grid(True)
    plt.savefig(os.path.join(analysis_dir, 'density_vs_time.png'))
    plt.close()
    
    # Plot 2: Combined density and box dimensions
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 1, 1)
    plt.plot(times, densities)
    plt.xlabel('Time (ps)')
    plt.ylabel('Density (g/cm³)')
    plt.title('Water Density vs Time')
    plt.grid(True)
    
    # Plot box dimensions
    plt.subplot(2, 1, 2)
    plt.plot(times, box_x, label='X')
    plt.plot(times, box_y, label='Y')
    plt.plot(times, box_z, label='Z')
    plt.xlabel('Time (ps)')
    plt.ylabel('Box dimension (nm)')
    plt.title('Box Dimensions vs Time')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_dir, 'density_and_box.png'))
    plt.close()
    
    print(f"Saved plots as {analysis_dir}/density_vs_time.png and {analysis_dir}/density_and_box.png")
    
    # Save numerical data
    np.savetxt(os.path.join(analysis_dir, 'density_vs_time.txt'),
               np.column_stack((times, densities, box_x, box_y, box_z)),
               header='Time(ps) Density(g/cm³) Box_X(nm) Box_Y(nm) Box_Z(nm)',
               delimiter='\t')
    print(f"Saved numerical data as {analysis_dir}/density_vs_time.txt")
    
    # Calculate and print statistics
    avg_density = np.mean(densities)
    std_density = np.std(densities)
    print(f"\nAverage density: {avg_density:.4f} ± {std_density:.4f} g/cm³")
    print(f"Box dimensions (average):")
    print(f"X: {np.mean(box_x):.4f} ± {np.std(box_x):.4f} nm")
    print(f"Y: {np.mean(box_y):.4f} ± {np.std(box_y):.4f} nm")
    print(f"Z: {np.mean(box_z):.4f} ± {np.std(box_z):.4f} nm")

if __name__ == "__main__":
    main() 