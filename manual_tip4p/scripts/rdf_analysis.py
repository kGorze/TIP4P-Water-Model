#!/usr/bin/env python3
"""
Radial Distribution Function (RDF) Analysis for TIP4P Water

This script calculates O-O, O-H, and H-H radial distribution functions
from GROMACS trajectory files and creates plots.

Usage:
    python rdf_analysis.py [tpr_file] [trajectory_file]

Default:
    Uses md.tpr and md.xtc in the current directory
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from utils import load_universe, save_plot

def calculate_and_plot_rdfs(universe, nbins=100, rmax=10.0):
    """
    Calculate and plot radial distribution functions for water
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory
    nbins : int
        Number of bins for RDF calculation
    rmax : float
        Maximum distance for RDF calculation (Å)
    
    Returns:
    --------
    dict
        Dictionary containing RDF data for different atom pairs
    """
    print("Calculating radial distribution functions...")
    
    # Select atoms
    try:
        oxygen = universe.select_atoms("name OW")
        hydrogen = universe.select_atoms("name HW*")
        
        if len(oxygen) == 0 or len(hydrogen) == 0:
            raise ValueError("Could not find water atoms with standard TIP4P names")
    except Exception as e:
        print(f"Error with standard TIP4P selection: {e}")
        print("Trying alternative selections...")
        
        # Try different water models
        for o_sel, h_sel in [
            ("name O", "name H*"),
            ("resname SOL and name O*", "resname SOL and name H*"),
            ("resname TIP4 and name O*", "resname TIP4 and name H*"),
            ("resname WAT and name O*", "resname WAT and name H*")
        ]:
            try:
                oxygen = universe.select_atoms(o_sel)
                hydrogen = universe.select_atoms(h_sel)
                if len(oxygen) > 0 and len(hydrogen) > 0:
                    print(f"Using alternative selections: {o_sel} and {h_sel}")
                    print(f"Found {len(oxygen)} oxygen atoms and {len(hydrogen)} hydrogen atoms")
                    break
            except:
                continue
        
        if len(oxygen) == 0 or len(hydrogen) == 0:
            raise ValueError("Could not find water atoms with any selection")
    
    print(f"Calculating RDFs for {len(oxygen)} oxygen atoms and {len(hydrogen)} hydrogen atoms")
    
    # Calculate RDFs
    from MDAnalysis.analysis import rdf
    
    # Define a minimum distance to avoid artifacts at very small distances
    # This is important because at very small distances, the RDF can have extremely large values
    # due to numerical issues or self-interactions
    rmin = 0.5  # Minimum distance in Angstroms
    
    # O-O RDF
    rdf_OO = rdf.InterRDF(oxygen, oxygen, nbins=nbins, range=(0.0, rmax))
    rdf_OO.run()
    
    # O-H RDF
    rdf_OH = rdf.InterRDF(oxygen, hydrogen, nbins=nbins, range=(0.0, rmax))
    rdf_OH.run()
    
    # H-H RDF
    rdf_HH = rdf.InterRDF(hydrogen, hydrogen, nbins=nbins, range=(0.0, rmax))
    rdf_HH.run()
    
    # Get bins and RDF data using the newer API
    # For MDAnalysis 2.0.0+, use results.bins and results.rdf
    # For older versions, fall back to bins and rdf attributes
    
    # Helper function to get bins and RDF data safely
    def get_rdf_data(rdf_obj):
        if hasattr(rdf_obj, 'results') and hasattr(rdf_obj.results, 'bins'):
            bins = rdf_obj.results.bins
            rdf_values = rdf_obj.results.rdf
            print("Using newer MDAnalysis API (results.bins and results.rdf)")
        else:
            bins = rdf_obj.bins
            rdf_values = rdf_obj.rdf
            print("Using older MDAnalysis API (bins and rdf attributes)")
        return bins, rdf_values
    
    # Get data for each RDF
    bins_OO, rdf_OO_values = get_rdf_data(rdf_OO)
    bins_OH, rdf_OH_values = get_rdf_data(rdf_OH)
    bins_HH, rdf_HH_values = get_rdf_data(rdf_HH)
    
    # Apply cutoff to avoid artifacts at very small distances
    # Create masks for valid regions (r > rmin)
    mask_OO = bins_OO >= rmin
    mask_OH = bins_OH >= rmin
    mask_HH = bins_HH >= rmin
    
    # Save data (full range)
    np.savetxt('../analysis/rdf_OO.dat', 
               np.column_stack((bins_OO, rdf_OO_values)),
               header='r (Å)\tg_OO(r)', 
               comments='# ')
    
    np.savetxt('../analysis/rdf_OH.dat', 
               np.column_stack((bins_OH, rdf_OH_values)),
               header='r (Å)\tg_OH(r)', 
               comments='# ')
    
    np.savetxt('../analysis/rdf_HH.dat', 
               np.column_stack((bins_HH, rdf_HH_values)),
               header='r (Å)\tg_HH(r)', 
               comments='# ')
    
    # Define a consistent x-axis limit for all plots
    x_limit = 6.0  # Focus on first two shells (0-6 Å)
    
    # Find reasonable y-limits by looking at the data in the valid range
    # This helps avoid the extreme values at very small distances
    max_OO = np.max(rdf_OO_values[mask_OO]) * 1.1  # Add 10% margin
    max_OH = np.max(rdf_OH_values[mask_OH]) * 1.1
    max_HH = np.max(rdf_HH_values[mask_HH]) * 1.1
    
    # Plot O-O RDF
    plt.figure(figsize=(8, 6))
    plt.plot(bins_OO[mask_OO], rdf_OO_values[mask_OO], linewidth=2)
    
    # Find and annotate key peaks for O-O RDF
    # First peak around 2.8 Å (nearest neighbor)
    first_peak_idx = np.argmax(rdf_OO_values[(bins_OO > 2.0) & (bins_OO < 3.5)])
    first_peak_r = bins_OO[first_peak_idx]
    first_peak_g = rdf_OO_values[first_peak_idx]
    
    # Second peak around 4.5 Å (second shell)
    second_peak_idx = np.argmax(rdf_OO_values[(bins_OO > 4.0) & (bins_OO < 5.5)])
    second_peak_r = bins_OO[second_peak_idx]
    second_peak_g = rdf_OO_values[second_peak_idx]
    
    # Add text annotations instead of arrows for cleaner look
    plt.text(first_peak_r - 0.3, first_peak_g + 0.2, f'First shell\n{first_peak_r:.2f} Å', 
             fontsize=9, ha='center', va='bottom', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    plt.text(second_peak_r - 0.3, second_peak_g + 0.2, f'Second shell\n{second_peak_r:.2f} Å', 
             fontsize=9, ha='center', va='bottom', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    plt.xlabel('r (Å)')
    plt.ylabel('g$_{O-O}$(r)')
    plt.title('Oxygen-Oxygen Radial Distribution Function')
    plt.grid(True, alpha=0.3)
    plt.xlim(rmin, x_limit)
    plt.ylim(0, max_OO)
    plt.savefig('../analysis/rdf_OO_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot O-H RDF
    plt.figure(figsize=(8, 6))
    plt.plot(bins_OH[mask_OH], rdf_OH_values[mask_OH], linewidth=2)
    
    # Find and annotate key peaks for O-H RDF
    # Intramolecular O-H peak around 1.0 Å
    intra_peak_idx = np.argmax(rdf_OH_values[(bins_OH > 0.9) & (bins_OH < 1.2)])
    intra_peak_r = bins_OH[intra_peak_idx]
    intra_peak_g = rdf_OH_values[intra_peak_idx]
    
    # Hydrogen bond peak around 1.8 Å
    hbond_peak_idx = np.argmax(rdf_OH_values[(bins_OH > 1.5) & (bins_OH < 2.0)])
    hbond_peak_r = bins_OH[hbond_peak_idx]
    hbond_peak_g = rdf_OH_values[hbond_peak_idx]
    
    # Add text annotations
    plt.text(intra_peak_r, intra_peak_g * 0.8, f'Intramolecular\nO-H bond\n{intra_peak_r:.2f} Å', 
             fontsize=9, ha='center', va='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    plt.text(hbond_peak_r, hbond_peak_g * 0.8, f'Hydrogen bond\n{hbond_peak_r:.2f} Å', 
             fontsize=9, ha='center', va='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    plt.xlabel('r (Å)')
    plt.ylabel('g$_{O-H}$(r)')
    plt.title('Oxygen-Hydrogen Radial Distribution Function')
    plt.grid(True, alpha=0.3)
    plt.xlim(rmin, x_limit)
    plt.ylim(0, max_OH)
    plt.savefig('../analysis/rdf_OH_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot H-H RDF
    plt.figure(figsize=(8, 6))
    plt.plot(bins_HH[mask_HH], rdf_HH_values[mask_HH], linewidth=2)
    
    # Find and annotate key peaks for H-H RDF
    # Intramolecular H-H peak around 1.6 Å
    intra_hh_peak_idx = np.argmax(rdf_HH_values[(bins_HH > 1.4) & (bins_HH < 1.8)])
    intra_hh_peak_r = bins_HH[intra_hh_peak_idx]
    intra_hh_peak_g = rdf_HH_values[intra_hh_peak_idx]
    
    # Add text annotation
    plt.text(intra_hh_peak_r, intra_hh_peak_g * 0.8, f'Intramolecular\nH-H distance\n{intra_hh_peak_r:.2f} Å', 
             fontsize=9, ha='center', va='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    plt.xlabel('r (Å)')
    plt.ylabel('g$_{H-H}$(r)')
    plt.title('Hydrogen-Hydrogen Radial Distribution Function')
    plt.grid(True, alpha=0.3)
    plt.xlim(rmin, x_limit)
    plt.ylim(0, max_HH)
    plt.savefig('../analysis/rdf_HH_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create a combined plot with improved styling
    plt.figure(figsize=(8, 6))
    plt.plot(bins_OO[mask_OO], rdf_OO_values[mask_OO], label='O-O', linewidth=2)
    plt.plot(bins_OH[mask_OH], rdf_OH_values[mask_OH], label='O-H', linewidth=2)
    plt.plot(bins_HH[mask_HH], rdf_HH_values[mask_HH], label='H-H', linewidth=2)
    
    # Add a legend
    plt.legend(loc='best')
    
    plt.xlabel('r (Å)')
    plt.ylabel('g(r)')
    plt.title('Combined RDF (TIP4P Water)')
    plt.grid(True, alpha=0.3)
    plt.xlim(rmin, x_limit)
    # Use a reasonable y-limit for the combined plot
    plt.ylim(0, max(max_OO, max_OH, max_HH))
    plt.savefig('../analysis/combined_rdf_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create a detailed subplot layout with consistent styling
    fig, axs = plt.subplots(3, 1, figsize=(8, 12), sharex=True)
    
    # O-O RDF
    axs[0].plot(bins_OO[mask_OO], rdf_OO_values[mask_OO], linewidth=2)
    axs[0].set_ylabel('g$_{O-O}$(r)')
    axs[0].set_title('Oxygen-Oxygen RDF')
    axs[0].grid(True, alpha=0.3)
    axs[0].set_ylim(0, max_OO)
    
    # Add text annotations for O-O
    axs[0].text(first_peak_r, first_peak_g, f'First shell: {first_peak_r:.2f} Å', 
               fontsize=9, ha='left', va='bottom', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    axs[0].text(second_peak_r, second_peak_g, f'Second shell: {second_peak_r:.2f} Å', 
               fontsize=9, ha='left', va='bottom', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    # O-H RDF
    axs[1].plot(bins_OH[mask_OH], rdf_OH_values[mask_OH], linewidth=2)
    axs[1].set_ylabel('g$_{O-H}$(r)')
    axs[1].set_title('Oxygen-Hydrogen RDF')
    axs[1].grid(True, alpha=0.3)
    axs[1].set_ylim(0, max_OH)
    
    # Add text annotations for O-H
    axs[1].text(intra_peak_r, intra_peak_g * 0.7, f'Intramolecular: {intra_peak_r:.2f} Å', 
               fontsize=9, ha='left', va='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    axs[1].text(hbond_peak_r, hbond_peak_g * 0.7, f'H-bond: {hbond_peak_r:.2f} Å', 
               fontsize=9, ha='left', va='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    # H-H RDF
    axs[2].plot(bins_HH[mask_HH], rdf_HH_values[mask_HH], linewidth=2)
    axs[2].set_xlabel('r (Å)')
    axs[2].set_ylabel('g$_{H-H}$(r)')
    axs[2].set_title('Hydrogen-Hydrogen RDF')
    axs[2].grid(True, alpha=0.3)
    axs[2].set_xlim(rmin, x_limit)
    axs[2].set_ylim(0, max_HH)
    
    # Add text annotation for H-H
    axs[2].text(intra_hh_peak_r, intra_hh_peak_g * 0.7, f'Intramolecular: {intra_hh_peak_r:.2f} Å', 
               fontsize=9, ha='left', va='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    plt.tight_layout()
    plt.savefig('../analysis/detailed_rdf_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create a zoomed-in version of the plots to show the important features
    # This is especially useful for the O-O RDF where the first and second peaks are important
    
    # O-O zoomed
    plt.figure(figsize=(8, 6))
    plt.plot(bins_OO, rdf_OO_values, linewidth=2)
    plt.xlabel('r (Å)')
    plt.ylabel('g$_{O-O}$(r)')
    plt.title('Oxygen-Oxygen RDF (Zoomed)')
    plt.grid(True, alpha=0.3)
    plt.xlim(2.0, 6.0)  # Focus on the first and second peaks
    plt.ylim(0, max_OO)
    plt.text(first_peak_r, first_peak_g, f'First shell: {first_peak_r:.2f} Å', 
            fontsize=9, ha='left', va='bottom', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    plt.text(second_peak_r, second_peak_g, f'Second shell: {second_peak_r:.2f} Å', 
            fontsize=9, ha='left', va='bottom', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    plt.savefig('../analysis/rdf_OO_zoomed_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # O-H zoomed
    plt.figure(figsize=(8, 6))
    plt.plot(bins_OH, rdf_OH_values, linewidth=2)
    plt.xlabel('r (Å)')
    plt.ylabel('g$_{O-H}$(r)')
    plt.title('Oxygen-Hydrogen RDF (Zoomed)')
    plt.grid(True, alpha=0.3)
    plt.xlim(0.8, 3.0)  # Focus on the intramolecular and H-bond peaks
    plt.ylim(0, max_OH)
    plt.text(intra_peak_r, intra_peak_g * 0.7, f'Intramolecular: {intra_peak_r:.2f} Å', 
            fontsize=9, ha='left', va='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    plt.text(hbond_peak_r, hbond_peak_g * 0.7, f'H-bond: {hbond_peak_r:.2f} Å', 
            fontsize=9, ha='left', va='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    plt.savefig('../analysis/rdf_OH_zoomed_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Return the RDF data
    return {
        'OO': {'bins': bins_OO, 'rdf': rdf_OO_values},
        'OH': {'bins': bins_OH, 'rdf': rdf_OH_values},
        'HH': {'bins': bins_HH, 'rdf': rdf_HH_values}
    }

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
    
    # Calculate and plot RDFs
    calculate_and_plot_rdfs(universe)

if __name__ == '__main__':
    main() 