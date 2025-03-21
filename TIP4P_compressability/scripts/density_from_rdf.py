#!/usr/bin/env python3
"""
Density Calculation from Radial Distribution Function (RDF) for TIP4P Water

This script calculates the local density and coordination numbers from RDF data
and provides insights into the structure of water.

Usage:
    python density_from_rdf.py [rdf_file_OO] [rdf_file_OH] [rdf_file_HH]

Default:
    Uses ../analysis/rdf_OO.dat, ../analysis/rdf_OH.dat, and ../analysis/rdf_HH.dat
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from scipy.signal import find_peaks
import argparse

def calculate_density_from_rdf(rdf_file, bulk_density=None, molecular_weight=18.01528, 
                               temperature=273, pressure=1.0, rmin_cutoff=1.0):
    """
    Calculate local density and coordination number from RDF data
    
    Parameters:
    -----------
    rdf_file : str
        Path to the RDF data file
    bulk_density : float, optional
        Bulk density of the system in g/cm³. If None, estimated for water at given T and P.
    molecular_weight : float, optional
        Molecular weight in g/mol (default: water)
    temperature : float, optional
        Temperature in K (default: 273K for TIP4P simulation)
    pressure : float, optional
        Pressure in bar (default: 1.0 bar for TIP4P simulation)
    rmin_cutoff : float, optional
        Minimum distance cutoff in Å to avoid artifacts at very small distances (default: 1.0 Å)
        
    Returns:
    --------
    tuple
        (r, g_r, local_density, coordination_number, mask)
        where mask is a boolean array indicating valid data points (r >= rmin_cutoff)
    """
    # Load RDF data
    try:
        data = np.loadtxt(rdf_file)
        r = data[:, 0]  # Distance in Å
        g_r = data[:, 1]  # RDF values
    except Exception as e:
        print(f"Error loading RDF data from {rdf_file}: {e}")
        return None, None, None, None, None
    
    # Create mask for valid data (r >= rmin_cutoff)
    mask = r >= rmin_cutoff
    
    # Estimate bulk density if not provided
    if bulk_density is None:
        # For TIP4P water at 273K (0°C) and 1 bar, density is approximately 0.9998 g/cm³
        # This is more accurate than the general formula for this specific case
        if abs(temperature - 273) < 1 and abs(pressure - 1.0) < 0.1:
            bulk_density = 0.9998  # g/cm³, specific for TIP4P at 273K and 1 bar
            print(f"Using TIP4P water density at 273K and 1 bar: {bulk_density:.4f} g/cm³")
        else:
            # Simple estimation for water density at different temperatures
            # More accurate models exist but this is a reasonable approximation
            # Density of water at 4°C is 1.0 g/cm³ and decreases with temperature
            bulk_density = 1.0 - 7e-5 * (temperature - 277)
            # Pressure correction (simplified)
            bulk_density += 5e-5 * (pressure - 1.0)  # Approximate compressibility effect
            print(f"Estimated bulk density of water at {temperature}K and {pressure} bar: {bulk_density:.4f} g/cm³")
    
    # Convert bulk density to particles per Å³
    # 1 g/cm³ = 1 g/cm³ * (1 mol / MW g) * (6.022e23 particles / 1 mol) * (1 cm³ / 1e24 Å³)
    # = 6.022e23 / (MW * 1e24) particles/Å³ = 6.022e-1 / MW particles/Å³
    avogadro = 6.022e23  # Avogadro's number
    bulk_number_density = bulk_density * avogadro / (molecular_weight * 1e24)  # particles/Å³
    print(f"Bulk number density: {bulk_number_density:.6e} particles/Å³")
    
    # Calculate local density: ρ(r) = ρ₀ * g(r)
    local_density = bulk_number_density * g_r
    
    # Calculate coordination number: n(r) = 4πρ₀∫₀ʳ g(r')r'²dr'
    # This gives the average number of particles within distance r of a central particle
    r_squared = r**2
    integrand = g_r * r_squared
    
    # Perform the integration using the trapezoidal rule
    integral = cumtrapz(integrand, r, initial=0)
    coordination_number = 4 * np.pi * bulk_number_density * integral
    
    return r, g_r, local_density, coordination_number, mask

def find_rdf_features(r, g_r, mask, rdf_type):
    """
    Find important features in the RDF curve (peaks and minima)
    
    Parameters:
    -----------
    r : numpy.ndarray
        Distance array
    g_r : numpy.ndarray
        RDF values
    mask : numpy.ndarray
        Boolean mask for valid data points
    rdf_type : str
        Type of RDF (O-O, O-H, or H-H)
        
    Returns:
    --------
    dict
        Dictionary with peak and minima information
    """
    # Apply mask to get valid data
    r_valid = r[mask]
    g_r_valid = g_r[mask]
    
    # Define search ranges based on RDF type
    if rdf_type == 'O-O':
        # For O-O, first peak is typically around 2.8 Å
        first_peak_range = (2.5, 3.1)
        first_min_range = (3.1, 3.8)
        second_peak_range = (4.0, 5.0)
    elif rdf_type == 'O-H':
        # For O-H, first peak is typically around 1.8 Å (hydrogen bond)
        first_peak_range = (1.6, 2.0)
        first_min_range = (2.0, 2.5)
        second_peak_range = (3.0, 4.0)
    elif rdf_type == 'H-H':
        # For H-H, peaks are less pronounced but first is around 2.4 Å
        first_peak_range = (2.0, 2.8)
        first_min_range = (2.8, 3.5)
        second_peak_range = (3.5, 4.5)
    else:
        # Default ranges if type is unknown
        first_peak_range = (2.0, 3.5)
        first_min_range = (3.5, 4.5)
        second_peak_range = (4.5, 6.0)
    
    # Find peaks within specified ranges
    features = {}
    
    # First peak
    mask_peak1 = (r_valid >= first_peak_range[0]) & (r_valid <= first_peak_range[1])
    if np.any(mask_peak1):
        idx_peak1 = np.argmax(g_r_valid[mask_peak1]) + np.where(mask_peak1)[0][0]
        features['first_peak'] = {
            'r': r_valid[idx_peak1],
            'g_r': g_r_valid[idx_peak1],
            'idx': idx_peak1
        }
    
    # First minimum (after first peak)
    if 'first_peak' in features:
        start_idx = features['first_peak']['idx']
        mask_min1 = (r_valid >= first_min_range[0]) & (r_valid <= first_min_range[1])
        if np.any(mask_min1):
            idx_min1 = np.argmin(g_r_valid[mask_min1]) + np.where(mask_min1)[0][0]
            features['first_min'] = {
                'r': r_valid[idx_min1],
                'g_r': g_r_valid[idx_min1],
                'idx': idx_min1
            }
    
    # Second peak
    mask_peak2 = (r_valid >= second_peak_range[0]) & (r_valid <= second_peak_range[1])
    if np.any(mask_peak2):
        idx_peak2 = np.argmax(g_r_valid[mask_peak2]) + np.where(mask_peak2)[0][0]
        features['second_peak'] = {
            'r': r_valid[idx_peak2],
            'g_r': g_r_valid[idx_peak2],
            'idx': idx_peak2
        }
    
    return features

def plot_density_and_coordination(r, g_r, local_density, coordination_number, mask, rdf_type, output_dir, features=None):
    """
    Create plots for RDF, local density, and coordination number
    
    Parameters:
    -----------
    r : numpy.ndarray
        Distance array
    g_r : numpy.ndarray
        RDF values
    local_density : numpy.ndarray
        Local density values
    coordination_number : numpy.ndarray
        Coordination number values
    mask : numpy.ndarray
        Boolean mask for valid data points
    rdf_type : str
        Type of RDF (O-O, O-H, or H-H)
    output_dir : str
        Directory to save the plots
    features : dict, optional
        Dictionary with peak and minima information
    """
    # Apply mask to get valid data for plotting
    r_valid = r[mask]
    g_r_valid = g_r[mask]
    local_density_valid = local_density[mask]
    coordination_number_valid = coordination_number[mask]
    
    # Create a figure with 3 subplots
    fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    
    # Plot RDF with proper y-axis limits
    axs[0].plot(r_valid, g_r_valid, 'b-', linewidth=2)
    axs[0].set_ylabel('g(r)')
    axs[0].set_title(f'{rdf_type} Radial Distribution Function (TIP4P, 273K)')
    axs[0].grid(True, alpha=0.3)
    
    # Set reasonable y-axis limits based on RDF type
    if rdf_type == 'O-O':
        max_g = min(5.0, np.percentile(g_r_valid, 99))
        axs[0].set_ylim(0, max_g * 1.1)
    elif rdf_type == 'O-H':
        max_g = min(25.0, np.percentile(g_r_valid, 99))
        axs[0].set_ylim(0, max_g * 1.1)
    elif rdf_type == 'H-H':
        max_g = min(5.0, np.percentile(g_r_valid, 99))
        axs[0].set_ylim(0, max_g * 1.1)
    
    # Plot local density with proper y-axis limits
    axs[1].plot(r_valid, local_density_valid, 'g-', linewidth=2)
    axs[1].set_ylabel('ρ(r) [particles/Å³]')
    axs[1].set_title(f'{rdf_type} Local Density')
    axs[1].grid(True, alpha=0.3)
    
    # Set reasonable y-axis limits for density
    max_density = np.percentile(local_density_valid, 99)
    axs[1].set_ylim(0, max_density * 1.1)
    
    # Plot coordination number
    axs[2].plot(r_valid, coordination_number_valid, 'r-', linewidth=2)
    axs[2].set_xlabel('r [Å]')
    axs[2].set_ylabel('n(r)')
    axs[2].set_title(f'{rdf_type} Coordination Number')
    axs[2].grid(True, alpha=0.3)
    
    # Set common x-axis limits
    for ax in axs:
        ax.set_xlim(1.0, 10.0)
    
    # Add annotations for important features if available
    if features:
        # First peak
        if 'first_peak' in features:
            peak = features['first_peak']
            axs[0].axvline(x=peak['r'], color='k', linestyle='--', alpha=0.5)
            axs[0].text(peak['r'], peak['g_r']*0.9, f'First peak: {peak["r"]:.2f} Å', 
                       ha='center', va='top', bbox=dict(facecolor='white', alpha=0.7))
        
        # First minimum (first coordination shell)
        if 'first_min' in features:
            min1 = features['first_min']
            cn_first_shell = coordination_number[np.where(r >= min1['r'])[0][0]]
            
            # Mark on RDF plot
            axs[0].axvline(x=min1['r'], color='k', linestyle=':', alpha=0.5)
            
            # Mark on coordination number plot
            axs[2].axvline(x=min1['r'], color='k', linestyle=':', alpha=0.5)
            axs[2].axhline(y=cn_first_shell, color='k', linestyle=':', alpha=0.5)
            axs[2].text(min1['r']+0.5, cn_first_shell, 
                       f'First shell: {cn_first_shell:.2f} molecules at {min1["r"]:.2f} Å', 
                       ha='left', va='center', bbox=dict(facecolor='white', alpha=0.7))
        
        # Second peak
        if 'second_peak' in features:
            peak2 = features['second_peak']
            axs[0].axvline(x=peak2['r'], color='k', linestyle='--', alpha=0.5)
            axs[0].text(peak2['r'], peak2['g_r']*0.9, f'Second peak: {peak2["r"]:.2f} Å', 
                       ha='center', va='top', bbox=dict(facecolor='white', alpha=0.7))
    
    # Add interpretation for water structure
    if rdf_type == 'O-O':
        axs[0].text(0.98, 0.98, 
                   "O-O RDF shows water molecule arrangement:\n"
                   "- First peak (~2.8 Å): nearest neighbor water molecules\n"
                   "- First minimum (~3.4 Å): boundary of first coordination shell\n"
                   "- Second peak (~4.5 Å): second coordination shell",
                   transform=axs[0].transAxes, ha='right', va='top',
                   bbox=dict(facecolor='white', alpha=0.7), fontsize=9)
    elif rdf_type == 'O-H':
        axs[0].text(0.98, 0.98, 
                   "O-H RDF shows hydrogen bonding:\n"
                   "- First peak (~1.8 Å): hydrogen bonds\n"
                   "- Second peak (~3.3 Å): H atoms in second coordination shell",
                   transform=axs[0].transAxes, ha='right', va='top',
                   bbox=dict(facecolor='white', alpha=0.7), fontsize=9)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'density_from_rdf_{rdf_type}.png'), dpi=300, bbox_inches='tight')
    plt.close()

def analyze_water_structure(rdf_files, output_dir='../analysis', rmin_cutoff=1.0):
    """
    Analyze water structure from RDF data
    
    Parameters:
    -----------
    rdf_files : dict
        Dictionary with RDF file paths for different atom pairs
    output_dir : str
        Directory to save the output files and plots
    rmin_cutoff : float
        Minimum distance cutoff in Å to avoid artifacts at very small distances
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each RDF file
    results = {}
    for rdf_type, rdf_file in rdf_files.items():
        print(f"\nAnalyzing {rdf_type} RDF from {rdf_file}...")
        
        # Calculate density and coordination number
        # Use TIP4P parameters: 273K and 1.0 bar
        r, g_r, local_density, coordination_number, mask = calculate_density_from_rdf(
            rdf_file, 
            bulk_density=None,  # Will use the TIP4P-specific value
            temperature=273,    # 273K (0°C) from simulation parameters
            pressure=1.0,       # 1.0 bar from simulation parameters
            rmin_cutoff=rmin_cutoff  # Minimum distance cutoff
        )
        
        if r is None:
            print(f"Skipping {rdf_type} analysis due to errors.")
            continue
        
        # Find important features in the RDF
        features = find_rdf_features(r, g_r, mask, rdf_type)
        print(f"Detected features for {rdf_type}:")
        for feature, data in features.items():
            print(f"  - {feature}: r = {data['r']:.2f} Å, g(r) = {data['g_r']:.2f}")
        
        # Store results
        results[rdf_type] = {
            'r': r,
            'g_r': g_r,
            'local_density': local_density,
            'coordination_number': coordination_number,
            'mask': mask,
            'features': features
        }
        
        # Create plots
        plot_density_and_coordination(r, g_r, local_density, coordination_number, mask, rdf_type, output_dir, features)
        
        # Save data to file (only valid data points)
        output_file = os.path.join(output_dir, f'density_from_rdf_{rdf_type}.dat')
        valid_data = np.column_stack((
            r[mask], 
            g_r[mask], 
            local_density[mask], 
            coordination_number[mask]
        ))
        np.savetxt(output_file, valid_data,
                  header='r [Å]\tg(r)\tρ(r) [particles/Å³]\tn(r)', 
                  comments='# ')
        print(f"Data saved to {output_file}")
    
    # Analyze water structure if O-O RDF is available
    if 'O-O' in results and 'first_min' in results['O-O']['features']:
        r = results['O-O']['r']
        g_r = results['O-O']['g_r']
        coordination = results['O-O']['coordination_number']
        features = results['O-O']['features']
        
        # Get first coordination shell information
        first_min_r = features['first_min']['r']
        idx_first_min = np.where(r >= first_min_r)[0][0]
        cn_first_shell = coordination[idx_first_min]
        
        # Try to find second coordination shell
        second_shell_info = "Not detected"
        cn_second_shell = None
        
        if 'second_peak' in features:
            # Look for minimum after second peak
            second_peak_r = features['second_peak']['r']
            idx_second_peak = np.where(r >= second_peak_r)[0][0]
            
            # Search for minimum after second peak
            search_range = slice(idx_second_peak, idx_second_peak + 20)
            if len(g_r[search_range]) > 0:
                idx_second_min = idx_second_peak + np.argmin(g_r[search_range])
                second_min_r = r[idx_second_min]
                cn_second_shell = coordination[idx_second_min]
                second_shell_info = f"Radius: {second_min_r:.2f} Å, Coordination: {cn_second_shell:.2f} molecules"
        
        # Write a summary of water structure
        with open(os.path.join(output_dir, 'water_structure_summary.txt'), 'w') as f:
            f.write("Water Structure Analysis from RDF (TIP4P, 273K, 1 bar)\n")
            f.write("==================================================\n\n")
            f.write(f"First coordination shell:\n")
            f.write(f"  - Radius: {first_min_r:.2f} Å\n")
            f.write(f"  - Coordination number: {cn_first_shell:.2f} water molecules\n\n")
            
            if cn_second_shell is not None:
                f.write(f"Second coordination shell:\n")
                f.write(f"  - Radius: {second_min_r:.2f} Å\n")
                f.write(f"  - Coordination number: {cn_second_shell:.2f} water molecules\n")
                f.write(f"  - Additional molecules in second shell: {cn_second_shell - cn_first_shell:.2f}\n\n")
            
            # Compare with literature values for TIP4P at 273K
            f.write("Comparison with literature values for TIP4P water at 273K:\n")
            f.write("  - Expected first shell coordination: 4.3-4.5 molecules\n")
            f.write("  - Expected second shell coordination: 15-20 molecules\n")
            f.write("  - Expected first peak position: ~2.75-2.80 Å\n")
            f.write("  - Expected first minimum position: ~3.2-3.4 Å\n\n")
            
            f.write("Note: TIP4P water model at 273K has a density of approximately 0.9998 g/cm³\n")
            f.write("      and slightly different coordination numbers compared to experimental water.\n")
        
        print(f"\nWater structure summary saved to {os.path.join(output_dir, 'water_structure_summary.txt')}")

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Calculate density from RDF data for TIP4P water at 273K')
    parser.add_argument('--oo', dest='rdf_oo', default='../analysis/rdf_OO.dat',
                        help='Path to O-O RDF data file (default: ../analysis/rdf_OO.dat)')
    parser.add_argument('--oh', dest='rdf_oh', default='../analysis/rdf_OH.dat',
                        help='Path to O-H RDF data file (default: ../analysis/rdf_OH.dat)')
    parser.add_argument('--hh', dest='rdf_hh', default='../analysis/rdf_HH.dat',
                        help='Path to H-H RDF data file (default: ../analysis/rdf_HH.dat)')
    parser.add_argument('--output', dest='output_dir', default='../analysis',
                        help='Output directory (default: ../analysis)')
    parser.add_argument('--density', dest='density', type=float, default=None,
                        help='Bulk density in g/cm³ (default: 0.9998 g/cm³ for TIP4P at 273K)')
    parser.add_argument('--temp', dest='temperature', type=float, default=273,
                        help='Temperature in K (default: 273K from simulation parameters)')
    parser.add_argument('--pressure', dest='pressure', type=float, default=1.0,
                        help='Pressure in bar (default: 1.0 bar from simulation parameters)')
    parser.add_argument('--rmin', dest='rmin_cutoff', type=float, default=1.0,
                        help='Minimum distance cutoff in Å (default: 1.0 Å)')
    
    args = parser.parse_args()
    
    # Check if RDF files exist
    rdf_files = {
        'O-O': args.rdf_oo,
        'O-H': args.rdf_oh,
        'H-H': args.rdf_hh
    }
    
    for rdf_type, rdf_file in rdf_files.items():
        if not os.path.exists(rdf_file):
            print(f"Warning: {rdf_type} RDF file {rdf_file} not found.")
    
    # Analyze water structure
    analyze_water_structure(rdf_files, args.output_dir, args.rmin_cutoff)

if __name__ == '__main__':
    main() 