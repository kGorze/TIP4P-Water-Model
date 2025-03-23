#!/usr/bin/env python3
"""
Extract Correct Temperature and Energy Data from GROMACS

This script uses multiple approaches to extract temperature and energy data from
GROMACS files, attempting to identify and correct unit conversion issues.

Usage:
    python extract_correct_data.py
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import glob
import re

def find_gromacs_files():
    """Find GROMACS files in the parent directory"""
    gromacs_files = {
        'topology': [],
        'trajectory': [],
        'energy': [],
        'structure': []
    }
    
    # Find potential topology files (.top, .tpr)
    top_files = glob.glob('../*.top') + glob.glob('../*.tpr')
    gromacs_files['topology'] = top_files
    
    # Find trajectory files (.trr, .xtc)
    traj_files = glob.glob('../*.trr') + glob.glob('../*.xtc')
    gromacs_files['trajectory'] = traj_files
    
    # Find energy files (.edr)
    energy_files = glob.glob('../*.edr')
    gromacs_files['energy'] = energy_files
    
    # Find structure files (.gro, .pdb)
    struct_files = glob.glob('../*.gro') + glob.glob('../*.pdb')
    gromacs_files['structure'] = struct_files
    
    return gromacs_files

def try_gmx_energy_directly(edr_file='../md.edr'):
    """
    Try to extract energy information using gmx energy directly from the command line
    """
    print(f"\nAttempting direct gmx energy extraction from {edr_file}")
    
    properties = {
        'Temperature': [],
        'Pressure': [],
        'Potential': [],
        'Kinetic': [],
        'Total-Energy': []
    }
    
    # First, let's see what's available in this file
    try:
        # List available energy terms
        cmd = f"echo 0 | gmx energy -f {edr_file} -o /dev/null"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Error running gmx energy to list available terms: {result.stderr}")
            return None
        
        # Extract terms from output
        output = result.stdout + result.stderr
        available_terms = []
        
        # Check for the line that lists the terms
        term_section = False
        for line in output.split('\n'):
            if "Select the terms you want from the following list" in line:
                term_section = True
                continue
            if term_section and len(line.strip()) > 0 and not line.startswith("Select"):
                # Extract the term name, which is typically after a number and spaces
                match = re.search(r'^\s*\d+\s+([\w-]+)', line)
                if match:
                    available_terms.append(match.group(1))
        
        print(f"Available energy terms: {available_terms}")
        
        # Now extract each property we're interested in
        for prop_name in properties.keys():
            if prop_name not in available_terms and prop_name.replace('-', ' ') not in available_terms:
                print(f"Warning: {prop_name} not found in available terms")
                continue
            
            # Create input file for gmx energy
            with open('prop_input.txt', 'w') as f:
                f.write(f"{prop_name}\n0\n")
            
            # Run gmx energy and extract the property
            cmd = f"gmx energy -f {edr_file} -o temp_prop.xvg < prop_input.txt"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"Error extracting {prop_name}: {result.stderr}")
                continue
            
            # Read the output XVG file
            if not os.path.exists('temp_prop.xvg'):
                print(f"Error: XVG file for {prop_name} was not created")
                continue
            
            # Extract header information to get units
            unit = "Unknown"
            with open('temp_prop.xvg', 'r') as f:
                for line in f:
                    if line.startswith("@ title"):
                        title = line.split('"')[1]
                        print(f"  {prop_name} title: {title}")
                    if line.startswith("@ yaxis label"):
                        label = line.split('"')[1]
                        print(f"  {prop_name} y-axis: {label}")
                        # Try to extract unit from label
                        if "(" in label and ")" in label:
                            unit = label.split("(")[1].split(")")[0]
                            print(f"  {prop_name} unit: {unit}")
                    if not line.startswith('@') and not line.startswith('#'):
                        break
            
            # Extract data values
            times = []
            values = []
            with open('temp_prop.xvg', 'r') as f:
                for line in f:
                    if not line.startswith('@') and not line.startswith('#'):
                        cols = line.strip().split()
                        if len(cols) >= 2:
                            times.append(float(cols[0]))
                            values.append(float(cols[1]))
            
            if len(values) > 0:
                print(f"  {prop_name} range: {min(values):.2f} to {max(values):.2f} {unit}")
                print(f"  {prop_name} mean: {np.mean(values):.2f} {unit}")
                
                properties[prop_name] = {
                    'times': np.array(times),
                    'values': np.array(values),
                    'unit': unit
                }
            
            # Clean up
            if os.path.exists('temp_prop.xvg'):
                os.remove('temp_prop.xvg')
            if os.path.exists('prop_input.txt'):
                os.remove('prop_input.txt')
        
        return properties
    
    except Exception as e:
        print(f"Error in gmx energy extraction: {str(e)}")
        return None

def analyze_temperature_data(properties):
    """
    Analyze temperature data to identify potential unit issues
    """
    if 'Temperature' not in properties or not properties['Temperature']:
        print("No temperature data available for analysis")
        return
    
    temp_data = properties['Temperature']
    temp_values = temp_data['values']
    temp_unit = temp_data['unit']
    
    print("\nTemperature Analysis:")
    print(f"Range: {min(temp_values):.2f} to {max(temp_values):.2f} {temp_unit}")
    print(f"Mean: {np.mean(temp_values):.2f} {temp_unit}")
    print(f"Standard deviation: {np.std(temp_values):.2f} {temp_unit}")
    
    # Check if temperature is in a reasonable range for water simulations
    if np.mean(temp_values) < 0:
        print("Warning: Negative temperature values detected (physically impossible)")
        
        # Possibilities: sign error, index instead of value, or unit conversion issue
        if -300 < np.mean(temp_values) < -200:
            print("Suggestion: Temperature values may have incorrect sign. Try multiplying by -1.")
            corrected = -temp_values
            print(f"  Corrected range: {min(corrected):.2f} to {max(corrected):.2f} K")
            print(f"  Corrected mean: {np.mean(corrected):.2f} K")
        elif np.mean(temp_values) < -10000:
            print("Warning: Temperature values are extremely negative")
            print("Suggestion: This might be a data structure issue or fundamental miscommunication from GROMACS")
    elif 0 < np.mean(temp_values) < 100 and temp_unit.lower() != 'c':
        print("Warning: Temperature values appear to be in Celsius, but should be in Kelvin")
        corrected = temp_values + 273.15
        print(f"  Corrected range: {min(corrected):.2f} to {max(corrected):.2f} K")
        print(f"  Corrected mean: {np.mean(corrected):.2f} K")
    elif 200 < np.mean(temp_values) < 400:
        print("Temperature values appear to be in a reasonable range for water simulations in Kelvin")

def analyze_energy_data(properties):
    """
    Analyze energy data to identify potential unit issues
    """
    energy_props = ['Potential', 'Kinetic', 'Total-Energy']
    
    for prop in energy_props:
        if prop in properties and properties[prop]:
            data = properties[prop]
            values = data['values']
            unit = data['unit']
            
            print(f"\n{prop} Energy Analysis:")
            print(f"Range: {min(values):.2f} to {max(values):.2f} {unit}")
            print(f"Mean: {np.mean(values):.2f} {unit}")
            
            # Check if energy values are in reasonable range
            mean_val = np.mean(values)
            if abs(mean_val) > 10000:
                print(f"Warning: {prop} values are very large in magnitude")
                print("Suggestion: Values may need scaling to per-molecule units")
                
                # TIP4P water systems typically have energy around -40 kJ/mol per molecule
                # Estimate number of molecules based on this
                est_molecules = abs(mean_val) / 40
                if est_molecules > 100:
                    print(f"  Estimated ~{est_molecules:.0f} molecules in system")
                    corrected = values / est_molecules
                    print(f"  Corrected per-molecule range: {min(corrected):.2f} to {max(corrected):.2f} {unit}/molecule")
                    print(f"  Corrected per-molecule mean: {np.mean(corrected):.2f} {unit}/molecule")

def plot_corrected_data(properties, output_dir='.'):
    """
    Generate plots of the data with appropriate corrections
    """
    print("\nGenerating corrected plots...")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Check if we have temperature data
    if 'Temperature' in properties and properties['Temperature']:
        temp_data = properties['Temperature']
        temp_times = temp_data['times']
        temp_values = temp_data['values']
        temp_unit = temp_data['unit']
        
        plt.figure(figsize=(10, 6))
        plt.plot(temp_times, temp_values, 'b-', label='Raw Temperature')
        
        # If temperature is negative, add corrected line
        if np.mean(temp_values) < 0:
            plt.plot(temp_times, -temp_values, 'r-', label='Sign-corrected Temperature')
        
        plt.xlabel('Time (ps)')
        plt.ylabel(f'Temperature ({temp_unit})')
        plt.title('Temperature vs Time')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.savefig(f"{output_dir}/corrected_temperature.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # Check if we have energy data
    energy_props = ['Potential', 'Kinetic', 'Total-Energy']
    for prop in energy_props:
        if prop in properties and properties[prop]:
            data = properties[prop]
            times = data['times']
            values = data['values']
            unit = data['unit']
            
            plt.figure(figsize=(10, 6))
            plt.plot(times, values, 'b-', label=f'Raw {prop}')
            
            # If energy values are very large, add scaled line
            if abs(np.mean(values)) > 10000:
                est_molecules = abs(np.mean(values)) / 40
                if est_molecules > 100:
                    plt.plot(times, values / est_molecules, 'r-', 
                            label=f'Scaled {prop} (per molecule)')
            
            plt.xlabel('Time (ps)')
            plt.ylabel(f'{prop} ({unit})')
            plt.title(f'{prop} vs Time')
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.savefig(f"{output_dir}/corrected_{prop.lower().replace('-', '_')}.png", dpi=300, bbox_inches='tight')
            plt.close()
    
    # Create combined energy plot
    has_energy = any(prop in properties and properties[prop] for prop in energy_props)
    if has_energy:
        plt.figure(figsize=(12, 8))
        
        for prop in energy_props:
            if prop in properties and properties[prop]:
                data = properties[prop]
                times = data['times']
                values = data['values']
                unit = data['unit']
                
                # For combined plot, always scale if very large
                if abs(np.mean(values)) > 10000:
                    est_molecules = abs(np.mean(values)) / 40
                    if est_molecules > 100:
                        plt.plot(times, values / est_molecules, label=f'{prop} (per molecule)')
                else:
                    plt.plot(times, values, label=prop)
        
        plt.xlabel('Time (ps)')
        plt.ylabel(f'Energy (kJ/mol)')
        plt.title('Energy Components vs Time')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.savefig(f"{output_dir}/corrected_combined_energy.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # Create a heat capacity plot using corrected values if we have temperature and energy
    if ('Temperature' in properties and properties['Temperature'] and 
        'Total-Energy' in properties and properties['Total-Energy']):
        
        temp_data = properties['Temperature']
        energy_data = properties['Total-Energy']
        temp_values = temp_data['values'] 
        energy_values = energy_data['values']
        
        # Apply corrections if needed
        if np.mean(temp_values) < 0:
            temp_values = -temp_values
        
        if abs(np.mean(energy_values)) > 10000:
            est_molecules = abs(np.mean(energy_values)) / 40
            if est_molecules > 100:
                energy_values = energy_values / est_molecules
        
        # Create temperature-energy scatter plot
        plt.figure(figsize=(10, 6))
        plt.scatter(temp_values, energy_values, alpha=0.5, s=10)
        plt.xlabel('Temperature (K)')
        plt.ylabel('Energy (kJ/mol)')
        plt.title('Energy vs Temperature')
        plt.grid(True, alpha=0.3)
        plt.savefig(f"{output_dir}/corrected_energy_vs_temp.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Calculate and plot heat capacity
        # Sort by temperature
        sorted_indices = np.argsort(temp_values)
        temp_sorted = temp_values[sorted_indices]
        energy_sorted = energy_values[sorted_indices]
        
        # Bin the data
        temp_min = np.floor(min(temp_sorted))
        temp_max = np.ceil(max(temp_sorted))
        temp_bins = np.arange(temp_min, temp_max + 0.5, 0.5)
        
        binned_temps = []
        binned_energies = []
        
        for i in range(len(temp_bins) - 1):
            bin_start = temp_bins[i]
            bin_end = temp_bins[i + 1]
            bin_indices = np.where((temp_sorted >= bin_start) & (temp_sorted < bin_end))[0]
            
            if len(bin_indices) > 0:
                binned_temps.append(np.mean(temp_sorted[bin_indices]))
                binned_energies.append(np.mean(energy_sorted[bin_indices]))
        
        # Calculate heat capacity (derivative of energy with respect to temperature)
        heat_capacity = []
        binned_temps_center = []
        
        for i in range(1, len(binned_temps)):
            dE = binned_energies[i] - binned_energies[i-1]
            dT = binned_temps[i] - binned_temps[i-1]
            
            if abs(dT) > 1e-6:  # Avoid division by very small values
                Cv = dE / dT
                heat_capacity.append(Cv)
                binned_temps_center.append((binned_temps[i] + binned_temps[i-1]) / 2)
        
        if len(heat_capacity) > 5:  # Need enough points for a meaningful plot
            plt.figure(figsize=(10, 6))
            plt.plot(binned_temps_center, heat_capacity, 'o-')
            plt.xlabel('Temperature (K)')
            plt.ylabel('Heat Capacity (kJ/mol·K)')
            plt.title('Heat Capacity vs Temperature')
            plt.grid(True, alpha=0.3)
            
            # Try to identify peaks
            for i in range(2, len(heat_capacity) - 2):
                if (heat_capacity[i] > heat_capacity[i-1] and 
                    heat_capacity[i] > heat_capacity[i-2] and
                    heat_capacity[i] > heat_capacity[i+1] and
                    heat_capacity[i] > heat_capacity[i+2] and
                    heat_capacity[i] > np.mean(heat_capacity) + np.std(heat_capacity)):
                    
                    plt.plot(binned_temps_center[i], heat_capacity[i], 'ro', markersize=10)
                    plt.annotate(f'{binned_temps_center[i]:.1f} K', 
                                xy=(binned_temps_center[i], heat_capacity[i]),
                                xytext=(binned_temps_center[i]+2, heat_capacity[i]*1.1),
                                arrowprops=dict(arrowstyle='->'))
            
            plt.savefig(f"{output_dir}/corrected_heat_capacity.png", dpi=300, bbox_inches='tight')
            plt.close()

def main():
    """Main function"""
    print("GROMACS Data Extraction and Correction")
    
    # Find GROMACS files
    print("\nSearching for GROMACS files...")
    gromacs_files = find_gromacs_files()
    
    for category, files in gromacs_files.items():
        print(f"Found {len(files)} {category} files:")
        for f in files:
            print(f"  {f}")
    
    # Try to extract energy data directly
    energy_data = None
    for edr_file in gromacs_files['energy']:
        print(f"\nExtracting data from {edr_file}...")
        energy_data = try_gmx_energy_directly(edr_file)
        if energy_data:
            break
    
    if not energy_data:
        print("\nCould not extract energy data from any of the files.")
        return
    
    # Analyze temperature data
    analyze_temperature_data(energy_data)
    
    # Analyze energy data
    analyze_energy_data(energy_data)
    
    # Generate corrected plots
    plot_corrected_data(energy_data)
    
    print("\nCorrection and analysis complete. Check visualization folder for plots.")

if __name__ == "__main__":
    main() 