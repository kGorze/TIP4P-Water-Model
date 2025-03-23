#!/usr/bin/env python3
"""
Inspect GROMACS Energy Data

This script directly examines GROMACS energy data to diagnose any unit conversion
issues or scaling problems. It extracts temperature, energy, and other properties
directly from the energy file without applying any automatic conversions.

Usage:
    python inspect_gromacs_data.py
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess

def extract_raw_gromacs_data(edr_file='../md.edr', tpr_file='../md.tpr'):
    """
    Extract raw thermodynamic data from GROMACS energy files without any conversions
    
    Parameters:
    -----------
    edr_file : str
        Path to GROMACS energy file
    tpr_file : str
        Path to GROMACS topology file
        
    Returns:
    --------
    dict
        Dictionary with extracted properties
    """
    properties = ["Temperature", "Pressure", "Potential", "Kinetic-En.", "Total-Energy", "Density"]
    property_data = {}
    
    for prop in properties:
        # Create gmx energy input file
        with open('prop_input.txt', 'w') as f:
            f.write(f"{prop}\n0\n")
        
        # Run gmx energy command to extract just this property
        cmd = f"gmx energy -f {edr_file} -s {tpr_file} -o {prop.lower().replace('-', '_')}.xvg < prop_input.txt"
        try:
            subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"Error running gmx energy for {prop}: {e}")
            continue
        finally:
            if os.path.exists('prop_input.txt'):
                os.remove('prop_input.txt')
        
        # Read the XVG file
        xvg_file = f"{prop.lower().replace('-', '_')}.xvg"
        if not os.path.exists(xvg_file):
            print(f"XVG file {xvg_file} was not created")
            continue
        
        try:
            # Extract data from XVG file, preserving header information
            header_lines = []
            data_lines = []
            
            with open(xvg_file, 'r') as f:
                for line in f:
                    if line.startswith('@') or line.startswith('#'):
                        header_lines.append(line)
                    else:
                        data_lines.append(line)
            
            # Parse units and labels from header
            x_label = "Time"
            y_label = prop
            x_unit = "ps"
            y_unit = "Unknown"
            
            for line in header_lines:
                if line.startswith("@ xaxis label"):
                    x_label = line.split('"')[1]
                elif line.startswith("@ yaxis label"):
                    y_label = line.split('"')[1]
                elif line.startswith("@ s0 legend"):
                    prop_name = line.split('"')[1]
                    y_label = prop_name
            
            # Extract units from labels if present
            if "(" in x_label and ")" in x_label:
                x_unit = x_label.split("(")[1].split(")")[0]
            if "(" in y_label and ")" in y_label:
                y_unit = y_label.split("(")[1].split(")")[0]
            
            # Parse data
            time = []
            values = []
            
            for line in data_lines:
                parts = line.strip().split()
                if len(parts) >= 2:
                    time.append(float(parts[0]))
                    values.append(float(parts[1]))
            
            # Store in dictionary
            property_data[prop] = {
                'time': np.array(time),
                'values': np.array(values),
                'x_unit': x_unit,
                'y_unit': y_unit,
                'x_label': x_label,
                'y_label': y_label,
                'header': header_lines
            }
            
            print(f"Extracted {prop} data: {len(time)} points, units: {y_unit}")
            
            # Calculate basic statistics
            if len(values) > 0:
                print(f"  Min: {np.min(values):.6g}, Max: {np.max(values):.6g}")
                print(f"  Mean: {np.mean(values):.6g}, Std: {np.std(values):.6g}")
            
        except Exception as e:
            print(f"Error reading {xvg_file}: {e}")
        finally:
            # Clean up XVG file
            if os.path.exists(xvg_file):
                os.remove(xvg_file)
    
    return property_data

def suggest_corrections(property_data):
    """
    Analyze the extracted data and suggest corrections for unit issues
    
    Parameters:
    -----------
    property_data : dict
        Dictionary with extracted properties
        
    Returns:
    --------
    dict
        Dictionary with suggested correction factors
    """
    corrections = {}
    
    # Check temperature
    if 'Temperature' in property_data:
        temp = property_data['Temperature']['values']
        temp_mean = np.mean(temp)
        
        if temp_mean < 0:
            print("\nTemperature Diagnostics:")
            print(f"Temperature values are negative (mean: {temp_mean:.6g}), which is physically impossible")
            
            # Check if it could be a sign issue
            if -350 < temp_mean < -250:
                corrections['Temperature'] = -1
                print(f"Suggestion: Temperature may have an incorrect sign. Try multiplying by -1")
            else:
                # Check if it could be an offset issue
                # For example, temperatures might be relative to a reference temperature
                for offset in [273.15, -273.15]:
                    adjusted = temp + offset
                    if 250 < np.mean(adjusted) < 350:
                        corrections['Temperature'] = {'offset': offset}
                        print(f"Suggestion: Temperature may need an offset of {offset}")
                        break
                
                # Check if it could be a scaling issue
                if 'Temperature' not in corrections:
                    possible_scales = [10**i for i in range(-6, 7)]
                    for scale in possible_scales:
                        adjusted = temp * scale
                        if 250 < np.mean(adjusted) < 350:
                            corrections['Temperature'] = {'scale': scale}
                            print(f"Suggestion: Temperature may need scaling by {scale}")
                            break
                
                # If still unable to correct, try more complex transformations
                if 'Temperature' not in corrections:
                    # Check if it could be an issue with indexing or missing values
                    print("Warning: Unable to automatically determine correction for temperature")
                    # Suggest different correction based on the magnitude
                    if abs(temp_mean) > 10000:
                        print("Values are very large. This could be numerical indices misinterpreted as temperatures.")
                    elif abs(temp_mean) < 0.01:
                        print("Values are very small. This could be temperatures stored in unusual units.")
    
    # Check energy values
    energy_properties = ['Potential', 'Kinetic-En.', 'Total-Energy']
    for prop in energy_properties:
        if prop in property_data:
            values = property_data[prop]['values']
            mean_val = np.mean(values)
            unit = property_data[prop]['y_unit']
            
            print(f"\n{prop} Diagnostics:")
            print(f"Mean value: {mean_val:.6g} {unit}")
            
            # Typical energy ranges for water simulations in different units
            # Values are approximate for a TIP4P water system
            if unit.lower() in ['kj/mol', 'kj/mole']:
                if abs(mean_val) < 1 or abs(mean_val) > 1000000:
                    print(f"Warning: {prop} values seem to be in unusual range for kJ/mol")
                    if abs(mean_val) > 10000:
                        corrections[prop] = {'scale': 0.001}
                        print(f"Suggestion: Values may need scaling to conventional units (×0.001)")
            elif unit.lower() in ['kcal/mol', 'kcal/mole']:
                if abs(mean_val) < 0.1 or abs(mean_val) > 100000:
                    print(f"Warning: {prop} values seem to be in unusual range for kcal/mol")
            else:
                print(f"Warning: Unknown energy unit: {unit}")
                print("Values may need unit conversion depending on analysis requirements")
    
    # Check pressure
    if 'Pressure' in property_data:
        pres = property_data['Pressure']['values']
        pres_mean = np.mean(pres)
        unit = property_data['Pressure']['y_unit']
        
        print("\nPressure Diagnostics:")
        print(f"Mean pressure: {pres_mean:.6g} {unit}")
        
        # Check if values are in expected range based on unit
        if unit.lower() == 'bar':
            if abs(pres_mean) > 10000:
                print("Warning: Pressure values very high for bar units")
                corrections['Pressure'] = {'scale': 0.001}
                print("Suggestion: Values may need scaling (×0.001)")
        elif unit.lower() == 'atm':
            if abs(pres_mean) > 10000:
                print("Warning: Pressure values very high for atm units")
                corrections['Pressure'] = {'scale': 0.001}
                print("Suggestion: Values may need scaling (×0.001)")
    
    return corrections

def generate_diagnostic_plots(property_data, corrections, output_dir='../visualization'):
    """
    Generate diagnostic plots for the raw data and with suggested corrections
    
    Parameters:
    -----------
    property_data : dict
        Dictionary with extracted properties
    corrections : dict
        Dictionary with suggested correction factors
    output_dir : str
        Output directory for plots
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Temperature plot
    if 'Temperature' in property_data:
        plt.figure(figsize=(10, 6))
        
        # Plot raw data
        time = property_data['Temperature']['time'] / 1000  # ps to ns
        temp = property_data['Temperature']['values']
        plt.plot(time, temp, 'b-', label='Raw Temperature')
        
        # Plot corrected data if corrections exist
        if 'Temperature' in corrections:
            if isinstance(corrections['Temperature'], dict):
                if 'scale' in corrections['Temperature']:
                    scale = corrections['Temperature']['scale']
                    plt.plot(time, temp * scale, 'r-', 
                            label=f'Corrected (×{scale})')
                elif 'offset' in corrections['Temperature']:
                    offset = corrections['Temperature']['offset']
                    plt.plot(time, temp + offset, 'r-', 
                            label=f'Corrected (+{offset})')
            else:
                # Simple scalar factor
                plt.plot(time, temp * corrections['Temperature'], 'r-', 
                        label=f'Corrected (×{corrections["Temperature"]})')
        
        plt.xlabel('Time (ns)')
        plt.ylabel('Temperature (K)')
        plt.title('Temperature Diagnostic Plot')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.savefig(f"{output_dir}/temperature_diagnostic.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # Energy plots
    energy_properties = ['Potential', 'Kinetic-En.', 'Total-Energy']
    for prop in energy_properties:
        if prop in property_data:
            plt.figure(figsize=(10, 6))
            
            # Plot raw data
            time = property_data[prop]['time'] / 1000  # ps to ns
            values = property_data[prop]['values']
            unit = property_data[prop]['y_unit']
            plt.plot(time, values, 'b-', label=f'Raw {prop}')
            
            # Plot corrected data if corrections exist
            if prop in corrections:
                if isinstance(corrections[prop], dict):
                    if 'scale' in corrections[prop]:
                        scale = corrections[prop]['scale']
                        plt.plot(time, values * scale, 'r-', 
                                label=f'Corrected (×{scale})')
            
            plt.xlabel('Time (ns)')
            plt.ylabel(f'{prop} ({unit})')
            plt.title(f'{prop} Diagnostic Plot')
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.savefig(f"{output_dir}/{prop.lower().replace('-', '_')}_diagnostic.png", dpi=300, bbox_inches='tight')
            plt.close()
    
    # Create a summary plot with all properties in separate subplots
    n_props = len(property_data)
    if n_props > 0:
        fig, axes = plt.subplots(n_props, 1, figsize=(12, 4 * n_props), sharex=True)
        
        if n_props == 1:
            axes = [axes]  # Make sure axes is a list for consistent indexing
        
        for i, (prop, data) in enumerate(property_data.items()):
            time = data['time'] / 1000  # ps to ns
            values = data['values']
            unit = data['y_unit']
            
            axes[i].plot(time, values, 'b-', label='Raw')
            
            # Add corrected values if available
            if prop in corrections:
                if isinstance(corrections[prop], dict):
                    if 'scale' in corrections[prop]:
                        scale = corrections[prop]['scale']
                        axes[i].plot(time, values * scale, 'r-', 
                                    label=f'Corrected (×{scale})')
                    elif 'offset' in corrections[prop]:
                        offset = corrections[prop]['offset']
                        axes[i].plot(time, values + offset, 'r-', 
                                    label=f'Corrected (+{offset})')
                else:
                    # Simple scalar factor
                    axes[i].plot(time, values * corrections[prop], 'r-', 
                                label=f'Corrected (×{corrections[prop]})')
            
            axes[i].set_ylabel(f'{prop} ({unit})')
            axes[i].set_title(f'{prop}')
            axes[i].grid(True, alpha=0.3)
            axes[i].legend()
        
        # Add shared x-axis label to bottom subplot
        axes[-1].set_xlabel('Time (ns)')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/gromacs_data_diagnostics.png", dpi=300, bbox_inches='tight')
        plt.close()

def main():
    """Main function to run the script"""
    print("GROMACS Energy Data Inspection")
    
    # Extract raw data from GROMACS files
    if os.path.exists('../md.edr') and os.path.exists('../md.tpr'):
        print("\nExtracting raw data from GROMACS energy files...")
        property_data = extract_raw_gromacs_data()
        
        if property_data:
            # Analyze data and suggest corrections
            print("\nAnalyzing data for potential unit/scaling issues...")
            corrections = suggest_corrections(property_data)
            
            # Generate diagnostic plots
            print("\nGenerating diagnostic plots...")
            generate_diagnostic_plots(property_data, corrections)
            
            print("\nDiagnostic complete. Check visualization folder for plots.")
        else:
            print("No data extracted from GROMACS files.")
    else:
        print("GROMACS files (md.edr, md.tpr) not found.")

if __name__ == "__main__":
    main() 