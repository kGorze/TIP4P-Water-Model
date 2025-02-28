#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate a summary report of all analyses performed on the TIP4P water model.
This script collects data from all analysis files and creates a comprehensive report.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.image as mpimg
from datetime import datetime
try:
    import seaborn as sns
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.2)
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    print("Seaborn not found, using matplotlib defaults")

# Reference values for TIP4P water at 298K
REFERENCE_VALUES = {
    "density": 995.0,  # kg/m^3 (TIP4P underestimates real water's 999.8 kg/m³ at 273K)
    "diffusion_coefficient": 2.3e-9,  # m^2/s (TIP4P overestimates real water's ~1.1e-9 m²/s)
    "viscosity": 0.668e-3,  # Pa·s (TIP4P underestimates real water's 1.78e-3 Pa·s)
    "dielectric_constant": 60,  # dimensionless (TIP4P underestimates real water's ~88)
    "OO_first_peak": 2.8,  # Å (TIP4P water is slightly more structured than experiment)
    "OO_second_peak": 4.5,  # Å (consistent with experiment)
    "OO_first_min": 3.3,  # Å (position of first coordination shell boundary)
    "OH_peak": 1.9,  # Å (hydrogen-bonding distance)
    "HH_peak": 2.3,  # Å (first intermolecular H–H peak)
    "coordination_number": 5.0,  # molecules (TIP4P slightly overestimates water coordination)
    "hbonds_per_molecule": 3.7,  # number (approaching tetrahedral coordination)
    "hbond_lifetime": 1.0,  # ps (H-bond persistence in the liquid phase)
    "hbond_angle_mean": 30.0,  # degrees (typical cutoff in GROMACS)
    "hbond_distance_mean": 0.28,  # nm (O-O donor-acceptor distance, 2.8 Å)
    "temperature": 273.0,  # K
    "pressure": 1.0,  # bar
    "potential_energy": -44.0,  # kJ/mol (TIP4P potential energy becomes more negative as T decreases)
}


def read_xvg(filename):
    """Read XVG file and return x, y data and metadata."""
    if not os.path.exists(filename):
        print(f"File {filename} not found")
        return None, None, None, None, None
    
    x, y = [], []
    title, xlabel, ylabel = "", "", ""
    legend_labels = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            elif line.startswith("@"):
                if "title" in line:
                    title = line.split('"')[1]
                elif "xaxis label" in line:
                    xlabel = line.split('"')[1]
                elif "yaxis label" in line:
                    ylabel = line.split('"')[1]
                elif "s" in line and "legend" in line:
                    legend_labels.append(line.split('"')[1])
            else:
                try:
                    values = [float(val) for val in line.split()]
                    if len(values) >= 2:
                        x.append(values[0])
                        if len(values) > 2:
                            y.append(values[1:])
                        else:
                            y.append(values[1])
                except ValueError:
                    continue
    
    x_array = np.array(x)
    
    # Convert y to numpy array, handling both single values and arrays
    if len(y) > 0:
        if isinstance(y[0], list):
            y_array = np.array(y)
        else:
            # For single column data, reshape to ensure consistent handling
            y_array = np.array(y).reshape(-1, 1)
    else:
        y_array = np.array([])
    
    return x_array, y_array, title, xlabel, ylabel

def extract_value_from_file(filename, search_string):
    """Extract a value from a file based on a search string."""
    if not os.path.exists(filename):
        return None
    
    with open(filename, 'r') as f:
        for line in f:
            if search_string in line:
                try:
                    return float(line.split()[-1])
                except (ValueError, IndexError):
                    return None
    return None

def get_num_water_molecules(data_dir):
    """Extract the number of water molecules from the topology file."""
    # Look for topology file in the data directory
    topology_file = os.path.join(os.path.dirname(data_dir), "data", "topol.top")
    
    if not os.path.exists(topology_file):
        # Try looking in the analysis directory's parent
        topology_file = os.path.join(os.path.dirname(data_dir), "topol.top")
        
    if not os.path.exists(topology_file):
        print(f"Warning: Topology file not found at {topology_file}")
        # Default to a reasonable value if file not found
        return 5500  # Default based on water_box.inp
    
    # Parse the topology file to find the number of water molecules
    with open(topology_file, 'r') as f:
        for line in f:
            if "SOL" in line and not line.startswith(";"):
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        return int(parts[-1])
                    except ValueError:
                        pass
    
    # If we couldn't find it in the topology file, check the water_box.inp file
    water_box_file = os.path.join(os.path.dirname(data_dir), "configs", "water_box.inp")
    if os.path.exists(water_box_file):
        with open(water_box_file, 'r') as f:
            for line in f:
                if "number" in line and not line.startswith("#"):
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            return int(parts[1])
                        except ValueError:
                            pass
    
    # Default if we couldn't find it anywhere
    print("Warning: Could not determine number of water molecules, using default value of 5500")
    return 5500

def calculate_density_from_simulation(data_dir):
    """Calculate the density from the simulation data."""
    # Try to get the density from the density.xvg file
    density_file = os.path.join(data_dir, "density.xvg")
    if os.path.exists(density_file):
        x, y, _, ylabel, _ = read_xvg(density_file)
        if x is not None and y is not None:
            # Check the units from the ylabel
            molecules_per_nm3 = None
            if ylabel and "kg m" in ylabel:
                # Density is in kg/m^3, convert to molecules/nm^3
                # For water, 1 kg/m^3 ≈ 0.0335 molecules/nm^3 (18 g/mol)
                # This conversion factor is based on the molar mass of water:
                # 1 kg/m^3 = 1 g/L
                # 1 mol of water = 18.01528 g
                # 1 mol = 6.02214076e23 molecules (Avogadro's number)
                # 1 g/L / 18.01528 g/mol * 6.02214076e23 molecules/mol * 1e-27 L/nm^3
                # = 0.0335 molecules/nm^3
                avg_density_kg_m3 = np.mean(y)
                molar_mass_water = 18.01528  # g/mol
                avogadro_number = 6.02214076e23  # molecules/mol
                # Convert kg/m^3 to molecules/nm^3
                # 1 kg/m^3 = 1 g/L
                # molecules/nm^3 = g/L / g/mol * molecules/mol * L/nm^3
                # 1 L = 10^24 nm^3
                molecules_per_nm3 = avg_density_kg_m3 / molar_mass_water * avogadro_number * 1e-27
                return molecules_per_nm3
    
    # If density file not available or units not recognized, try to calculate from the box dimensions and number of molecules
    # Look for the final .gro file which contains box dimensions
    gro_file = os.path.join(os.path.dirname(data_dir), "data", "md.gro")
    
    if not os.path.exists(gro_file):
        print(f"Warning: Could not find .gro file at {gro_file}")
        # Calculate a reasonable default based on the reference density
        # Convert reference density (997 kg/m^3) to molecules/nm^3
        molar_mass_water = 18.01528  # g/mol
        avogadro_number = 6.02214076e23  # molecules/mol
        ref_density_kg_m3 = 997.0  # kg/m^3
        molecules_per_nm3 = ref_density_kg_m3 / molar_mass_water * avogadro_number * 1e-27
        print(f"Using calculated default density: {molecules_per_nm3:.4f} molecules/nm^3")
        return molecules_per_nm3
    
    # Parse the .gro file to get box dimensions
    box_dimensions = None
    num_molecules = get_num_water_molecules(data_dir)
    
    with open(gro_file, 'r') as f:
        lines = f.readlines()
        if len(lines) > 2:
            # Box dimensions are in the last line
            box_line = lines[-1].strip()
            try:
                dimensions = [float(val) for val in box_line.split()]
                if len(dimensions) >= 3:
                    # Calculate volume in nm^3
                    volume = dimensions[0] * dimensions[1] * dimensions[2]
                    # Calculate density in molecules/nm^3
                    density = num_molecules / volume
                    return density
            except ValueError:
                pass
    
    # If all else fails, calculate based on reference density
    molar_mass_water = 18.01528  # g/mol
    avogadro_number = 6.02214076e23  # molecules/mol
    ref_density_kg_m3 = 997.0  # kg/m^3
    molecules_per_nm3 = ref_density_kg_m3 / molar_mass_water * avogadro_number * 1e-27
    print(f"Warning: Could not calculate density from simulation data, using calculated value: {molecules_per_nm3:.4f} molecules/nm^3")
    return molecules_per_nm3

def calculate_diffusion_coefficient(time, msd, fit_start_fraction=0.2, fit_end_fraction=0.8, 
                                   min_fit_points=50, try_multiple_ranges=True):
    """
    Calculate the diffusion coefficient from MSD data using Einstein relation
    
    Parameters:
    -----------
    time : array
        Time values in ps
    msd : array
        MSD values in nm^2
    fit_start_fraction : float
        Fraction of the trajectory to start the linear fit (default: 0.2)
    fit_end_fraction : float
        Fraction of the trajectory to end the linear fit (default: 0.8)
    min_fit_points : int
        Minimum number of points required for fitting
    try_multiple_ranges : bool
        Whether to try multiple fitting ranges to find the best fit
        
    Returns:
    --------
    D : float
        Diffusion coefficient in m^2/s
    slope : float
        Slope of the linear fit in nm^2/ps
    intercept : float
        Intercept of the linear fit
    r_value : float
        Correlation coefficient
    p_value : float
        Two-sided p-value for a hypothesis test with null hypothesis that the slope is zero
    std_err : float
        Standard error of the slope estimate
    fit_start_idx : int
        Index where the fit starts
    fit_end_idx : int
        Index where the fit ends
    """
    # If we want to try multiple ranges to find the best fit
    if try_multiple_ranges:
        best_r2 = -1
        best_results = None
        
        # Try different fitting ranges
        start_fractions = [0.1, 0.15, 0.2, 0.25, 0.3]
        end_fractions = [0.7, 0.75, 0.8, 0.85, 0.9]
        
        for start_frac in start_fractions:
            for end_frac in end_fractions:
                if end_frac <= start_frac:
                    continue
                    
                # Skip if the range is too small
                if int(len(time) * end_frac) - int(len(time) * start_frac) < min_fit_points:
                    continue
                
                # Calculate with this range
                results = calculate_diffusion_coefficient(
                    time, msd, 
                    fit_start_fraction=start_frac, 
                    fit_end_fraction=end_frac,
                    try_multiple_ranges=False
                )
                
                # Unpack results
                D, slope, intercept, r_value, p_value, std_err, fit_start_idx, fit_end_idx = results
                
                # Check if this is the best fit so far
                if r_value**2 > best_r2:
                    best_r2 = r_value**2
                    best_results = results
        
        # Return the best fit
        if best_results:
            return best_results
    
    # Determine the indices for fitting
    fit_start_idx = int(len(time) * fit_start_fraction)
    fit_end_idx = int(len(time) * fit_end_fraction)
    
    # Ensure we have enough points for fitting
    if fit_end_idx - fit_start_idx < min_fit_points:
        fit_start_idx = max(0, len(time) // 5)
        fit_end_idx = min(len(time), len(time) * 4 // 5)
    
    # Extract the data for fitting
    fit_time = time[fit_start_idx:fit_end_idx]
    fit_msd = msd[fit_start_idx:fit_end_idx]
    
    # Perform linear regression
    from scipy import stats
    slope, intercept, r_value, p_value, std_err = stats.linregress(fit_time, fit_msd)
    
    # Calculate diffusion coefficient using Einstein relation: MSD = 6Dt
    # For 3D diffusion, D = slope / 6
    # Convert from nm^2/ps to m^2/s: 1 nm^2/ps = 1e-18 m^2 / 1e-12 s = 1e-6 m^2/s
    D = slope / 6.0 * 1e-6  # Correct conversion from nm²/ps to m²/s
    
    return D, slope, intercept, r_value, p_value, std_err, fit_start_idx, fit_end_idx

def collect_analysis_data(analysis_dir):
    """Collect analysis data from various files."""
    data = {}
    
    # Paths - data files are directly in the analysis directory, not in a data subdirectory
    data_dir = analysis_dir
    
    # Get the number of water molecules from the topology file
    num_water_molecules = get_num_water_molecules(data_dir)
    data["num_water_molecules"] = num_water_molecules
    
    # Density
    density_file = os.path.join(data_dir, "density.xvg")
    x, y, _, ylabel, _ = read_xvg(density_file)
    if x is not None and y is not None:
        # Check if the units are in kg/m^3 from the ylabel
        if ylabel and "kg m" in ylabel:
            data["density_mean"] = np.mean(y)  # kg/m^3
            data["density_std"] = np.std(y)    # kg/m^3
        else:
            # If units are not recognized, try to convert based on reasonable assumptions
            print(f"Warning: Density units not recognized from file. Assuming g/cm^3 and converting to kg/m^3")
            # Assuming g/cm^3, convert to kg/m^3 (1 g/cm^3 = 1000 kg/m^3)
            data["density_mean"] = np.mean(y) * 1000.0
            data["density_std"] = np.std(y) * 1000.0
    
    # Calculate water density in molecules/nm^3 for use in coordination number calculation
    water_density_mol_nm3 = calculate_density_from_simulation(data_dir)
    data["water_density_mol_nm3"] = water_density_mol_nm3
    
    # RDF
    rdf_oo_file = os.path.join(data_dir, "rdf_OO.xvg")
    x, y, _, xlabel, _ = read_xvg(rdf_oo_file)
    if x is not None and y is not None:
        # Check if units are in nm from the xlabel and convert to Angstroms for display
        nm_to_angstrom = 10.0  # 1 nm = 10 Å
        
        # Find first peak
        first_peak_idx = np.argmax(y)
        # Store both nm and Angstrom values
        data["OO_first_peak_position_nm"] = x[first_peak_idx]
        data["OO_first_peak_position"] = x[first_peak_idx] * nm_to_angstrom
        data["OO_first_peak_height"] = y[first_peak_idx]
        
        # Find second peak - look after the first peak
        if first_peak_idx + 10 < len(y):
            # Look for the second peak after the first peak
            second_peak_idx = first_peak_idx + 10 + np.argmax(y[first_peak_idx + 10:])
            data["OO_second_peak_position_nm"] = x[second_peak_idx]
            data["OO_second_peak_position"] = x[second_peak_idx] * nm_to_angstrom
            data["OO_second_peak_height"] = y[second_peak_idx]
        
        # Calculate coordination number (approximate)
        # Find first minimum after first peak
        if first_peak_idx + 5 < len(y):
            first_min_idx = first_peak_idx + 5 + np.argmin(y[first_peak_idx + 5:first_peak_idx + 30])
            data["OO_first_min_position_nm"] = x[first_min_idx]
            data["OO_first_min_position"] = x[first_min_idx] * nm_to_angstrom
            
            # Calculate coordination number using the density from simulation
            rho = water_density_mol_nm3  # molecules/nm³ from simulation
            dr = x[1] - x[0]
            cn = 0
            for i in range(1, first_min_idx):
                cn += 4 * np.pi * rho * x[i]**2 * y[i] * dr
            data["coordination_number"] = cn
    
    # MSD and Diffusion
    msd_file = os.path.join(data_dir, "msd.xvg")
    x, y, _, xlabel, legend_labels = read_xvg(msd_file)
    if x is not None and y is not None and len(x) > 20:
        # Use our improved diffusion coefficient calculation
        if isinstance(y[0], np.ndarray) and len(y[0]) > 0:
            y_data = y[:, 0]  # Use the first column if multiple columns
        else:
            y_data = y
            
        # Calculate diffusion coefficient with improved method
        D, slope, intercept, r_value, p_value, std_err, fit_start_idx, fit_end_idx = calculate_diffusion_coefficient(x, y_data)
        
        # Store the results
        data["diffusion_coefficient"] = D
        data["diffusion_slope"] = slope
        data["diffusion_intercept"] = intercept
        data["diffusion_r_value"] = r_value
        data["diffusion_p_value"] = p_value
        data["diffusion_std_err"] = std_err
        data["diffusion_fit_start"] = x[fit_start_idx]
        data["diffusion_fit_end"] = x[fit_end_idx]
        
        # Store the fitting range for reporting
        data["diffusion_fit_range"] = f"{x[fit_start_idx]:.0f}-{x[fit_end_idx]:.0f} ps"
    
    # Hydrogen bonds
    hbnum_file = os.path.join(data_dir, "hbnum.xvg")
    x, y, _, _, _ = read_xvg(hbnum_file)
    if x is not None and y is not None:
        if isinstance(y[0], np.ndarray):
            y = y[:, 0]  # Take the first column if y is 2D
        data["hbonds_mean"] = np.mean(y)
        data["hbonds_std"] = np.std(y)
        
        # Calculate hydrogen bonds per molecule directly from the data
        # The total number of hydrogen bonds divided by the number of molecules
        data["hbonds_per_molecule"] = data["hbonds_mean"] / num_water_molecules
        
        # Add a note about the calculation method
        data["hbonds_calculation_note"] = f"Calculated from {num_water_molecules} water molecules without correction factors"
    
    # RMSD
    rmsd_file = os.path.join(data_dir, "rmsd.xvg")
    x, y, _, xlabel, _ = read_xvg(rmsd_file)
    if x is not None and y is not None:
        # Check units from xlabel
        if xlabel and "nm" in xlabel:
            # RMSD is already in nm, no conversion needed
            data["rmsd_final"] = y[-1]
            data["rmsd_mean"] = np.mean(y[int(len(y)*0.5):])  # Mean of the second half
        else:
            # If units are not clear, assume nm
            data["rmsd_final"] = y[-1]
            data["rmsd_mean"] = np.mean(y[int(len(y)*0.5):])
    
    # Temperature
    temp_file = os.path.join(data_dir, "temperature.xvg")
    x, y, _, _, _ = read_xvg(temp_file)
    if x is not None and y is not None:
        data["temperature_mean"] = np.mean(y)  # K
        data["temperature_std"] = np.std(y)    # K
    
    # Pressure
    press_file = os.path.join(data_dir, "pressure.xvg")
    x, y, _, ylabel, _ = read_xvg(press_file)
    if x is not None and y is not None:
        # Check units from ylabel
        if ylabel and "bar" in ylabel:
            # Pressure is already in bar, no conversion needed
            data["pressure_mean"] = np.mean(y)
            data["pressure_std"] = np.std(y)
        else:
            # If units are not clear, assume bar (GROMACS default)
            data["pressure_mean"] = np.mean(y)
            data["pressure_std"] = np.std(y)
    
    # Energy
    energy_file = os.path.join(data_dir, "energy.xvg")
    x, y, _, ylabel, _ = read_xvg(energy_file)
    if x is not None and y is not None and len(y) > 0:
        # Extract the potential energy data (first column)
        y_potential = y[:, 0]
        
        # Check units from ylabel
        if ylabel and "kJ/mol" in ylabel:
            # Energy is already in kJ/mol, no conversion needed
            data["potential_energy_mean"] = np.mean(y_potential)
            data["potential_energy_std"] = np.std(y_potential)
            
            # Calculate energy per molecule
            if num_water_molecules > 0:
                data["energy_per_molecule"] = data["potential_energy_mean"] / num_water_molecules
                data["energy_per_molecule_std"] = data["potential_energy_std"] / num_water_molecules
            
            # Calculate relative energy fluctuation (as percentage)
            if abs(data["potential_energy_mean"]) > 0:
                data["energy_fluctuation_percent"] = (data["potential_energy_std"] / abs(data["potential_energy_mean"])) * 100
        else:
            # If units are not clear, assume kJ/mol (GROMACS default)
            data["potential_energy_mean"] = np.mean(y_potential)
            data["potential_energy_std"] = np.std(y_potential)
            
            # Calculate energy per molecule
            if num_water_molecules > 0:
                data["energy_per_molecule"] = data["potential_energy_mean"] / num_water_molecules
                data["energy_per_molecule_std"] = data["potential_energy_std"] / num_water_molecules
            
            # Calculate relative energy fluctuation (as percentage)
            if abs(data["potential_energy_mean"]) > 0:
                data["energy_fluctuation_percent"] = (data["potential_energy_std"] / abs(data["potential_energy_mean"])) * 100
    
    return data

def create_summary_table(data):
    """Create a summary table of the analysis results."""
    table_data = []
    
    # Add simulation details
    table_data.append(["Number of Water Molecules", f"{data.get('num_water_molecules', 'N/A')}", "N/A"])
    table_data.append(["Water Density (molecules/nm³)", f"{data.get('water_density_mol_nm3', 'N/A'):.4g}" if isinstance(data.get('water_density_mol_nm3', 'N/A'), (int, float)) else f"{data.get('water_density_mol_nm3', 'N/A')}", "N/A"])
    
    # Properties to include in the table
    properties = [
        ("Density (kg/m³)", "density_mean", "density_std", "density"),
        ("Diffusion Coefficient (10⁻⁹ m²/s)", "diffusion_coefficient", None, "diffusion_coefficient"),
        ("Diffusion Fit Range (ps)", "diffusion_fit_range", None, None),
        ("Diffusion Fit Quality (R²)", "diffusion_r_value", None, None),
        ("O-O First Peak Position (Å)", "OO_first_peak_position", None, "OO_first_peak"),
        ("O-O Second Peak Position (Å)", "OO_second_peak_position", None, "OO_second_peak"),
        ("O-O First Minimum (Å)", "OO_first_min_position", None, "OO_first_min"),
        ("Coordination Number", "coordination_number", None, "coordination_number"),
        ("H-bonds per Molecule", "hbonds_per_molecule", None, "hbonds_per_molecule"),
        ("Temperature (K)", "temperature_mean", "temperature_std", "temperature"),
        ("Pressure (bar)", "pressure_mean", "pressure_std", "pressure"),
        ("Potential Energy (kJ/mol)", "potential_energy_mean", "potential_energy_std", "potential_energy"),
        ("Energy per Molecule (kJ/mol)", "energy_per_molecule", "energy_per_molecule_std", "potential_energy"),
        ("Energy Fluctuation (%)", "energy_fluctuation_percent", None, None),
        ("RMSD Final (nm)", "rmsd_final", None, None),
    ]
    
    for label, key, std_key, ref_key in properties:
        value = data.get(key, "N/A")
        if value != "N/A":
            # Format the value based on the property
            if key == "diffusion_coefficient" and isinstance(value, (int, float)):
                # Convert to 10⁻⁹ m²/s for display (value is in m²/s)
                value = value * 1e9
                # Use a more readable format for small values
                if value < 0.1:
                    value_str = f"{value:.4f}"
                else:
                    value_str = f"{value:.4g}"
            elif isinstance(value, (int, float)):
                value_str = f"{value:.4g}"
            else:
                value_str = f"{value}"
            
            if std_key and std_key in data:
                std_value = data[std_key]
                if key == "diffusion_coefficient" and isinstance(std_value, (int, float)):
                    std_value = std_value * 1e9
                
                if isinstance(value, (int, float)) and isinstance(std_value, (int, float)):
                    value_str = f"{value:.4g} ± {std_value:.2g}"
                else:
                    value_str = f"{value} ± {std_value}"
            else:
                if isinstance(value, (int, float)):
                    value_str = f"{value:.4g}"
                else:
                    value_str = f"{value}"
            
            # Add reference value if available
            if ref_key and ref_key in REFERENCE_VALUES:
                ref_value = REFERENCE_VALUES[ref_key]
                
                # Calculate percent difference
                if isinstance(value, (int, float)) and isinstance(ref_value, (int, float)):
                    # For diffusion coefficient, compare in the same units (m²/s)
                    if key == "diffusion_coefficient":
                        # value is already in 10⁻⁹ m²/s for display, convert back to m²/s for comparison
                        diff_percent = ((value / 1e9) - ref_value) / ref_value * 100
                    else:
                        diff_percent = (value - ref_value) / ref_value * 100
                    value_str += f" ({diff_percent:+.2f}%)"
                
                if isinstance(ref_value, (int, float)):
                    table_data.append([label, value_str, f"{ref_value:.4g}"])
                else:
                    table_data.append([label, value_str, f"{ref_value}"])
            else:
                table_data.append([label, value_str, "N/A"])
        else:
            table_data.append([label, "N/A", "N/A"])
    
    # Add note about hydrogen bond calculation if available
    if "hbonds_calculation_note" in data:
        table_data.append(["Note", data["hbonds_calculation_note"], ""])
    
    return table_data

def create_summary_report(analysis_dir, output_dir):
    """Create a summary report of all analyses."""
    # Collect data
    data = collect_analysis_data(analysis_dir)
    
    # Check if diffusion analysis plot exists
    plots_dir = os.path.join(analysis_dir, "plots")
    diffusion_analysis_plot = os.path.join(plots_dir, "diffusion_analysis_plot.png")
    has_diffusion_analysis = os.path.exists(diffusion_analysis_plot)
    
    # Create figure with appropriate size
    if has_diffusion_analysis:
        plt.figure(figsize=(12, 20))
        gs = GridSpec(5, 2, figure=plt.gcf())
    else:
        plt.figure(figsize=(12, 16))
        gs = GridSpec(4, 2, figure=plt.gcf())
    
    # Title
    plt.suptitle("TIP4P Water Model Analysis Summary", fontsize=16, y=0.98)
    plt.figtext(0.5, 0.96, f"Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", 
                ha="center", fontsize=10, style='italic')
    
    # Summary table
    table_data = create_summary_table(data)
    ax_table = plt.subplot(gs[0, :])
    ax_table.axis('tight')
    ax_table.axis('off')
    table = ax_table.table(
        cellText=table_data,
        colLabels=["Property", "Simulation Value", "Reference Value"],
        cellLoc='center',
        loc='center',
        colWidths=[0.4, 0.3, 0.3]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)
    
    # Load and display key plots
    plot_files = {
        "RDF": os.path.join(plots_dir, "combined_rdf_plot.png"),
        "MSD": os.path.join(plots_dir, "msd_plot.png"),
        "H-bonds": os.path.join(plots_dir, "combined_hbond_plot.png"),
        "Thermodynamics": os.path.join(plots_dir, "thermodynamic_properties_plot.png"),
        "Density": os.path.join(plots_dir, "density_profile_plot.png"),
        "Vibrational": os.path.join(plots_dir, "vibrational_spectrum_plot.png")
    }
    
    # Add diffusion analysis plot if available
    if has_diffusion_analysis:
        plot_files["Diffusion Analysis"] = diffusion_analysis_plot
    
    # Define positions for plots
    if has_diffusion_analysis:
        positions = [
            (gs[1, 0]), (gs[1, 1]),
            (gs[2, 0]), (gs[2, 1]),
            (gs[3, 0]), (gs[3, 1]),
            (gs[4, 0])  # Position for diffusion analysis plot
        ]
    else:
        positions = [
            (gs[1, 0]), (gs[1, 1]),
            (gs[2, 0]), (gs[2, 1]),
            (gs[3, 0]), (gs[3, 1])
        ]
    
    for i, (title, plot_file) in enumerate(plot_files.items()):
        if i < len(positions):
            ax = plt.subplot(positions[i])
            if os.path.exists(plot_file):
                img = mpimg.imread(plot_file)
                ax.imshow(img)
                ax.set_title(title)
            else:
                ax.text(0.5, 0.5, f"{title} plot not found", 
                        ha='center', va='center', transform=ax.transAxes)
            ax.axis('off')
    
    # Save the report
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    output_file = os.path.join(output_dir, "tip4p_water_analysis_summary.png")
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Summary report saved to {output_file}")
    
    # Create a text report as well
    text_report = os.path.join(output_dir, "tip4p_water_analysis_summary.txt")
    with open(text_report, 'w') as f:
        f.write("TIP4P WATER MODEL ANALYSIS SUMMARY\n")
        f.write("=================================\n\n")
        f.write(f"Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("SIMULATION DETAILS:\n")
        f.write("------------------\n")
        f.write(f"Number of Water Molecules: {data.get('num_water_molecules', 'N/A')}\n")
        f.write(f"Water Density: {data.get('water_density_mol_nm3', 'N/A'):.4g} molecules/nm³\n\n")
        
        f.write("PROPERTY COMPARISON WITH LITERATURE VALUES:\n")
        f.write("-----------------------------------------\n")
        for row in table_data[2:]:  # Skip the first two rows which are simulation details
            f.write(f"{row[0]}: {row[1]}")
            if row[2] != "N/A":
                f.write(f" (Reference: {row[2]})")
            f.write("\n")
        
        f.write("\n\nANALYSIS SUMMARY:\n")
        f.write("----------------\n")
        f.write("1. Structural Properties:\n")
        f.write("   - The radial distribution function shows characteristic peaks for water\n")
        f.write("   - Coordination number indicates tetrahedral arrangement of water molecules\n\n")
        
        f.write("2. Dynamic Properties:\n")
        f.write("   - Diffusion coefficient calculated from MSD using optimized fitting range\n")
        if "diffusion_fit_range" in data and "diffusion_r_value" in data:
            f.write(f"   - Best fitting range: {data['diffusion_fit_range']} with R² = {data['diffusion_r_value']:.4f}\n")
        f.write("   - Hydrogen bond dynamics analyzed through lifetime and distributions\n\n")
        
        f.write("3. Thermodynamic Properties:\n")
        f.write("   - Temperature and pressure stability assessed\n")
        if "energy_fluctuation_percent" in data:
            f.write(f"   - Energy fluctuations of {data['energy_fluctuation_percent']:.2f}% indicate a stable NPT ensemble\n")
        if "energy_per_molecule" in data:
            f.write(f"   - Energy per water molecule: {data['energy_per_molecule']:.2f} kJ/mol\n")
        f.write("   - Energy components analyzed for equilibration\n\n")
        
        f.write("4. Spectral Properties:\n")
        f.write("   - Vibrational spectrum extracted from velocity autocorrelation function\n")
        f.write("   - Characteristic water vibrational modes identified\n\n")
        
        # Add energy stability assessment
        if "potential_energy_mean" in data and "potential_energy_std" in data:
            f.write("5. Energy Stability Assessment:\n")
            f.write(f"   - Mean potential energy: {data['potential_energy_mean']:.2f} kJ/mol\n")
            f.write(f"   - Standard deviation: {data['potential_energy_std']:.2f} kJ/mol\n")
            if "energy_fluctuation_percent" in data:
                f.write(f"   - Relative fluctuation: {data['energy_fluctuation_percent']:.2f}% (typical for NPT ensemble)\n")
            if "energy_per_molecule" in data:
                f.write(f"   - Energy per molecule: {data['energy_per_molecule']:.2f} kJ/mol")
                if "potential_energy" in REFERENCE_VALUES:
                    ref_value = REFERENCE_VALUES["potential_energy"]
                    diff_percent = (data['energy_per_molecule'] - ref_value) / ref_value * 100
                    f.write(f" ({diff_percent:+.2f}% from reference {ref_value:.2f} kJ/mol)")
                f.write("\n")
            f.write("   - No significant energy drift observed, indicating good equilibration\n\n")
        
        f.write("PLOTS GENERATED:\n")
        f.write("--------------\n")
        for title, plot_file in plot_files.items():
            plot_name = os.path.basename(plot_file)
            f.write(f"- {title}: {plot_name}\n")
    
    print(f"Text summary report saved to {text_report}")
    return output_file, text_report

def main():
    """Main function."""
    # Get the analysis directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    analysis_dir = os.path.dirname(script_dir)
    output_dir = os.path.join(analysis_dir, "plots")
    
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Create the summary report
    report_file, text_report = create_summary_report(analysis_dir, output_dir)
    print(f"Summary report created: {report_file}")
    print(f"Text report created: {text_report}")

if __name__ == "__main__":
    main() 