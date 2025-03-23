#!/usr/bin/env python3
"""
Phase Transition Markers Visualization

This script creates visualizations for the quantitative markers of phase transitions,
including:
1. First derivatives of free energy with respect to temperature and pressure
2. Heat capacity peaks at transition points
3. Density changes across temperature range
4. Tetrahedral order parameter distribution at different phases

Usage:
    python phase_transition_markers.py
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter
from scipy import stats
import matplotlib.gridspec as gridspec

# Constants
k_B_kcal = 0.0019872041  # kcal/(mol·K)

# Function to load data from energy files
def load_energy_data(edr_file='../md.edr', tpr_file='../md.tpr'):
    """Extract temperature, pressure and energy data from a GROMACS energy file"""
    import subprocess
    
    # Create input for gmx energy
    with open('temp_input.txt', 'w') as f:
        f.write("Temperature\nPressure\nPotential\nKinetic-En.\nTotal-Energy\n0\n")
    
    # Run gmx energy command
    cmd = f"gmx energy -f {edr_file} -s {tpr_file} -o energy_data.xvg < temp_input.txt"
    try:
        process = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running gmx energy: {e}")
        print(f"Error output: {e.stderr}")
        sys.exit(1)
    finally:
        if os.path.exists('temp_input.txt'):
            os.remove('temp_input.txt')
    
    # Read the generated xvg file
    if not os.path.exists('energy_data.xvg'):
        raise FileNotFoundError("energy_data.xvg was not created!")
    
    try:
        # Skip comment lines (those starting with # or @)
        with open('energy_data.xvg', 'r') as f:
            lines = f.readlines()
        data_lines = [l for l in lines if not (l.startswith('#') or l.startswith('@'))]
        
        # Convert to DataFrame
        data = pd.DataFrame([line.split() for line in data_lines], 
                          columns=['Time', 'Temperature', 'Pressure', 'Potential', 'Kinetic', 'Total'],
                          dtype=float)
        
        # Convert pressure from bar to MPa and time to ns
        data['Pressure'] = data['Pressure'] * 0.1  # bar to MPa
        data['Time'] = data['Time'] / 1000  # ps to ns
        
        return data
    except Exception as e:
        print(f"Error processing energy_data.xvg: {e}")
        sys.exit(1)
    finally:
        if os.path.exists('energy_data.xvg') and os.path.exists('energy_data.xvg'):
            os.remove('energy_data.xvg')

# Generate synthetic data for demonstration if real data is not available
def generate_synthetic_data():
    """Generate synthetic data for phase transition demonstration"""
    print("Generating synthetic data for demonstration...")
    
    # Temperature range covering ice, liquid, and gas phases
    temp = np.linspace(250, 350, 500)  # K
    
    # Free energy data with phase transitions
    # Simulate transitions at around 273K and 373K (water freezing and boiling)
    f_energy = 10 * np.exp(-0.01 * (temp - 273)**2) + 5 * np.exp(-0.02 * (temp - 373)**2)
    
    # Add some noise
    f_energy += np.random.normal(0, 0.2, size=len(temp))
    
    # Compute free energy derivative with respect to temperature
    dF_dT = np.gradient(f_energy, temp)
    
    # Pressure data (MPa)
    pressure = np.linspace(0.1, 5, 500)
    
    # Free energy with respect to pressure
    p_f_energy = 8 * np.exp(-1 * (pressure - 0.5)**2) + 4 * np.exp(-2 * (pressure - 3)**2)
    p_f_energy += np.random.normal(0, 0.1, size=len(pressure))
    
    # Derivative with respect to pressure
    dF_dP = np.gradient(p_f_energy, pressure)
    
    # Heat capacity data
    # Sharp peak at ice-liquid transition (273K)
    heat_capacity = 40 + 45 * np.exp(-0.5 * (temp - 273)**2 / 4)
    heat_capacity += np.random.normal(0, 1, size=len(temp))
    
    # Density data
    # Ice (less dense) to liquid water (more dense) transition
    density = np.zeros_like(temp)
    ice_mask = temp < 273
    liquid_mask = (temp >= 273) & (temp < 373)
    gas_mask = temp >= 373
    
    density[ice_mask] = 0.92 - 0.0001 * (273 - temp[ice_mask])
    density[liquid_mask] = 0.997 - 0.0003 * (temp[liquid_mask] - 273)
    density[gas_mask] = 0.1 + 0.001 * (373 - temp[gas_mask])
    
    # Add transition noise
    transition_noise = np.zeros_like(temp)
    transition_mask = (temp > 270) & (temp < 276)
    transition_noise[transition_mask] = np.random.normal(0, 0.01, size=np.sum(transition_mask))
    density += transition_noise
    
    # Tetrahedral order parameter
    # High for ice, medium for liquid, low for gas
    order_param = np.zeros_like(temp)
    order_param[ice_mask] = 0.85 - 0.05 * np.random.random(size=np.sum(ice_mask))
    order_param[liquid_mask] = 0.6 - 0.1 * np.random.random(size=np.sum(liquid_mask))
    order_param[gas_mask] = 0.2 - 0.2 * np.random.random(size=np.sum(gas_mask))
    
    # Add transition noise
    order_transition_noise = np.zeros_like(temp)
    order_transition_noise[transition_mask] = np.random.normal(0, 0.05, size=np.sum(transition_mask))
    order_param += order_transition_noise
    
    return {
        'temperature': temp,
        'pressure': pressure,
        'free_energy_T': f_energy,
        'free_energy_P': p_f_energy,
        'dF_dT': dF_dT,
        'dF_dP': dF_dP,
        'heat_capacity': heat_capacity,
        'density': density,
        'order_parameter': order_param
    }

def calculate_tetrahedral_order(positions, k=4):
    """
    Calculates tetrahedral parameter for each oxygen atom.
    
    Parameters:
    -----------
    positions : numpy.ndarray
        Array of atom positions
    k : int
        Number of nearest neighbors to consider (typically 4 for water)
        
    Returns:
    --------
    numpy.ndarray
        Array of tetrahedral order parameters
    """
    from scipy.spatial import cKDTree
    
    n_atoms = len(positions)
    q_values = np.zeros(n_atoms)
    
    kdtree = cKDTree(positions)
    
    for i in range(n_atoms):
        # k+1 nearest neighbors (first is the atom itself)
        distances, indices = kdtree.query(positions[i], k=k+1)
        neighbors = indices[1:k+1]
        
        vectors = positions[neighbors] - positions[i]
        norms = np.linalg.norm(vectors, axis=1)
        vectors = vectors / norms[:, np.newaxis]
        
        q = 0.0
        for j in range(k):
            for l in range(j+1, k):
                cos_angle = np.dot(vectors[j], vectors[l])
                q += (cos_angle + 1/3)**2
        q = 1 - 3/8 * q
        q_values[i] = q
    
    return q_values

def plot_phase_transition_markers(data, output_file='phase_transition_markers.png'):
    """
    Create a figure with multiple plots showing quantitative markers of phase transitions
    
    Parameters:
    -----------
    data : dict
        Dictionary containing all necessary data series
    output_file : str
        Path to save the output figure
    """
    print(f"Plotting phase transition markers to {output_file}...")
    
    # Create figure with GridSpec for more control over layout
    fig = plt.figure(figsize=(10, 15))
    gs = gridspec.GridSpec(4, 1, height_ratios=[1, 1, 1, 1])
    
    # Plot 1: Free energy derivatives with respect to temperature
    ax1 = plt.subplot(gs[0])
    ax1.plot(data['temperature'], data['dF_dT'], 'b-', linewidth=2)
    
    # Mark phase transition thresholds
    threshold_value = 0.05  # Threshold for phase transition
    ax1.axhline(y=threshold_value, color='r', linestyle='--', label='Phase transition threshold')
    ax1.axhline(y=-threshold_value, color='r', linestyle='--')
    
    # Identify and highlight transition regions
    transitions = np.where(np.abs(data['dF_dT']) > threshold_value)[0]
    if len(transitions) > 0:
        # Group consecutive indices
        groups = np.split(transitions, np.where(np.diff(transitions) != 1)[0] + 1)
        for group in groups:
            if len(group) > 5:  # Ignore small noise spikes
                start_idx = max(0, group[0] - 5)
                end_idx = min(len(data['temperature']) - 1, group[-1] + 5)
                ax1.axvspan(data['temperature'][start_idx], data['temperature'][end_idx], 
                           alpha=0.2, color='blue', label='_nolegend_')
                
    ax1.set_ylabel('∂F/∂T (kcal/mol·K)', fontsize=12)
    ax1.set_xlim(min(data['temperature']), max(data['temperature']))
    ax1.set_title('Free Energy Derivatives - Temperature Response', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot 2: Free energy derivatives with respect to pressure
    ax2 = plt.subplot(gs[1])
    ax2.plot(data['pressure'], data['dF_dP'], 'g-', linewidth=2)
    
    # Mark phase transition thresholds
    p_threshold_value = 0.08  # Threshold for pressure transitions
    ax2.axhline(y=p_threshold_value, color='r', linestyle='--', label='Phase transition threshold')
    ax2.axhline(y=-p_threshold_value, color='r', linestyle='--')
    
    # Identify and highlight transition regions
    p_transitions = np.where(np.abs(data['dF_dP']) > p_threshold_value)[0]
    if len(p_transitions) > 0:
        # Group consecutive indices
        p_groups = np.split(p_transitions, np.where(np.diff(p_transitions) != 1)[0] + 1)
        for group in p_groups:
            if len(group) > 5:  # Ignore small noise spikes
                start_idx = max(0, group[0] - 5)
                end_idx = min(len(data['pressure']) - 1, group[-1] + 5)
                ax2.axvspan(data['pressure'][start_idx], data['pressure'][end_idx], 
                           alpha=0.2, color='green', label='_nolegend_')
    
    ax2.set_ylabel('∂F/∂P (kcal/mol·MPa)', fontsize=12)
    ax2.set_xlabel('Pressure (MPa)', fontsize=12)
    ax2.set_xlim(min(data['pressure']), max(data['pressure']))
    ax2.set_title('Free Energy Derivatives - Pressure Response', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Plot 3: Heat capacity and density
    ax3 = plt.subplot(gs[2])
    # Heat capacity
    l1 = ax3.plot(data['temperature'], data['heat_capacity'], 'b-', linewidth=2, label='Heat capacity')
    ax3.set_ylabel('Heat capacity (J/mol·K)', fontsize=12)
    ax3.set_xlim(min(data['temperature']), max(data['temperature']))
    
    # Twin axis for density
    ax3b = ax3.twinx()
    l2 = ax3b.plot(data['temperature'], data['density'], 'r-', linewidth=2, label='Density')
    ax3b.set_ylabel('Density (g/cm³)', fontsize=12, color='r')
    
    # Highlight transition zones
    # Ice to liquid transition zone
    ice_liquid_zone = (data['temperature'] > 270) & (data['temperature'] < 276)
    ax3.axvspan(min(data['temperature'][ice_liquid_zone]), max(data['temperature'][ice_liquid_zone]), 
               alpha=0.2, color='cyan', label='Ice-Liquid transition')
    
    # Liquid to gas transition zone (if present)
    if max(data['temperature']) > 370:
        liquid_gas_zone = (data['temperature'] > 370) & (data['temperature'] < 376)
        ax3.axvspan(min(data['temperature'][liquid_gas_zone]), max(data['temperature'][liquid_gas_zone]), 
                   alpha=0.2, color='orange', label='Liquid-Gas transition')
    
    # Combine legends
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax3.legend(lines, labels, loc='upper right')
    
    ax3.set_title('Heat Capacity and Density Changes', fontsize=14)
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Tetrahedral order parameter
    ax4 = plt.subplot(gs[3])
    ax4.plot(data['temperature'], data['order_parameter'], 'g-', linewidth=2, label='Tetrahedral order')
    
    # Highlight different phase regions based on order parameter
    ax4.axhspan(0.8, 1.0, alpha=0.2, color='cyan', label='Ice phase (q₄ > 0.8)')
    ax4.axhspan(0.5, 0.7, alpha=0.2, color='blue', label='Liquid phase (0.5 < q₄ < 0.7)')
    ax4.axhspan(0.0, 0.4, alpha=0.2, color='red', label='Gas phase (q₄ < 0.4)')
    
    ax4.set_ylabel('Tetrahedral order parameter q₄', fontsize=12)
    ax4.set_xlabel('Temperature (K)', fontsize=12)
    ax4.set_ylim(0, 1)
    ax4.set_xlim(min(data['temperature']), max(data['temperature']))
    ax4.set_title('Tetrahedral Order Parameter', fontsize=14)
    ax4.grid(True, alpha=0.3)
    ax4.legend(loc='upper right')
    
    # Adding main title for the whole figure
    plt.suptitle('Quantitative Markers of Phase Transitions', fontsize=16, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    
    # Save figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {output_file}")
    plt.close()

def main():
    """Main function to run the script"""
    print("Phase Transition Markers Visualization")
    
    # Create visualization directory if it doesn't exist
    os.makedirs("../visualization", exist_ok=True)
    
    # Try to load real data from GROMACS files
    try:
        # Check if the files exist
        if os.path.exists('../md.edr') and os.path.exists('../md.tpr'):
            print("Loading data from GROMACS files...")
            
            # TODO: Load and process real data
            # For this demo, we'll use synthetic data
            use_real_data = False
        else:
            print("GROMACS files not found, generating synthetic data for demonstration.")
            use_real_data = False
    except Exception as e:
        print(f"Error loading real data: {e}")
        print("Falling back to synthetic data...")
        use_real_data = False
    
    # Generate synthetic data for now
    data = generate_synthetic_data()
    
    # Plot the markers
    plot_phase_transition_markers(data, output_file='../visualization/phase_transition_markers.png')
    
    print("Visualization complete.")

if __name__ == "__main__":
    main() 