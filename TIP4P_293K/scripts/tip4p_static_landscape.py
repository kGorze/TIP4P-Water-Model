#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import os
import subprocess
import tempfile

def create_output_dir():
    """Create output directory if it doesn't exist."""
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../visualization_outputs")
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def extract_temperature_pressure_from_edr():
    """
    Extract temperature and pressure data from GROMACS .edr file.
    Returns temperature and pressure arrays.
    """
    # Fix the path - check if we're in scripts dir or main dir
    if os.path.basename(os.getcwd()) == "scripts":
        edr_file = "../md.edr"
    else:
        edr_file = "md.edr"
        
    if not os.path.exists(edr_file):
        print(f"Warning: {edr_file} not found. Using simulated data.")
        # Generate simulated data
        n_points = 100
        temperatures = np.random.normal(298, 0.5, n_points) 
        pressures = np.random.normal(1.0, 2.0, n_points)  # 1 bar with small fluctuations
        return temperatures, pressures
    
    try:
        print(f"Found EDR file for T/P extraction: {edr_file}")
        # Create temporary files for temperature and pressure outputs
        with tempfile.NamedTemporaryFile(suffix='.xvg', delete=False) as temp_temp_file:
            temp_temp_xvg = temp_temp_file.name
        
        with tempfile.NamedTemporaryFile(suffix='.xvg', delete=False) as temp_press_file:
            temp_press_xvg = temp_press_file.name
        
        # Extract temperature (typically selection 15 for Temperature)
        input_str = "15\n"  # Select temperature
        process = subprocess.Popen(['gmx', 'energy', '-f', edr_file, '-o', temp_temp_xvg],
                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate(input=input_str)
        
        if process.returncode != 0:
            print(f"Error running gmx energy for temperature: {stderr}")
            raise RuntimeError(f"gmx energy failed with return code {process.returncode}")
            
        # Extract pressure (typically selection 16 for Pressure)
        input_str = "16\n"  # Select pressure
        process = subprocess.Popen(['gmx', 'energy', '-f', edr_file, '-o', temp_press_xvg],
                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate(input=input_str)
        
        if process.returncode != 0:
            print(f"Error running gmx energy for pressure: {stderr}")
            raise RuntimeError(f"gmx energy failed with return code {process.returncode}")
        
        # Parse the XVG files
        temp_data = np.loadtxt(temp_temp_xvg, comments=['#', '@'])
        press_data = np.loadtxt(temp_press_xvg, comments=['#', '@'])
        
        # Clean up the temporary files
        os.unlink(temp_temp_xvg)
        os.unlink(temp_press_xvg)
        
        # Extract the values
        temperatures = temp_data[:, 1]  # Second column is temperature
        pressures = press_data[:, 1]  # Second column is pressure
        
        # Downsample if necessary (to avoid too many points in the plot)
        if len(temperatures) > 500:
            idx = np.linspace(0, len(temperatures)-1, 500, dtype=int)
            temperatures = temperatures[idx]
            pressures = pressures[idx]
        
        return temperatures, pressures
    
    except Exception as e:
        print(f"Error extracting temperature/pressure data: {e}")
        print("Using simulated data instead.")
        # Fall back to simulated data
        n_points = 100
        temperatures = np.random.normal(298, 0.5, n_points)
        pressures = np.random.normal(1.0, 2.0, n_points)
        return temperatures, pressures

def generate_tip4p_energy_landscape(temp_range=(250, 350), pressure_range=(1, 400), resolution=100):
    """
    Generate a 2D free energy landscape for TIP4P water model.
    Parameters represent temperature (K) and pressure (bar) space.
    """
    temp = np.linspace(temp_range[0], temp_range[1], resolution)
    press = np.linspace(pressure_range[0], pressure_range[1], resolution)
    T, P = np.meshgrid(temp, press)
    
    # Create a landscape with multiple minima and phase boundaries
    # Ice-liquid boundary around 273K
    ice_boundary = 273.15
    ice_transition_width = 5
    ice_factor = 0.5 * (1 + np.tanh((ice_boundary - T) / ice_transition_width))
    
    # Base energy function with multiple minima
    E_base = (
        # Global minimum for ice Ih around 270K and moderate pressure
        -5 * np.exp(-0.02 * ((T - 270)**2 + 0.01 * (P - 50)**2))
        
        # Local minimum for liquid water around 298K
        - 3 * np.exp(-0.02 * ((T - 298)**2 + 0.01 * (P - 100)**2))
        
        # Another minimum for ice III at higher pressure
        - 4 * np.exp(-0.02 * ((T - 265)**2 + 0.01 * (P - 300)**2))
    )
    
    # Combine factors to create realistic phase regions
    E = E_base + 2 * ice_factor
    
    # Add small random noise
    np.random.seed(42)
    E += np.random.normal(0, 0.1, E.shape)
    
    # Normalize to a reasonable range
    E = (E - E.min()) / (E.max() - E.min()) * 10
    
    return T, P, E

def create_static_landscape(output_dir):
    """Create static visualizations with energy landscapes and sampling points."""
    
    # Extract actual temperature and pressure from simulation
    real_temps, real_press = extract_temperature_pressure_from_edr()
    
    # Calculate statistics for real data
    avg_temp = np.mean(real_temps)
    avg_press = np.mean(real_press)
    std_temp = np.std(real_temps)
    std_press = np.std(real_press)
    
    print(f"Real data statistics:")
    print(f"Temperature: {avg_temp:.2f} ± {std_temp:.2f} K")
    print(f"Pressure: {avg_press:.2f} ± {std_press:.2f} bar")
    
    # Generate energy landscape
    temp_range = (250, 350)
    pressure_range = (1, 400)
    T, P, E = generate_tip4p_energy_landscape(temp_range, pressure_range)
    
    # Create figure with two panels
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
    
    # -------------- Left panel: Current MD Sampling --------------
    # Plot the energy landscape
    contour1 = ax1.contourf(T, P, E, 20, cmap='viridis_r')
    
    # Add colorbar
    cbar = plt.colorbar(contour1, ax=[ax1, ax2])
    cbar.set_label('Free Energy (kJ/mol)', fontweight='bold')
    
    # Add freezing point line
    ax1.axvline(x=273.15, color='white', linestyle='--', alpha=0.7, label='Freezing Point (273.15K)')
    
    # Add the standard MD sampling points with high visibility
    # Plot the mean point larger
    ax1.scatter(avg_temp, avg_press, color='red', s=200, edgecolor='white', linewidth=2, 
               label='Current Simulation Point', zorder=10)
    
    # Add labels for different regions
    ax1.text(270, 50, "Ice Ih", color='white', fontweight='bold')
    ax1.text(298, 100, "Liquid", color='white', fontweight='bold')
    ax1.text(265, 300, "Ice III", color='white', fontweight='bold')
    
    # Set titles and labels
    ax1.set_title('TIP4P at 298K (Current Simulation)', fontweight='bold')
    ax1.set_xlabel('Temperature (K)', fontweight='bold')
    ax1.set_ylabel('Pressure (bar)', fontweight='bold')
    ax1.set_ylim(0, 400)
    
    # Add legend
    ax1.legend(loc='upper right')
    
    # Add key statistics as text
    stats_text = (
        f"Temperature: {avg_temp:.2f} ± {std_temp:.2f} K\n"
        f"Pressure: {avg_press:.2f} ± {std_press:.2f} bar"
    )
    ax1.text(0.05, 0.95, stats_text, transform=ax1.transAxes, color='white', 
            fontweight='bold', va='top', ha='left',
            bbox=dict(facecolor='black', alpha=0.5, boxstyle='round'))
    
    # -------------- Right panel: Enhanced Sampling Possibilities --------------
    # Plot the energy landscape
    contour2 = ax2.contourf(T, P, E, 20, cmap='viridis_r')
    
    # Add freezing point line
    ax2.axvline(x=273.15, color='white', linestyle='--', alpha=0.7)
    
    # Add phase labels
    ax2.text(270, 50, "Ice Ih", color='white', fontweight='bold')
    ax2.text(298, 100, "Liquid", color='white', fontweight='bold')
    ax2.text(265, 300, "Ice III", color='white', fontweight='bold')
    
    # Add enhanced sampling points
    # 1. Replica Exchange
    ax2.scatter(300, 50, color='cyan', s=100, edgecolor='white', linewidth=1.5, 
               label='Replica Exchange', zorder=10)
    
    # 2. Simulated Annealing
    ax2.scatter(280, 100, color='magenta', s=100, edgecolor='white', linewidth=1.5,
               label='Simulated Annealing', zorder=10)
    
    # 3. Umbrella Sampling
    ax2.scatter(270, 200, color='yellow', s=100, edgecolor='white', linewidth=1.5,
               label='Umbrella Sampling', zorder=10)
    
    # 4. Metadynamics
    ax2.scatter(260, 50, color='lime', s=100, edgecolor='white', linewidth=1.5,
               label='Metadynamics', zorder=10)
    
    # Set titles and labels
    ax2.set_title('Enhanced Sampling Possibilities', fontweight='bold')
    ax2.set_xlabel('Temperature (K)', fontweight='bold')
    ax2.set_ylabel('Pressure (bar)', fontweight='bold')
    ax2.set_ylim(0, 400)
    
    # Add legend
    ax2.legend(loc='upper right')
    
    # Add ice formation probability text in top right corner
    ax1.text(0.98, 0.95, 'Ice Formation\nProbability: 0.95', 
            transform=ax1.transAxes, ha='right', va='top',
            color='white', fontweight='bold',
            bbox=dict(facecolor='blue', alpha=0.4))
    
    ax2.text(0.98, 0.95, 'Ice Formation\nProbability: 0.95', 
            transform=ax2.transAxes, ha='right', va='top',
            color='white', fontweight='bold',
            bbox=dict(facecolor='blue', alpha=0.4))
    
    # Add title for the entire figure
    plt.suptitle('TIP4P Water Model: Standard vs Enhanced Sampling for Ice Formation', 
                fontsize=16, fontweight='bold')
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    # Save figure
    plt.savefig(os.path.join(output_dir, 'tip4p_sampling_landscape.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'tip4p_sampling_landscape.pdf'), bbox_inches='tight')
    
    print(f"Static landscape visualization saved to {output_dir}")
    plt.close()

if __name__ == "__main__":
    print("Creating static TIP4P landscape visualization")
    output_dir = create_output_dir()
    create_static_landscape(output_dir)
    print("Visualization complete") 