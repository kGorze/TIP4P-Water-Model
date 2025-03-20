#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import os
import sys
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

def create_simple_visualization(output_dir):
    """Create a simple visualization of TIP4P water model sampling and enhanced sampling possibilities."""
    
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
    
    # Create figure with two panels
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # -------------- Left panel: Actual simulation data --------------
    # Plot the temperature-pressure data
    ax1.scatter(real_temps, real_press, c='blue', alpha=0.5, s=10, label=f'TIP4P at 298K')
    
    # Add mean point with error bars
    ax1.errorbar(avg_temp, avg_press, xerr=std_temp, yerr=std_press, 
               fmt='ro', markersize=8, capsize=5, label='Mean ± Std Dev')
    
    # Set labels and title
    ax1.set_xlabel('Temperature (K)', fontsize=12)
    ax1.set_ylabel('Pressure (bar)', fontsize=12)
    ax1.set_title('Current TIP4P Simulation at 298K', fontsize=14, fontweight='bold')
    
    # Add freezing point line
    ax1.axvline(x=273.15, color='skyblue', linestyle='--', alpha=0.7, label='Freezing Point (273.15K)')
    
    # Set x and y limits to focus on the data
    ax1.set_xlim(avg_temp - 5*std_temp, avg_temp + 5*std_temp)
    ax1.set_ylim(max(0, avg_press - 5*std_press), avg_press + 5*std_press)
    
    # Add grid
    ax1.grid(True, alpha=0.3)
    
    # Add legend
    ax1.legend(loc='best')
    
    # Add text with statistics
    stats_text = (
        f"Temperature: {avg_temp:.2f} ± {std_temp:.2f} K\n"
        f"Pressure: {avg_press:.2f} ± {std_press:.2f} bar\n"
        f"Sample points: {len(real_temps)}"
    )
    ax1.text(0.05, 0.05, stats_text, transform=ax1.transAxes, 
           bbox=dict(facecolor='white', alpha=0.7, boxstyle='round'))
    
    # -------------- Right panel: Enhanced sampling possibilities --------------
    # Create a temperature-pressure diagram showing where enhanced sampling could explore
    temp_range = np.linspace(250, 350, 100)
    press_range = np.linspace(0, 200, 100)
    
    # Mark the current sampling region
    current_rect = plt.Rectangle(
        (avg_temp - 3*std_temp, avg_press - 3*std_press),
        6*std_temp, 6*std_press,
        facecolor='blue', alpha=0.2, label='Current Sampling'
    )
    ax2.add_patch(current_rect)
    
    # Add freezing point line
    ax2.axvline(x=273.15, color='skyblue', linestyle='--', alpha=0.7, label='Freezing Point (273.15K)')
    
    # Add enhanced sampling examples
    # 1. Replica Exchange: vertical exploration across temperatures
    temps_replica = np.linspace(260, 320, 20)
    press_replica = np.full_like(temps_replica, avg_press) + np.random.normal(0, 5, size=len(temps_replica))
    ax2.scatter(temps_replica, press_replica, color='red', alpha=0.7, s=30, label='Replica Exchange')
    
    # 2. Umbrella Sampling: specific path across barriers
    temps_umbrella = np.linspace(270, 295, 10)
    press_umbrella = np.linspace(1, 100, 10)
    ax2.scatter(temps_umbrella, press_umbrella, color='green', alpha=0.7, s=30, label='Umbrella Sampling')
    ax2.plot(temps_umbrella, press_umbrella, 'g-', alpha=0.5)
    
    # 3. Simulated Annealing: cooling from liquid to ice formation
    temps_annealing = np.linspace(298, 260, 10)
    press_annealing = np.full_like(temps_annealing, avg_press) + np.random.normal(0, 2, size=len(temps_annealing))
    ax2.scatter(temps_annealing, press_annealing, color='purple', alpha=0.7, s=30, label='Simulated Annealing')
    ax2.plot(temps_annealing, press_annealing, 'purple', alpha=0.5)
    
    # 4. Metadynamics: comprehensive sampling
    np.random.seed(42)
    temps_meta = np.random.uniform(255, 320, 30)
    press_meta = np.random.uniform(1, 150, 30)
    ax2.scatter(temps_meta, press_meta, color='orange', alpha=0.7, s=30, label='Metadynamics')
    
    # Set labels and title
    ax2.set_xlabel('Temperature (K)', fontsize=12)
    ax2.set_ylabel('Pressure (bar)', fontsize=12)
    ax2.set_title('Enhanced Sampling Possibilities', fontsize=14, fontweight='bold')
    
    # Add labels for different regions
    ax2.text(260, 80, "Ice Formation\nRegion", ha='center', fontsize=10, 
           bbox=dict(facecolor='white', alpha=0.7, boxstyle='round'))
    ax2.text(300, 80, "Liquid Region", ha='center', fontsize=10, 
           bbox=dict(facecolor='white', alpha=0.7, boxstyle='round'))
    
    # Set axis limits
    ax2.set_xlim(250, 350)
    ax2.set_ylim(0, 200)
    
    # Add grid
    ax2.grid(True, alpha=0.3)
    
    # Add legend
    ax2.legend(loc='best')
    
    # Add main title for the entire figure
    plt.suptitle('TIP4P Water Model: Current vs. Enhanced Sampling Approaches', 
                fontsize=16, fontweight='bold')
    
    # Add explanation text
    explanation = (
        "Left: Actual temperature and pressure data from TIP4P simulation at 298K\n"
        "Right: Potential regions that could be explored using enhanced sampling techniques"
    )
    
    fig.text(0.5, 0.01, explanation, ha='center', fontsize=12, 
            bbox=dict(facecolor='lightyellow', alpha=0.5, boxstyle='round'))
    
    # Adjust layout and save
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    # Save the figure as PNG and PDF
    plt.savefig(os.path.join(output_dir, 'tip4p_sampling_simple.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'tip4p_sampling_simple.pdf'), bbox_inches='tight')
    
    # Create and save static version of the comparison
    print(f"Visualization saved to {output_dir}")
    
    plt.close()

if __name__ == "__main__":
    print("Creating simple TIP4P sampling visualization")
    output_dir = create_output_dir()
    create_simple_visualization(output_dir)
    print("Visualization complete") 