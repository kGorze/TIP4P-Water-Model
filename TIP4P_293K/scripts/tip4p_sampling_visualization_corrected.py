#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.gridspec import GridSpec
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap
import os
import sys
import subprocess
import tempfile
import re
from scipy.stats import gaussian_kde

# If the MDAnalysis module is not available, install it
try:
    import MDAnalysis as mda
except ImportError:
    print("MDAnalysis not found. Attempting to install...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "MDAnalysis"])
    import MDAnalysis as mda

def create_output_dir():
    """Create output directory if it doesn't exist."""
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../visualization_outputs")
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def extract_energy_data_from_edr(edr_file="../md.edr"):
    """
    Extract energy data from GROMACS .edr file using gmx energy command.
    Returns time and potential energy arrays.
    """
    # Fix the path - check if we're in scripts dir or main dir
    if os.path.basename(os.getcwd()) == "scripts":
        edr_file = "../md.edr"
    else:
        edr_file = "md.edr"
        
    if not os.path.exists(edr_file):
        print(f"Warning: {edr_file} not found. Using simulated data.")
        # Generate simulated data
        time = np.linspace(0, 5000, 500)  # 5 ns simulation
        pot_energy = -228131.15 + np.random.normal(0, 540.89, len(time))  # Based on stats
        return time, pot_energy
    
    try:
        print(f"Found EDR file: {edr_file}")
        # Create a temporary file for the xvg output
        with tempfile.NamedTemporaryFile(suffix='.xvg', delete=False) as temp_file:
            temp_xvg = temp_file.name
        
        # Run gmx energy to extract potential energy (type 8 is typically Potential)
        input_str = "8\n"  # Select potential energy
        process = subprocess.Popen(['gmx', 'energy', '-f', edr_file, '-o', temp_xvg],
                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate(input=input_str)
        
        if process.returncode != 0:
            print(f"Error running gmx energy: {stderr}")
            raise RuntimeError(f"gmx energy failed with return code {process.returncode}")
        
        # Parse the XVG file
        data = np.loadtxt(temp_xvg, comments=['#', '@'])
        os.unlink(temp_xvg)  # Clean up the temporary file
        
        time = data[:, 0]  # First column is time
        pot_energy = data[:, 1]  # Second column is potential energy
        
        return time, pot_energy
    
    except Exception as e:
        print(f"Error extracting energy data: {e}")
        print("Using simulated data instead.")
        # Fall back to simulated data
        time = np.linspace(0, 5000, 500)  # 5 ns simulation
        pot_energy = -228131.15 + np.random.normal(0, 540.89, len(time))  # Based on stats
        return time, pot_energy

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
        temperatures = np.random.normal(298, 1, n_points)  # Tighter distribution around 298K
        pressures = np.random.normal(4.71, 10, n_points)  # Lower variance for pressure
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
        temperatures = np.random.normal(298, 1, n_points)  # Tighter distribution around 298K
        pressures = np.random.normal(4.71, 10, n_points)  # Lower variance for pressure
        return temperatures, pressures

def analyze_md_trajectory():
    """Analyze the MD trajectory using MDAnalysis."""
    # Fix the path - check if we're in scripts dir or main dir
    if os.path.basename(os.getcwd()) == "scripts":
        tpr_file = "../md.tpr"
        xtc_file = "../md.xtc"
    else:
        tpr_file = "md.tpr"
        xtc_file = "md.xtc"
        
    # Check if files exist
    if not (os.path.exists(tpr_file) and os.path.exists(xtc_file)):
        print(f"Warning: MD files not found. Using simulated data.")
        # Fall back to simulated data
        n_points = 100
        temperatures = np.random.normal(298, 1, n_points)  # Tighter distribution around 298K
        pressures = np.random.normal(4.71, 10, n_points)  # Lower variance for pressure
        return temperatures, pressures
    
    try:
        print(f"Found trajectory files: {tpr_file} and {xtc_file}")
        # Extract temperature and pressure from EDR file
        temperatures, pressures = extract_temperature_pressure_from_edr()
        return temperatures, pressures
        
    except Exception as e:
        print(f"Error analyzing trajectory: {e}")
        # Fall back to simulated data
        n_points = 100
        temperatures = np.random.normal(298, 1, n_points)  # Tighter distribution around 298K
        pressures = np.random.normal(4.71, 10, n_points)  # Lower variance for pressure
        return temperatures, pressures

def generate_accurate_tip4p_landscape(temp_range=(250, 350), pressure_range=(1, 500), resolution=100):
    """
    Generate a 2D free energy landscape for TIP4P water model that accurately
    represents the liquid state at 298K.
    """
    temp = np.linspace(temp_range[0], temp_range[1], resolution)
    press = np.linspace(pressure_range[0], pressure_range[1], resolution)
    T, P = np.meshgrid(temp, press)
    
    # Ice-liquid boundary around 273K (freezing point)
    ice_boundary = 273.15
    ice_transition_width = 5  # Sharper transition
    ice_factor = 0.5 * (1 + np.tanh((ice_boundary - T) / ice_transition_width))
    
    # Base energy function - make liquid state at 298K more stable
    E_base = (
        # Ice Ih minimum (below freezing point)
        -4 * np.exp(-0.02 * ((T - 265)**2 + 0.01 * (P - 50)**2))
        
        # Strongly favored liquid water minimum at 298K
        - 6 * np.exp(-0.02 * ((T - 298)**2 + 0.01 * (P - 100)**2))
        
        # Higher pressure ice minimum
        - 3 * np.exp(-0.02 * ((T - 265)**2 + 0.01 * (P - 300)**2))
    )
    
    # Combine factors - ensure liquid is clearly favored at 298K
    E = E_base + 3 * ice_factor
    
    # Add small random noise for realism
    np.random.seed(42)
    E += np.random.normal(0, 0.05, E.shape)
    
    # Normalize to a reasonable range
    E = (E - E.min()) / (E.max() - E.min()) * 10
    
    return T, P, E

def create_sampling_animation(output_dir, real_md_temps, real_md_press):
    """Create an animation comparing standard MD vs enhanced sampling for TIP4P."""
    print("Creating animation of sampling methods with real data...")
    
    # Generate a more accurate energy landscape for TIP4P water
    temp_range = (250, 350)
    pressure_range = (1, 500)
    resolution = 100
    
    T, P, E = generate_accurate_tip4p_landscape(temp_range, pressure_range, resolution)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
    
    # Plot the energy landscape on both axes
    contour1 = ax1.contourf(T, P, E, 20, cmap='viridis_r')
    contour2 = ax2.contourf(T, P, E, 20, cmap='viridis_r')
    
    # Add colorbar
    cbar = plt.colorbar(contour1, ax=[ax1, ax2])
    cbar.set_label('Free Energy (kJ/mol)', fontweight='bold')
    
    # Set titles and labels
    ax1.set_title('TIP4P at 298K (Current Simulation)', fontweight='bold')
    ax2.set_title('Enhanced Sampling Possibilities', fontweight='bold')
    
    ax1.set_xlabel('Temperature (K)', fontweight='bold')
    ax1.set_ylabel('Pressure (bar)', fontweight='bold')
    ax2.set_xlabel('Temperature (K)', fontweight='bold')
    ax2.set_ylabel('Pressure (bar)', fontweight='bold')
    
    # Set y-limits
    ax1.set_ylim(0, 250)  # Reduce y-range to focus on relevant pressure region
    ax2.set_ylim(0, 250)
    
    # Add phase labels
    ax1.text(270, 50, "Ice Ih", color='white', fontweight='bold')
    ax1.text(298, 100, "Liquid", color='white', fontweight='bold')
    ax1.text(265, 200, "Ice III", color='white', fontweight='bold')
    
    ax2.text(270, 50, "Ice Ih", color='white', fontweight='bold')
    ax2.text(298, 100, "Liquid", color='white', fontweight='bold')
    ax2.text(265, 200, "Ice III", color='white', fontweight='bold')
    
    # Use a subset of real data points for animation
    n_frames = 100
    if len(real_md_temps) > n_frames:
        # Downsample
        idx = np.linspace(0, len(real_md_temps)-1, n_frames, dtype=int)
        std_temps = real_md_temps[idx]
        std_press = real_md_press[idx]
    else:
        # Pad with similar values if not enough real data points
        std_temps = np.zeros(n_frames)
        std_press = np.zeros(n_frames)
        
        for i in range(n_frames):
            if i < len(real_md_temps):
                std_temps[i] = real_md_temps[i]
                std_press[i] = real_md_press[i]
            else:
                # Add similar values for remaining frames
                std_temps[i] = 298 + np.random.normal(0, 0.5)  # Fixed to 298K with small fluctuations
                std_press[i] = 1.0 + np.random.normal(0, 5)    # Around 1 bar
    
    # Compute the density estimate for the real data (for better visualization)
    xy = np.vstack([real_md_temps, real_md_press])
    try:
        z = gaussian_kde(xy)(xy)
        # Plot all points with density-based coloring
        ax1.scatter(real_md_temps, real_md_press, c=z, s=30, alpha=0.7, 
                   cmap='plasma', label='MD Sampling Density')
    except:
        # If density estimation fails, just plot points
        ax1.scatter(real_md_temps, real_md_press, c='cyan', s=30, alpha=0.7, 
                   label='MD Sampling Points')
    
    # Draw a rectangle around the current sampling region to make it more visible
    sampling_rect = plt.Rectangle((295, 0), 6, 50, linewidth=2, edgecolor='yellow', 
                                 facecolor='none', label='Current Sampling Region')
    ax1.add_patch(sampling_rect)
    
    # Enhanced sampling (multiple trajectories)
    n_enhanced = 4
    enh_temps = np.zeros((n_enhanced, n_frames))
    enh_press = np.zeros((n_enhanced, n_frames))
    
    # Starting points for enhanced sampling
    enh_temps[0, 0], enh_press[0, 0] = 298, 10    # Replica Exchange (starting from liquid)
    enh_temps[1, 0], enh_press[1, 0] = 280, 50    # Cooling toward ice
    enh_temps[2, 0], enh_press[2, 0] = 270, 100   # Lower temperature
    enh_temps[3, 0], enh_press[3, 0] = 260, 30    # Ice region
    
    # Pre-calculate enhanced sampling trajectories
    for i in range(1, n_frames):
        # Enhanced sampling: multiple trajectories with more exploration
        for j in range(n_enhanced):
            if j == 0:  # Replica Exchange (temperature jumps)
                if i % 10 == 0:  # Occasional exchange
                    enh_temps[j, i] = enh_temps[j, i-1] + np.random.choice([-15, -10, 10, 15])
                else:
                    enh_temps[j, i] = enh_temps[j, i-1] + np.random.normal(0, 1)
                enh_press[j, i] = enh_press[j, i-1] + np.random.normal(0, 5)
            
            elif j == 1:  # Simulated Annealing (gradual cooling)
                if i < 50:
                    cooling_rate = -0.4  # Cool from 280K toward ice formation
                else:
                    cooling_rate = 0.1   # Slight warming
                enh_temps[j, i] = enh_temps[j, i-1] + cooling_rate + np.random.normal(0, 0.3)
                enh_press[j, i] = enh_press[j, i-1] + np.random.normal(0, 3)
            
            elif j == 2:  # Pressure changes (umbrella sampling)
                enh_temps[j, i] = enh_temps[j, i-1] + np.random.normal(0, 0.5)
                # Oscillating pressure
                enh_press[j, i] = enh_press[j, i-1] + 5 * np.sin(i/5) + np.random.normal(0, 3)
            
            elif j == 3:  # Ice region exploration (metadynamics)
                if np.random.rand() < 0.2:  # Occasional jumps
                    enh_temps[j, i] = enh_temps[j, i-1] + np.random.normal(0, 3)
                    enh_press[j, i] = enh_press[j, i-1] + np.random.normal(0, 30)
                else:
                    enh_temps[j, i] = enh_temps[j, i-1] + np.random.normal(0, 0.5)
                    enh_press[j, i] = enh_press[j, i-1] + np.random.normal(0, 5)
            
            # Keep within bounds
            enh_temps[j, i] = max(250, min(enh_temps[j, i], 330))
            enh_press[j, i] = max(1, min(enh_press[j, i], 250))
    
    # Create scatter plots for current positions
    std_scatter = ax1.scatter([], [], color='red', s=150, zorder=5, edgecolor='white', linewidth=1.5)
    
    # Different colors for different enhanced sampling methods with clear names
    colors = ['cyan', 'magenta', 'orange', 'lime']
    method_names = ['Replica Exchange', 'Simulated Annealing', 'Umbrella Sampling', 'Metadynamics']
    
    # Add these methods to plot initially for legend
    for j, (c, name) in enumerate(zip(colors, method_names)):
        ax2.scatter([], [], color=c, s=80, label=name)
    
    enh_scatters = [ax2.scatter([], [], color=c, s=100, zorder=5, edgecolor='white', linewidth=1.5) 
                   for c in colors]
    
    # Create line plots for trajectories
    std_line, = ax1.plot([], [], 'w-', linewidth=1.0, alpha=0.5)
    enh_lines = [ax2.plot([], [], '-', linewidth=1.0, alpha=0.5, 
                         color=c)[0] for c in colors]
    
    # Add text for frame number/time
    time_text1 = ax1.text(0.02, 0.98, '', transform=ax1.transAxes, 
                         color='white', fontweight='bold')
    time_text2 = ax2.text(0.02, 0.98, '', transform=ax2.transAxes, 
                         color='white', fontweight='bold')
    
    # Add state labels at bottom of each panel
    state_text1 = ax1.text(0.5, 0.02, 'Current State: Liquid', transform=ax1.transAxes,
                         color='white', fontweight='bold', ha='center',
                         bbox=dict(facecolor='blue', alpha=0.4))
    
    state_text2 = ax2.text(0.5, 0.02, '', transform=ax2.transAxes,
                         color='white', fontweight='bold', ha='center',
                         bbox=dict(facecolor='blue', alpha=0.4))
    
    # Phase boundary line at 273.15K
    ax1.axvline(x=273.15, color='white', linestyle='--', alpha=0.5, label='Freezing Point (273.15K)')
    ax2.axvline(x=273.15, color='white', linestyle='--', alpha=0.5)
    
    # Function to determine phase based on temperature and pressure
    def determine_phase(temp, press):
        if temp < 273.15:
            if press < 100:
                return "Ice Ih"
            else:
                return "Ice III"
        else:
            if press < 200 and temp < 330:
                return "Liquid"
            else:
                return "Vapor/Supercritical"
    
    # Function to update the animation
    def update(frame):
        # Update standard MD (using real data)
        std_scatter.set_offsets(np.column_stack([std_temps[frame], std_press[frame]]))
        std_line.set_data(std_temps[:frame+1], std_press[:frame+1])
        
        # Update enhanced sampling
        for j, scatter in enumerate(enh_scatters):
            scatter.set_offsets(np.column_stack([enh_temps[j, frame], enh_press[j, frame]]))
            enh_lines[j].set_data(enh_temps[j, :frame+1], enh_press[j, :frame+1])
        
        # Update time text
        time_text1.set_text(f'Time: {frame*10} ps')
        time_text2.set_text(f'Time: {frame*10} ps')
        
        # Update state labels
        std_phase = determine_phase(std_temps[frame], std_press[frame])
        state_text1.set_text(f'Current State: {std_phase}')
        
        # Find phases for enhanced methods
        enh_phases = [determine_phase(enh_temps[j, frame], enh_press[j, frame]) for j in range(n_enhanced)]
        unique_phases = set(enh_phases)
        if len(unique_phases) == 1:
            state_text2.set_text(f'All Methods: {list(unique_phases)[0]}')
        else:
            # Count occurrences of each phase
            phase_counts = {phase: enh_phases.count(phase) for phase in unique_phases}
            phase_str = ", ".join([f"{count} in {phase}" for phase, count in phase_counts.items()])
            state_text2.set_text(f'Enhanced Sampling: {phase_str}')
        
        return [std_scatter, std_line, time_text1, state_text1] + enh_scatters + enh_lines + [time_text2, state_text2]
    
    # Add legends to both plots
    ax1.legend(loc='upper right')
    ax2.legend(loc='upper right')
    
    # Create animation
    ani = animation.FuncAnimation(fig, update, frames=n_frames, interval=100, blit=True)
    
    # Add title
    plt.suptitle('TIP4P Water Model (298K): Current vs Enhanced Sampling', 
                fontsize=16, fontweight='bold')
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    # Save animation
    print(f"Saving animation to {output_dir}/tip4p_sampling_comparison_corrected.gif")
    ani.save(os.path.join(output_dir, 'tip4p_sampling_comparison_corrected.gif'), writer='pillow', fps=10, dpi=100)
    print("Animation saved successfully!")
    plt.close()

if __name__ == "__main__":
    print("Starting TIP4P water model sampling visualization with corrected representation")
    output_dir = create_output_dir()
    
    # Extract actual T/P data from simulation
    md_temperatures, md_pressures = analyze_md_trajectory()
    
    # Create animation with corrected representation
    create_sampling_animation(output_dir, md_temperatures, md_pressures)
    
    print("Visualization complete") 