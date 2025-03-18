#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec
import os
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
from scipy import stats  # Add scipy.stats for maxwell distribution

def create_output_dir():
    """Create output directory if it doesn't exist."""
    output_dir = "visualization_outputs"
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def plot_nvt_temperature_equilibration():
    """Plot a visualization of temperature equilibration during NVT."""
    output_dir = create_output_dir()
    
    # Simulation parameters
    target_temp = 298.0  # K (room temperature)
    simulation_time = 500  # ps
    
    # Create time points
    time_points = np.linspace(0, simulation_time, 500)
    
    # Create temperature data with initial randomness that converges to target
    # Initial random temperatures between 0-200K (cold start)
    initial_temp = 100.0
    temp_noise_amplitude = 50.0  # Initial noise amplitude
    
    # Generate temperature trajectory
    temperatures = []
    for t in time_points:
        # Exponential decay of noise amplitude over time
        decay_factor = np.exp(-t / 100)
        noise = (np.random.random() - 0.5) * temp_noise_amplitude * decay_factor
        
        # Temperature converges to target temperature
        temp = target_temp + (initial_temp - target_temp) * np.exp(-t / 50) + noise
        temperatures.append(temp)
    
    temperatures = np.array(temperatures)
    
    # Create figure
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[2, 1, 1])
    
    # Temperature plot
    ax_temp = fig.add_subplot(gs[0, :])
    ax_temp.plot(time_points, temperatures, 'b-', linewidth=1.5, alpha=0.7)
    ax_temp.axhline(y=target_temp, color='r', linestyle='--', linewidth=1.5, 
                    label=f'Target Temperature ({target_temp} K)')
    
    # Add shaded region for equilibrated portion
    equil_start = 200  # ps
    equil_idx = np.where(time_points >= equil_start)[0][0]
    ax_temp.axvspan(equil_start, simulation_time, alpha=0.2, color='green', 
                   label='Equilibrated Region')
    
    # Calculate statistics for equilibrated region
    equil_temps = temperatures[equil_idx:]
    equil_mean = np.mean(equil_temps)
    equil_std = np.std(equil_temps)
    
    # Add statistics to plot
    stats_text = f"Equilibrated Temperature: {equil_mean:.2f} ± {equil_std:.2f} K"
    ax_temp.text(0.02, 0.05, stats_text, transform=ax_temp.transAxes, 
                fontsize=10, bbox=dict(facecolor='white', alpha=0.7))
    
    # Set labels and title
    ax_temp.set_xlabel('Simulation Time (ps)')
    ax_temp.set_ylabel('Temperature (K)')
    ax_temp.set_title('NVT Equilibration: Temperature vs. Time')
    ax_temp.legend(loc='upper right')
    ax_temp.grid(True, alpha=0.3)
    
    # --- System energy fluctuations plot ---
    ax_energy = fig.add_subplot(gs[1, :])
    
    # Simulate energy fluctuations that become more regular over time
    energy_noise_scale = 500  # Scale of energy fluctuations
    energies = []
    
    for t in time_points:
        # Energy fluctuations become more regular as temperature stabilizes
        decay_factor = np.exp(-t / 100)
        irregular_noise = np.random.normal(0, energy_noise_scale * decay_factor)
        
        # Regular thermal fluctuations around constant energy
        regular_amp = 200  # Amplitude of regular fluctuations
        regular_freq = 0.2  # Frequency of regular fluctuations
        regular_fluct = regular_amp * np.sin(regular_freq * t)
        
        # Blend from irregular to regular fluctuations
        blend_factor = 1 - np.exp(-t / 150)
        
        # Total energy is sum of regular and irregular components
        energy = -200000 + regular_fluct * blend_factor + irregular_noise * (1 - blend_factor)
        energies.append(energy)
    
    energies = np.array(energies)
    
    # Plot energy
    ax_energy.plot(time_points, energies, 'g-', linewidth=1.5, alpha=0.7)
    ax_energy.axvspan(equil_start, simulation_time, alpha=0.2, color='green')
    
    # Set labels
    ax_energy.set_xlabel('Simulation Time (ps)')
    ax_energy.set_ylabel('Total Energy (kJ/mol)')
    ax_energy.set_title('System Energy During NVT Equilibration')
    ax_energy.grid(True, alpha=0.3)
    
    # --- Visualization of particles with velocities ---
    ax_visual = fig.add_subplot(gs[2, :])
    ax_visual.set_xlim(0, 10)
    ax_visual.set_ylim(0, 5)
    ax_visual.set_aspect('equal')
    ax_visual.set_xlabel('X (nm)')
    ax_visual.set_ylabel('Y (nm)')
    ax_visual.set_title('Velocity Distribution During NVT Equilibration')
    
    # Generate particles at three time points (start, middle, equilibrated)
    time_indices = [0, len(time_points) // 4, len(time_points) // 2, -1]
    time_labels = ['Initial', 'Early Equilibration', 'Mid Equilibration', 'Fully Equilibrated']
    colors = ['blue', 'purple', 'orange', 'red']
    
    # Generate random positions for particles
    np.random.seed(42)  # For reproducibility
    num_particles = 40
    positions = np.random.rand(num_particles, 2) * np.array([10, 5])
    
    # Draw particles for each time point
    legend_patches = []
    for i, (idx, label, color) in enumerate(zip(time_indices, time_labels, colors)):
        # Calculate temperature scale for velocity vectors
        temp_scale = temperatures[idx] / target_temp
        
        # Generate velocities based on temperature (direction is random, magnitude scales with sqrt(T))
        angles = np.random.rand(num_particles) * 2 * np.pi
        vel_magnitudes = stats.maxwell.rvs(scale=temp_scale, size=num_particles) / 2
        
        # Convert to velocity vectors
        velocities = np.column_stack([
            vel_magnitudes * np.cos(angles),
            vel_magnitudes * np.sin(angles)
        ])
        
        # Plot particles and velocity vectors in specific region
        x_start = 0 + i * 2.5
        x_end = 2.5 + i * 2.5
        
        # Only include particles in this region
        mask = (positions[:, 0] >= x_start) & (positions[:, 0] < x_end)
        pos_subset = positions[mask]
        vel_subset = velocities[mask]
        
        # Plot particles
        ax_visual.scatter(
            pos_subset[:, 0], pos_subset[:, 1], 
            color=color, s=50, alpha=0.7, 
        )
        
        # Plot velocity vectors
        ax_visual.quiver(
            pos_subset[:, 0], pos_subset[:, 1],
            vel_subset[:, 0], vel_subset[:, 1],
            color=color, scale=3, width=0.005, 
            scale_units='xy'
        )
        
        # Add vertical separator
        if i < len(time_indices) - 1:
            ax_visual.axvline(x=x_end, color='gray', linestyle='--', alpha=0.5)
        
        # Add to legend
        legend_patches.append(mpatches.Patch(color=color, label=f"{label} ({temperatures[idx]:.1f} K)"))
    
    ax_visual.legend(handles=legend_patches, loc='upper center', ncol=4, fontsize=8)
    
    # --- Add annotations explaining the NVT process ---
    explanation = """
    NVT Equilibration in Molecular Dynamics:
    
    1. Purpose: To bring the system to the target temperature before production runs
    
    2. Key Characteristics:
       • N (Number of particles): Fixed
       • V (Volume): Constant
       • T (Temperature): Controlled by a thermostat
    
    3. Process Steps:
       • Initial velocities are assigned from a Maxwell-Boltzmann distribution
       • Thermostat (e.g., V-rescale or Nosé-Hoover) adjusts velocities to reach target temperature
       • Temperature equilibrates with natural fluctuations around the target value
       • Energy fluctuations become more regular as the system reaches thermal equilibrium
       
    4. Thermostats:
       • V-rescale: Velocity rescaling with a stochastic term (efficient temperature control)
       • Nosé-Hoover: More accurately samples the canonical ensemble, but slower to equilibrate
       • Berendsen: Simple but distorts the ensemble; useful for initial equilibration
       
    5. Next Steps after NVT:
       • NPT equilibration to stabilize pressure and density
       • Production MD with proper ensemble sampling
    """
    
    # Add textbox with explanation
    fig.text(0.5, 0.01, explanation, ha='center', va='bottom', fontsize=10,
            bbox=dict(facecolor='whitesmoke', alpha=0.9, boxstyle='round,pad=1'))
    
    plt.tight_layout(rect=[0, 0.15, 1, 1])
    
    # Save figure
    plt.savefig(os.path.join(output_dir, "nvt_equilibration.png"), dpi=300)
    plt.savefig(os.path.join(output_dir, "nvt_equilibration.pdf"))
    
    return fig

def animate_temperature_equilibration():
    """Create an animation showing how temperature equilibrates during NVT simulation."""
    output_dir = create_output_dir()
    
    # Simulation parameters
    target_temp = 298.0  # K (room temperature)
    simulation_time = 500  # ps
    
    # Create time points
    time_points = np.linspace(0, simulation_time, 500)
    
    # Create temperature data with initial randomness that converges to target
    # Initial random temperatures between 0-200K (cold start)
    initial_temp = 100.0
    temp_noise_amplitude = 50.0  # Initial noise amplitude
    
    # Generate temperature trajectory
    temperatures = []
    for t in time_points:
        # Exponential decay of noise amplitude over time
        decay_factor = np.exp(-t / 100)
        noise = (np.random.random() - 0.5) * temp_noise_amplitude * decay_factor
        
        # Temperature converges to target temperature
        temp = target_temp + (initial_temp - target_temp) * np.exp(-t / 50) + noise
        temperatures.append(temp)
    
    temperatures = np.array(temperatures)
    
    # Create figure with two subplots - temperature vs time and energy histogram
    fig, (ax_temp, ax_energy) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [2, 1]})
    
    # Temperature line
    line, = ax_temp.plot([], [], 'b-', lw=2)
    target_line = ax_temp.axhline(y=target_temp, color='r', linestyle='--', linewidth=1.5, 
                      label=f'Target Temperature ({target_temp} K)')
    
    # Add shaded region for equilibrated portion
    equil_start = 200  # ps
    ax_temp.axvspan(equil_start, simulation_time, alpha=0.2, color='green', 
                   label='Equilibrated Region')
    
    # Set up temperature plot
    ax_temp.set_xlim(0, simulation_time)
    ax_temp.set_ylim(min(temperatures)*0.9, max(temperatures)*1.1)
    ax_temp.set_xlabel('Simulation Time (ps)')
    ax_temp.set_ylabel('Temperature (K)')
    ax_temp.set_title('NVT Equilibration: Temperature vs. Time')
    ax_temp.grid(True, alpha=0.3)
    ax_temp.legend(loc='upper right')
    
    # Text for stats
    stats_text = ax_temp.text(0.02, 0.05, '', transform=ax_temp.transAxes, 
                fontsize=10, bbox=dict(facecolor='white', alpha=0.7))
    
    # Set up energy histogram
    # Generate synthetic energy data
    energy_mean = -200000
    energy_data = np.zeros(100)  # Initial empty histogram
    ax_energy.set_xlim(energy_mean-1000, energy_mean+1000)
    ax_energy.set_ylim(0, 100)
    ax_energy.set_xlabel('Total Energy (kJ/mol)')
    ax_energy.set_ylabel('Frequency')
    ax_energy.set_title('Energy Distribution')
    
    # Energy histogram plot
    energy_bars = ax_energy.bar(
        np.linspace(energy_mean-1000, energy_mean+1000, 100),
        energy_data,
        width=20, alpha=0.7, color='g'
    )
    
    time_text = ax_temp.text(0.98, 0.95, '', transform=ax_temp.transAxes, 
                            ha='right', fontsize=12, 
                            bbox=dict(facecolor='white', alpha=0.7))
    
    plt.tight_layout()
    
    def init():
        line.set_data([], [])
        time_text.set_text('')
        stats_text.set_text('')
        for rect in energy_bars:
            rect.set_height(0)
        return [line, time_text, stats_text] + list(energy_bars)
    
    def update(frame):
        # Update temperature line
        current_time_points = time_points[:frame+1]
        current_temps = temperatures[:frame+1]
        line.set_data(current_time_points, current_temps)
        
        # Update time text
        time_text.set_text(f'Time: {time_points[frame]:.1f} ps')
        
        # Update stats if we're in equilibrated region
        if frame * simulation_time / len(time_points) >= equil_start:
            equil_idx = np.where(time_points >= equil_start)[0][0]
            equil_temps = temperatures[equil_idx:frame+1]
            equil_mean = np.mean(equil_temps)
            equil_std = np.std(equil_temps)
            stats_text.set_text(f"Equilibrated Temp: {equil_mean:.2f} ± {equil_std:.2f} K")
        
        # Update energy histogram - simulate energy fluctuations
        # As system equilibrates, the distribution becomes narrower and more gaussian
        if frame > 0:
            # Width of distribution narrows with equilibration
            energy_std = 500 * (1.0 - 0.5 * min(frame / (len(time_points) * 0.7), 1.0))
            # Center approaches equilibrium value
            temp_factor = temperatures[frame] / target_temp
            energy_shift = energy_mean + 500 * (temp_factor - 1.0)
            
            # Generate energy samples
            energy_samples = np.random.normal(energy_shift, energy_std, 1000)
            
            # Update histogram
            counts, bin_edges = np.histogram(
                energy_samples, 
                bins=100, 
                range=(energy_mean-1000, energy_mean+1000)
            )
            
            # Normalize counts 
            counts = counts * 100.0 / max(counts.max(), 1)
            
            # Update bar heights
            for i, rect in enumerate(energy_bars):
                rect.set_height(counts[i])
        
        return [line, time_text, stats_text] + list(energy_bars)
    
    # Create animation
    frames = range(0, len(time_points), 3)  # Use every 3rd frame for smoother animation
    ani = animation.FuncAnimation(fig, update, frames=frames,
                                init_func=init, blit=True, interval=30)
    
    # Save animation
    writer = animation.PillowWriter(fps=20)
    ani.save(os.path.join(output_dir, "nvt_temperature_animation.gif"), writer=writer, dpi=120)
    
    plt.close(fig)
    
    return "Temperature equilibration animation created"

def animate_velocity_distribution():
    """Create an animation showing how velocity distribution evolves during NVT equilibration."""
    output_dir = create_output_dir()
    
    # Simulation parameters
    target_temp = 298.0  # K (room temperature)
    num_particles = 1000  # Number of particles for distribution
    num_frames = 60      # Number of frames for animation
    
    # Create figure with two subplots
    fig, (ax_vel, ax_mb) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Set up scatter plot for velocities
    ax_vel.set_xlim(-2, 2)
    ax_vel.set_ylim(-2, 2)
    ax_vel.set_xlabel('Velocity X (nm/ps)')
    ax_vel.set_ylabel('Velocity Y (nm/ps)')
    ax_vel.set_title('2D Velocity Distribution')
    
    # Set up histogram for Maxwell-Boltzmann distribution
    ax_mb.set_xlim(0, 3)
    ax_mb.set_ylim(0, 1.5)
    ax_mb.set_xlabel('Velocity Magnitude (nm/ps)')
    ax_mb.set_ylabel('Probability Density')
    ax_mb.set_title('Maxwell-Boltzmann Distribution')
    
    # Generate reference Maxwell-Boltzmann distributions
    vel_range = np.linspace(0, 3, 1000)
    ideal_mb = stats.maxwell.pdf(vel_range, scale=np.sqrt(target_temp/298.0) * 0.5)
    ideal_line, = ax_mb.plot(vel_range, ideal_mb, 'r-', lw=2, alpha=0.7, 
                           label=f'Theoretical ({target_temp} K)')
    
    # Initial distribution with low temperature
    initial_temp = 50.0  # K (cold start)
    
    # Setup histogram
    n_bins = 30
    hist_vals = np.zeros(n_bins)
    bin_centers = np.linspace(0, 3, n_bins)
    bin_width = bin_centers[1] - bin_centers[0]
    bars = ax_mb.bar(bin_centers, hist_vals, width=bin_width, alpha=0.7, color='b', label='Simulation')
    
    # Scatter plot for 2D velocities
    scatter = ax_vel.scatter([], [], c=[], cmap='plasma', s=10, alpha=0.7)
    
    # Text labels
    temp_text = ax_vel.text(0.05, 0.95, '', transform=ax_vel.transAxes, 
                          fontsize=10, bbox=dict(facecolor='white', alpha=0.7))
    ax_mb.legend()
    
    # Temperature transition
    temperatures = initial_temp + (target_temp - initial_temp) * (1 - np.exp(-np.linspace(0, 5, num_frames)))
    
    def init():
        scatter.set_offsets(np.empty((0, 2)))
        scatter.set_array(np.array([]))
        temp_text.set_text('')
        for rect in bars:
            rect.set_height(0)
        return [scatter, temp_text] + list(bars)
    
    def update(frame):
        # Current temperature based on frame
        current_temp = temperatures[frame]
        
        # Generate velocities based on current temperature
        # Scale by sqrt(T) for Maxwell-Boltzmann distribution
        scale_factor = np.sqrt(current_temp/298.0) * 0.5
        
        # Generate random velocities for X and Y
        vx = np.random.normal(0, scale_factor, num_particles)
        vy = np.random.normal(0, scale_factor, num_particles)
        v_data = np.column_stack([vx, vy])
        
        # Update scatter plot
        scatter.set_offsets(v_data)
        
        # Calculate velocity magnitudes
        v_mags = np.sqrt(vx**2 + vy**2)
        
        # Update colormap based on velocity magnitude
        scatter.set_array(v_mags)
        
        # Update velocity histogram
        counts, _ = np.histogram(v_mags, bins=n_bins, range=(0, 3), density=True)
        for i, rect in enumerate(bars):
            rect.set_height(counts[i])
        
        # Update temperature text
        temp_text.set_text(f'Temperature: {current_temp:.1f} K')
        
        return [scatter, temp_text] + list(bars)
    
    # Create animation
    ani = animation.FuncAnimation(fig, update, frames=range(num_frames),
                                init_func=init, blit=True, interval=200)
    
    # Save animation
    writer = animation.PillowWriter(fps=10)
    ani.save(os.path.join(output_dir, "nvt_velocity_distribution.gif"), writer=writer, dpi=120)
    
    plt.close(fig)
    
    return "Velocity distribution animation created"

def animate_thermostat_action():
    """Create an animation showing how a thermostat works in NVT."""
    output_dir = create_output_dir()
    
    try:
        print("Starting thermostat animation with simplified approach...")
        
        # Simulation parameters
        num_particles = 50
        target_temp = 298.0  # K
        num_frames = 20  # Even fewer frames to ensure success
        box_size = 10.0  # nm
        
        # Initialize particle positions randomly in the box
        np.random.seed(42)  # For reproducibility
        positions = np.random.rand(num_particles, 2) * box_size
        
        # Initialize particle velocities with random directions
        angles = np.random.rand(num_particles) * 2 * np.pi
        
        # Initial temperature is very low
        initial_temp = 100.0  # K
        
        # Velocity magnitudes follow Maxwell-Boltzmann distribution
        vel_scale = np.sqrt(initial_temp / target_temp)  # Scale factor based on temperature ratio
        vel_magnitudes = stats.maxwell.rvs(scale=vel_scale, size=num_particles) / 5
        
        # Convert to velocity vectors
        velocities = np.column_stack([
            vel_magnitudes * np.cos(angles),
            vel_magnitudes * np.sin(angles)
        ])
        
        # Calculate temperature from velocities (simplified)
        def calculate_temperature(vels):
            # T ∝ 1/N * sum(m*v²)
            # Since mass is constant for all particles, we can simplify
            return np.mean(np.sum(vels**2, axis=1)) * 100  # Scale factor
        
        # Function to update velocities with thermostat
        def apply_thermostat(vels, current_temp, target, coupling=0.1):
            # Simple velocity rescaling (V-rescale like)
            # lambda = sqrt(T_target/T_current)
            if current_temp > 0:  # Avoid division by zero
                scale_factor = np.sqrt(target / current_temp)
                
                # Apply partial coupling
                effective_scale = 1.0 + coupling * (scale_factor - 1.0)
                return vels * effective_scale
            else:
                return vels
        
        # Precompute all frame data for smooth animation
        print("Precomputing frames...")
        
        # Store data for each frame
        all_positions = []
        all_temp_values = []
        
        # Start with initial positions and velocities
        current_positions = positions.copy()
        current_velocities = velocities.copy()
        current_temp = calculate_temperature(current_velocities)
        
        # Store initial state
        all_positions.append(current_positions.copy())
        all_temp_values.append(current_temp)
        
        # Generate all subsequent frames
        for i in range(1, num_frames):
            # Apply thermostat to velocities
            coupling_factor = min(0.1 + (i-1) * 0.01, 0.5)
            current_velocities = apply_thermostat(
                current_velocities, 
                all_temp_values[-1],
                target_temp, 
                coupling_factor
            )
            
            # Update positions
            current_positions = current_positions + current_velocities * 0.1
            current_positions = current_positions % box_size  # Periodic boundary
            
            # Calculate new temperature
            current_temp = calculate_temperature(current_velocities)
            
            # Store this frame's data
            all_positions.append(current_positions.copy())
            all_temp_values.append(current_temp)
        
        print(f"Generated {len(all_positions)} frames")
        
        # Create figure with two subplots - top for temperature, bottom for particles
        fig, (ax_temp, ax_particles) = plt.subplots(2, 1, figsize=(8, 10), 
                                                  gridspec_kw={'height_ratios': [1, 2]})
        
        # Setup temperature plot
        time_points = np.arange(num_frames)
        temp_line, = ax_temp.plot([], [], 'r-', lw=2)
        ax_temp.axhline(y=target_temp, color='k', linestyle='--', 
                      label=f'Target: {target_temp} K')
        ax_temp.set_xlim(0, num_frames)
        ax_temp.set_ylim(0, max(max(all_temp_values) * 1.1, target_temp * 1.5))
        ax_temp.set_xlabel('Time Step')
        ax_temp.set_ylabel('Temperature (K)')
        ax_temp.set_title('Temperature Evolution during NVT')
        ax_temp.grid(True, alpha=0.3)
        ax_temp.legend()
        
        # Setup particle visualization
        ax_particles.set_xlim(0, box_size)
        ax_particles.set_ylim(0, box_size)
        ax_particles.set_aspect('equal')
        ax_particles.set_xlabel('X position (nm)')
        ax_particles.set_ylabel('Y position (nm)')
        
        # Initial scatter plot (empty)
        scatter = ax_particles.scatter([], [], s=80, c=[], cmap='plasma', 
                                     vmin=0, vmax=1, alpha=0.7)
        
        # Information text
        time_text = ax_particles.text(0.02, 0.95, '', transform=ax_particles.transAxes, 
                                    fontsize=10, bbox=dict(facecolor='white', alpha=0.7))
        
        plt.tight_layout()
        
        def init():
            # Initialize with empty data
            temp_line.set_data([], [])
            scatter.set_offsets(np.zeros((0, 2)))
            scatter.set_array(np.array([]))
            time_text.set_text('')
            return temp_line, scatter, time_text
        
        def update(frame):
            # Update temperature plot
            temp_line.set_data(time_points[:frame+1], all_temp_values[:frame+1])
            
            # Update particle positions
            pos = all_positions[frame]
            
            # Calculate color based on frame progress (temperature increases)
            # Use a consistent color scheme where color intensity indicates temperature
            temp_scale = all_temp_values[frame] / target_temp  # Normalized temperature
            particle_colors = np.ones(num_particles) * temp_scale
            
            # Update particle visualization
            scatter.set_offsets(pos)
            scatter.set_array(particle_colors)
            
            # Update text
            time_text.set_text(f'Time Step: {frame}\n'
                             f'Temperature: {all_temp_values[frame]:.1f} K\n'
                             f'Target: {target_temp:.1f} K')
            
            # Update title to indicate stage of equilibration
            progress = frame / (num_frames - 1)
            if progress < 0.25:
                stage = "Initial Heating"
            elif progress < 0.75:
                stage = "Approaching Equilibrium"
            else:
                stage = "Equilibrated"
                
            ax_particles.set_title(f'NVT Thermostat: {stage}')
            
            return temp_line, scatter, time_text
        
        print("Creating animation...")
        
        # Create animation
        ani = animation.FuncAnimation(
            fig, 
            update, 
            frames=range(num_frames),
            init_func=init,
            blit=True,
            interval=250  # slower for better visibility
        )
        
        print("Saving animation...")
        
        # Save as GIF with simpler settings
        output_path = os.path.join(output_dir, "nvt_thermostat_animation.gif")
        writer = animation.PillowWriter(fps=4)
        ani.save(output_path, writer=writer, dpi=80)
        
        print(f"Animation saved to {output_path}")
        
        # Also save the final frame as a static image
        plt.savefig(os.path.join(output_dir, "thermostat_final_state.png"), dpi=300)
        
        plt.close(fig)
        
        return "Thermostat animation created successfully"
        
    except Exception as e:
        import traceback
        print(f"Error in animate_thermostat_action: {e}")
        traceback.print_exc()
        return f"Error: {e}"

def main():
    """Run all NVT visualization functions."""
    print("Creating NVT equilibration visualizations...")
    output_dir = create_output_dir()
    
    # Create main NVT temperature equilibration plot
    print("Creating NVT temperature equilibration plot...")
    plot_nvt_temperature_equilibration()
    
    # Create animations
    print("Creating NVT temperature animation...")
    try:
        animate_temperature_equilibration()
        print("Temperature animation created successfully.")
    except Exception as e:
        print(f"Warning: Could not create temperature animation: {e}")
    
    print("Creating velocity distribution animation...")
    try:
        animate_velocity_distribution()
        print("Velocity distribution animation created successfully.")
    except Exception as e:
        print(f"Warning: Could not create velocity distribution animation: {e}")
    
    # Create thermostat action animation
    print("Creating thermostat action animation...")
    try:
        animate_thermostat_action()
        print("Thermostat animation created successfully.")
    except Exception as e:
        print(f"Warning: Could not create thermostat animation: {e}")
    
    print(f"All NVT visualizations saved to {output_dir}/")

if __name__ == "__main__":
    main() 