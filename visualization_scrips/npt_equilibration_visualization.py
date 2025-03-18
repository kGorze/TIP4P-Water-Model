#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec
import os
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
from scipy import stats

def create_output_dir():
    """Create output directory if it doesn't exist."""
    output_dir = "visualization_outputs"
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def plot_npt_equilibration():
    """Plot a visualization of pressure and volume equilibration during NPT."""
    output_dir = create_output_dir()
    
    # Simulation parameters
    target_pressure = 1.0  # bar (standard pressure)
    target_temp = 298.0    # K (room temperature)
    simulation_time = 500  # ps
    
    # Create time points
    time_points = np.linspace(0, simulation_time, 500)
    
    # Generate pressure trajectory
    # Start with high pressure due to initial configuration
    initial_pressure = 500.0  # bar (very high)
    pressure_noise_amplitude = 50.0  # bar
    
    pressures = []
    for t in time_points:
        # Exponential decay of pressure to target with noise
        decay_factor = np.exp(-t / 100)
        noise_level = pressure_noise_amplitude * decay_factor
        
        # Add both decaying noise and persistent equilibrium fluctuations
        noise = np.random.normal(0, noise_level)
        equil_fluct = np.random.normal(0, 5.0)  # Equilibrium pressure fluctuations
        
        # Pressure converges to target
        pressure = target_pressure + (initial_pressure - target_pressure) * np.exp(-t / 80) + noise + equil_fluct
        pressures.append(pressure)
    
    pressures = np.array(pressures)
    
    # Generate volume trajectory (inversely related to pressure)
    # Starting with lower volume and expanding
    initial_volume = 90.0  # nm^3
    target_volume = 100.0  # nm^3
    volume_noise_amplitude = 2.0  # nm^3
    
    volumes = []
    for t in time_points:
        # Volume relaxation with noise
        decay_factor = np.exp(-t / 120)  # Slower equilibration than pressure
        noise = np.random.normal(0, volume_noise_amplitude * decay_factor)
        equil_fluct = np.random.normal(0, 0.2)  # Small equilibrium volume fluctuations
        
        # Volume converges to target
        volume = target_volume + (initial_volume - target_volume) * np.exp(-t / 100) + noise + equil_fluct
        volumes.append(volume)
    
    volumes = np.array(volumes)
    
    # Generate density trajectory (derived from volume)
    # For simplicity, assuming a fixed number of particles
    initial_density = 1000.0  # kg/m^3 (example typical value for water)
    densities = initial_density * initial_volume / volumes
    
    # Create figure
    fig = plt.figure(figsize=(14, 12))
    gs = GridSpec(5, 2, figure=fig, height_ratios=[2, 2, 2, 2, 3])
    
    # --- Pressure Plot ---
    ax_pressure = fig.add_subplot(gs[0, :])
    ax_pressure.plot(time_points, pressures, 'r-', linewidth=1.5, alpha=0.7)
    ax_pressure.axhline(y=target_pressure, color='k', linestyle='--', linewidth=1.5, 
                       label=f'Target Pressure ({target_pressure} bar)')
    
    # Add shaded region for equilibrated portion
    equil_start = 250  # ps
    equil_idx = np.where(time_points >= equil_start)[0][0]
    ax_pressure.axvspan(equil_start, simulation_time, alpha=0.2, color='green', 
                       label='Equilibrated Region')
    
    # Calculate statistics for equilibrated pressure
    equil_press = pressures[equil_idx:]
    press_mean = np.mean(equil_press)
    press_std = np.std(equil_press)
    
    # Add statistics to plot
    stats_text = f"Equilibrated Pressure: {press_mean:.2f} ± {press_std:.2f} bar"
    ax_pressure.text(0.02, 0.05, stats_text, transform=ax_pressure.transAxes, 
                    fontsize=10, bbox=dict(facecolor='white', alpha=0.7))
    
    # Set labels and title
    ax_pressure.set_xlabel('Simulation Time (ps)')
    ax_pressure.set_ylabel('Pressure (bar)')
    ax_pressure.set_title('NPT Equilibration: Pressure vs. Time')
    ax_pressure.legend(loc='upper right')
    ax_pressure.grid(True, alpha=0.3)
    
    # --- Volume Plot ---
    ax_volume = fig.add_subplot(gs[1, :])
    ax_volume.plot(time_points, volumes, 'b-', linewidth=1.5, alpha=0.7)
    ax_volume.axhline(y=target_volume, color='k', linestyle='--', linewidth=1.5, 
                     label=f'Equilibrium Volume ({target_volume} nm³)')
    
    # Add shaded region for equilibrated portion
    ax_volume.axvspan(equil_start, simulation_time, alpha=0.2, color='green')
    
    # Calculate statistics for equilibrated volume
    equil_vol = volumes[equil_idx:]
    vol_mean = np.mean(equil_vol)
    vol_std = np.std(equil_vol)
    
    # Add statistics to plot
    stats_text = f"Equilibrated Volume: {vol_mean:.2f} ± {vol_std:.2f} nm³"
    ax_volume.text(0.02, 0.05, stats_text, transform=ax_volume.transAxes, 
                  fontsize=10, bbox=dict(facecolor='white', alpha=0.7))
    
    # Set labels and title
    ax_volume.set_xlabel('Simulation Time (ps)')
    ax_volume.set_ylabel('System Volume (nm³)')
    ax_volume.set_title('NPT Equilibration: Volume vs. Time')
    ax_volume.legend(loc='lower right')
    ax_volume.grid(True, alpha=0.3)
    
    # --- Density Plot ---
    ax_density = fig.add_subplot(gs[2, :])
    ax_density.plot(time_points, densities, 'g-', linewidth=1.5, alpha=0.7)
    
    # Calculate equilibrium density and add reference line
    equil_density = np.mean(densities[equil_idx:])
    ax_density.axhline(y=equil_density, color='k', linestyle='--', linewidth=1.5, 
                      label=f'Equilibrium Density ({equil_density:.1f} kg/m³)')
    
    # Add shaded region for equilibrated portion
    ax_density.axvspan(equil_start, simulation_time, alpha=0.2, color='green')
    
    # Calculate statistics for equilibrated density
    equil_dens = densities[equil_idx:]
    dens_mean = np.mean(equil_dens)
    dens_std = np.std(equil_dens)
    
    # Add statistics to plot
    stats_text = f"Equilibrated Density: {dens_mean:.2f} ± {dens_std:.2f} kg/m³"
    ax_density.text(0.02, 0.05, stats_text, transform=ax_density.transAxes, 
                   fontsize=10, bbox=dict(facecolor='white', alpha=0.7))
    
    # Set labels and title
    ax_density.set_xlabel('Simulation Time (ps)')
    ax_density.set_ylabel('Density (kg/m³)')
    ax_density.set_title('NPT Equilibration: Density vs. Time')
    ax_density.legend(loc='upper right')
    ax_density.grid(True, alpha=0.3)
    
    # --- Box Scaling Visualization ---
    ax_box = fig.add_subplot(gs[3, :])
    
    # Create a visualization of box size changes
    num_boxes = 8
    box_times = np.linspace(0, simulation_time, num_boxes)
    box_volumes = np.interp(box_times, time_points, volumes)
    
    # Calculate box dimensions (assuming cubic box for simplicity)
    # Normalize box sizes for visualization
    box_sizes = np.power(box_volumes / target_volume, 1/3)
    
    # Plot boxes at different time points
    box_positions = np.linspace(0, 10, num_boxes)
    box_colors = plt.cm.viridis(np.linspace(0, 1, num_boxes))
    
    for i, (pos, size, vol, t, color) in enumerate(zip(box_positions, box_sizes, box_volumes, box_times, box_colors)):
        # Draw box
        rect = plt.Rectangle((pos - size/2, 0), size, size, color=color, alpha=0.7)
        ax_box.add_patch(rect)
        
        # Add time and volume labels
        ax_box.text(pos, size + 0.05, f"{t:.0f} ps", ha='center', va='bottom', fontsize=8)
        ax_box.text(pos, -0.05, f"{vol:.1f} nm³", ha='center', va='top', fontsize=8)
    
    # Add a color bar to indicate time progression
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(0, simulation_time))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax_box, orientation='horizontal', pad=0.2)
    cbar.set_label('Simulation Time (ps)')
    
    ax_box.set_xlim(-1, 11)
    ax_box.set_ylim(-0.2, 1.5)
    ax_box.set_aspect('equal')
    ax_box.set_title('NPT Equilibration: Simulation Box Size Evolution')
    ax_box.set_xticks([])
    ax_box.set_yticks([])
    
    # --- PV and Compressibility Relationship ---
    ax_pv = fig.add_subplot(gs[4, :])
    
    # Create scatter plot of pressure vs volume
    ax_pv.scatter(volumes, pressures, c=time_points, cmap='viridis', alpha=0.7)
    
    # Add arrow to show direction of time
    mid_idx = len(time_points) // 4
    ax_pv.annotate('', xy=(volumes[mid_idx*3], pressures[mid_idx*3]), 
                  xytext=(volumes[mid_idx], pressures[mid_idx]),
                  arrowprops=dict(arrowstyle='->', lw=2, color='red'))
    ax_pv.text(volumes[mid_idx*2], pressures[mid_idx*2], 'Time', 
              color='red', fontsize=12, ha='center', va='bottom')
    
    # Calculate and plot isotherm (PV=constant line)
    # For an ideal system at constant temperature
    v_range = np.linspace(np.min(volumes)*0.95, np.max(volumes)*1.05, 100)
    equil_p_v_product = press_mean * vol_mean
    p_ideal = equil_p_v_product / v_range
    
    ax_pv.plot(v_range, p_ideal, 'k--', label='Ideal Gas Isotherm (PV=const)')
    
    # Mark equilibrium point
    ax_pv.plot(vol_mean, press_mean, 'ro', markersize=8, label='Equilibrium State')
    
    # Add colorbar for time
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(0, simulation_time))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax_pv)
    cbar.set_label('Simulation Time (ps)')
    
    ax_pv.set_xlabel('Volume (nm³)')
    ax_pv.set_ylabel('Pressure (bar)')
    ax_pv.set_title('Pressure-Volume Relationship During NPT Equilibration')
    ax_pv.legend(loc='upper right')
    ax_pv.grid(True, alpha=0.3)
    
    # Adjust layout to fit all elements
    plt.tight_layout()
    
    # Save figure
    plt.savefig(os.path.join(output_dir, "npt_equilibration.png"), dpi=300)
    plt.savefig(os.path.join(output_dir, "npt_equilibration.pdf"))
    
    return fig, pressures, volumes, densities, time_points, equil_start

def visualize_box_fluctuations():
    """Create a visualization of box size fluctuations during NPT."""
    output_dir = create_output_dir()
    
    # Create figure for visualization
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Simulation parameters
    simulation_steps = 200
    
    # Initial box dimensions (nm)
    box_x = 5.0
    box_y = 5.0
    box_z = 5.0
    
    # Target pressure and compressibility
    target_pressure = 1.0  # bar
    compressibility = 4.5e-5  # bar^-1 (typical for water)
    
    # Arrays to store box dimensions and pressure values
    steps = np.arange(simulation_steps)
    pressures = []
    box_dimensions = []
    
    # Starting pressure - high pressure to see box adjustments
    current_pressure = 200.0  # bar
    
    # Coupling constant for pressure adjustment (like tau_p in MD settings)
    pressure_coupling = 0.1
    
    # Pressure adjustment damping factor (increases over time to simulate equilibration)
    damping_factor = 0.2
    
    # Generate pressure and box size trajectory
    for step in range(simulation_steps):
        # Store current values
        pressures.append(current_pressure)
        box_dimensions.append((box_x, box_y, box_z))
        
        # Calculate pressure adjustment
        pressure_error = current_pressure - target_pressure
        
        # Adjust box based on pressure error (simplified barostat)
        # Box scales with cube root of pressure change and compressibility
        # For isotropic system: dV/V = -compressibility * dP
        volume_scale_factor = 1.0 - compressibility * pressure_error * pressure_coupling * damping_factor
        
        # Convert to linear scaling factor (cube root for isotropic scaling)
        linear_scale = volume_scale_factor**(1/3)
        
        # Apply scaling to box dimensions
        box_x *= linear_scale
        box_y *= linear_scale
        box_z *= linear_scale
        
        # Update pressure based on new volume (simple inverse relationship)
        # P1V1 = P2V2 is approximate for small changes
        # Add some random noise to simulate fluctuations
        old_volume = box_x * box_y * box_z / (linear_scale**3)
        new_volume = box_x * box_y * box_z
        current_pressure = current_pressure * old_volume / new_volume
        
        # Add random pressure fluctuations (decreasing over time)
        fluctuation_scale = 10.0 * np.exp(-step / 50)
        current_pressure += np.random.normal(0, fluctuation_scale)
        
        # Ensure pressure stays positive
        current_pressure = max(current_pressure, 0.1)
        
        # Increase damping factor over time to approach equilibrium more smoothly
        damping_factor = min(0.2 + step * 0.005, 1.0)
    
    # Convert to arrays
    pressures = np.array(pressures)
    box_dimensions = np.array(box_dimensions)
    volumes = box_dimensions[:, 0] * box_dimensions[:, 1] * box_dimensions[:, 2]
    
    # Create figure
    fig = plt.figure(figsize=(14, 12))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[2, 2, 3])
    
    # --- Pressure Plot ---
    ax_pressure = fig.add_subplot(gs[0, :])
    ax_pressure.plot(steps, pressures, 'r-', linewidth=1.5, alpha=0.7)
    ax_pressure.axhline(y=target_pressure, color='k', linestyle='--', 
                      label=f'Target Pressure ({target_pressure} bar)')
    
    # Mark equilibration point
    equil_point = simulation_steps // 2
    ax_pressure.axvspan(equil_point, simulation_steps, alpha=0.2, color='green', 
                      label='Equilibrated Region')
    
    ax_pressure.set_xlabel('Simulation Step')
    ax_pressure.set_ylabel('Pressure (bar)')
    ax_pressure.set_title('Barostat Action: Pressure Equilibration')
    ax_pressure.legend()
    ax_pressure.grid(True, alpha=0.3)
    
    # --- Volume Plot ---
    ax_volume = fig.add_subplot(gs[1, :])
    ax_volume.plot(steps, volumes, 'b-', linewidth=1.5, alpha=0.7)
    
    # Calculate equilibrium volume
    equil_volume = np.mean(volumes[equil_point:])
    ax_volume.axhline(y=equil_volume, color='k', linestyle='--', 
                    label=f'Equilibrium Volume ({equil_volume:.2f} nm³)')
    
    ax_volume.axvspan(equil_point, simulation_steps, alpha=0.2, color='green')
    ax_volume.set_xlabel('Simulation Step')
    ax_volume.set_ylabel('Volume (nm³)')
    ax_volume.set_title('Barostat Action: Volume Adjustment')
    ax_volume.legend()
    ax_volume.grid(True, alpha=0.3)
    
    # --- Box Visualization ---
    ax_box = fig.add_subplot(gs[2, :], projection='3d')
    
    # Choose time points to visualize
    selected_steps = [0, 10, 25, 50, 100, 199]
    colors = plt.cm.viridis(np.linspace(0, 1, len(selected_steps)))
    
    # Plot each box state
    for i, step in enumerate(selected_steps):
        # Get dimensions from stored values
        x, y, z = box_dimensions[step]
        color = colors[i]
        
        # Create corner points for the cuboid
        corners = np.array([
            [0, 0, 0], [x, 0, 0], [x, y, 0], [0, y, 0],
            [0, 0, z], [x, 0, z], [x, y, z], [0, y, z]
        ])
        
        # Define faces using corner indices
        faces = [
            [0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4],
            [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]
        ]
        
        # Set position to offset boxes from each other
        offset = i * x * 0.3
        corners[:, 0] += offset
        
        # Plot each face
        for face in faces:
            # Get coordinates of this face
            x_face = [corners[idx, 0] for idx in face]
            y_face = [corners[idx, 1] for idx in face]
            z_face = [corners[idx, 2] for idx in face]
            
            # Close the loop for the face
            x_face.append(x_face[0])
            y_face.append(y_face[0])
            z_face.append(z_face[0])
            
            # Plot the face
            ax_box.plot3D(x_face, y_face, z_face, color=color, alpha=0.7)
        
        # Add label with step and pressure
        ax_box.text(offset + x/2, y/2, z*1.1, 
                  f"Step {step}\nP: {pressures[step]:.1f} bar\nV: {volumes[step]:.1f} nm³", 
                  ha='center', fontsize=9)
    
    # Add color bar for time progression
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(0, simulation_steps))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax_box)
    cbar.set_label('Simulation Time')
    
    # Set equal aspect ratio for 3D plot (compatible approach)
    ax_box.set_xlabel('X (nm)')
    ax_box.set_ylabel('Y (nm)') 
    ax_box.set_zlabel('Z (nm)')
    # Calculate the plot bounds
    max_range = np.array([
        corners[:, 0].max() - corners[:, 0].min(),
        corners[:, 1].max() - corners[:, 1].min(),
        corners[:, 2].max() - corners[:, 2].min()
    ]).max() / 2.0
    
    # Set axes limits to ensure equal scaling
    mid_x = np.mean([corners[:, 0].min(), corners[:, 0].max()])
    mid_y = np.mean([corners[:, 1].min(), corners[:, 1].max()])
    mid_z = np.mean([corners[:, 2].min(), corners[:, 2].max()])
    
    ax_box.set_xlim(mid_x - max_range * 3, mid_x + max_range * 3)  # Wider in x to show sequence
    ax_box.set_ylim(mid_y - max_range, mid_y + max_range)
    ax_box.set_zlim(mid_z - max_range, mid_z + max_range)

    ax_box.set_title('Simulation Box Evolution During NPT Equilibration')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    plt.savefig(os.path.join(output_dir, "npt_barostat_visualization.png"), dpi=300)
    plt.savefig(os.path.join(output_dir, "npt_barostat_visualization.pdf"))
    
    return fig, pressures, volumes, box_dimensions

def create_npt_animation():
    """Create an animation showing NPT equilibration process."""
    output_dir = create_output_dir()
    
    # Simulation parameters
    simulation_time = 500  # ps
    time_points = np.linspace(0, simulation_time, 500)
    
    # Generate data using the same method as plot_npt_equilibration
    target_pressure = 1.0  # bar
    target_volume = 100.0  # nm^3
    initial_pressure = 500.0  # bar
    initial_volume = 90.0  # nm^3
    pressure_noise_amplitude = 50.0  # bar
    volume_noise_amplitude = 2.0  # nm^3
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Generate pressure trajectory
    pressures = []
    for t in time_points:
        decay_factor = np.exp(-t / 100)
        noise_level = pressure_noise_amplitude * decay_factor
        noise = np.random.normal(0, noise_level)
        equil_fluct = np.random.normal(0, 5.0)
        pressure = target_pressure + (initial_pressure - target_pressure) * np.exp(-t / 80) + noise + equil_fluct
        pressures.append(pressure)
    
    pressures = np.array(pressures)
    
    # Generate volume trajectory
    volumes = []
    for t in time_points:
        decay_factor = np.exp(-t / 120)
        noise = np.random.normal(0, volume_noise_amplitude * decay_factor)
        equil_fluct = np.random.normal(0, 0.2)
        volume = target_volume + (initial_volume - target_volume) * np.exp(-t / 100) + noise + equil_fluct
        volumes.append(volume)
    
    volumes = np.array(volumes)
    
    # Generate density trajectory
    initial_density = 1000.0  # kg/m^3
    densities = initial_density * initial_volume / volumes
    
    # Create animation
    fig = plt.figure(figsize=(12, 9))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[1, 1, 1])
    
    # Pressure subplot
    ax_pressure = fig.add_subplot(gs[0, :])
    pressure_line, = ax_pressure.plot([], [], 'r-', linewidth=1.5)
    ax_pressure.axhline(y=target_pressure, color='k', linestyle='--', linewidth=1.5, 
                       label=f'Target Pressure ({target_pressure} bar)')
    ax_pressure.set_xlim(0, simulation_time)
    ax_pressure.set_ylim(0, max(pressures) * 1.1)
    ax_pressure.set_xlabel('Simulation Time (ps)')
    ax_pressure.set_ylabel('Pressure (bar)')
    ax_pressure.set_title('NPT Equilibration: Pressure vs. Time')
    ax_pressure.legend(loc='upper right')
    ax_pressure.grid(True, alpha=0.3)
    
    # Volume subplot
    ax_volume = fig.add_subplot(gs[1, :])
    volume_line, = ax_volume.plot([], [], 'b-', linewidth=1.5)
    ax_volume.axhline(y=target_volume, color='k', linestyle='--', linewidth=1.5, 
                     label=f'Target Volume ({target_volume} nm³)')
    ax_volume.set_xlim(0, simulation_time)
    ax_volume.set_ylim(min(volumes) * 0.95, max(volumes) * 1.05)
    ax_volume.set_xlabel('Simulation Time (ps)')
    ax_volume.set_ylabel('Volume (nm³)')
    ax_volume.set_title('NPT Equilibration: Volume vs. Time')
    ax_volume.legend(loc='lower right')
    ax_volume.grid(True, alpha=0.3)
    
    # Visualization of the box
    ax_box = fig.add_subplot(gs[2, :])
    ax_box.set_xlim(-1, 2)
    ax_box.set_ylim(-1, 2)
    ax_box.set_aspect('equal')
    ax_box.set_title('Simulation Box Size')
    ax_box.set_xticks([])
    ax_box.set_yticks([])
    
    # Create the box patch
    box_patch = plt.Rectangle((-0.5, -0.5), 1, 1, color='b', alpha=0.5)
    ax_box.add_patch(box_patch)
    
    # Text for displaying current values
    time_text = ax_box.text(0.02, 0.95, '', transform=ax_box.transAxes)
    pressure_text = ax_box.text(0.02, 0.90, '', transform=ax_box.transAxes)
    volume_text = ax_box.text(0.02, 0.85, '', transform=ax_box.transAxes)
    density_text = ax_box.text(0.02, 0.80, '', transform=ax_box.transAxes)
    
    def init():
        pressure_line.set_data([], [])
        volume_line.set_data([], [])
        time_text.set_text('')
        pressure_text.set_text('')
        volume_text.set_text('')
        density_text.set_text('')
        box_patch.set_width(1)
        box_patch.set_height(1)
        box_patch.set_xy((-0.5, -0.5))
        return pressure_line, volume_line, time_text, pressure_text, volume_text, density_text, box_patch
    
    def update(frame):
        # Update the lines
        pressure_line.set_data(time_points[:frame], pressures[:frame])
        volume_line.set_data(time_points[:frame], volumes[:frame])
        
        # Update the texts
        time_text.set_text(f'Time: {time_points[frame]:.1f} ps')
        pressure_text.set_text(f'Pressure: {pressures[frame]:.1f} bar')
        volume_text.set_text(f'Volume: {volumes[frame]:.1f} nm³')
        density_text.set_text(f'Density: {densities[frame]:.1f} kg/m³')
        
        # Update the box size
        # Calculate box dimensions (assuming cubic box for simplicity)
        box_scale = (volumes[frame] / target_volume) ** (1/3)
        box_patch.set_width(box_scale)
        box_patch.set_height(box_scale)
        box_patch.set_xy((-box_scale/2, -box_scale/2))
        
        # Color based on time progression
        color = plt.cm.viridis(frame / len(time_points))
        box_patch.set_color(color)
        
        return pressure_line, volume_line, time_text, pressure_text, volume_text, density_text, box_patch
    
    # Create animation
    ani = animation.FuncAnimation(fig, update, frames=range(0, len(time_points), 5),
                                  init_func=init, blit=True, interval=50)
    
    # Save animation
    writer = animation.PillowWriter(fps=20)
    ani.save(os.path.join(output_dir, "npt_equilibration_animation.gif"), writer=writer, dpi=100)
    
    plt.close(fig)
    
    # Create a 3D animation of box evolution
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Calculate box dimensions for each frame
    box_sizes = np.power(volumes / target_volume, 1/3)
    
    def init_3d():
        ax.clear()
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_zlim(-1, 1)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('NPT Equilibration: Box Size Evolution')
        return []
    
    def update_3d(frame):
        ax.clear()
        
        # Set limits
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_zlim(-1, 1)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        # Scale box dimensions
        box_size = box_sizes[frame*5]  # Using every 5th frame
        
        # Create corner points for the cube
        r = box_size / 2  # half-width of the box
        corners = np.array([
            [-r, -r, -r], [r, -r, -r], [r, r, -r], [-r, r, -r],
            [-r, -r, r], [r, -r, r], [r, r, r], [-r, r, r]
        ])
        
        # Define faces using corner indices
        faces = [
            [0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4],
            [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]
        ]
        
        # Color based on time progression
        color = plt.cm.viridis(frame / (len(time_points) // 5))
        
        # Plot each face
        for face in faces:
            x = [corners[i, 0] for i in face]
            y = [corners[i, 1] for i in face]
            z = [corners[i, 2] for i in face]
            # Close the polygon
            x.append(x[0])
            y.append(y[0])
            z.append(z[0])
            ax.plot3D(x, y, z, color=color, alpha=0.7)
        
        # Show pressure and volume info
        time_val = time_points[frame*5]
        pressure_val = pressures[frame*5]
        volume_val = volumes[frame*5]
        
        ax.set_title(f'Time: {time_val:.1f} ps, P: {pressure_val:.1f} bar, V: {volume_val:.1f} nm³')
        
        return []
    
    # Create 3D animation
    ani_3d = animation.FuncAnimation(fig, update_3d, frames=range(len(time_points) // 5),
                                     init_func=init_3d, blit=False, interval=100)
    
    # Save 3D animation
    ani_3d.save(os.path.join(output_dir, "npt_box_evolution_3d.gif"), writer=writer, dpi=100)
    
    plt.close(fig)
    
    return "Animations created"

def create_pv_curve_animation():
    """Create an animation showing the pressure-volume curve during NPT equilibration."""
    output_dir = create_output_dir()
    
    # Generate data - reusing the same method as plot_npt_equilibration
    simulation_time = 500  # ps
    time_points = np.linspace(0, simulation_time, 500)
    
    # Set parameters
    target_pressure = 1.0  # bar
    target_volume = 100.0  # nm^3
    initial_pressure = 500.0  # bar
    initial_volume = 90.0  # nm^3
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Generate pressure and volume trajectories
    pressures = []
    volumes = []
    
    for t in time_points:
        # Pressure
        decay_factor_p = np.exp(-t / 100)
        noise_p = np.random.normal(0, 50.0 * decay_factor_p)
        pressure = target_pressure + (initial_pressure - target_pressure) * np.exp(-t / 80) + noise_p
        pressures.append(pressure)
        
        # Volume
        decay_factor_v = np.exp(-t / 120)
        noise_v = np.random.normal(0, 2.0 * decay_factor_v)
        volume = target_volume + (initial_volume - target_volume) * np.exp(-t / 100) + noise_v
        volumes.append(volume)
    
    pressures = np.array(pressures)
    volumes = np.array(volumes)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Set up the plot
    ax.set_xlim(min(volumes) * 0.95, max(volumes) * 1.05)
    ax.set_ylim(0, max(pressures) * 1.1)
    ax.set_xlabel('Volume (nm³)')
    ax.set_ylabel('Pressure (bar)')
    ax.set_title('NPT Equilibration: Pressure-Volume Relationship')
    ax.grid(True, alpha=0.3)
    
    # Calculate equilibrium point
    equil_start = 250  # ps
    equil_idx = np.where(time_points >= equil_start)[0][0]
    press_mean = np.mean(pressures[equil_idx:])
    vol_mean = np.mean(volumes[equil_idx:])
    
    # Plot isotherm (PV=constant line)
    v_range = np.linspace(min(volumes) * 0.95, max(volumes) * 1.05, 100)
    equil_p_v_product = press_mean * vol_mean
    p_ideal = equil_p_v_product / v_range
    ax.plot(v_range, p_ideal, 'k--', label='Ideal Gas Isotherm (PV=const)')
    
    # Highlight equilibrium point
    ax.plot(vol_mean, press_mean, 'ro', markersize=8, label='Equilibrium State')
    
    # Initial scatter point
    scatter = ax.scatter([], [], c=[], cmap='viridis', s=50, alpha=0.8)
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
    
    # Create colorbar
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(0, simulation_time))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Simulation Time (ps)')
    
    # Animation functions
    def init():
        scatter.set_offsets(np.empty((0, 2)))
        scatter.set_array(np.array([]))
        time_text.set_text('')
        return scatter, time_text
    
    def update(frame):
        # Create data for scatter plot - show points up to current frame
        data = np.column_stack([volumes[:frame+1], pressures[:frame+1]])
        scatter.set_offsets(data)
        
        # Update colors
        scatter.set_array(time_points[:frame+1])
        
        # Update time text
        time_text.set_text(f'Time: {time_points[frame]:.1f} ps')
        
        return scatter, time_text
    
    # Create animation
    ani = animation.FuncAnimation(fig, update, frames=range(0, len(time_points), 3),
                                  init_func=init, blit=True, interval=30)
    
    # Add legend
    ax.legend()
    
    # Save animation
    writer = animation.PillowWriter(fps=20)
    ani.save(os.path.join(output_dir, "npt_pv_curve_animation.gif"), writer=writer, dpi=120)
    
    plt.close(fig)
    
    return "PV curve animation created"

def main():
    """Run all NPT visualization functions."""
    print("Creating NPT equilibration visualizations...")
    output_dir = create_output_dir()
    
    # Create main NPT equilibration plots
    print("Creating NPT equilibration plots...")
    plot_npt_equilibration()
    
    # Create barostat visualization
    print("Creating box fluctuation visualization...")
    try:
        visualize_box_fluctuations()
    except Exception as e:
        print(f"Warning: Could not create box fluctuation visualization: {e}")
        print("Main NPT plots were still created successfully.")
    
    # Create animations
    print("Creating NPT equilibration animations...")
    try:
        create_npt_animation()
        create_pv_curve_animation()
        print("NPT animations created successfully.")
    except Exception as e:
        print(f"Warning: Could not create animations: {e}")
    
    print(f"All NPT visualizations saved to {output_dir}/")

if __name__ == "__main__":
    main() 