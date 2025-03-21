#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec
import os
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
from scipy import stats
from itertools import product, combinations  # Added import for 3D box drawing

def create_output_dir():
    """Create output directory if it doesn't exist."""
    output_dir = "visualization_outputs"
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def plot_md_production():
    """Plot static visualizations of the regular MD production run."""
    output_dir = create_output_dir()
    
    # Simulation parameters
    simulation_time = 1000  # ps (1 ns)
    time_points = np.linspace(0, simulation_time, 1000)
    
    # Create a figure with multiple subplots
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 2, figure=fig)
    
    # 1. Temperature plot
    ax1 = fig.add_subplot(gs[0, 0])
    # Generate simulated temperature data with small fluctuations around 298K
    target_temp = 298.0  # K
    temp_noise = np.random.normal(0, 0.5, len(time_points))  # Small fluctuations
    temperatures = target_temp + temp_noise
    
    ax1.plot(time_points, temperatures, 'r-', linewidth=1.5)
    ax1.axhline(y=target_temp, color='k', linestyle='--', alpha=0.5)
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Temperature (K)')
    ax1.set_title('Temperature Stability')
    ax1.set_ylim(target_temp - 5, target_temp + 5)
    ax1.grid(True, alpha=0.3)
    
    # 2. Pressure plot
    ax2 = fig.add_subplot(gs[0, 1])
    # Generate simulated pressure data with fluctuations around 1 bar
    target_pressure = 1.0  # bar
    pressure_noise = np.random.normal(0, 5, len(time_points))  # Larger fluctuations for pressure
    pressures = target_pressure + pressure_noise
    
    ax2.plot(time_points, pressures, 'b-', linewidth=1.5)
    ax2.axhline(y=target_pressure, color='k', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Time (ps)')
    ax2.set_ylabel('Pressure (bar)')
    ax2.set_title('Pressure Stability')
    ax2.grid(True, alpha=0.3)
    
    # 3. Potential Energy plot
    ax3 = fig.add_subplot(gs[1, 0])
    # Generate simulated potential energy data
    potential_energy_mean = -50000  # Typical value for water systems
    potential_noise = np.random.normal(0, 100, len(time_points))
    potential_energies = potential_energy_mean + potential_noise
    
    ax3.plot(time_points, potential_energies, 'g-', linewidth=1.5)
    ax3.set_xlabel('Time (ps)')
    ax3.set_ylabel('Potential Energy (kJ/mol)')
    ax3.set_title('Potential Energy')
    ax3.grid(True, alpha=0.3)
    
    # 4. Kinetic Energy plot
    ax4 = fig.add_subplot(gs[1, 1])
    # Generate simulated kinetic energy data related to temperature
    # For 5500 water molecules = 16500 atoms, kinetic energy at 298K is around 40000 kJ/mol
    kinetic_energy_mean = 40000  
    kinetic_noise = np.random.normal(0, 80, len(time_points))
    kinetic_energies = kinetic_energy_mean + (temperatures - target_temp) * 150 + kinetic_noise
    
    ax4.plot(time_points, kinetic_energies, 'orange', linewidth=1.5)
    ax4.set_xlabel('Time (ps)')
    ax4.set_ylabel('Kinetic Energy (kJ/mol)')
    ax4.set_title('Kinetic Energy')
    ax4.grid(True, alpha=0.3)
    
    # 5. Total Energy plot
    ax5 = fig.add_subplot(gs[2, 0])
    # Calculate total energy
    total_energies = potential_energies + kinetic_energies
    
    ax5.plot(time_points, total_energies, 'purple', linewidth=1.5)
    ax5.set_xlabel('Time (ps)')
    ax5.set_ylabel('Total Energy (kJ/mol)')
    ax5.set_title('Total Energy Conservation')
    ax5.grid(True, alpha=0.3)
    
    # 6. Instantaneous RMSD plot
    ax6 = fig.add_subplot(gs[2, 1])
    # Generate simulated RMSD data
    rmsd_time = np.linspace(0, simulation_time, 100)  # Fewer points for RMSD
    rmsd_values = 0.05 * np.sqrt(rmsd_time) + np.random.normal(0, 0.01, len(rmsd_time))
    
    ax6.plot(rmsd_time, rmsd_values, 'brown', linewidth=1.5)
    ax6.set_xlabel('Time (ps)')
    ax6.set_ylabel('RMSD (nm)')
    ax6.set_title('Root Mean Square Deviation')
    ax6.grid(True, alpha=0.3)
    
    # Final layout adjustments
    plt.tight_layout()
    plt.suptitle('Regular MD Production Run Analysis', fontsize=16, y=1.02)
    
    # Save the figure
    plt.savefig(f"{output_dir}/regular_md_static_plots.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/regular_md_static_plots.pdf", bbox_inches='tight')
    
    plt.close()
    print(f"Generated static plots for regular MD production run in {output_dir}/regular_md_static_plots.png/pdf")

def animate_water_diffusion():
    """Create an animation showing water diffusion during the MD simulation."""
    output_dir = create_output_dir()
    
    # Parameters for better visualization
    n_particles = 49  # Square number for perfect grid
    n_frames = 50  # Number of frames in animation
    box_size = 4.0  # nm, size of the simulation box

    # Create a perfect grid pattern for initial positions
    particles_per_side = int(np.sqrt(n_particles))  # Should be 7 for 49 particles
    spacing = box_size / (particles_per_side + 1)
    
    positions = []
    for i in range(particles_per_side):
        for j in range(particles_per_side):
            positions.append([spacing * (i + 1), spacing * (j + 1)])
    
    # Convert to numpy array
    positions = np.array(positions)
    
    # Create frames with smooth, correlated motion
    all_positions = [positions.copy()]
    
    # Initialize velocity vectors (initial random direction, but consistent over time)
    np.random.seed(42)  # For reproducibility
    velocities = np.random.normal(0, 0.05, size=positions.shape)
    
    # Normalize velocities to have consistent initial speed
    speeds = np.sqrt(np.sum(velocities**2, axis=1)).reshape(-1, 1)
    velocities = velocities / speeds * 0.04  # Initial speed of 0.04 nm per frame
    
    # Generate all frames with smooth motion
    for frame in range(1, n_frames):
        # Start with previous positions
        new_positions = all_positions[-1].copy()
        
        # Add small random fluctuations to velocities (Brownian motion component)
        # The fluctuation gets stronger over time to show diffusion increasing
        brownian_strength = 0.002 * np.sqrt(frame)  # Gradually increase randomness
        velocity_changes = np.random.normal(0, brownian_strength, size=velocities.shape)
        velocities += velocity_changes
        
        # Apply the velocities to move particles
        new_positions += velocities
        
        # Apply periodic boundary conditions
        new_positions = new_positions % box_size
        
        # When particles cross boundaries, maintain their velocity direction
        # This prevents sudden direction changes at boundaries
        
        all_positions.append(new_positions)
    
    # Create animation with a single figure
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_xlim(0, box_size)
    ax.set_ylim(0, box_size)
    ax.set_xlabel('x (nm)')
    ax.set_ylabel('y (nm)')
    ax.set_title('Water Molecule Diffusion')
    
    # Create initial scatter plot with larger markers for better visibility
    scatter = ax.scatter(all_positions[0][:, 0], all_positions[0][:, 1], s=40, alpha=0.7, c='blue', marker='o')
    time_text = ax.text(0.05, 0.95, f'Time: 0 ps', transform=ax.transAxes, fontsize=12)
    msd_text = ax.text(0.05, 0.90, f'MSD: 0.000 nm²', transform=ax.transAxes, fontsize=12)
    
    # Draw box boundary
    box = plt.Rectangle((0, 0), box_size, box_size, fill=False, edgecolor='black')
    ax.add_patch(box)
    
    # Animation update function
    def update(frame):
        # Update scatter positions
        scatter.set_offsets(all_positions[frame])
        
        # Calculate time - assuming 10 ps per frame
        time = frame * 10
        time_text.set_text(f'Time: {time:.0f} ps')
        
        # Calculate actual MSD
        initial_pos = all_positions[0]
        current_pos = all_positions[frame]
        
        # Calculate displacements considering periodic boundary conditions
        dx = current_pos[:, 0] - initial_pos[:, 0]
        dy = current_pos[:, 1] - initial_pos[:, 1]
        
        # Correct for boundary crossings (particles that crossed the boundary)
        dx = np.where(dx > box_size/2, dx - box_size, dx)
        dx = np.where(dx < -box_size/2, dx + box_size, dx)
        dy = np.where(dy > box_size/2, dy - box_size, dy)
        dy = np.where(dy < -box_size/2, dy + box_size, dy)
        
        # Calculate MSD properly
        squared_displacement = dx**2 + dy**2
        msd = np.mean(squared_displacement)
        msd_text.set_text(f'MSD: {msd:.3f} nm²')
        
        return scatter, time_text, msd_text
    
    # Generate animation with smoother frame rate
    ani = animation.FuncAnimation(fig, update, frames=n_frames, blit=True, interval=200)
    
    # Save animation with higher quality
    ani.save(f"{output_dir}/water_diffusion.gif", writer='pillow', fps=8, dpi=120)
    
    plt.close()
    print(f"Generated water diffusion animation in {output_dir}/water_diffusion.gif")

def animate_trajectory_evolution():
    """Create an animation showing the evolution of a water trajectory from different angles."""
    output_dir = create_output_dir()
    
    # Parameters - reduced for faster animation
    box_size = 4.0  # nm
    n_molecules = 50  # Reduced number of water molecules
    n_frames = 30  # Reduced number of frames
    rotation_degrees = 360  # Total rotation degrees
    
    # Create figure
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Generate initial positions
    np.random.seed(42)
    positions = np.random.uniform(0, box_size, size=(n_molecules, 3))
    
    # Create trajectory data - pre-calculate all frames
    trajectories = []
    current_positions = positions.copy()
    
    for frame in range(n_frames):
        # Add random motion (diffusion)
        disp = np.random.normal(0, 0.05, size=(n_molecules, 3))  # Reduced displacement size
        current_positions += disp
        
        # Apply periodic boundary conditions
        current_positions = current_positions % box_size
        
        trajectories.append(current_positions.copy())
    
    # Set up the 3D box once
    ax.set_xlim(0, box_size)
    ax.set_ylim(0, box_size)
    ax.set_zlim(0, box_size)
    ax.set_xlabel('X (nm)')
    ax.set_ylabel('Y (nm)')
    ax.set_zlabel('Z (nm)')
    ax.set_title('Water Molecule Trajectories in 3D')
    
    # Draw box
    r = [0, box_size]
    for s, e in combinations(np.array(list(product(r, r, r))), 2):
        if np.sum(np.abs(s-e)) == r[1]-r[0]:
            ax.plot3D(*zip(s, e), color="gray", alpha=0.2)
            
    # Create initial scatter plot
    scatter = ax.scatter(
        trajectories[0][:, 0],
        trajectories[0][:, 1],
        trajectories[0][:, 2],
        c='blue', s=20, alpha=0.6
    )
    
    # Create text for timestep
    time_text = ax.text2D(0.05, 0.95, f'Time: 0 ps', transform=ax.transAxes, fontsize=12)
    
    # Simplified update function
    def update(frame):
        # Update scatter positions
        scatter._offsets3d = (
            trajectories[frame][:, 0],
            trajectories[frame][:, 1],
            trajectories[frame][:, 2]
        )
        
        # Update view angle
        angle = (frame / n_frames) * rotation_degrees
        ax.view_init(elev=30, azim=angle)
        
        # Update time text
        time_text.set_text(f'Time: {frame*10:.0f} ps')
        
        return scatter, time_text
    
    # Generate animation with simpler settings
    ani = animation.FuncAnimation(fig, update, frames=n_frames, blit=False)
    
    # Save animation with lower quality for faster generation
    ani.save(f"{output_dir}/trajectory_3d.gif", writer='pillow', fps=6, dpi=80)
    
    plt.close()
    print(f"Generated 3D trajectory animation in {output_dir}/trajectory_3d.gif")

def main():
    """Run all regular MD visualizations."""
    print("Generating Regular MD Visualizations...")
    
    # Create static plots
    plot_md_production()
    
    # Create animations
    animate_water_diffusion()
    animate_trajectory_evolution()
    
    print("All Regular MD visualizations complete!")

if __name__ == "__main__":
    main() 