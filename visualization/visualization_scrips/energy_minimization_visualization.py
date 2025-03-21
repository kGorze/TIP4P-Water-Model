#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import os
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec

def create_output_dir():
    """Create output directory if it doesn't exist."""
    output_dir = "visualization_outputs"
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def generate_initial_configuration(num_particles=100, box_size=10, displacement=2.0):
    """Generate an initial non-equilibrium configuration for minimization."""
    # Create a grid-like arrangement but with some randomization
    positions = []
    grid_points = int(np.ceil(num_particles**(1/3)))
    spacing = box_size / grid_points
    
    for i in range(grid_points):
        for j in range(grid_points):
            for k in range(grid_points):
                if len(positions) < num_particles:
                    # Add random displacement to create high-energy starting configuration
                    pos = [
                        i * spacing + np.random.uniform(-displacement, displacement),
                        j * spacing + np.random.uniform(-displacement, displacement),
                        k * spacing + np.random.uniform(-displacement, displacement)
                    ]
                    positions.append(pos)
    
    return np.array(positions)

def lennard_jones_potential(r, epsilon=1.0, sigma=1.0):
    """Calculate the Lennard-Jones potential."""
    if r == 0:
        return float('inf')
    sr6 = (sigma / r)**6
    sr12 = sr6**2
    return 4 * epsilon * (sr12 - sr6)

def calculate_system_energy(positions, box_size, cutoff=2.5):
    """Calculate the total potential energy of the system."""
    num_particles = len(positions)
    total_energy = 0.0
    
    # Apply minimum image convention for periodic boundary conditions
    for i in range(num_particles):
        for j in range(i+1, num_particles):
            r_ij = positions[i] - positions[j]
            
            # Minimum image convention (periodic boundary)
            r_ij = r_ij - box_size * np.round(r_ij / box_size)
            
            # Calculate distance
            distance = np.sqrt(np.sum(r_ij**2))
            
            # Apply cutoff
            if distance < cutoff:
                # Calculate LJ potential
                pair_energy = lennard_jones_potential(distance)
                total_energy += pair_energy
    
    return total_energy

def calculate_forces(positions, box_size, cutoff=2.5):
    """Calculate forces on all particles."""
    num_particles = len(positions)
    forces = np.zeros_like(positions)
    
    # Apply minimum image convention and calculate forces
    for i in range(num_particles):
        for j in range(num_particles):
            if i != j:
                r_ij = positions[i] - positions[j]
                
                # Minimum image convention (periodic boundary)
                r_ij = r_ij - box_size * np.round(r_ij / box_size)
                
                # Calculate distance
                distance = np.sqrt(np.sum(r_ij**2))
                
                # Apply cutoff
                if distance < cutoff:
                    # Calculate LJ force (negative gradient of potential)
                    # F = -dU/dr
                    epsilon = 1.0
                    sigma = 1.0
                    sr6 = (sigma / distance)**6
                    sr12 = sr6**2
                    magnitude = 24 * epsilon * (2 * sr12 - sr6) / distance**2
                    
                    # Direction is along r_ij
                    forces[i] += magnitude * r_ij / distance
    
    return forces

def steepest_descent_step(positions, forces, step_size, box_size):
    """Perform one step of steepest descent minimization."""
    # Calculate the maximum force magnitude
    force_magnitudes = np.sqrt(np.sum(forces**2, axis=1))
    max_force = np.max(force_magnitudes)
    
    # Normalize step size by maximum force
    if max_force > 0:
        normalized_forces = forces / max_force
    else:
        normalized_forces = forces
    
    # Update positions
    new_positions = positions + step_size * normalized_forces
    
    # Apply periodic boundary conditions
    new_positions = new_positions % box_size
    
    return new_positions, max_force

def energy_minimization(init_positions, box_size, max_steps=100, tolerance=0.1, init_step_size=0.01):
    """Perform energy minimization using steepest descent."""
    positions = init_positions.copy()
    step_size = init_step_size
    
    # Store energy trajectory
    energy_trajectory = []
    position_trajectory = []
    max_force_trajectory = []
    
    # Initial energy
    energy = calculate_system_energy(positions, box_size)
    energy_trajectory.append(energy)
    position_trajectory.append(positions.copy())
    
    for step in range(max_steps):
        # Calculate forces
        forces = calculate_forces(positions, box_size)
        
        # Store maximum force
        force_magnitudes = np.sqrt(np.sum(forces**2, axis=1))
        max_force = np.max(force_magnitudes)
        max_force_trajectory.append(max_force)
        
        # Check convergence
        if max_force < tolerance:
            print(f"Converged after {step} steps with max force {max_force:.6f}")
            break
            
        # Perform steepest descent step
        positions, _ = steepest_descent_step(positions, forces, step_size, box_size)
        
        # Calculate new energy
        energy = calculate_system_energy(positions, box_size)
        energy_trajectory.append(energy)
        position_trajectory.append(positions.copy())
        
        # Adaptive step size (simplified)
        if step > 0 and energy < energy_trajectory[-2]:
            step_size *= 1.2  # Increase step size if energy decreased
        else:
            step_size *= 0.5  # Decrease step size if energy increased
            
        # Limit step size
        step_size = min(step_size, 0.1)
    
    return energy_trajectory, position_trajectory, max_force_trajectory

def plot_energy_minimization(energy_trajectory, max_force_trajectory):
    """Create plots showing energy and force during minimization."""
    output_dir = create_output_dir()
    
    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Plot potential energy vs step
    steps = np.arange(len(energy_trajectory))
    ax1.plot(steps, energy_trajectory, 'b-', linewidth=2)
    ax1.set_ylabel('Potential Energy')
    ax1.set_title('Energy Minimization in Molecular Dynamics')
    ax1.grid(True, alpha=0.3)
    
    # Plot max force vs step
    ax2.plot(steps[:-1], max_force_trajectory, 'r-', linewidth=2)  # One fewer force points
    ax2.set_ylabel('Maximum Force')
    ax2.set_xlabel('Minimization Step')
    ax2.set_yscale('log')  # Log scale for forces
    ax2.grid(True, alpha=0.3)
    
    # Highlight regions of rapid change
    for i in range(1, len(energy_trajectory)):
        if i < len(energy_trajectory) - 1:
            # Check for significant drops in energy
            if (energy_trajectory[i] - energy_trajectory[i+1]) / abs(energy_trajectory[0]) > 0.05:
                ax1.axvspan(i, i+1, alpha=0.2, color='green')
                ax2.axvspan(i, i+1, alpha=0.2, color='green')
    
    plt.tight_layout()
    
    # Save figure
    plt.savefig(os.path.join(output_dir, "energy_minimization.png"), dpi=300)
    plt.savefig(os.path.join(output_dir, "energy_minimization.pdf"))
    
    # Create a combined plot with annotations
    plt.figure(figsize=(12, 10))
    gs = GridSpec(3, 3, figure=plt.gcf())
    
    # Energy plot
    ax_energy = plt.subplot(gs[0, :])
    ax_energy.plot(steps, energy_trajectory, 'b-', linewidth=2)
    ax_energy.set_ylabel('Potential Energy')
    ax_energy.set_title('Energy Minimization Process in Molecular Dynamics')
    ax_energy.grid(True, alpha=0.3)
    
    # Force plot
    ax_force = plt.subplot(gs[1, :])
    ax_force.plot(steps[:-1], max_force_trajectory, 'r-', linewidth=2)
    ax_force.set_ylabel('Maximum Force')
    ax_force.set_yscale('log')
    ax_force.grid(True, alpha=0.3)
    
    # Annotations explaining the minimization process
    ax_text = plt.subplot(gs[2, :])
    ax_text.axis('off')
    explanation = """
    Energy Minimization in Molecular Dynamics:
    
    1. Initial configuration: System starts in a high-energy, non-equilibrium state 
       with overlapping particles or unfavorable interactions.
       
    2. Steepest Descent Method: Forces on each atom are calculated from the gradient
       of the potential energy function. Atoms are moved in the direction of the force.
       
    3. Process continues iteratively until:
       • Maximum force falls below tolerance (convergence criterion)
       • Maximum number of steps is reached
       
    4. Step size adaptation: 
       • Increased when energy decreases (faster convergence)
       • Decreased when energy increases (finer adjustments)
       
    Applications:
    • Removing steric clashes in initial structures
    • Relaxing systems before dynamics
    • Preparing systems for equilibration and production runs
    """
    ax_text.text(0.02, 0.98, explanation, fontsize=11, va='top', ha='left', 
                transform=ax_text.transAxes, bbox=dict(facecolor='whitesmoke', alpha=0.7))
    
    plt.tight_layout()
    
    # Save figure
    plt.savefig(os.path.join(output_dir, "energy_minimization_explained.png"), dpi=300)
    plt.savefig(os.path.join(output_dir, "energy_minimization_explained.pdf"))
    
    return fig

def visualize_3d_system(position_trajectory, box_size, energy_trajectory):
    """Create an animation of the minimization process in 3D."""
    output_dir = create_output_dir()
    
    # Create figure for 3D visualization
    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[3, 1])
    
    # 3D plot of particles
    ax_3d = fig.add_subplot(gs[0, :], projection='3d')
    
    # Energy plot
    ax_energy = fig.add_subplot(gs[1, :])
    
    # Normalized energies for visualization (initially at 100%)
    normalized_energy = np.array(energy_trajectory) / energy_trajectory[0] * 100
    steps = np.arange(len(energy_trajectory))
    
    # Empty line for energy plot
    line_energy, = ax_energy.plot([], [], 'b-', linewidth=2)
    marker_energy, = ax_energy.plot([], [], 'ro', markersize=8)
    
    # Set up energy plot
    ax_energy.set_xlim(0, len(energy_trajectory))
    ax_energy.set_ylim(min(normalized_energy) - 5, 105)
    ax_energy.set_xlabel('Minimization Step')
    ax_energy.set_ylabel('Relative Energy (%)')
    ax_energy.grid(True, alpha=0.3)
    
    # Create a colormap for particles (blue for low energy, red for high energy)
    cmap = LinearSegmentedColormap.from_list('energy_cmap', ['blue', 'green', 'red'])
    
    # Animation initialization function
    def init():
        ax_3d.clear()
        ax_3d.set_xlim(0, box_size)
        ax_3d.set_ylim(0, box_size)
        ax_3d.set_zlim(0, box_size)
        ax_3d.set_xlabel('X')
        ax_3d.set_ylabel('Y')
        ax_3d.set_zlabel('Z')
        ax_3d.set_title('Energy Minimization Visualization')
        
        # Draw box boundaries
        draw_box(ax_3d, box_size)
        
        line_energy.set_data([], [])
        marker_energy.set_data([], [])
        
        return []
    
    # Animation update function
    def update(frame):
        ax_3d.clear()
        ax_3d.set_xlim(0, box_size)
        ax_3d.set_ylim(0, box_size)
        ax_3d.set_zlim(0, box_size)
        ax_3d.set_xlabel('X')
        ax_3d.set_ylabel('Y')
        ax_3d.set_zlabel('Z')
        
        # Draw box
        draw_box(ax_3d, box_size)
        
        # Get positions for current frame
        positions = position_trajectory[frame]
        
        # Color based on energy (redder means higher energy)
        energy_percentage = normalized_energy[frame]
        
        # Plot particles with color indicating energy state
        color_val = energy_percentage / 100  # Normalize to 0-1 for colormap
        ax_3d.scatter(
            positions[:, 0], positions[:, 1], positions[:, 2],
            c=[cmap(color_val)], s=50, alpha=0.8
        )
        
        # Update energy plot
        line_energy.set_data(steps[:frame+1], normalized_energy[:frame+1])
        marker_energy.set_data([frame], [normalized_energy[frame]])
        
        ax_3d.set_title(f'Energy Minimization - Step {frame}, Energy: {energy_percentage:.1f}%')
        
        return []
    
    # Draw box function
    def draw_box(ax, size):
        # Create a list of 8 vertices for a cube
        vertices = np.array([
            [0, 0, 0], [size, 0, 0], [size, size, 0], [0, size, 0],
            [0, 0, size], [size, 0, size], [size, size, size], [0, size, size]
        ])
        
        # Define the edges connecting the vertices
        edges = [
            (0, 1), (1, 2), (2, 3), (3, 0),  # Bottom face
            (4, 5), (5, 6), (6, 7), (7, 4),  # Top face
            (0, 4), (1, 5), (2, 6), (3, 7)   # Connecting edges
        ]
        
        # Plot the edges
        for edge in edges:
            ax.plot3D(
                [vertices[edge[0], 0], vertices[edge[1], 0]],
                [vertices[edge[0], 1], vertices[edge[1], 1]],
                [vertices[edge[0], 2], vertices[edge[1], 2]],
                'gray', alpha=0.5, linestyle='--'
            )
    
    # Create animation
    anim = animation.FuncAnimation(
        fig, update, frames=len(position_trajectory),
        init_func=init, blit=False, interval=200
    )
    
    # Save animation
    anim.save(os.path.join(output_dir, 'energy_minimization_3d.gif'), 
              writer='pillow', fps=5, dpi=100)
    
    return fig

def plot_energy_landscape():
    """Create a 2D visualization of energy landscape during minimization."""
    output_dir = create_output_dir()
    
    # Create a 2D grid for visualization
    x = np.linspace(-2, 2, 100)
    y = np.linspace(-2, 2, 100)
    X, Y = np.meshgrid(x, y)
    
    # Simple 2D energy landscape function (combination of Lennard-Jones-like potentials)
    def energy_function(x, y):
        # Center positions for some "atoms" in our 2D landscape
        centers = [
            (0.5, 0.5), (-0.5, -0.5), (0.8, -0.8), (-0.8, 0.8),
            (0, 0), (1.2, 0), (0, 1.2), (-1.2, 0), (0, -1.2)
        ]
        
        energy = 0
        for cx, cy in centers:
            r = np.sqrt((x - cx)**2 + (y - cy)**2)
            # Avoid division by zero
            r = np.maximum(r, 0.01)
            energy += 0.5 * ((1/r)**12 - (1/r)**6)
        
        return energy
    
    # Calculate energy on grid
    Z = energy_function(X, Y)
    
    # Clip extreme values for better visualization
    Z = np.clip(Z, -5, 5)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot energy landscape as contour
    contour = ax.contourf(X, Y, Z, 50, cmap='viridis')
    ax.contour(X, Y, Z, 20, colors='k', alpha=0.3, linewidths=0.5)
    plt.colorbar(contour, ax=ax, label='Potential Energy')
    
    # Simulated minimization path
    # Start positions (high energy regions)
    start_points = [
        (1.5, 1.5), (-1.5, -1.5), (1.5, -1.5), (-1.5, 1.5),
        (0.9, 0.1), (-0.9, -0.1), (0.1, 0.9), (-0.1, -0.9)
    ]
    
    # Simplified minimization path function
    def simple_minimization_path(start_x, start_y, steps=15):
        path_x = [start_x]
        path_y = [start_y]
        
        x, y = start_x, start_y
        
        for _ in range(steps):
            # Calculate gradient (forces)
            eps = 0.01
            grad_x = (energy_function(x+eps, y) - energy_function(x-eps, y)) / (2*eps)
            grad_y = (energy_function(x, y+eps) - energy_function(x, y-eps)) / (2*eps)
            
            # Normalize gradient
            grad_mag = np.sqrt(grad_x**2 + grad_y**2)
            if grad_mag > 1e-10:
                grad_x /= grad_mag
                grad_y /= grad_mag
            
            # Step in the negative gradient direction
            step_size = 0.1
            x -= step_size * grad_x
            y -= step_size * grad_y
            
            # Add to path
            path_x.append(x)
            path_y.append(y)
            
            # Break if gradient is very small (converged)
            if grad_mag < 0.01:
                break
        
        return path_x, path_y
    
    # Plot multiple minimization paths
    for i, (start_x, start_y) in enumerate(start_points):
        path_x, path_y = simple_minimization_path(start_x, start_y)
        
        # Use different colors for each path
        color = plt.cm.tab10(i % 10)
        
        # Plot path
        ax.plot(path_x, path_y, '-', color=color, linewidth=2, alpha=0.7)
        
        # Mark start and end
        ax.plot(start_x, start_y, 'o', color=color, markersize=8, alpha=0.8)
        ax.plot(path_x[-1], path_y[-1], '*', color=color, markersize=10)
    
    # Labels and title
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    ax.set_title('Energy Landscape and Minimization Paths')
    
    # Add annotations explaining minimization
    textbox = """
    Energy Minimization Process:
    
    • Circles: Starting structures (high energy)
    • Stars: Minimized structures (energy minima)
    • Lines: Minimization pathways
    
    Each minimization follows the direction of
    steepest descent on the energy landscape.
    """
    
    ax.text(1.4, -1.4, textbox, fontsize=10, va='bottom', ha='right',
            bbox=dict(facecolor='white', alpha=0.7))
    
    plt.tight_layout()
    
    # Save figure
    plt.savefig(os.path.join(output_dir, "energy_landscape.png"), dpi=300)
    plt.savefig(os.path.join(output_dir, "energy_landscape.pdf"))
    
    return fig

def main():
    """Run all energy minimization visualizations."""
    print("Creating energy minimization visualizations...")
    output_dir = create_output_dir()
    
    # Settings for the simulation
    box_size = 10.0
    num_particles = 50
    max_steps = 50
    
    # Generate initial configuration
    print("Generating initial particle configuration...")
    initial_positions = generate_initial_configuration(num_particles, box_size)
    
    # Run energy minimization
    print("Running energy minimization simulation...")
    energy_trajectory, position_trajectory, max_force_trajectory = energy_minimization(
        initial_positions, box_size, max_steps
    )
    
    # Create static plots
    print("Creating energy minimization plots...")
    plot_energy_minimization(energy_trajectory, max_force_trajectory)
    
    # Create 2D energy landscape visualization
    print("Creating energy landscape visualization...")
    plot_energy_landscape()
    
    try:
        # Create 3D visualization/animation
        print("Creating 3D animation of minimization process...")
        visualize_3d_system(position_trajectory, box_size, energy_trajectory)
    except Exception as e:
        print(f"Warning: Could not create 3D animation: {e}")
        print("Static visualizations were still created successfully.")
    
    print(f"All visualizations saved to {output_dir}/")

if __name__ == "__main__":
    main() 