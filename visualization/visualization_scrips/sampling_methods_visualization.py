#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Ellipse, Rectangle
from matplotlib.gridspec import GridSpec
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap
import os

def create_output_dir():
    """Create output directory if it doesn't exist."""
    output_dir = "visualization_outputs"
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def generate_2d_energy_landscape(x_range=(-3, 3), y_range=(-3, 3), resolution=100):
    """Generate a 2D free energy landscape with multiple minima."""
    x = np.linspace(x_range[0], x_range[1], resolution)
    y = np.linspace(y_range[0], y_range[1], resolution)
    X, Y = np.meshgrid(x, y)
    
    # Create a landscape with multiple minima (representing different water configurations)
    Z = (
        3 * (1 - X)**2 * np.exp(-X**2 - (Y + 1)**2) 
        - 10 * (X/5 - X**3 - Y**5) * np.exp(-X**2 - Y**2)
        - 1/3 * np.exp(-(X + 1)**2 - Y**2)
        + 5 * np.exp(-(X - 1.5)**2 - (Y - 1.5)**2)
        + 3 * np.exp(-(X + 1.5)**2 - (Y - 1.5)**2)
    )
    
    # Normalize to a reasonable range for visualization
    Z = (Z - Z.min()) / (Z.max() - Z.min()) * 10
    
    return X, Y, Z

def visualize_sampling_methods():
    """Create visualizations for different sampling methods in MD simulations."""
    output_dir = create_output_dir()
    
    # Create a 2D free energy landscape
    X, Y, Z = generate_2d_energy_landscape()
    
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(15, 12))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 1])
    
    # 1. Free Energy Landscape Overview
    ax1 = fig.add_subplot(gs[0, :])
    
    # Plot the energy landscape as a contour plot
    contour = ax1.contourf(X, Y, Z, 20, cmap='viridis_r')
    ax1.contour(X, Y, Z, 10, colors='white', alpha=0.3, linewidths=0.5)
    
    # Add colorbar
    cbar = plt.colorbar(contour, ax=ax1)
    cbar.set_label('Free Energy (kJ/mol)', fontweight='bold')
    
    # Mark the global minimum and some local minima
    minima_coords = [
        (0.5, 0.5),    # Global minimum
        (-1.5, -1.0),  # Local minimum 1
        (1.5, 1.5),    # Local minimum 2
        (-1.5, 1.5)    # Local minimum 3
    ]
    
    for i, (x, y) in enumerate(minima_coords):
        if i == 0:
            ax1.scatter(x, y, color='red', s=100, edgecolor='white', linewidth=1.5, 
                       label='Global Minimum')
        else:
            if i == 1:
                ax1.scatter(x, y, color='orange', s=80, edgecolor='white', linewidth=1.5,
                           label='Local Minima')
            else:
                ax1.scatter(x, y, color='orange', s=80, edgecolor='white', linewidth=1.5)
    
    # Set labels and title
    ax1.set_xlabel('Reaction Coordinate 1 (e.g., O-H distance)', fontweight='bold')
    ax1.set_ylabel('Reaction Coordinate 2 (e.g., H-O-H angle)', fontweight='bold')
    ax1.set_title('Free Energy Landscape of Water Configurations', fontweight='bold')
    
    # Add legend
    ax1.legend(loc='upper right')
    
    # Add annotations for different regions
    ax1.annotate('Ice-like\nStructures', xy=(-2, -2), xytext=(-2.5, -2.5),
                color='white', fontweight='bold', ha='center',
                arrowprops=dict(facecolor='white', shrink=0.05, width=1.5, alpha=0.7))
    
    ax1.annotate('Liquid Water\nConfigurations', xy=(0.5, 0.5), xytext=(1.5, 0),
                color='white', fontweight='bold', ha='center',
                arrowprops=dict(facecolor='white', shrink=0.05, width=1.5, alpha=0.7))
    
    ax1.annotate('Gas-like\nStates', xy=(2.5, 2.5), xytext=(2.5, 2),
                color='white', fontweight='bold', ha='center',
                arrowprops=dict(facecolor='white', shrink=0.05, width=1.5, alpha=0.7))
    
    # 2. Standard MD Sampling
    ax2 = fig.add_subplot(gs[1, 0])
    
    # Plot the energy landscape
    ax2.contourf(X, Y, Z, 20, cmap='viridis_r')
    
    # Simulate a standard MD trajectory (gets trapped in local minimum)
    np.random.seed(42)
    n_steps = 100
    
    # Start near a local minimum
    x_start, y_start = -1.5, -1.0
    
    # Generate trajectory with small random steps (easily trapped in local minimum)
    x_traj = np.zeros(n_steps)
    y_traj = np.zeros(n_steps)
    x_traj[0], y_traj[0] = x_start, y_start
    
    for i in range(1, n_steps):
        # Calculate gradient (force) at current position
        # We'll approximate it using the energy landscape
        ix, iy = int((x_traj[i-1] - X.min()) / (X.max() - X.min()) * (X.shape[0] - 1)), \
                 int((y_traj[i-1] - Y.min()) / (Y.max() - Y.min()) * (Y.shape[1] - 1))
        
        # Limit indices to valid range
        ix = max(0, min(ix, X.shape[0] - 2))
        iy = max(0, min(iy, Y.shape[1] - 2))
        
        # Approximate gradient
        grad_x = Z[iy, ix+1] - Z[iy, ix]
        grad_y = Z[iy+1, ix] - Z[iy, ix]
        
        # Move downhill with some random noise (thermal fluctuations)
        step_size = 0.05
        random_factor = 0.03
        
        x_traj[i] = x_traj[i-1] - step_size * grad_x + random_factor * np.random.randn()
        y_traj[i] = y_traj[i-1] - step_size * grad_y + random_factor * np.random.randn()
        
        # Keep within bounds
        x_traj[i] = max(X.min(), min(x_traj[i], X.max()))
        y_traj[i] = max(Y.min(), min(y_traj[i], Y.max()))
    
    # Plot the trajectory
    ax2.plot(x_traj, y_traj, 'w-', alpha=0.7, linewidth=1)
    ax2.plot(x_traj, y_traj, 'wo', markersize=2, alpha=0.5)
    
    # Mark start and end points
    ax2.plot(x_traj[0], y_traj[0], 'go', markersize=8, label='Start')
    ax2.plot(x_traj[-1], y_traj[-1], 'ro', markersize=8, label='End')
    
    # Set labels and title
    ax2.set_xlabel('Reaction Coordinate 1', fontweight='bold')
    ax2.set_ylabel('Reaction Coordinate 2', fontweight='bold')
    ax2.set_title('Standard MD Sampling (Limited Exploration)', fontweight='bold')
    
    # Add legend
    ax2.legend(loc='upper right')
    
    # Add annotation about limitations
    ax2.text(0.5, 0.05, 
            "• Gets trapped in local minima\n• Limited by energy barriers\n• Requires long simulation times",
            transform=ax2.transAxes, fontsize=10, ha='center', va='bottom',
            bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.5'))
    
    # 3. Enhanced Sampling Methods
    ax3 = fig.add_subplot(gs[1, 1])
    
    # Plot the energy landscape
    ax3.contourf(X, Y, Z, 20, cmap='viridis_r')
    
    # Simulate enhanced sampling trajectories
    np.random.seed(42)
    
    # Generate multiple short trajectories from different starting points
    n_trajectories = 5
    n_steps_per_traj = 30
    
    # Starting points distributed across the landscape
    start_points = [
        (-1.5, -1.0),  # Near local minimum 1
        (0.5, 0.5),    # Near global minimum
        (1.5, 1.5),    # Near local minimum 2
        (-1.5, 1.5),   # Near local minimum 3
        (0.0, -2.0)    # Another region
    ]
    
    # Different colors for different trajectories
    colors = ['cyan', 'magenta', 'yellow', 'lime', 'white']
    
    for t in range(n_trajectories):
        x_traj = np.zeros(n_steps_per_traj)
        y_traj = np.zeros(n_steps_per_traj)
        x_traj[0], y_traj[0] = start_points[t]
        
        for i in range(1, n_steps_per_traj):
            # Enhanced sampling: occasionally make large jumps
            if np.random.rand() < 0.1:  # 10% chance of a large jump
                x_traj[i] = x_traj[i-1] + np.random.randn() * 0.5
                y_traj[i] = y_traj[i-1] + np.random.randn() * 0.5
            else:
                # Regular MD step with bias to explore higher energy regions
                ix, iy = int((x_traj[i-1] - X.min()) / (X.max() - X.min()) * (X.shape[0] - 1)), \
                         int((y_traj[i-1] - Y.min()) / (Y.max() - Y.min()) * (Y.shape[1] - 1))
                
                # Limit indices to valid range
                ix = max(0, min(ix, X.shape[0] - 2))
                iy = max(0, min(iy, Y.shape[1] - 2))
                
                # Approximate gradient with bias
                grad_x = Z[iy, ix+1] - Z[iy, ix]
                grad_y = Z[iy+1, ix] - Z[iy, ix]
                
                # Biased move with random component
                bias_factor = 0.5  # Reduces the effect of the gradient
                step_size = 0.05
                random_factor = 0.1  # Larger random component
                
                x_traj[i] = x_traj[i-1] - bias_factor * step_size * grad_x + random_factor * np.random.randn()
                y_traj[i] = y_traj[i-1] - bias_factor * step_size * grad_y + random_factor * np.random.randn()
            
            # Keep within bounds
            x_traj[i] = max(X.min(), min(x_traj[i], X.max()))
            y_traj[i] = max(Y.min(), min(y_traj[i], Y.max()))
        
        # Plot the trajectory
        ax3.plot(x_traj, y_traj, '-', color=colors[t], alpha=0.7, linewidth=1)
        ax3.plot(x_traj, y_traj, 'o', color=colors[t], markersize=2, alpha=0.5)
        
        # Mark start points
        if t == 0:
            ax3.plot(x_traj[0], y_traj[0], 'go', markersize=8, label='Starting Points')
        else:
            ax3.plot(x_traj[0], y_traj[0], 'go', markersize=8)
    
    # Set labels and title
    ax3.set_xlabel('Reaction Coordinate 1', fontweight='bold')
    ax3.set_ylabel('Reaction Coordinate 2', fontweight='bold')
    ax3.set_title('Enhanced Sampling Methods', fontweight='bold')
    
    # Add legend
    ax3.legend(loc='upper right')
    
    # Add annotation about enhanced sampling methods
    methods_text = (
        "Enhanced Sampling Techniques:\n"
        "• Replica Exchange MD\n"
        "• Umbrella Sampling\n"
        "• Metadynamics\n"
        "• Simulated Annealing\n"
        "• Parallel Tempering"
    )
    
    ax3.text(0.5, 0.05, methods_text,
            transform=ax3.transAxes, fontsize=10, ha='center', va='bottom',
            bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.5'))
    
    # Add explanation text for the entire figure
    explanation = (
        "Sampling Methods in Water Model Simulations:\n\n"
        "• Standard MD: Limited by energy barriers, may get trapped in local minima, requires long simulation times\n"
        "• Enhanced Sampling: Uses various techniques to overcome energy barriers and explore the free energy landscape more efficiently\n"
        "• For uncorrelated samples: Multiple short simulations from different starting points can be more efficient than a single long one"
    )
    
    plt.figtext(0.5, 0.01, explanation, ha='center', fontsize=12, 
               bbox=dict(facecolor='lightyellow', alpha=0.5, boxstyle='round,pad=0.5'))
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.08, 1, 0.98])
    
    # Save figure
    plt.savefig(os.path.join(output_dir, 'sampling_methods.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'sampling_methods.pdf'), bbox_inches='tight')
    
    # Create an animation showing sampling over time
    create_sampling_animation(output_dir)
    
    print(f"Sampling methods visualization saved to {output_dir}")
    plt.close()

def create_sampling_animation(output_dir):
    """Create an animation comparing standard MD vs enhanced sampling."""
    # Create a simpler energy landscape for animation
    x = np.linspace(-3, 3, 100)
    y = np.linspace(-3, 3, 100)
    X, Y = np.meshgrid(x, y)
    
    Z = (
        3 * (1 - X)**2 * np.exp(-X**2 - (Y + 1)**2) 
        - 10 * (X/5 - X**3 - Y**5) * np.exp(-X**2 - Y**2)
        - 1/3 * np.exp(-(X + 1)**2 - Y**2)
    )
    
    # Normalize
    Z = (Z - Z.min()) / (Z.max() - Z.min()) * 10
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Plot the energy landscape on both axes
    contour1 = ax1.contourf(X, Y, Z, 20, cmap='viridis_r')
    contour2 = ax2.contourf(X, Y, Z, 20, cmap='viridis_r')
    
    # Set titles and labels
    ax1.set_title('Standard MD Sampling', fontweight='bold')
    ax2.set_title('Enhanced Sampling', fontweight='bold')
    
    ax1.set_xlabel('Reaction Coordinate 1', fontweight='bold')
    ax1.set_ylabel('Reaction Coordinate 2', fontweight='bold')
    ax2.set_xlabel('Reaction Coordinate 1', fontweight='bold')
    ax2.set_ylabel('Reaction Coordinate 2', fontweight='bold')
    
    # Initialize trajectories
    n_frames = 100
    
    # Standard MD (single trajectory)
    x_std = np.zeros(n_frames)
    y_std = np.zeros(n_frames)
    x_std[0], y_std[0] = -1.5, -1.0  # Start near a local minimum
    
    # Enhanced sampling (multiple trajectories)
    n_enhanced = 3
    x_enh = np.zeros((n_enhanced, n_frames))
    y_enh = np.zeros((n_enhanced, n_frames))
    
    # Starting points for enhanced sampling
    x_enh[0, 0], y_enh[0, 0] = -1.5, -1.0  # Same as standard MD
    x_enh[1, 0], y_enh[1, 0] = 0.5, 0.5    # Near global minimum
    x_enh[2, 0], y_enh[2, 0] = 1.5, 1.5    # Another region
    
    # Pre-calculate trajectories
    for i in range(1, n_frames):
        # Standard MD: small steps, easily trapped
        ix, iy = int((x_std[i-1] - x.min()) / (x.max() - x.min()) * (len(x) - 1)), \
                 int((y_std[i-1] - y.min()) / (y.max() - y.min()) * (len(y) - 1))
        
        # Limit indices
        ix = max(0, min(ix, len(x) - 2))
        iy = max(0, min(iy, len(y) - 2))
        
        # Gradient
        grad_x = Z[iy, ix+1] - Z[iy, ix]
        grad_y = Z[iy+1, ix] - Z[iy, ix]
        
        # Move downhill with small random component
        step_size = 0.05
        random_factor = 0.03
        
        x_std[i] = x_std[i-1] - step_size * grad_x + random_factor * np.random.randn()
        y_std[i] = y_std[i-1] - step_size * grad_y + random_factor * np.random.randn()
        
        # Keep within bounds
        x_std[i] = max(x.min(), min(x_std[i], x.max()))
        y_std[i] = max(y.min(), min(y_std[i], y.max()))
        
        # Enhanced sampling: multiple trajectories with occasional jumps
        for j in range(n_enhanced):
            if np.random.rand() < 0.1:  # 10% chance of a large jump
                x_enh[j, i] = x_enh[j, i-1] + np.random.randn() * 0.5
                y_enh[j, i] = y_enh[j, i-1] + np.random.randn() * 0.5
            else:
                # Regular step with bias
                ix, iy = int((x_enh[j, i-1] - x.min()) / (x.max() - x.min()) * (len(x) - 1)), \
                         int((y_enh[j, i-1] - y.min()) / (y.max() - y.min()) * (len(y) - 1))
                
                # Limit indices
                ix = max(0, min(ix, len(x) - 2))
                iy = max(0, min(iy, len(y) - 2))
                
                # Gradient with bias
                grad_x = Z[iy, ix+1] - Z[iy, ix]
                grad_y = Z[iy+1, ix] - Z[iy, ix]
                
                bias_factor = 0.5
                step_size = 0.05
                random_factor = 0.1
                
                x_enh[j, i] = x_enh[j, i-1] - bias_factor * step_size * grad_x + random_factor * np.random.randn()
                y_enh[j, i] = y_enh[j, i-1] - bias_factor * step_size * grad_y + random_factor * np.random.randn()
            
            # Keep within bounds
            x_enh[j, i] = max(x.min(), min(x_enh[j, i], x.max()))
            y_enh[j, i] = max(y.min(), min(y_enh[j, i], y.max()))
    
    # Create scatter plots for current positions
    std_scatter = ax1.scatter([], [], color='red', s=100, zorder=5)
    enh_scatters = [ax2.scatter([], [], color=c, s=80, zorder=5) 
                   for c in ['red', 'cyan', 'yellow']]
    
    # Create line plots for trajectories
    std_line, = ax1.plot([], [], 'w-', linewidth=1.5, alpha=0.7)
    enh_lines = [ax2.plot([], [], '-', linewidth=1.5, alpha=0.7, 
                         color=c)[0] for c in ['red', 'cyan', 'yellow']]
    
    # Add text for frame number/time
    time_text1 = ax1.text(0.02, 0.98, '', transform=ax1.transAxes, 
                         color='white', fontweight='bold')
    time_text2 = ax2.text(0.02, 0.98, '', transform=ax2.transAxes, 
                         color='white', fontweight='bold')
    
    # Add text for visited states
    visited_text1 = ax1.text(0.02, 0.90, '', transform=ax1.transAxes,
                            color='white', fontweight='bold')
    visited_text2 = ax2.text(0.02, 0.90, '', transform=ax2.transAxes,
                            color='white', fontweight='bold')
    
    # Function to update the animation
    def update(frame):
        # Update standard MD
        std_scatter.set_offsets(np.column_stack([x_std[frame], y_std[frame]]))
        std_line.set_data(x_std[:frame+1], y_std[:frame+1])
        
        # Update enhanced sampling
        for j, scatter in enumerate(enh_scatters):
            scatter.set_offsets(np.column_stack([x_enh[j, frame], y_enh[j, frame]]))
            enh_lines[j].set_data(x_enh[j, :frame+1], y_enh[j, :frame+1])
        
        # Update time text
        time_text1.set_text(f'Time: {frame} fs')
        time_text2.set_text(f'Time: {frame} fs')
        
        # Calculate visited states (discretized)
        grid_size = 0.5
        x_bins = np.arange(x.min(), x.max() + grid_size, grid_size)
        y_bins = np.arange(y.min(), y.max() + grid_size, grid_size)
        
        # Count unique bins visited
        std_states = set()
        for i in range(frame + 1):
            x_bin = int((x_std[i] - x.min()) / grid_size)
            y_bin = int((y_std[i] - y.min()) / grid_size)
            std_states.add((x_bin, y_bin))
        
        enh_states = set()
        for j in range(n_enhanced):
            for i in range(frame + 1):
                x_bin = int((x_enh[j, i] - x.min()) / grid_size)
                y_bin = int((y_enh[j, i] - y.min()) / grid_size)
                enh_states.add((x_bin, y_bin))
        
        # Update visited states text
        visited_text1.set_text(f'Visited States: {len(std_states)}')
        visited_text2.set_text(f'Visited States: {len(enh_states)}')
        
        return [std_scatter, std_line, time_text1, visited_text1] + enh_scatters + enh_lines + [time_text2, visited_text2]
    
    # Create animation
    ani = animation.FuncAnimation(fig, update, frames=n_frames, interval=100, blit=True)
    
    # Add title
    plt.suptitle('Sampling Methods Comparison', fontsize=16, fontweight='bold')
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    # Save animation
    ani.save(os.path.join(output_dir, 'sampling_comparison.gif'), writer='pillow', fps=10, dpi=100)
    plt.close()

if __name__ == "__main__":
    visualize_sampling_methods() 