#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec
import os

def create_output_dir():
    """Create output directory if it doesn't exist."""
    output_dir = "visualization_outputs"
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def harmonic_oscillator(t, freq=100, amplitude=1.0, phase=0.0):
    """Simple harmonic oscillator function."""
    return amplitude * np.cos(2 * np.pi * freq * t + phase)

def visualize_timesteps():
    """Create visualizations comparing different integration time steps."""
    output_dir = create_output_dir()
    
    # Create figure with multiple subplots using GridSpec for more control
    fig = plt.figure(figsize=(15, 12))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[1, 1, 1.2])
    
    # 1. Comparison of different time steps for bond vibration
    ax1 = fig.add_subplot(gs[0, :])
    
    # Time range (ps)
    t = np.linspace(0, 0.1, 10000)  # 100 fs total with high resolution
    
    # Bond vibration frequency (O-H stretch ~3500 cm^-1 = ~100 THz)
    freq = 100  # THz
    
    # Calculate bond vibration
    bond_vibration = harmonic_oscillator(t, freq=freq)
    
    # Define different time steps
    dt_small = 0.5e-3  # 0.5 fs
    dt_medium = 1.0e-3  # 1.0 fs
    dt_large = 2.0e-3  # 2.0 fs (with constraints)
    
    # Sample the vibration at different time steps
    t_small = np.arange(0, 0.1, dt_small)
    t_medium = np.arange(0, 0.1, dt_medium)
    t_large = np.arange(0, 0.1, dt_large)
    
    vibration_small = harmonic_oscillator(t_small, freq=freq)
    vibration_medium = harmonic_oscillator(t_medium, freq=freq)
    vibration_large = harmonic_oscillator(t_large, freq=freq)
    
    # Plot continuous vibration
    ax1.plot(t, bond_vibration, 'k-', linewidth=1, alpha=0.3, label='Actual Bond Vibration')
    
    # Plot sampled points
    ax1.plot(t_small, vibration_small, 'bo', markersize=4, label=f'0.5 fs timestep ({len(t_small)} steps)')
    ax1.plot(t_medium, vibration_medium, 'go', markersize=6, label=f'1.0 fs timestep ({len(t_medium)} steps)')
    ax1.plot(t_large, vibration_large, 'ro', markersize=8, label=f'2.0 fs timestep ({len(t_large)} steps)')
    
    # Add Nyquist frequency limit for 2 fs timestep
    nyquist_freq = 1 / (2 * dt_large)  # Nyquist frequency in THz
    ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    
    # Add annotation for Nyquist limit
    ax1.text(0.05, 0.9, f"Nyquist limit for 2 fs timestep: {nyquist_freq:.1f} THz\nO-H vibration: ~100 THz", 
            transform=ax1.transAxes, fontsize=10, va='top', ha='left',
            bbox=dict(facecolor='lightyellow', alpha=0.8, boxstyle='round,pad=0.5'))
    
    # Set labels and title
    ax1.set_xlabel('Time (ps)', fontweight='bold')
    ax1.set_ylabel('Bond Length Deviation', fontweight='bold')
    ax1.set_title('Sampling of Bond Vibration with Different Time Steps', fontweight='bold')
    
    # Add grid and legend
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend(loc='upper right')
    
    # Zoom in to show a smaller time window
    ax1.set_xlim(0, 0.02)  # Show only first 20 fs
    
    # 2. Timeline visualization of integration steps
    ax2 = fig.add_subplot(gs[1, :])
    
    # Create timeline for different time steps
    timeline_length = 0.05  # 50 fs
    
    # Create markers for different time steps
    t_small_timeline = np.arange(0, timeline_length, dt_small)
    t_medium_timeline = np.arange(0, timeline_length, dt_medium)
    t_large_timeline = np.arange(0, timeline_length, dt_large)
    
    # Plot timelines
    ax2.scatter(t_small_timeline, np.ones_like(t_small_timeline)*3, color='blue', s=30, marker='|')
    ax2.scatter(t_medium_timeline, np.ones_like(t_medium_timeline)*2, color='green', s=50, marker='|')
    ax2.scatter(t_large_timeline, np.ones_like(t_large_timeline)*1, color='red', s=70, marker='|')
    
    # Add labels for each timeline
    ax2.text(-0.005, 3, "0.5 fs timestep", ha='right', va='center', fontweight='bold', color='blue')
    ax2.text(-0.005, 2, "1.0 fs timestep", ha='right', va='center', fontweight='bold', color='green')
    ax2.text(-0.005, 1, "2.0 fs timestep\n(with constraints)", ha='right', va='center', fontweight='bold', color='red')
    
    # Add computational cost comparison
    steps_small = len(t_small_timeline)
    steps_medium = len(t_medium_timeline)
    steps_large = len(t_large_timeline)
    
    ax2.text(timeline_length + 0.002, 3, f"{steps_small} steps", ha='left', va='center', fontweight='bold', color='blue')
    ax2.text(timeline_length + 0.002, 2, f"{steps_medium} steps\n({steps_small/steps_medium:.1f}× faster)", ha='left', va='center', fontweight='bold', color='green')
    ax2.text(timeline_length + 0.002, 1, f"{steps_large} steps\n({steps_small/steps_large:.1f}× faster)", ha='left', va='center', fontweight='bold', color='red')
    
    # Set labels and title
    ax2.set_xlabel('Time (ps)', fontweight='bold')
    ax2.set_title('Integration Steps Comparison for 50 fs of Simulation', fontweight='bold')
    
    # Remove y-axis ticks and labels
    ax2.set_yticks([])
    ax2.set_yticklabels([])
    
    # Set x-axis limits
    ax2.set_xlim(-0.01, timeline_length + 0.02)
    ax2.set_ylim(0.5, 3.5)
    
    # Add grid
    ax2.grid(True, axis='x', linestyle='--', alpha=0.7)
    
    # 3. Constraint vs. No Constraint Comparison
    ax3 = fig.add_subplot(gs[2, 0])
    ax4 = fig.add_subplot(gs[2, 1])
    
    # Time for animation
    t_anim = np.linspace(0, 2*np.pi, 100)
    
    # No constraints - bond vibrates
    bond_length = 0.1  # nm
    vibration_amplitude = 0.02  # nm
    
    # Create water molecule coordinates
    o_pos = np.array([0, 0])
    angle = np.deg2rad(104.52)
    
    # Function to update water molecule with vibrating bonds
    def update_vibrating_water(i, line1, line2, h1_scatter, h2_scatter, time_text):
        # Vibrating bond lengths
        bond1_length = bond_length + vibration_amplitude * np.sin(5*t_anim[i])
        bond2_length = bond_length + vibration_amplitude * np.sin(5*t_anim[i] + np.pi/3)
        
        # Update hydrogen positions
        h1_pos = o_pos + bond1_length * np.array([np.sin(angle/2), np.cos(angle/2)])
        h2_pos = o_pos + bond2_length * np.array([-np.sin(angle/2), np.cos(angle/2)])
        
        # Update lines
        line1.set_data([o_pos[0], h1_pos[0]], [o_pos[1], h1_pos[1]])
        line2.set_data([o_pos[0], h2_pos[0]], [o_pos[1], h2_pos[1]])
        
        # Update scatter points
        h1_scatter.set_offsets([h1_pos[0], h1_pos[1]])
        h2_scatter.set_offsets([h2_pos[0], h2_pos[1]])
        
        # Update time text
        time_text.set_text(f'Time: {i*0.5:.1f} fs')
        
        return line1, line2, h1_scatter, h2_scatter, time_text
    
    # Function to update water molecule with constrained bonds
    def update_constrained_water(i, line1, line2, h1_scatter, h2_scatter, time_text):
        # Fixed bond lengths
        bond1_length = bond_length
        bond2_length = bond_length
        
        # Update hydrogen positions - only the angle changes slightly
        angle_variation = np.deg2rad(2) * np.sin(t_anim[i])
        current_angle = angle + angle_variation
        
        h1_pos = o_pos + bond1_length * np.array([np.sin(current_angle/2), np.cos(current_angle/2)])
        h2_pos = o_pos + bond2_length * np.array([-np.sin(current_angle/2), np.cos(current_angle/2)])
        
        # Update lines
        line1.set_data([o_pos[0], h1_pos[0]], [o_pos[1], h1_pos[1]])
        line2.set_data([o_pos[0], h2_pos[0]], [o_pos[1], h2_pos[1]])
        
        # Update scatter points
        h1_scatter.set_offsets([h1_pos[0], h1_pos[1]])
        h2_scatter.set_offsets([h2_pos[0], h2_pos[1]])
        
        # Update time text
        time_text.set_text(f'Time: {i*2.0:.1f} fs')
        
        return line1, line2, h1_scatter, h2_scatter, time_text
    
    # Setup for unconstrained water
    ax3.set_xlim(-0.2, 0.2)
    ax3.set_ylim(-0.05, 0.25)
    ax3.set_aspect('equal')
    ax3.set_title('Without Constraints (0.5 fs timestep)', fontweight='bold')
    ax3.set_xlabel('X (nm)', fontweight='bold')
    ax3.set_ylabel('Y (nm)', fontweight='bold')
    
    # Initial positions for unconstrained
    h1_pos_init = o_pos + bond_length * np.array([np.sin(angle/2), np.cos(angle/2)])
    h2_pos_init = o_pos + bond_length * np.array([-np.sin(angle/2), np.cos(angle/2)])
    
    # Create initial plot elements for unconstrained
    line1_unconstrained, = ax3.plot([o_pos[0], h1_pos_init[0]], [o_pos[1], h1_pos_init[1]], 'k-', linewidth=2)
    line2_unconstrained, = ax3.plot([o_pos[0], h2_pos_init[0]], [o_pos[1], h2_pos_init[1]], 'k-', linewidth=2)
    o_scatter = ax3.scatter(*o_pos, color='red', s=200, zorder=3)
    h1_scatter_unconstrained = ax3.scatter(*h1_pos_init, color='white', s=100, edgecolors='black', zorder=3)
    h2_scatter_unconstrained = ax3.scatter(*h2_pos_init, color='white', s=100, edgecolors='black', zorder=3)
    time_text_unconstrained = ax3.text(0.05, 0.95, '', transform=ax3.transAxes)
    
    # Add annotation for high frequency vibrations
    ax3.text(0.5, 0.1, "High-frequency\nbond vibrations", ha='center', va='center',
            transform=ax3.transAxes, fontsize=10,
            bbox=dict(facecolor='lightyellow', alpha=0.8, boxstyle='round,pad=0.5'))
    
    # Setup for constrained water
    ax4.set_xlim(-0.2, 0.2)
    ax4.set_ylim(-0.05, 0.25)
    ax4.set_aspect('equal')
    ax4.set_title('With Constraints (2.0 fs timestep)', fontweight='bold')
    ax4.set_xlabel('X (nm)', fontweight='bold')
    ax4.set_ylabel('Y (nm)', fontweight='bold')
    
    # Create initial plot elements for constrained
    line1_constrained, = ax4.plot([o_pos[0], h1_pos_init[0]], [o_pos[1], h1_pos_init[1]], 'k-', linewidth=2)
    line2_constrained, = ax4.plot([o_pos[0], h2_pos_init[0]], [o_pos[1], h2_pos_init[1]], 'k-', linewidth=2)
    ax4.scatter(*o_pos, color='red', s=200, zorder=3)
    h1_scatter_constrained = ax4.scatter(*h1_pos_init, color='white', s=100, edgecolors='black', zorder=3)
    h2_scatter_constrained = ax4.scatter(*h2_pos_init, color='white', s=100, edgecolors='black', zorder=3)
    time_text_constrained = ax4.text(0.05, 0.95, '', transform=ax4.transAxes)
    
    # Add annotation for constrained bonds
    ax4.text(0.5, 0.1, "Fixed bond lengths\nallow larger timesteps", ha='center', va='center',
            transform=ax4.transAxes, fontsize=10,
            bbox=dict(facecolor='lightyellow', alpha=0.8, boxstyle='round,pad=0.5'))
    
    # Add constraint indicators
    constraint_line1 = ax4.plot([o_pos[0], h1_pos_init[0]], [o_pos[1], h1_pos_init[1]], 'r--', linewidth=1, alpha=0.7)[0]
    constraint_line2 = ax4.plot([o_pos[0], h2_pos_init[0]], [o_pos[1], h2_pos_init[1]], 'r--', linewidth=1, alpha=0.7)[0]
    
    # Create animations
    ani1 = animation.FuncAnimation(fig, update_vibrating_water, frames=len(t_anim),
                                  fargs=(line1_unconstrained, line2_unconstrained, 
                                         h1_scatter_unconstrained, h2_scatter_unconstrained,
                                         time_text_unconstrained),
                                  interval=50, blit=True)
    
    ani2 = animation.FuncAnimation(fig, update_constrained_water, frames=len(t_anim),
                                  fargs=(line1_constrained, line2_constrained, 
                                         h1_scatter_constrained, h2_scatter_constrained,
                                         time_text_constrained),
                                  interval=50, blit=True)
    
    # Add explanation text
    explanation = (
        "Integration Time Step in MD Simulations:\n\n"
        "• Without constraints: Bond vibrations limit time step to ~0.5-1.0 fs\n"
        "• With constraints (SHAKE/LINCS): Fixed bond lengths eliminate high-frequency vibrations\n"
        "• Larger time step (2.0 fs) possible with constraints, making simulations 2-4× faster\n"
        "• Constraints preserve the essential physics while improving computational efficiency"
    )
    
    plt.figtext(0.5, 0.01, explanation, ha='center', fontsize=12, 
               bbox=dict(facecolor='lightyellow', alpha=0.5, boxstyle='round,pad=0.5'))
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.08, 1, 0.98])
    
    # Save animations as GIFs
    ani1.save(os.path.join(output_dir, 'unconstrained_water.gif'), writer='pillow', fps=15, dpi=100)
    ani2.save(os.path.join(output_dir, 'constrained_water.gif'), writer='pillow', fps=15, dpi=100)
    
    # Save static figure
    plt.savefig(os.path.join(output_dir, 'timestep_comparison.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'timestep_comparison.pdf'), bbox_inches='tight')
    
    print(f"Time step visualization saved to {output_dir}")
    plt.close()

if __name__ == "__main__":
    visualize_timesteps() 