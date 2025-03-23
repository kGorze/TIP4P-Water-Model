#!/usr/bin/env python3
"""
Boltzmann Distribution Visualization

This script creates a visualization of the Boltzmann distribution to illustrate
how energy states are populated according to their energy differences.
Lower energy states are exponentially more probable than higher energy states.

Figure 5: Boltzmann distribution illustrating how states are populated according 
to their energy differences. Lower energy states (left) are exponentially more probable 
than higher energy states (right). This fundamental principle underlies our ability 
to convert probability distributions observed in simulations into meaningful free energy landscapes.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec

# Constants
k_B = 8.3145e-3  # Boltzmann constant in kJ/(mol·K)
T_room = 298.15  # Room temperature in K
T_water = 273.15  # Water freezing point in K
T_ice = 250.0   # Cold ice temperature in K

def boltzmann_probability(energy, temperature):
    """
    Calculate the Boltzmann probability for a given energy and temperature
    
    Parameters:
    -----------
    energy : float or array-like
        Energy level(s) in kJ/mol
    temperature : float
        Temperature in Kelvin
        
    Returns:
    --------
    float or array-like
        Probability according to Boltzmann distribution
    """
    # Calculate the exponential factor: exp(-E/kT)
    exp_factor = np.exp(-energy / (k_B * temperature))
    
    # For visualization purposes, we'll normalize to the maximum
    # In a real system, we would divide by the partition function
    return exp_factor / np.max(exp_factor)

def energy_levels_diagram():
    """
    Create a visualization of energy levels and their populations based on the Boltzmann distribution
    
    Returns:
    --------
    matplotlib.figure.Figure
        The complete figure with energy levels and Boltzmann distribution
    """
    # Create figure with better layout
    fig = plt.figure(figsize=(12, 10))
    
    # Create GridSpec for flexible layout
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1.5], width_ratios=[1.5, 1])
    
    # Energy levels and distribution plot
    ax1 = plt.subplot(gs[0, 0])  # Top left - Energy levels
    ax2 = plt.subplot(gs[1, 0])  # Bottom left - Boltzmann distribution
    ax3 = plt.subplot(gs[:, 1])  # Right column - Free energy landscape
    
    # Define a range of energies in kJ/mol
    energies = np.linspace(0, 15, 1000)
    
    # Calculate Boltzmann distributions at different temperatures
    p_room = boltzmann_probability(energies, T_room)
    p_water = boltzmann_probability(energies, T_water)
    p_ice = boltzmann_probability(energies, T_ice)
    
    # Define some discrete energy levels for illustration
    energy_levels = [0, 2, 5, 8, 12]
    level_labels = ["Ground\nState", "State 1", "State 2", "State 3", "State 4"]
    populations = boltzmann_probability(np.array(energy_levels), T_room)
    
    # Colors for different temperatures
    colors = {
        'room': '#FF6347',  # Tomato red for room temperature
        'water': '#1E90FF', # Dodger blue for water
        'ice': '#4682B4',   # Steel blue for ice
        'levels': '#2E8B57' # Sea green for energy levels
    }
    
    # PLOT 1: Energy Levels Diagram
    # -----------------------------
    # Draw horizontal lines for energy levels
    for i, (e, label, pop) in enumerate(zip(energy_levels, level_labels, populations)):
        # Line width proportional to population
        lw = max(1, 10 * pop)
        ax1.plot([0, 1], [e, e], linewidth=lw, color=colors['levels'])
        
        # Add label
        ax1.text(-0.1, e, label, verticalalignment='center', horizontalalignment='right', fontsize=12)
        
        # Add population percentage
        ax1.text(1.1, e, f"{pop*100:.1f}%", verticalalignment='center', horizontalalignment='left', fontsize=10)
        
        # Add small circles at ends to represent states
        ax1.plot(0, e, 'o', markersize=10, color=colors['levels'], alpha=0.7)
        ax1.plot(1, e, 'o', markersize=10, color=colors['levels'], alpha=0.7)
    
    # Add labels and title
    ax1.set_title("Energy Level Populations", fontsize=14)
    ax1.set_xlabel("States", fontsize=12)
    ax1.set_ylabel("Energy (kJ/mol)", fontsize=12)
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(["", ""])
    ax1.set_xlim(-0.2, 1.2)
    ax1.set_ylim(-1, energy_levels[-1] + 2)
    
    # Add annotations for energy gaps
    for i in range(len(energy_levels)-1):
        gap = energy_levels[i+1] - energy_levels[i]
        mid_point = (energy_levels[i] + energy_levels[i+1]) / 2
        ax1.annotate(f"ΔE = {gap:.1f} kJ/mol", 
                    xy=(0.5, mid_point),
                    xytext=(0.5, mid_point),
                    textcoords='data', 
                    ha='center', va='center',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.3),
                    fontsize=10)
    
    # PLOT 2: Boltzmann Distribution
    # -----------------------------
    # Plot the Boltzmann distributions
    ax2.plot(energies, p_room, '-', color=colors['room'], linewidth=2.5, label=f'T = {T_room}K (Room Temp)')
    ax2.plot(energies, p_water, '-', color=colors['water'], linewidth=2.5, label=f'T = {T_water}K (Water Freezing)')
    ax2.plot(energies, p_ice, '-', color=colors['ice'], linewidth=2.5, label=f'T = {T_ice}K (Ice)')
    
    # Add markers for the specific energy levels
    for e, pop in zip(energy_levels, populations):
        ax2.plot(e, boltzmann_probability(e, T_room), 'o', markersize=8, color=colors['room'])
        # Add dashed vertical lines from energy levels to distribution
        ax2.plot([e, e], [0, boltzmann_probability(e, T_room)], 'k--', alpha=0.3)
    
    # Add labels and title
    ax2.set_title("Boltzmann Distribution", fontsize=14)
    ax2.set_xlabel("Energy (kJ/mol)", fontsize=12)
    ax2.set_ylabel("Relative Probability", fontsize=12)
    ax2.set_xlim(0, energies[-1])
    ax2.set_ylim(0, 1.05)
    ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    # Add key equations
    eqn_text = r"$p_i \propto e^{-\frac{E_i}{k_B T}}$"
    ax2.text(0.1, 0.92, eqn_text, transform=ax2.transAxes, fontsize=14, 
             bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))
    
    # Add explanatory text
    explanation = "Lower energy states have\nexponentially higher probability"
    ax2.text(0.7, 0.5, explanation, transform=ax2.transAxes, fontsize=10,
            bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.5'))
    
    # Add temperature effect explanation
    temp_explanation = "Higher temperatures\nflatten the distribution"
    ax2.text(0.7, 0.3, temp_explanation, transform=ax2.transAxes, fontsize=10,
            bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.5'))
    
    # PLOT 3: Free Energy Landscape
    # -----------------------------
    # Create a simple 2D potential energy surface
    x = np.linspace(-3, 3, 100)
    y = np.linspace(-3, 3, 100)
    X, Y = np.meshgrid(x, y)
    
    # Simple double-well potential with two minima
    Z = 5*(X**2 - 1)**2 + Y**2
    
    # Calculate Boltzmann probabilities
    P = np.exp(-Z / (k_B * T_room))
    # Convert to free energy: F = -kT ln(P)
    F = -k_B * T_room * np.log(P / np.max(P))
    
    # Plot free energy landscape
    contour = ax3.contourf(X, Y, F, 20, cmap='viridis')
    ax3.contour(X, Y, F, 10, colors='k', alpha=0.3)
    
    # Add a colorbar
    cbar = plt.colorbar(contour, ax=ax3)
    cbar.set_label('Free Energy (kJ/mol)', fontsize=10)
    
    # Mark minima and transition state
    min1_x, min1_y = -1, 0
    min2_x, min2_y = 1, 0
    trans_x, trans_y = 0, 0
    
    ax3.plot(min1_x, min1_y, 'o', color='white', markersize=10)
    ax3.plot(min2_x, min2_y, 'o', color='white', markersize=10)
    ax3.plot(trans_x, trans_y, 'D', color='red', markersize=8)
    
    # Add labels
    ax3.text(min1_x-0.5, min1_y-0.5, "State A\n(Low Energy)", color='white', 
             fontsize=10, ha='center', va='center',
             bbox=dict(facecolor='black', alpha=0.5, boxstyle='round,pad=0.3'))
    ax3.text(min2_x+0.5, min2_y-0.5, "State B\n(Higher Energy)", color='white', 
             fontsize=10, ha='center', va='center',
             bbox=dict(facecolor='black', alpha=0.5, boxstyle='round,pad=0.3'))
    ax3.text(trans_x, trans_y+0.8, "Transition\nState", color='black', 
             fontsize=10, ha='center', va='center',
             bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
    
    # Draw a reaction pathway with an arrow
    path = FancyArrowPatch((min1_x, min1_y), (min2_x, min2_y), 
                          connectionstyle="arc3,rad=0.3", 
                          color='white', linewidth=2, arrowstyle='-|>', 
                          mutation_scale=20)
    ax3.add_patch(path)
    
    # Add title and labels
    ax3.set_title("Free Energy Landscape", fontsize=14)
    ax3.set_xlabel("Reaction Coordinate X", fontsize=12)
    ax3.set_ylabel("Reaction Coordinate Y", fontsize=12)
    
    # Add explanation of connection to Boltzmann
    connection = (
        "Free Energy (F) relates to probability (p)\n"
        "via the Boltzmann distribution:\n"
        r"$F = -k_B T \ln(p)$"
    )
    ax3.text(0.03, 0.03, connection, transform=ax3.transAxes, fontsize=10,
            bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.5'))
    
    # Add overall title to the figure
    fig.suptitle("Boltzmann Distribution and Free Energy Landscapes", fontsize=16, y=0.98)
    
    # Add explanatory footer
    footer_text = (
        "The Boltzmann distribution describes how states are populated based on their energy differences. "
        "Lower energy states (left) are exponentially more probable than higher energy states (right). "
        "This principle allows us to convert probability distributions from simulations into free energy landscapes."
    )
    fig.text(0.5, 0.01, footer_text, ha='center', fontsize=10, 
             bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.5'))
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.02, 1, 0.95])
    
    return fig

def main():
    """Main function to run the script"""
    # Create the visualization
    fig = energy_levels_diagram()
    
    # Save the figure
    output_file = 'boltzmann_distribution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {output_file}")
    
    # Show the figure
    plt.close()

if __name__ == "__main__":
    main() 