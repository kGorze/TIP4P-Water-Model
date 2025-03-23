#!/usr/bin/env python3
"""
Umbrella Sampling Visualization

This script creates a visualization of the umbrella sampling technique, illustrating
how multiple biased simulations with overlapping windows allow a system to overcome
energy barriers and how the complete free energy profile is reconstructed.

Figure 6: Schematic representation of umbrella sampling. The unbiased free energy landscape 
(solid line) shows high barriers between states A and B. Multiple biased simulations (colored curves) 
with overlapping sampling windows allow the system to overcome these barriers. The complete free energy 
profile is then reconstructed by combining these biased simulations using methods like the 
Weighted Histogram Analysis Method (WHAM).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
import matplotlib.gridspec as gridspec
from matplotlib.collections import PatchCollection
import matplotlib.colors as mcolors

# Constants
k_B = 8.3145e-3  # Boltzmann constant in kJ/(mol·K)
T = 300  # Temperature in Kelvin

def create_double_well_potential(x):
    """
    Create a double well potential with a barrier
    
    Parameters:
    -----------
    x : array-like
        Reaction coordinate values
        
    Returns:
    --------
    array-like
        Energy values for the potential
    """
    # Double well potential: A*(x^4 - B*x^2)
    A = 1.0
    B = 2.0
    return A * (x**4 - B * x**2) + 5

def biasing_potential(x, center, k=10.0):
    """
    Harmonic biasing potential centered at a specific value
    
    Parameters:
    -----------
    x : array-like
        Reaction coordinate values
    center : float
        Center of the harmonic potential
    k : float
        Force constant for the harmonic potential
        
    Returns:
    --------
    array-like
        Energy values for the biasing potential
    """
    return 0.5 * k * (x - center)**2

def boltzmann_probability(energy, normalize=True):
    """
    Calculate the Boltzmann probability for a given energy
    
    Parameters:
    -----------
    energy : array-like
        Energy values
    normalize : bool
        Whether to normalize the probabilities
        
    Returns:
    --------
    array-like
        Probability values
    """
    # Calculate probabilities: p ∝ exp(-E/kT)
    p = np.exp(-energy / (k_B * T))
    
    if normalize:
        # Normalize probabilities
        return p / np.sum(p)
    else:
        return p

def get_free_energy(probability, shift_to_min=True):
    """
    Convert probability to free energy
    
    Parameters:
    -----------
    probability : array-like
        Probability values
    shift_to_min : bool
        Whether to shift the free energy so the minimum is at zero
        
    Returns:
    --------
    array-like
        Free energy values
    """
    # Small constant to avoid log(0)
    epsilon = 1e-10
    
    # Calculate free energy: F = -kT*ln(p)
    free_energy = -k_B * T * np.log(probability + epsilon)
    
    if shift_to_min:
        free_energy -= np.min(free_energy)
        
    return free_energy

def create_umbrella_sampling_visualization():
    """
    Create a visualization of umbrella sampling method
    
    Returns:
    --------
    matplotlib.figure.Figure
        The complete figure showing umbrella sampling
    """
    # Create figure
    fig = plt.figure(figsize=(14, 12))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1.5, 1, 1.5])
    
    # Create subplots
    ax1 = fig.add_subplot(gs[0])  # Potential energy with umbrella potentials
    ax2 = fig.add_subplot(gs[1])  # Histograms from biased simulations
    ax3 = fig.add_subplot(gs[2])  # Reconstructed free energy
    
    # Reaction coordinate
    x = np.linspace(-2.5, 2.5, 1000)
    
    # Define the unbiased potential (double well)
    unbiased_potential = create_double_well_potential(x)
    
    # Define centers for umbrella potentials
    umbrella_centers = np.linspace(-2.0, 2.0, 9)
    
    # Define colors for the umbrellas - use a color map for nice gradient
    colors = plt.cm.viridis(np.linspace(0, 1, len(umbrella_centers)))
    
    # PLOT 1: Potential Energy with Umbrella Potentials
    # -------------------------------------------------
    # Plot the unbiased potential
    ax1.plot(x, unbiased_potential, 'k-', linewidth=3, label='Unbiased Potential')
    
    # Plot umbrella potentials and corresponding biased potentials
    umbrella_biased_potentials = []
    for i, center in enumerate(umbrella_centers):
        # Calculate the biasing potential
        bias = biasing_potential(x, center, k=8.0)
        
        # Calculate the biased potential (sum of unbiased and bias)
        biased_potential = unbiased_potential + bias
        umbrella_biased_potentials.append(biased_potential)
        
        # Plot the biasing potential
        ax1.plot(x, bias, '--', color=colors[i], alpha=0.5, linewidth=1.5)
        
        # Plot the biased potential
        ax1.plot(x, biased_potential, '-', color=colors[i], alpha=0.7, linewidth=2, 
                label=f'Window {i+1}' if i % 2 == 0 else "")
        
        # Mark the center of the umbrella
        ax1.plot(center, bias[np.argmin(np.abs(x - center))], 'o', 
                markersize=6, color=colors[i])
        
        # Add shaded area to show the sampling region
        sampling_width = 0.6  # Width of effective sampling
        mask = (x >= center - sampling_width) & (x <= center + sampling_width)
        ax1.fill_between(x[mask], 
                        np.min(unbiased_potential) * np.ones_like(x[mask]), 
                        biased_potential[mask], 
                        color=colors[i], alpha=0.2)
    
    # Add state labels at minima
    minima_idx = [np.argmin(np.abs(x - (-1))), np.argmin(np.abs(x - 1))]
    state_labels = ['State A', 'State B']
    
    for i, idx in enumerate(minima_idx):
        ax1.text(x[idx], unbiased_potential[idx] - 1, state_labels[i], 
                ha='center', va='top', fontsize=12, 
                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
    
    # Add barrier label
    barrier_idx = np.argmin(np.abs(x - 0))
    ax1.text(x[barrier_idx], unbiased_potential[barrier_idx] + 1, 'Energy\nBarrier', 
            ha='center', va='bottom', fontsize=12, 
            bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
    
    # Add umbrella label
    ax1.text(1.8, 10, 'Biasing\nPotentials', 
            ha='center', va='center', fontsize=12, 
            bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
    
    # Add equation for umbrella potential
    eq_text = r"$V_{bias}(x) = \frac{1}{2}k(x-x_0)^2$"
    ax1.text(0.02, 0.95, eq_text, transform=ax1.transAxes, fontsize=12,
            bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
    
    # Add title and labels
    ax1.set_title('Umbrella Sampling: Biased Potentials', fontsize=14)
    ax1.set_xlabel('Reaction Coordinate', fontsize=12)
    ax1.set_ylabel('Energy (kJ/mol)', fontsize=12)
    ax1.legend(loc='upper right', fontsize=10, framealpha=0.8)
    ax1.set_ylim(np.min(unbiased_potential) - 1, 20)
    
    # PLOT 2: Histograms from Biased Simulations
    # ------------------------------------------
    # Create histograms for each umbrella window
    for i, center in enumerate(umbrella_centers):
        # Calculate Boltzmann probabilities for the biased potential
        biased_potential = umbrella_biased_potentials[i]
        biased_prob = boltzmann_probability(biased_potential, normalize=False)
        
        # Apply normalization to simulate histogram for this window
        norm_factor = 1.0 / np.sum(biased_prob)
        hist_data = biased_prob * norm_factor * 100  # Scale for visibility
        
        # Plot histogram
        ax2.fill_between(x, 0, hist_data, color=colors[i], alpha=0.5,
                        label=f'Window {i+1}' if i % 3 == 0 else "")
        
        # Add outline to the histogram for clarity
        ax2.plot(x, hist_data, '-', color=colors[i], alpha=0.8, linewidth=1.5)
    
    # Add overlap annotation
    overlap_positions = [(umbrella_centers[i] + umbrella_centers[i+1])/2 for i in range(len(umbrella_centers)-1)]
    for pos in overlap_positions[::2]:  # Add only every other annotation to avoid crowding
        idx = np.argmin(np.abs(x - pos))
        ax2.annotate('Overlap', xy=(pos, 1.5), xytext=(pos, 5),
                    arrowprops=dict(arrowstyle='->', color='black', lw=1.5),
                    ha='center', va='bottom', fontsize=10,
                    bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.8))
    
    # Add title and labels
    ax2.set_title('Histograms from Biased Simulations', fontsize=14)
    ax2.set_xlabel('Reaction Coordinate', fontsize=12)
    ax2.set_ylabel('Probability Density', fontsize=12)
    ax2.legend(loc='upper right', fontsize=10, framealpha=0.8)
    ax2.set_ylim(0, 15)
    
    # PLOT 3: Reconstructed Free Energy
    # ---------------------------------
    # Calculate the true free energy from the unbiased potential
    true_prob = boltzmann_probability(unbiased_potential)
    true_free_energy = get_free_energy(true_prob)
    
    # Create a simple "reweighted" free energy for illustration
    # In reality, this would use WHAM or another method
    reconstructed_free_energy = true_free_energy + np.random.normal(0, 0.2, size=len(x))
    
    # Plot the true free energy profile
    ax3.plot(x, true_free_energy, 'k-', linewidth=3, label='True Free Energy')
    
    # Plot the reconstructed free energy with uncertainty
    ax3.plot(x, reconstructed_free_energy, 'r--', linewidth=2.5, label='Reconstructed (WHAM)')
    
    # Add uncertainty band
    uncertainty = 0.3 + 0.2 * np.sin(x * 5)  # Varying uncertainty for illustration
    ax3.fill_between(x, reconstructed_free_energy - uncertainty, 
                    reconstructed_free_energy + uncertainty,
                    color='red', alpha=0.2, label='Uncertainty')
    
    # Highlight regions where each window contributes most
    for i, center in enumerate(umbrella_centers):
        # Find region where this window contributes significantly
        width = 0.5
        mask = (x >= center - width) & (x <= center + width)
        
        # Draw rectangle to show the contribution region
        rect = Rectangle((center - width, ax3.get_ylim()[0]), 
                        2 * width, ax3.get_ylim()[1] - ax3.get_ylim()[0],
                        alpha=0.1, color=colors[i])
        ax3.add_patch(rect)
    
    # Mark states on the free energy profile
    for i, idx in enumerate(minima_idx):
        ax3.plot(x[idx], true_free_energy[idx], 'o', markersize=8, color='darkblue')
        ax3.text(x[idx], true_free_energy[idx] - 0.5, state_labels[i], 
                ha='center', va='top', fontsize=12,
                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
    
    # Mark the barrier on the free energy profile
    barrier_idx = np.argmin(np.abs(x - 0))
    ax3.plot(x[barrier_idx], true_free_energy[barrier_idx], 'D', markersize=8, color='darkred')
    
    # Add explanation for WHAM
    wham_text = (
        "Weighted Histogram Analysis Method (WHAM):\n"
        "Combines data from all windows to\n"
        "reconstruct the unbiased free energy profile"
    )
    ax3.text(0.02, 0.95, wham_text, transform=ax3.transAxes, fontsize=10,
            bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
    
    # Add equation for free energy calculation
    eq_text = r"$\Delta G = -k_B T \ln \frac{P(x_B)}{P(x_A)}$"
    ax3.text(0.65, 0.15, eq_text, transform=ax3.transAxes, fontsize=12,
            bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
    
    # Add title and labels
    ax3.set_title('Reconstructed Free Energy Profile', fontsize=14)
    ax3.set_xlabel('Reaction Coordinate', fontsize=12)
    ax3.set_ylabel('Free Energy (kJ/mol)', fontsize=12)
    ax3.legend(loc='upper right', fontsize=10, framealpha=0.8)
    
    # Set y-axis range to focus on the energy wells
    ax3.set_ylim(0, 10)
    
    # Add overall title
    fig.suptitle('Umbrella Sampling Method', fontsize=16, y=0.98)
    
    # Add explanatory footer
    footer_text = (
        "Umbrella sampling uses multiple biased simulations with overlapping windows to overcome energy barriers. "
        "Each simulation is biased to sample a specific region of the reaction coordinate. "
        "The Weighted Histogram Analysis Method (WHAM) combines these biased simulations to reconstruct the complete free energy profile."
    )
    fig.text(0.5, 0.01, footer_text, ha='center', fontsize=10, 
            bbox=dict(facecolor='white', alpha=0.9, boxstyle='round,pad=0.3'))
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    return fig

def main():
    """Main function to run the script"""
    # Create the visualization
    fig = create_umbrella_sampling_visualization()
    
    # Save the figure
    output_file = 'umbrella_sampling.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {output_file}")
    
    # Close the figure
    plt.close()

if __name__ == "__main__":
    main() 