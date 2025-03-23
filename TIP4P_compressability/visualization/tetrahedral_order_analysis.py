#!/usr/bin/env python3
"""
Tetrahedral Order Parameter Analysis

This script analyzes the tetrahedral order parameter (q₄) distributions in different phases
of water (ice, liquid, gas) based on simulation data.

Usage:
    python tetrahedral_order_analysis.py
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.ndimage import gaussian_filter
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde

# Check if MDAnalysis is available
try:
    import MDAnalysis as mda
    HAS_MDA = True
except ImportError:
    print("Warning: MDAnalysis not available, some functionality will be limited")
    HAS_MDA = False

def calculate_tetrahedral_order(positions, k=4):
    """
    Calculates tetrahedral parameter for each oxygen atom.
    
    Parameters:
    -----------
    positions : numpy.ndarray
        Array of atom positions
    k : int
        Number of nearest neighbors to consider (typically 4 for water)
        
    Returns:
    --------
    numpy.ndarray
        Array of tetrahedral order parameters
    """
    from scipy.spatial import cKDTree
    
    n_atoms = len(positions)
    q_values = np.zeros(n_atoms)
    
    kdtree = cKDTree(positions)
    
    for i in range(n_atoms):
        # k+1 nearest neighbors (first is the atom itself)
        distances, indices = kdtree.query(positions[i], k=k+1)
        neighbors = indices[1:k+1]
        
        vectors = positions[neighbors] - positions[i]
        norms = np.linalg.norm(vectors, axis=1)
        vectors = vectors / norms[:, np.newaxis]
        
        q = 0.0
        for j in range(k):
            for l in range(j+1, k):
                cos_angle = np.dot(vectors[j], vectors[l])
                q += (cos_angle + 1/3)**2
        q = 1 - 3/8 * q
        q_values[i] = q
    
    return q_values

def analyze_trajectory(tpr_file='../md.tpr', xtc_file='../md.xtc', 
                       temperature_file='../analysis/temperature_stats.dat',
                       max_frames=10):
    """
    Analyze the trajectory and calculate tetrahedral order parameters
    for water molecules at different temperatures
    
    Parameters:
    -----------
    tpr_file : str
        Path to TPR file
    xtc_file : str
        Path to trajectory file
    temperature_file : str
        File with temperature data
    max_frames : int
        Maximum number of frames to analyze
        
    Returns:
    --------
    dict
        Dictionary with analysis results
    """
    if not HAS_MDA:
        print("MDAnalysis is required for trajectory analysis")
        return generate_synthetic_data()
    
    print(f"Loading trajectory from {tpr_file} and {xtc_file}...")
    
    try:
        # Load the universe
        u = mda.Universe(tpr_file, xtc_file)
        
        # Select oxygen atoms (for TIP4P: OW)
        oxygen_atoms = u.select_atoms("name OW")
        print(f"Found {len(oxygen_atoms)} oxygen atoms")
        
        # Read temperature data if available
        temperature = None
        if os.path.exists(temperature_file):
            with open(temperature_file, 'r') as f:
                for line in f:
                    if "Mean temperature:" in line:
                        temperature = float(line.split(":")[1].split()[0])
                        break
        
        if temperature is None:
            temperature = 273  # Default to water freezing point
        
        print(f"Using temperature: {temperature} K")
        
        # Initialize results
        all_q_values = []
        frame_temps = []
        
        # Process frames
        for i, ts in enumerate(u.trajectory[:max_frames]):
            if i % 2 == 0:
                print(f"Processing frame {i+1}/{min(max_frames, len(u.trajectory))}")
            
            # Extract oxygen positions
            positions = oxygen_atoms.positions
            
            # Calculate tetrahedral order
            q_values = calculate_tetrahedral_order(positions)
            
            # Store results
            all_q_values.append(q_values)
            frame_temps.append(temperature)
        
        # Combine all values
        combined_q = np.concatenate(all_q_values)
        
        # Create synthetic phase data for demonstration
        # In a real analysis, this would be based on actual phase detection
        ice_mask = combined_q >= 0.8
        liquid_mask = (combined_q >= 0.5) & (combined_q < 0.8)
        gas_mask = combined_q < 0.5
        
        phase_distributions = {
            'ice': combined_q[ice_mask] if np.any(ice_mask) else np.array([]),
            'liquid': combined_q[liquid_mask] if np.any(liquid_mask) else np.array([]),
            'gas': combined_q[gas_mask] if np.any(gas_mask) else np.array([])
        }
        
        return {
            'all_q': combined_q,
            'temperature': temperature,
            'phase_distributions': phase_distributions
        }
        
    except Exception as e:
        print(f"Error analyzing trajectory: {e}")
        return generate_synthetic_data()

def generate_synthetic_data():
    """Generate synthetic data for demonstration"""
    print("Generating synthetic tetrahedral order data for demonstration...")
    
    # Generate synthetic q values for different phases
    np.random.seed(42)  # For reproducibility
    
    n_samples = 5000
    
    # Ice phase: high q values (0.8-1.0)
    ice_q = 0.9 + 0.1 * np.random.randn(n_samples)
    ice_q = np.clip(ice_q, 0.8, 1.0)
    
    # Liquid phase: medium q values (0.5-0.8)
    liquid_q = 0.65 + 0.15 * np.random.randn(n_samples)
    liquid_q = np.clip(liquid_q, 0.5, 0.8)
    
    # Gas phase: low q values (0.0-0.4)
    gas_q = 0.2 + 0.2 * np.random.randn(n_samples)
    gas_q = np.clip(gas_q, 0.0, 0.5)
    
    # Combine all for distribution
    all_q = np.concatenate([ice_q, liquid_q, gas_q])
    
    phase_distributions = {
        'ice': ice_q,
        'liquid': liquid_q,
        'gas': gas_q
    }
    
    return {
        'all_q': all_q,
        'temperature': 273,  # Simulated average temperature
        'phase_distributions': phase_distributions
    }

def plot_tetrahedral_distributions(data, output_file='../visualization/tetrahedral_order_distributions.png'):
    """
    Plot tetrahedral order parameter distributions for different phases
    
    Parameters:
    -----------
    data : dict
        Dictionary with analysis results
    output_file : str
        Path to save the output figure
    """
    print(f"Plotting tetrahedral order distributions to {output_file}...")
    
    # Create figure with GridSpec for better layout control
    fig = plt.figure(figsize=(10, 12))
    gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 1])
    
    # Extract data
    all_q = data['all_q']
    phase_distributions = data['phase_distributions']
    
    # Plot 1: Overall distribution
    ax1 = plt.subplot(gs[0, :])
    
    # Calculate KDE for smoother histogram
    kde = gaussian_kde(all_q)
    x_eval = np.linspace(0, 1, 1000)
    y_eval = kde(x_eval)
    
    # Plot KDE
    ax1.plot(x_eval, y_eval, 'k-', linewidth=2)
    
    # Mark phase regions
    ax1.axvspan(0.8, 1.0, alpha=0.2, color='cyan', label='Ice phase (q₄ > 0.8)')
    ax1.axvspan(0.5, 0.8, alpha=0.2, color='blue', label='Liquid phase (0.5 ≤ q₄ < 0.8)')
    ax1.axvspan(0.0, 0.5, alpha=0.2, color='red', label='Gas phase (q₄ < 0.5)')
    
    # Add statistical information
    mean_q = np.mean(all_q)
    median_q = np.median(all_q)
    std_q = np.std(all_q)
    
    stats_text = f"Mean q₄: {mean_q:.4f}\nMedian q₄: {median_q:.4f}\nStd Dev: {std_q:.4f}"
    ax1.text(0.75, 0.7 * max(y_eval), stats_text, 
            bbox=dict(facecolor='white', alpha=0.8, boxstyle='round'),
            transform=ax1.transData)
    
    ax1.set_xlabel('Tetrahedral Order Parameter (q₄)', fontsize=12)
    ax1.set_ylabel('Probability Density', fontsize=12)
    ax1.set_title(f'Overall Distribution of Tetrahedral Order Parameter (T ≈ {data["temperature"]:.1f} K)', 
                 fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper left')
    
    # Plot 2: Ice phase distribution
    ax2 = plt.subplot(gs[1, 0])
    if len(phase_distributions['ice']) > 0:
        ice_kde = gaussian_kde(phase_distributions['ice'])
        x_ice = np.linspace(0.8, 1.0, 500)
        y_ice = ice_kde(x_ice)
        ax2.plot(x_ice, y_ice, 'c-', linewidth=2)
        ax2.fill_between(x_ice, y_ice, alpha=0.4, color='cyan')
        
        # Add statistics
        ice_mean = np.mean(phase_distributions['ice'])
        ice_std = np.std(phase_distributions['ice'])
        ice_text = f"Mean: {ice_mean:.4f}\nStd Dev: {ice_std:.4f}\nSample Size: {len(phase_distributions['ice'])}"
        ax2.text(0.85, 0.7 * max(y_ice), ice_text, 
                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round'),
                transform=ax2.transData)
    else:
        ax2.text(0.9, 0.5, "No ice phase detected", 
                ha='center', va='center', 
                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round'))
    
    ax2.set_xlabel('Tetrahedral Order Parameter (q₄)', fontsize=12)
    ax2.set_ylabel('Probability Density', fontsize=12)
    ax2.set_title('Ice Phase (q₄ > 0.8)', fontsize=14)
    ax2.set_xlim(0.8, 1.0)
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Liquid phase distribution
    ax3 = plt.subplot(gs[1, 1])
    if len(phase_distributions['liquid']) > 0:
        liquid_kde = gaussian_kde(phase_distributions['liquid'])
        x_liquid = np.linspace(0.5, 0.8, 500)
        y_liquid = liquid_kde(x_liquid)
        ax3.plot(x_liquid, y_liquid, 'b-', linewidth=2)
        ax3.fill_between(x_liquid, y_liquid, alpha=0.4, color='blue')
        
        # Add statistics
        liquid_mean = np.mean(phase_distributions['liquid'])
        liquid_std = np.std(phase_distributions['liquid'])
        liquid_text = f"Mean: {liquid_mean:.4f}\nStd Dev: {liquid_std:.4f}\nSample Size: {len(phase_distributions['liquid'])}"
        ax3.text(0.7, 0.7 * max(y_liquid), liquid_text, 
                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round'),
                transform=ax3.transData)
    else:
        ax3.text(0.65, 0.5, "No liquid phase detected", 
                ha='center', va='center', 
                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round'))
    
    ax3.set_xlabel('Tetrahedral Order Parameter (q₄)', fontsize=12)
    ax3.set_ylabel('Probability Density', fontsize=12)
    ax3.set_title('Liquid Phase (0.5 ≤ q₄ < 0.8)', fontsize=14)
    ax3.set_xlim(0.5, 0.8)
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Gas phase distribution
    ax4 = plt.subplot(gs[2, 0])
    if len(phase_distributions['gas']) > 0:
        gas_kde = gaussian_kde(phase_distributions['gas'])
        x_gas = np.linspace(0.0, 0.5, 500)
        y_gas = gas_kde(x_gas)
        ax4.plot(x_gas, y_gas, 'r-', linewidth=2)
        ax4.fill_between(x_gas, y_gas, alpha=0.4, color='red')
        
        # Add statistics
        gas_mean = np.mean(phase_distributions['gas'])
        gas_std = np.std(phase_distributions['gas'])
        gas_text = f"Mean: {gas_mean:.4f}\nStd Dev: {gas_std:.4f}\nSample Size: {len(phase_distributions['gas'])}"
        ax4.text(0.35, 0.7 * max(y_gas), gas_text, 
                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round'),
                transform=ax4.transData)
    else:
        ax4.text(0.25, 0.5, "No gas phase detected", 
                ha='center', va='center', 
                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round'))
    
    ax4.set_xlabel('Tetrahedral Order Parameter (q₄)', fontsize=12)
    ax4.set_ylabel('Probability Density', fontsize=12)
    ax4.set_title('Gas Phase (q₄ < 0.5)', fontsize=14)
    ax4.set_xlim(0.0, 0.5)
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: Phase comparison (all on same plot)
    ax5 = plt.subplot(gs[2, 1])
    
    # Define colors and labels
    colors = ['cyan', 'blue', 'red']
    labels = ['Ice', 'Liquid', 'Gas']
    phase_keys = ['ice', 'liquid', 'gas']
    
    # Plot KDEs for each phase on same axes
    for i, phase in enumerate(phase_keys):
        if len(phase_distributions[phase]) > 0:
            phase_kde = gaussian_kde(phase_distributions[phase])
            x_phase = np.linspace(0, 1, 1000)
            y_phase = phase_kde(x_phase)
            ax5.plot(x_phase, y_phase, color=colors[i], linewidth=2, label=labels[i])
    
    ax5.set_xlabel('Tetrahedral Order Parameter (q₄)', fontsize=12)
    ax5.set_ylabel('Probability Density', fontsize=12)
    ax5.set_title('Phase Comparison', fontsize=14)
    ax5.set_xlim(0, 1)
    ax5.grid(True, alpha=0.3)
    ax5.legend()
    
    # Add overall title
    plt.suptitle('Tetrahedral Order Parameter Analysis for Water Phases', fontsize=16, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    
    # Save figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {output_file}")
    plt.close()

def main():
    """Main function to run the script"""
    print("Tetrahedral Order Parameter Analysis")
    
    # Create visualization directory if it doesn't exist
    os.makedirs("../visualization", exist_ok=True)
    
    # Try to analyze trajectory
    try:
        # Check if trajectory files exist
        if os.path.exists('../md.tpr') and os.path.exists('../md.xtc'):
            data = analyze_trajectory(max_frames=5)  # Limit to 5 frames for faster processing
        else:
            print("Trajectory files not found, using synthetic data")
            data = generate_synthetic_data()
    except Exception as e:
        print(f"Error in trajectory analysis: {e}")
        print("Using synthetic data for demonstration")
        data = generate_synthetic_data()
    
    # Plot the distributions
    plot_tetrahedral_distributions(data)
    
    print("Analysis complete.")

if __name__ == "__main__":
    main() 