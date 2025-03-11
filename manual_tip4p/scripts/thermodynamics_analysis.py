#!/usr/bin/env python3
"""
Thermodynamic Properties Analysis for TIP4P Water

This script analyzes thermodynamic properties from GROMACS energy files (.edr),
including temperature, pressure, energy components, etc.

Usage:
    python thermodynamics_analysis.py [edr_file]

Default:
    Uses md.edr in the current directory
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
import pandas as pd
from scipy import stats, signal
from utils import parse_edr_file, save_plot

def analyze_temperature(edr_file='md.edr'):
    """
    Analyze temperature data from GROMACS energy file
    
    Parameters:
    -----------
    edr_file : str
        Path to the .edr file
    """
    print("Analyzing temperature...")
    
    # Parse temperature data
    temp_data = parse_edr_file(edr_file, properties=['Temperature'])
    
    if 'Temperature' not in temp_data or temp_data['Temperature'] is None:
        print("Could not extract temperature data from EDR file")
        return None
    
    df = temp_data['Temperature']
    
    # Calculate statistics
    temp_mean = df['Temperature'].mean()
    temp_std = df['Temperature'].std()
    
    # Save statistics
    with open('../analysis/temperature_stats.dat', 'w') as f:
        f.write(f"# Temperature statistics from {edr_file}\n")
        f.write(f"# Mean temperature: {temp_mean:.4f} K\n")
        f.write(f"# Std deviation: {temp_std:.4f} K\n")
    
    # Plot temperature vs time
    plt.figure(figsize=(10, 6))
    plt.plot(df['Time (ps)'], df['Temperature'], 'b-')
    plt.axhline(temp_mean, color='r', linestyle='--', 
                label=f'Mean: {temp_mean:.2f} K')
    plt.axhline(temp_mean + temp_std, color='g', linestyle=':', 
                label=f'±1 Std Dev: {temp_std:.2f} K')
    plt.axhline(temp_mean - temp_std, color='g', linestyle=':')
    plt.xlabel('Time (ps)')
    plt.ylabel('Temperature (K)')
    plt.title('Temperature vs Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    save_plot(plt, 'temperature_plot.png')
    
    # Plot temperature histogram
    plt.figure(figsize=(10, 6))
    plt.hist(df['Temperature'], bins=30, alpha=0.7, color='blue', edgecolor='black')
    plt.axvline(temp_mean, color='r', linestyle='--', 
                label=f'Mean: {temp_mean:.2f} K')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Frequency')
    plt.title('Temperature Distribution')
    plt.legend()
    plt.grid(True, alpha=0.3)
    save_plot(plt, 'temperature_histogram.png')
    
    return df

def analyze_pressure(edr_file='md.edr'):
    """
    Analyze pressure data from GROMACS energy file
    
    Parameters:
    -----------
    edr_file : str
        Path to the .edr file
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing pressure data
    """
    print("Analyzing pressure...")
    
    # Parse pressure data
    pressure_data = parse_edr_file(edr_file, properties=['Pressure'])
    
    if 'Pressure' not in pressure_data or pressure_data['Pressure'] is None:
        print("Could not extract pressure data from EDR file")
        return None
    
    df = pressure_data['Pressure']
    
    # Calculate statistics
    pressure_mean = df['Pressure'].mean()
    pressure_std = df['Pressure'].std()
    
    # Save statistics
    with open('../analysis/pressure_stats.dat', 'w') as f:
        f.write(f"# Pressure statistics\n")
        f.write(f"# Mean pressure: {pressure_mean:.2f} bar\n")
        f.write(f"# Standard deviation: {pressure_std:.2f} bar\n")
        f.write(f"# Coefficient of variation: {pressure_std/abs(pressure_mean):.4f}\n")
        f.write(f"# Target pressure: 1.0 bar\n")
        f.write(f"# Deviation from target: {pressure_mean - 1.0:.2f} bar\n")
    
    # Save data
    df.to_csv('../analysis/pressure_data.dat', sep='\t', index=False, 
              float_format='%.6f', header=True)
    
    # Plot pressure
    plt.figure(figsize=(12, 8))
    
    # Plot the pressure data
    plt.plot(df['Time'], df['Pressure'], 'b-', alpha=0.7, label='Pressure')
    
    # Add horizontal line for the mean
    plt.axhline(y=pressure_mean, color='r', linestyle='--',
               label=f'Mean: {pressure_mean:.2f} bar')
    
    # Add horizontal lines for ±1 standard deviation
    plt.axhline(y=pressure_mean + pressure_std, color='g', linestyle=':',
               label=f'±1σ: {pressure_std:.2f} bar')
    plt.axhline(y=pressure_mean - pressure_std, color='g', linestyle=':', alpha=0.7)
    
    # Add target pressure line
    plt.axhline(y=1.0, color='k', linestyle='-', alpha=0.5,
               label='Target: 1.0 bar')
    
    # Add shaded region for ±1 standard deviation
    plt.fill_between(df['Time'], 
                    pressure_mean - pressure_std, 
                    pressure_mean + pressure_std, 
                    color='b', alpha=0.1)
    
    plt.xlabel('Time (ps)')
    plt.ylabel('Pressure (bar)')
    plt.title('Pressure vs. Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add text box with context about pressure fluctuations
    textstr = '\n'.join((
        f'Mean: {pressure_mean:.2f} bar',
        f'Std Dev: {pressure_std:.2f} bar',
        f'Target: 1.0 bar',
        f'Deviation: {pressure_mean - 1.0:.2f} bar'))
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0.02, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    # Add annotation about pressure fluctuations
    plt.annotate('Large pressure fluctuations are normal\nin molecular simulations with small systems', 
                 xy=(0.5, 0.05), xycoords='axes fraction',
                 bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="orange", alpha=0.8),
                 ha='center')
    
    plt.tight_layout()
    plt.savefig('../analysis/pressure_plot.png', dpi=300, bbox_inches='tight')
    
    # Create a more detailed analysis plot
    plt.figure(figsize=(12, 12))
    
    # Plot 1: Pressure vs. time
    plt.subplot(2, 1, 1)
    plt.plot(df['Time'], df['Pressure'], 'b-', alpha=0.7)
    plt.axhline(y=pressure_mean, color='r', linestyle='--',
               label=f'Mean: {pressure_mean:.2f} bar')
    plt.axhline(y=1.0, color='k', linestyle='-', alpha=0.5,
               label='Target: 1.0 bar')
    plt.fill_between(df['Time'], 
                    pressure_mean - pressure_std, 
                    pressure_mean + pressure_std, 
                    color='b', alpha=0.1,
                    label=f'±1σ: {pressure_std:.2f} bar')
    plt.xlabel('Time (ps)')
    plt.ylabel('Pressure (bar)')
    plt.title('Pressure vs. Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot 2: Pressure distribution
    plt.subplot(2, 1, 2)
    
    # Create histogram with appropriate bins
    n, bins, patches = plt.hist(df['Pressure'], bins=30, alpha=0.7, 
                               color='blue', edgecolor='black', density=True)
    
    # Add vertical lines for statistics
    plt.axvline(pressure_mean, color='r', linestyle='--',
               label=f'Mean: {pressure_mean:.2f} bar')
    plt.axvline(1.0, color='k', linestyle='-', alpha=0.5,
               label='Target: 1.0 bar')
    
    # Add a fitted normal distribution
    from scipy.stats import norm
    x = np.linspace(pressure_mean - 4*pressure_std, pressure_mean + 4*pressure_std, 1000)
    plt.plot(x, norm.pdf(x, pressure_mean, pressure_std), 'r-', 
             label='Normal distribution')
    
    plt.xlabel('Pressure (bar)')
    plt.ylabel('Probability density')
    plt.title('Pressure Distribution')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add text about pressure fluctuations in NPT simulations
    textstr = '\n'.join((
        'Pressure fluctuations in NPT simulations:',
        '• Fluctuations scale as 1/√N (N = number of particles)',
        '• Typical std dev: thousands of bar for small systems',
        '• Only the average should be close to the target pressure'))
    
    props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5)
    plt.text(0.02, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig('../analysis/pressure_analysis_plot.png', dpi=300, bbox_inches='tight')
    
    # Create a running average plot to show convergence
    plt.figure(figsize=(12, 6))
    
    # Calculate running average
    running_avg = np.cumsum(df['Pressure']) / np.arange(1, len(df['Pressure'])+1)
    
    # Plot running average
    plt.plot(df['Time'], running_avg, 'g-', 
             label='Running average')
    
    # Add target and final average lines
    plt.axhline(y=1.0, color='k', linestyle='-', alpha=0.5,
               label='Target: 1.0 bar')
    plt.axhline(y=pressure_mean, color='r', linestyle='--',
               label=f'Final average: {pressure_mean:.2f} bar')
    
    plt.xlabel('Time (ps)')
    plt.ylabel('Pressure (bar)')
    plt.title('Pressure Running Average vs. Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add annotation about convergence
    if abs(pressure_mean - 1.0) < 5.0:
        convergence = "Good convergence to target pressure"
    elif abs(pressure_mean - 1.0) < 10.0:
        convergence = "Reasonable convergence to target pressure"
    else:
        convergence = "Poor convergence to target pressure"
    
    plt.annotate(convergence, 
                 xy=(0.5, 0.05), xycoords='axes fraction',
                 bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="orange", alpha=0.8),
                 ha='center')
    
    plt.tight_layout()
    plt.savefig('../analysis/pressure_convergence_plot.png', dpi=300, bbox_inches='tight')
    
    return df

def analyze_energy_components(edr_file='md.edr'):
    """
    Analyze energy components from GROMACS energy file
    
    Parameters:
    -----------
    edr_file : str
        Path to the .edr file
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing energy data
    """
    print("Analyzing energy components...")
    
    # Parse energy data
    energy_data = parse_edr_file(edr_file, 
                                properties=['Potential', 'Kinetic', 'Total-Energy'])
    
    if not all(key in energy_data for key in ['Potential', 'Kinetic', 'Total-Energy']):
        print("Could not extract all energy components from EDR file")
        return None
    
    # Combine data into a single DataFrame
    df = pd.DataFrame({
        'Time': energy_data['Potential']['Time'],
        'Potential': energy_data['Potential']['Potential'],
        'Kinetic': energy_data['Kinetic']['Kinetic'],
        'Total': energy_data['Total-Energy']['Total-Energy']
    })
    
    # Calculate statistics
    potential_mean = df['Potential'].mean()
    potential_std = df['Potential'].std()
    kinetic_mean = df['Kinetic'].mean()
    kinetic_std = df['Kinetic'].std()
    total_mean = df['Total'].mean()
    total_std = df['Total'].std()
    
    # Calculate drift in total energy (in kJ/mol/ns)
    time_ns = df['Time'] / 1000  # Convert ps to ns
    total_energy = df['Total']
    
    # Linear regression to find the slope (drift)
    slope, intercept, r_value, p_value, std_err = stats.linregress(time_ns, total_energy)
    drift_rate = slope  # kJ/mol/ns
    
    # Save statistics
    with open('../analysis/energy_stats.dat', 'w') as f:
        f.write(f"# Energy statistics\n")
        f.write(f"# Potential energy: {potential_mean:.2f} ± {potential_std:.2f} kJ/mol\n")
        f.write(f"# Kinetic energy: {kinetic_mean:.2f} ± {kinetic_std:.2f} kJ/mol\n")
        f.write(f"# Total energy: {total_mean:.2f} ± {total_std:.2f} kJ/mol\n")
        f.write(f"# Energy drift: {drift_rate:.4f} kJ/mol/ns (R² = {r_value**2:.4f})\n")
        f.write(f"# Energy conservation: {abs(drift_rate/total_mean)*100:.6f}% drift per ns\n")
    
    # Save data
    df.to_csv('../analysis/energy_data.dat', sep='\t', index=False, 
              float_format='%.6f', header=True)
    
    # Plot energy components
    plt.figure(figsize=(12, 8))
    plt.plot(df['Time'], df['Potential'], 'b-', label=f'Potential: {potential_mean:.1f} ± {potential_std:.1f} kJ/mol')
    plt.plot(df['Time'], df['Kinetic'], 'r-', label=f'Kinetic: {kinetic_mean:.1f} ± {kinetic_std:.1f} kJ/mol')
    plt.plot(df['Time'], df['Total'], 'g-', label=f'Total: {total_mean:.1f} ± {total_std:.1f} kJ/mol')
    
    # Add horizontal lines for the means
    plt.axhline(y=potential_mean, color='b', linestyle='--', alpha=0.5)
    plt.axhline(y=kinetic_mean, color='r', linestyle='--', alpha=0.5)
    plt.axhline(y=total_mean, color='g', linestyle='--', alpha=0.5)
    
    # Add drift line for total energy
    drift_line = intercept + slope * time_ns
    plt.plot(df['Time'], drift_line, 'k--', alpha=0.7, 
             label=f'Energy drift: {drift_rate:.4f} kJ/mol/ns')
    
    plt.xlabel('Time (ps)')
    plt.ylabel('Energy (kJ/mol)')
    plt.title('Energy Components vs. Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add text box with energy conservation info
    conservation_percent = abs(drift_rate/total_mean)*100
    textstr = '\n'.join((
        f'Energy drift: {drift_rate:.4f} kJ/mol/ns',
        f'Relative drift: {conservation_percent:.6f}% per ns',
        f'R² of drift fit: {r_value**2:.4f}'))
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0.02, 0.05, textstr, transform=plt.gca().transAxes, fontsize=10,
            verticalalignment='bottom', bbox=props)
    
    plt.tight_layout()
    plt.savefig('../analysis/energy_plot.png', dpi=300, bbox_inches='tight')
    
    # Create a more detailed plot with subplots
    fig, axs = plt.subplots(3, 1, figsize=(12, 15), sharex=True)
    
    # Plot 1: Potential energy
    axs[0].plot(df['Time'], df['Potential'], 'b-')
    axs[0].axhline(y=potential_mean, color='b', linestyle='--', alpha=0.5)
    axs[0].fill_between(df['Time'], 
                       potential_mean - potential_std, 
                       potential_mean + potential_std, 
                       color='b', alpha=0.2)
    axs[0].set_ylabel('Potential Energy (kJ/mol)')
    axs[0].set_title(f'Potential Energy: {potential_mean:.1f} ± {potential_std:.1f} kJ/mol')
    axs[0].grid(True, alpha=0.3)
    
    # Plot 2: Kinetic energy
    axs[1].plot(df['Time'], df['Kinetic'], 'r-')
    axs[1].axhline(y=kinetic_mean, color='r', linestyle='--', alpha=0.5)
    axs[1].fill_between(df['Time'], 
                       kinetic_mean - kinetic_std, 
                       kinetic_mean + kinetic_std, 
                       color='r', alpha=0.2)
    axs[1].set_ylabel('Kinetic Energy (kJ/mol)')
    axs[1].set_title(f'Kinetic Energy: {kinetic_mean:.1f} ± {kinetic_std:.1f} kJ/mol')
    axs[1].grid(True, alpha=0.3)
    
    # Plot 3: Total energy with drift analysis
    axs[2].plot(df['Time'], df['Total'], 'g-')
    axs[2].plot(df['Time'], drift_line, 'k--', alpha=0.7, 
               label=f'Drift: {drift_rate:.4f} kJ/mol/ns')
    axs[2].axhline(y=total_mean, color='g', linestyle='--', alpha=0.5)
    axs[2].fill_between(df['Time'], 
                       total_mean - total_std, 
                       total_mean + total_std, 
                       color='g', alpha=0.2)
    axs[2].set_xlabel('Time (ps)')
    axs[2].set_ylabel('Total Energy (kJ/mol)')
    axs[2].set_title(f'Total Energy: {total_mean:.1f} ± {total_std:.1f} kJ/mol')
    axs[2].legend()
    axs[2].grid(True, alpha=0.3)
    
    # Add annotation about energy conservation
    if abs(conservation_percent) < 0.01:
        quality = "Excellent"
    elif abs(conservation_percent) < 0.1:
        quality = "Very good"
    elif abs(conservation_percent) < 0.5:
        quality = "Good"
    elif abs(conservation_percent) < 1.0:
        quality = "Acceptable"
    else:
        quality = "Poor"
    
    axs[2].annotate(f'Energy conservation: {quality}\n({conservation_percent:.6f}% drift per ns)', 
                   xy=(0.5, 0.05), xycoords='axes fraction',
                   bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="orange", alpha=0.8),
                   ha='center')
    
    plt.tight_layout()
    plt.savefig('../analysis/detailed_energy_plot.png', dpi=300, bbox_inches='tight')
    
    # Plot energy fluctuations (detrended)
    plt.figure(figsize=(12, 6))
    
    # Detrend the total energy
    detrended_total = df['Total'] - drift_line
    
    # Plot the detrended total energy
    plt.plot(df['Time'], detrended_total, 'g-', label='Detrended Total Energy')
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    # Add standard deviation bands
    detrended_std = np.std(detrended_total)
    plt.fill_between(df['Time'], 
                    -detrended_std, 
                    detrended_std, 
                    color='g', alpha=0.2,
                    label=f'±1σ ({detrended_std:.2f} kJ/mol)')
    
    plt.xlabel('Time (ps)')
    plt.ylabel('Detrended Energy (kJ/mol)')
    plt.title('Total Energy Fluctuations (Drift Removed)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add text about energy fluctuations
    textstr = '\n'.join((
        f'Std Dev: {detrended_std:.2f} kJ/mol',
        f'Relative fluctuation: {(detrended_std/abs(total_mean))*100:.4f}%'))
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0.02, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig('../analysis/energy_fluctuations_plot.png', dpi=300, bbox_inches='tight')
    
    return df

def calculate_energy_autocorrelation(energy_df, max_lag=1000):
    """
    Calculate autocorrelation function for energy components
    
    Parameters:
    -----------
    energy_df : pandas.DataFrame
        DataFrame containing energy data with 'Time (ps)' column
    max_lag : int
        Maximum lag time for autocorrelation
    """
    print("Calculating energy autocorrelation...")
    
    # Get time step
    time = energy_df['Time (ps)'].values
    dt = time[1] - time[0] if len(time) > 1 else 1.0
    
    # Calculate autocorrelation for each energy component
    acf_data = {}
    
    for column in energy_df.columns:
        if column == 'Time (ps)':
            continue
        
        # Get data and normalize
        data = energy_df[column].values
        data_norm = (data - np.mean(data)) / np.std(data)
        
        # Calculate autocorrelation
        acf = np.zeros(max_lag)
        for lag in range(max_lag):
            if lag == 0:
                acf[lag] = 1.0  # Autocorrelation at zero lag is 1 by definition
            else:
                # Manual calculation to avoid issues with missing data
                valid_indices = np.arange(len(data_norm) - lag)
                if len(valid_indices) > 0:
                    acf[lag] = np.mean(data_norm[valid_indices] * data_norm[valid_indices + lag])
                else:
                    acf[lag] = np.nan
        
        acf_data[column] = acf
    
    # Create lag time array
    lag_time = np.arange(max_lag) * dt
    
    # Save autocorrelation data
    acf_df = pd.DataFrame({'Time (ps)': lag_time})
    for column, acf in acf_data.items():
        acf_df[column] = acf
    
    acf_df.to_csv('../analysis/energy_acf_data.dat', sep='\t', float_format='%.6f', index=False)
    
    # Plot autocorrelation functions
    plt.figure(figsize=(12, 8))
    for column, acf in acf_data.items():
        plt.plot(lag_time, acf, label=column)
    
    plt.xlabel('Lag Time (ps)')
    plt.ylabel('Autocorrelation')
    plt.title('Energy Autocorrelation Functions')
    plt.legend()
    plt.grid(True, alpha=0.3)
    save_plot(plt, 'energy_acf.png')
    
    # Plot autocorrelation on log scale
    plt.figure(figsize=(12, 8))
    for column, acf in acf_data.items():
        # Convert to absolute value for log scale
        log_acf = np.abs(acf)
        plt.semilogy(lag_time, log_acf, label=column)
    
    plt.xlabel('Lag Time (ps)')
    plt.ylabel('|Autocorrelation|')
    plt.title('Energy Autocorrelation Functions (Log Scale)')
    plt.legend()
    plt.grid(True, alpha=0.3, which='both')
    save_plot(plt, 'energy_acf_log.png')
    
    return acf_df

def create_enhanced_analysis_plots(energy_df, temp_df=None, pressure_df=None):
    """
    Create enhanced analysis plots combining multiple thermodynamic properties
    
    Parameters:
    -----------
    energy_df : pandas.DataFrame
        DataFrame containing energy data
    temp_df : pandas.DataFrame or None
        DataFrame containing temperature data
    pressure_df : pandas.DataFrame or None
        DataFrame containing pressure data
    """
    print("Creating enhanced analysis plots...")
    
    # Create comprehensive energy analysis figure
    if 'Potential' in energy_df.columns and 'Kinetic' in energy_df.columns:
        plt.figure(figsize=(15, 12))
        gs = gridspec.GridSpec(3, 2, height_ratios=[2, 1, 1])
        
        # Energy vs time plot
        ax1 = plt.subplot(gs[0, :])
        ax1.plot(energy_df['Time (ps)'], energy_df['Potential'], 'r-', label='Potential')
        ax1.plot(energy_df['Time (ps)'], energy_df['Kinetic'], 'b-', label='Kinetic')
        
        if 'Total-Energy' in energy_df.columns:
            ax1.plot(energy_df['Time (ps)'], energy_df['Total-Energy'], 'g-', label='Total')
        
        ax1.set_xlabel('Time (ps)')
        ax1.set_ylabel('Energy (kJ/mol)')
        ax1.set_title('Energy Components vs Time')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Potential energy distribution
        ax2 = plt.subplot(gs[1, 0])
        ax2.hist(energy_df['Potential'], bins=30, alpha=0.7, color='red', edgecolor='black')
        ax2.axvline(energy_df['Potential'].mean(), color='k', linestyle='--')
        ax2.set_xlabel('Potential Energy (kJ/mol)')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Potential Energy Distribution')
        ax2.grid(True, alpha=0.3)
        
        # Kinetic energy distribution
        ax3 = plt.subplot(gs[1, 1])
        ax3.hist(energy_df['Kinetic'], bins=30, alpha=0.7, color='blue', edgecolor='black')
        ax3.axvline(energy_df['Kinetic'].mean(), color='k', linestyle='--')
        ax3.set_xlabel('Kinetic Energy (kJ/mol)')
        ax3.set_ylabel('Frequency')
        ax3.set_title('Kinetic Energy Distribution')
        ax3.grid(True, alpha=0.3)
        
        # Energy drift (Total - (Initial Total))
        if 'Total-Energy' in energy_df.columns:
            ax4 = plt.subplot(gs[2, 0])
            energy_drift = energy_df['Total-Energy'] - energy_df['Total-Energy'].iloc[0]
            ax4.plot(energy_df['Time (ps)'], energy_drift, 'g-')
            ax4.set_xlabel('Time (ps)')
            ax4.set_ylabel('Energy Drift (kJ/mol)')
            ax4.set_title('Total Energy Drift')
            ax4.grid(True, alpha=0.3)
            
            # Energy conservation
            ax5 = plt.subplot(gs[2, 1])
            energy_conservation = (
                np.abs(energy_df['Total-Energy'] - energy_df['Total-Energy'].iloc[0]) / 
                np.abs(energy_df['Total-Energy'].iloc[0])
            ) * 100  # In percent
            ax5.plot(energy_df['Time (ps)'], energy_conservation, 'g-')
            ax5.set_xlabel('Time (ps)')
            ax5.set_ylabel('Energy Drift (%)')
            ax5.set_title('Relative Energy Conservation')
            ax5.grid(True, alpha=0.3)
        
        plt.tight_layout()
        save_plot(plt, 'energy_enhanced_analysis.png')
    
    # Examine non-bonded interaction terms if available
    interaction_terms = [c for c in ['Coulomb-SR', 'LJ-SR', 'Coul.-recip.'] 
                         if c in energy_df.columns]
    
    if len(interaction_terms) >= 2:
        plt.figure(figsize=(15, 12))
        gs = gridspec.GridSpec(2, 2)
        
        # Terms vs time
        ax1 = plt.subplot(gs[0, :])
        for term in interaction_terms:
            ax1.plot(energy_df['Time (ps)'], energy_df[term], label=term)
        
        ax1.set_xlabel('Time (ps)')
        ax1.set_ylabel('Energy (kJ/mol)')
        ax1.set_title('Non-bonded Interaction Terms')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Individual histograms
        for i, term in enumerate(interaction_terms[:3]):  # At most 3 terms
            ax = plt.subplot(gs[1, i])
            ax.hist(energy_df[term], bins=30, alpha=0.7, edgecolor='black')
            ax.axvline(energy_df[term].mean(), color='r', linestyle='--')
            ax.set_xlabel(f'{term} (kJ/mol)')
            ax.set_ylabel('Frequency')
            ax.set_title(f'{term} Distribution')
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        save_plot(plt, 'energy_terms_enhanced_analysis.png')
    
    # Combined thermodynamic properties plot
    if temp_df is not None and pressure_df is not None:
        # Resample all dataframes to a common time base if necessary
        # For simplicity, we'll use energy_df's time points
        
        plt.figure(figsize=(15, 12))
        gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])
        
        # Top plot: Energy
        ax1 = plt.subplot(gs[0])
        if 'Potential' in energy_df.columns:
            ax1.plot(energy_df['Time (ps)'], energy_df['Potential'], 'r-', label='Potential Energy')
        if 'Kinetic' in energy_df.columns:
            ax1.plot(energy_df['Time (ps)'], energy_df['Kinetic'], 'b-', label='Kinetic Energy')
        
        ax1.set_xlabel('Time (ps)')
        ax1.set_ylabel('Energy (kJ/mol)')
        ax1.set_title('Energy Components')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Middle plot: Temperature
        ax2 = plt.subplot(gs[1], sharex=ax1)
        ax2.plot(temp_df['Time (ps)'], temp_df['Temperature'], 'g-')
        ax2.axhline(temp_df['Temperature'].mean(), color='k', linestyle='--',
                    label=f'Mean: {temp_df["Temperature"].mean():.2f} K')
        ax2.set_xlabel('Time (ps)')
        ax2.set_ylabel('Temperature (K)')
        ax2.set_title('Temperature')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Bottom plot: Pressure
        ax3 = plt.subplot(gs[2], sharex=ax1)
        ax3.plot(pressure_df['Time (ps)'], pressure_df['Pressure'], 'm-')
        ax3.axhline(pressure_df['Pressure'].mean(), color='k', linestyle='--',
                    label=f'Mean: {pressure_df["Pressure"].mean():.2f} bar')
        ax3.set_xlabel('Time (ps)')
        ax3.set_ylabel('Pressure (bar)')
        ax3.set_title('Pressure')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        plt.tight_layout()
        save_plot(plt, 'thermodynamic_properties_plot.png')

def main():
    # Get command line arguments if provided
    if len(sys.argv) > 1:
        edr_file = sys.argv[1]
    else:
        edr_file = 'md.edr'
    
    # Check if the EDR file exists
    parent_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
    
    # Check if file exists in current directory
    if not os.path.exists(edr_file):
        # Try in parent directory
        parent_edr = os.path.join(parent_dir, os.path.basename(edr_file))
        if os.path.exists(parent_edr):
            print(f"File {edr_file} not found in current directory, using {parent_edr} instead")
            edr_file = parent_edr
        else:
            print(f"Error: EDR file {edr_file} not found in current or parent directory")
            return
    
    # Analyze temperature
    temp_df = analyze_temperature(edr_file)
    
    # Analyze pressure
    pressure_df = analyze_pressure(edr_file)
    
    # Analyze energy components
    energy_df = analyze_energy_components(edr_file)
    
    if energy_df is not None:
        # Calculate energy autocorrelation
        acf_df = calculate_energy_autocorrelation(energy_df, max_lag=500)
        
        # Create enhanced analysis plots
        create_enhanced_analysis_plots(energy_df, temp_df, pressure_df)
    
    print("Thermodynamic analysis complete!")

if __name__ == '__main__':
    main() 