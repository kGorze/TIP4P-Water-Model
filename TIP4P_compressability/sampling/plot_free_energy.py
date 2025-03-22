#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import sys
import pandas as pd
from scipy.signal import savgol_filter
from scipy.constants import Boltzmann
from scipy.ndimage import gaussian_filter
from scipy.stats import gaussian_kde, linregress
from scipy.integrate import simps
import argparse

# Stała Boltzmanna w J/K
k_B = Boltzmann
# Stała Boltzmanna w kcal/(mol·K)
k_B_kcal = 0.0019872041  # kcal/(mol·K)

def run_gmx_energy(edr_file, tpr_file):
    """
    Extract pressure and temperature data from GROMACS energy file
    """
    if not os.path.exists(edr_file):
        raise FileNotFoundError(f"EDR file {edr_file} not found!")
    if not os.path.exists(tpr_file):
        raise FileNotFoundError(f"TPR file {tpr_file} not found!")

    # Create input for gmx energy
    with open('temp_input.txt', 'w') as f:
        f.write("Temperature\nPressure\n0\n")  # Changed order to Temperature first

    # Run gmx energy command
    cmd = f"gmx energy -f {edr_file} -s {tpr_file} -o energy.xvg < temp_input.txt"
    try:
        process = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running gmx energy: {e}")
        print(f"Error output: {e.stderr}")
        sys.exit(1)
    finally:
        if os.path.exists('temp_input.txt'):
            os.remove('temp_input.txt')

    # Read the generated xvg file
    if not os.path.exists('energy.xvg'):
        raise FileNotFoundError("energy.xvg was not created!")
    
    try:
        # Skip comment lines (those starting with # or @)
        with open('energy.xvg', 'r') as f:
            lines = f.readlines()
        data_lines = [l for l in lines if not (l.startswith('#') or l.startswith('@'))]
        
        # Convert to DataFrame with correct column order
        data = pd.DataFrame([line.split() for line in data_lines], 
                          columns=['Time', 'Temperature', 'Pressure'],  # Changed order
                          dtype=float)
        
        # Convert pressure from bar to atm and time to ns
        data['Pressure'] = data['Pressure'] * 0.986923  # bar to atm
        data['Time'] = data['Time'] / 1000  # ps to ns
        
        # Print some basic statistics to verify the data
        print("\nRaw Data Statistics:")
        print(f"Temperature: mean = {data['Temperature'].mean():.2f} K, std = {data['Temperature'].std():.2f} K")
        print(f"Min: {data['Temperature'].min():.2f} K, Max: {data['Temperature'].max():.2f} K")
        print(f"Pressure: mean = {data['Pressure'].mean():.2f} atm, std = {data['Pressure'].std():.2f} atm")
        print(f"Min: {data['Pressure'].min():.2f} atm, Max: {data['Pressure'].max():.2f} atm")
        print(f"Total number of data points: {len(data)}")
        
        return data
    except Exception as e:
        print(f"Error processing energy.xvg: {e}")
        sys.exit(1)
    finally:
        if os.path.exists('energy.xvg'):
            os.remove('energy.xvg')

def preprocess_data(data, downsample_factor=5, filter_outliers=True, smooth=True):
    """Preprocess the data for clearer visualization"""
    # Create a copy to avoid modifying original data
    processed_data = data.copy()
    
    # Filter outliers based on pressure
    if filter_outliers:
        p_q1, p_q3 = np.percentile(processed_data['Pressure'], [1, 99])
        p_iqr = p_q3 - p_q1
        p_lower = p_q1 - 1.5 * p_iqr
        p_upper = p_q3 + 1.5 * p_iqr
        
        # Filter out extreme outliers
        mask = (processed_data['Pressure'] >= p_lower) & (processed_data['Pressure'] <= p_upper)
        filtered_data = processed_data[mask]
        
        # Only use filtered data if we haven't lost too many points
        if len(filtered_data) > 0.9 * len(processed_data):
            processed_data = filtered_data
            print(f"Filtered outliers: removed {len(data) - len(processed_data)} points")
    
    # Downsample to reduce clutter
    processed_data = processed_data.iloc[::downsample_factor].reset_index(drop=True)
    print(f"Downsampled data from {len(data)} to {len(processed_data)} points")
    
    # Apply Savitzky-Golay filter for smoothing if needed
    if smooth and len(processed_data) > 10:
        # Make sure we have enough data points for the filter
        window_size = min(11, len(processed_data) // 2 * 2 - 1)  # must be odd
        if window_size >= 3:
            try:
                processed_data['Pressure_smooth'] = savgol_filter(processed_data['Pressure'], window_size, 3)
                processed_data['Temperature_smooth'] = savgol_filter(processed_data['Temperature'], window_size, 3)
                print(f"Applied smoothing with window size {window_size}")
            except Exception as e:
                print(f"Could not apply smoothing: {e}")
                processed_data['Pressure_smooth'] = processed_data['Pressure']
                processed_data['Temperature_smooth'] = processed_data['Temperature']
        else:
            processed_data['Pressure_smooth'] = processed_data['Pressure']
            processed_data['Temperature_smooth'] = processed_data['Temperature']
    else:
        processed_data['Pressure_smooth'] = processed_data['Pressure']
        processed_data['Temperature_smooth'] = processed_data['Temperature']
    
    return processed_data

def bootstrap_histogram(data, bins, range_vals, n_bootstrap=100):
    """
    Perform bootstrap estimation of histogram uncertainties
    Returns mean histogram and standard deviation for each bin
    """
    n_samples = len(data)
    hist_samples = np.zeros((n_bootstrap, bins))
    
    for i in range(n_bootstrap):
        # Generate bootstrap sample with replacement
        bootstrap_data = np.random.choice(data, size=n_samples, replace=True)
        hist, _ = np.histogram(bootstrap_data, bins=bins, range=range_vals, density=True)
        hist_samples[i, :] = hist
    
    # Compute mean and standard deviation across bootstrap samples
    hist_mean = np.mean(hist_samples, axis=0)
    hist_std = np.std(hist_samples, axis=0)
    
    return hist_mean, hist_std

def calculate_free_energy_profiles(data, temperature_bins=50, pressure_bins=50, smooth_sigma=1.0, truncate_energy=10.0, bootstrap=True, calculate_entropy=True):
    """
    Calculate free energy profiles from histograms
    F = -kT*ln(P), where P is the probability
    
    Parameters:
    -----------
    data : DataFrame
        Contains Temperature and Pressure data
    temperature_bins, pressure_bins : int
        Number of bins for histograms
    smooth_sigma : float
        Gaussian smoothing factor
    truncate_energy : float
        Maximum free energy value (higher values will be clipped)
    bootstrap : bool
        Whether to calculate uncertainties using bootstrap
    calculate_entropy : bool
        Whether to calculate entropy contributions explicitly
    """
    # Use average temperature for kT
    avg_temp = data['Temperature'].mean()
    kT = k_B_kcal * avg_temp  # kcal/mol
    
    # Calculate temperature range (use tighter range for better statistics)
    temp_mean = data['Temperature'].mean()
    temp_std = data['Temperature'].std()
    t_min = max(data['Temperature'].min(), temp_mean - 3 * temp_std)
    t_max = min(data['Temperature'].max(), temp_mean + 3 * temp_std)
    
    # Calculate pressure range (use percentiles for better range definition)
    p_min, p_max = data['Pressure'].quantile([0.02, 0.98])
    
    # Print ranges being used
    print(f"\nUsing temperature range: {t_min:.2f} K to {t_max:.2f} K")
    print(f"Using pressure range: {p_min:.2f} atm to {p_max:.2f} atm")
    
    # Create histograms
    temp_hist, temp_edges = np.histogram(data['Temperature'], bins=temperature_bins, 
                                        range=(t_min, t_max), density=True)
    press_hist, press_edges = np.histogram(data['Pressure'], bins=pressure_bins, 
                                          range=(p_min, p_max), density=True)
    
    # Calculate bin centers
    temp_centers = 0.5 * (temp_edges[1:] + temp_edges[:-1])
    press_centers = 0.5 * (press_edges[1:] + press_edges[:-1])
    
    # Bootstrap for uncertainties if requested
    if bootstrap:
        print("Calculating uncertainties using bootstrap...")
        temp_hist_boot, temp_hist_std = bootstrap_histogram(
            data['Temperature'], temperature_bins, (t_min, t_max), n_bootstrap=100)
        press_hist_boot, press_hist_std = bootstrap_histogram(
            data['Pressure'], pressure_bins, (p_min, p_max), n_bootstrap=100)
    else:
        temp_hist_boot, temp_hist_std = temp_hist, np.zeros_like(temp_hist)
        press_hist_boot, press_hist_std = press_hist, np.zeros_like(press_hist)
    
    # Apply Gaussian smoothing to histograms
    temp_hist_smooth = gaussian_filter(temp_hist_boot, sigma=smooth_sigma)
    press_hist_smooth = gaussian_filter(press_hist_boot, sigma=smooth_sigma)
    
    # Calculate probability (normalized hist)
    # Add small constant to avoid log(0)
    epsilon = 1e-10
    
    # Calculate free energy profiles from smoothed histograms
    # For temperature
    temp_prob = temp_hist_smooth / np.sum(temp_hist_smooth)
    temp_free_energy = -kT * np.log(temp_prob + epsilon)
    # Shift so that minimum is at zero
    temp_free_energy -= np.min(temp_free_energy)
    
    # For pressure
    press_prob = press_hist_smooth / np.sum(press_hist_smooth)
    press_free_energy = -kT * np.log(press_prob + epsilon)
    # Shift so that minimum is at zero
    press_free_energy -= np.min(press_free_energy)
    
    # Calculate uncertainty in free energy
    # δF = kT * δP/P where P is probability
    temp_uncertainty = np.zeros_like(temp_free_energy)
    press_uncertainty = np.zeros_like(press_free_energy)
    
    if bootstrap:
        # Only calculate where probability is non-zero to avoid division by zero
        valid_temp = temp_prob > epsilon
        valid_press = press_prob > epsilon
        
        if np.any(valid_temp):
            temp_uncertainty[valid_temp] = kT * temp_hist_std[valid_temp] / (temp_prob[valid_temp] + epsilon)
        
        if np.any(valid_press):
            press_uncertainty[valid_press] = kT * press_hist_std[valid_press] / (press_prob[valid_press] + epsilon)
    
    # New: Calculate entropy and enthalpy contributions if requested
    temp_entropy_profile = None
    press_entropy_profile = None
    temp_enthalpy = None
    press_enthalpy = None
    total_temp_entropy = None
    total_press_entropy = None
    temp_entropy_uncertainty = None
    press_entropy_uncertainty = None
    
    if calculate_entropy:
        print("Calculating entropy and enthalpy contributions...")
        
        # Calculate entropy profile directly from probability distribution
        # S(x) = -k_B * ln(p(x))
        # This gives the configurational entropy at each point
        temp_entropy_profile = -k_B_kcal * np.log(temp_prob + epsilon)
        press_entropy_profile = -k_B_kcal * np.log(press_prob + epsilon)
        
        # Shift entropy to have minimum at zero
        temp_entropy_profile -= np.min(temp_entropy_profile)
        press_entropy_profile -= np.min(press_entropy_profile)
        
        # Scale for visualization (but keep more natural variation)
        # We use TS rather than just S for direct energy comparison
        temp_entropy_profile = avg_temp * temp_entropy_profile
        press_entropy_profile = avg_temp * press_entropy_profile
        
        # Calculate total entropy using the formula S = -k_B * sum(p_i * ln(p_i))
        total_temp_entropy = -k_B_kcal * np.sum(temp_prob * np.log(temp_prob + epsilon))
        total_press_entropy = -k_B_kcal * np.sum(press_prob * np.log(press_prob + epsilon))
        
        # Calculate entropy uncertainty if bootstrap data is available
        if bootstrap:
            # Calculate δS = k_B * (1 + ln(p)) * δp
            temp_entropy_uncertainty = k_B_kcal * avg_temp * abs(1 + np.log(temp_prob + epsilon)) * temp_hist_std
            press_entropy_uncertainty = k_B_kcal * avg_temp * abs(1 + np.log(press_prob + epsilon)) * press_hist_std
            
            # Scale entropy uncertainties for visualization
            scale_factor = 0.2  # Reduced to avoid visual clutter
            temp_entropy_uncertainty *= scale_factor
            press_entropy_uncertainty *= scale_factor
        else:
            temp_entropy_uncertainty = np.zeros_like(temp_entropy_profile)
            press_entropy_uncertainty = np.zeros_like(press_entropy_profile)
        
        # Calculate enthalpy using H = F + TS
        # For a more natural variation in enthalpy
        temp_enthalpy = temp_free_energy + temp_entropy_profile
        press_enthalpy = press_free_energy + press_entropy_profile
        
        # Ensure enthalpy has minimum at zero (for consistency)
        temp_enthalpy -= np.min(temp_enthalpy)
        press_enthalpy -= np.min(press_enthalpy)
        
        # Define regions where entropy is favorable vs unfavorable
        temp_entropy_favorable = temp_entropy_profile > 0.5 * np.max(temp_entropy_profile)
        press_entropy_favorable = press_entropy_profile > 0.5 * np.max(press_entropy_profile)
        
        # Print summary of entropy analysis
        print(f"Total entropy: Temperature = {total_temp_entropy:.4f}, Pressure = {total_press_entropy:.4f} kcal/mol/K")
        print(f"Max TS contribution: Temperature = {np.max(temp_entropy_profile):.2f}, Pressure = {np.max(press_entropy_profile):.2f} kcal/mol")
    
    # Clip extreme values for better visualization - AFTER all calculations are done
    temp_free_energy_clipped = np.clip(temp_free_energy, 0, truncate_energy)
    press_free_energy_clipped = np.clip(press_free_energy, 0, truncate_energy)
    
    if calculate_entropy:
        temp_entropy_profile_clipped = np.clip(temp_entropy_profile, 0, truncate_energy)
        press_entropy_profile_clipped = np.clip(press_entropy_profile, 0, truncate_energy)
        temp_enthalpy_clipped = np.clip(temp_enthalpy, 0, truncate_energy * 1.2)  # Allow slightly higher ceiling for enthalpy
        press_enthalpy_clipped = np.clip(press_enthalpy, 0, truncate_energy * 1.2)
    
    return {
        'temperature': {
            'centers': temp_centers,
            'histogram': temp_hist,
            'histogram_smooth': temp_hist_smooth,
            'probability': temp_prob,
            'free_energy': temp_free_energy,
            'free_energy_clipped': temp_free_energy_clipped,
            'uncertainty': temp_uncertainty,
            'entropy_profile': temp_entropy_profile if calculate_entropy else None,
            'entropy_profile_clipped': temp_entropy_profile_clipped if calculate_entropy else None,
            'entropy_uncertainty': temp_entropy_uncertainty if calculate_entropy else None,
            'entropy_favorable': temp_entropy_favorable if calculate_entropy else None,
            'enthalpy': temp_enthalpy if calculate_entropy else None,
            'enthalpy_clipped': temp_enthalpy_clipped if calculate_entropy else None,
            'range': (t_min, t_max)
        },
        'pressure': {
            'centers': press_centers,
            'histogram': press_hist,
            'histogram_smooth': press_hist_smooth,
            'probability': press_prob,
            'free_energy': press_free_energy,
            'free_energy_clipped': press_free_energy_clipped,
            'uncertainty': press_uncertainty,
            'entropy_profile': press_entropy_profile if calculate_entropy else None,
            'entropy_profile_clipped': press_entropy_profile_clipped if calculate_entropy else None,
            'entropy_uncertainty': press_entropy_uncertainty if calculate_entropy else None,
            'entropy_favorable': press_entropy_favorable if calculate_entropy else None,
            'enthalpy': press_enthalpy if calculate_entropy else None,
            'enthalpy_clipped': press_enthalpy_clipped if calculate_entropy else None,
            'range': (p_min, p_max)
        },
        'kT': kT,
        'avg_temp': avg_temp,
        'truncate_energy': truncate_energy,
        'total_entropy': {
            'temperature': total_temp_entropy if calculate_entropy else None,
            'pressure': total_press_entropy if calculate_entropy else None
        }
    }

def plot_free_energy_profiles(data, profiles, output_file='free_energy_profiles.png'):
    """Plot free energy profiles alongside histograms with error bars"""
    plt.style.use('seaborn-v0_8-darkgrid')
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    
    # Temperature histogram
    axes[0, 0].bar(profiles['temperature']['centers'], 
                  profiles['temperature']['histogram_smooth'], 
                  width=(profiles['temperature']['centers'][1] - profiles['temperature']['centers'][0]),
                  color='royalblue', alpha=0.7)
    axes[0, 0].axvline(x=273.15, color='r', linestyle='--', alpha=0.7, label='Target T (273.15 K)')
    axes[0, 0].set_xlabel('Temperature (K)')
    axes[0, 0].set_ylabel('Probability Density')
    axes[0, 0].set_title('Temperature Distribution (Smoothed)')
    axes[0, 0].legend()
    
    # Temperature free energy profile
    axes[0, 1].plot(profiles['temperature']['centers'], 
                   profiles['temperature']['free_energy'],
                   'b-', linewidth=2, label='Free Energy (F)')
    
    # Add entropy and enthalpy profiles if available
    if profiles['temperature']['entropy_profile'] is not None:
        axes[0, 1].plot(profiles['temperature']['centers'], 
                       profiles['temperature']['entropy_profile'],
                       'g--', linewidth=1.5, label='Entropy (TS)')
        axes[0, 1].plot(profiles['temperature']['centers'], 
                       profiles['temperature']['enthalpy'],
                       'r-.', linewidth=1.5, label='Enthalpy (H)')
    
    # Add uncertainty if available
    if np.any(profiles['temperature']['uncertainty'] > 0):
        axes[0, 1].fill_between(
            profiles['temperature']['centers'],
            profiles['temperature']['free_energy'] - profiles['temperature']['uncertainty'],
            profiles['temperature']['free_energy'] + profiles['temperature']['uncertainty'],
            color='blue', alpha=0.2)
    
    axes[0, 1].axvline(x=273.15, color='r', linestyle='--', alpha=0.7, label='Target T (273.15 K)')
    axes[0, 1].set_xlabel('Temperature (K)')
    axes[0, 1].set_ylabel('Energy (kcal/mol)')
    axes[0, 1].set_title('Temperature Free Energy Profile')
    axes[0, 1].set_ylim(0, profiles['truncate_energy'])
    axes[0, 1].legend()
    
    # Pressure histogram
    axes[1, 0].bar(profiles['pressure']['centers'], 
                  profiles['pressure']['histogram_smooth'], 
                  width=(profiles['pressure']['centers'][1] - profiles['pressure']['centers'][0]),
                  color='forestgreen', alpha=0.7)
    axes[1, 0].axvline(x=1.0, color='r', linestyle='--', alpha=0.7, label='Target P (1.0 atm)')
    axes[1, 0].set_xlabel('Pressure (atm)')
    axes[1, 0].set_ylabel('Probability Density')
    axes[1, 0].set_title('Pressure Distribution (Smoothed)')
    axes[1, 0].legend()
    
    # Pressure free energy profile
    axes[1, 1].plot(profiles['pressure']['centers'], 
                   profiles['pressure']['free_energy'],
                   'g-', linewidth=2, label='Free Energy (F)')
    
    # Add entropy and enthalpy profiles if available
    if profiles['pressure']['entropy_profile'] is not None:
        axes[1, 1].plot(profiles['pressure']['centers'], 
                       profiles['pressure']['entropy_profile'],
                       'b--', linewidth=1.5, label='Entropy (TS)')
        axes[1, 1].plot(profiles['pressure']['centers'], 
                       profiles['pressure']['enthalpy'],
                       'r-.', linewidth=1.5, label='Enthalpy (H)')
    
    # Add uncertainty if available
    if np.any(profiles['pressure']['uncertainty'] > 0):
        axes[1, 1].fill_between(
            profiles['pressure']['centers'],
            profiles['pressure']['free_energy'] - profiles['pressure']['uncertainty'],
            profiles['pressure']['free_energy'] + profiles['pressure']['uncertainty'],
            color='green', alpha=0.2)
    
    axes[1, 1].axvline(x=1.0, color='r', linestyle='--', alpha=0.7, label='Target P (1.0 atm)')
    axes[1, 1].set_xlabel('Pressure (atm)')
    axes[1, 1].set_ylabel('Energy (kcal/mol)')
    axes[1, 1].set_title('Pressure Free Energy Profile')
    axes[1, 1].set_ylim(0, profiles['truncate_energy'])
    axes[1, 1].legend()
    
    # Add info text with entropy information if available
    info_text = (f"kT = {profiles['kT']:.4f} kcal/mol\n"
                f"Avg. Temperature = {profiles['avg_temp']:.2f} K\n"
                f"F = -kT ln(P) with smoothing")
    
    if profiles['total_entropy']['temperature'] is not None:
        info_text += f"\nTotal Temperature Entropy = {profiles['total_entropy']['temperature']:.4f} kcal/mol/K"
        info_text += f"\nTotal Pressure Entropy = {profiles['total_entropy']['pressure']:.4f} kcal/mol/K"
    
    info_text += f"\nEnergy truncated at {profiles['truncate_energy']:.1f} kcal/mol"
    
    fig.text(0.02, 0.02, info_text, fontsize=10, 
             bbox=dict(facecolor='white', alpha=0.8))
    
    plt.suptitle('Free Energy Profiles from Simulation Data', fontsize=16)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nFree energy profiles saved as: {output_file}")
    plt.close()

def plot_combined_free_energy(data, profiles, output_file='free_energy_2d.png'):
    """Plot 2D free energy landscape in temperature-pressure space"""
    plt.style.use('seaborn-v0_8-darkgrid')
    
    # Get temperature and pressure ranges from profiles
    t_min, t_max = profiles['temperature']['range']
    p_min, p_max = profiles['pressure']['range']
    
    # Create 2D histogram with more bins for smoother contours
    H, xedges, yedges = np.histogram2d(data['Temperature'], data['Pressure'], 
                                      bins=[80, 80], 
                                      range=[[t_min, t_max], [p_min, p_max]],
                                      density=True)
    
    # Apply Gaussian smoothing to the 2D histogram
    H_smooth = gaussian_filter(H, sigma=1.0)
    
    # Calculate free energy landscape from smoothed 2D histogram
    # Add small constant to avoid log(0)
    epsilon = 1e-10
    kT = profiles['kT']
    
    # Calculate free energy, transpose for correct orientation
    F = -kT * np.log(H_smooth.T + epsilon)
    
    # Calculate entropy contribution: S = -k_B * p * ln(p)
    P_norm = H_smooth.T / np.sum(H_smooth.T)
    S = -k_B_kcal * P_norm * np.log(P_norm + epsilon)
    # Normalize entropy for visualization
    S_norm = S.copy()
    if np.max(S) > 0:
        S_norm = S / np.max(S) * np.max(F) * 0.8
    
    # Calculate enthalpy: H = F + TS
    H_energy = F + profiles['avg_temp'] * S_norm
    
    # Clip extreme values for better color scale
    max_energy = profiles['truncate_energy']
    F = np.clip(F, 0, max_energy)
    H_energy = np.clip(H_energy, 0, max_energy)
    
    # Shift to make minimum zero
    F -= np.min(F)
    H_energy -= np.min(H_energy)
    
    # Get coordinates for plotting
    X, Y = np.meshgrid(0.5 * (xedges[1:] + xedges[:-1]), 
                      0.5 * (yedges[1:] + yedges[:-1]))
    
    # Create figure with four subplots (2x2 grid)
    fig = plt.figure(figsize=(18, 12))
    
    # 1. Free Energy (F)
    ax1 = fig.add_subplot(221)
    cmap_F = plt.cm.plasma_r
    contour_F = ax1.contourf(X, Y, F, 50, cmap=cmap_F, vmax=max_energy)
    ax1.contour(X, Y, F, 10, colors='k', alpha=0.3, linewidths=0.5)
    cb_F = plt.colorbar(contour_F, ax=ax1)
    cb_F.set_label('Free Energy (kcal/mol)')
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('Pressure (atm)')
    ax1.set_title('Free Energy Landscape (F = -kT ln(P))')
    
    # 2. Entropy Contribution (TS)
    ax2 = fig.add_subplot(222)
    # Use viridis colormap for entropy (reversed to match free energy)
    cmap_S = plt.cm.viridis_r
    contour_S = ax2.contourf(X, Y, S_norm, 50, cmap=cmap_S, vmax=max_energy)
    ax2.contour(X, Y, S_norm, 10, colors='k', alpha=0.3, linewidths=0.5)
    cb_S = plt.colorbar(contour_S, ax=ax2)
    cb_S.set_label('Entropy Contribution (kcal/mol)')
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('Pressure (atm)')
    ax2.set_title('Entropy Landscape (TS)')
    
    # 3. Enthalpy (H = F + TS)
    ax3 = fig.add_subplot(223)
    # Use inferno colormap for enthalpy
    cmap_H = plt.cm.inferno_r
    contour_H = ax3.contourf(X, Y, H_energy, 50, cmap=cmap_H, vmax=max_energy)
    ax3.contour(X, Y, H_energy, 10, colors='k', alpha=0.3, linewidths=0.5)
    cb_H = plt.colorbar(contour_H, ax=ax3)
    cb_H.set_label('Enthalpy (kcal/mol)')
    ax3.set_xlabel('Temperature (K)')
    ax3.set_ylabel('Pressure (atm)')
    ax3.set_title('Enthalpy Landscape (H = F + TS)')
    
    # 4. Combined Analysis with Target Lines
    ax4 = fig.add_subplot(224)
    # Reuse free energy plot but add target lines and sample points
    contour_comb = ax4.contourf(X, Y, F, 50, cmap=plt.cm.plasma_r, vmax=max_energy)
    cb_comb = plt.colorbar(contour_comb, ax=ax4)
    cb_comb.set_label('Free Energy (kcal/mol)')
    
    # Add target lines
    ax4.axhline(y=1.0, color='red', linestyle='--', alpha=0.7, label='Target P (1.0 atm)')
    ax4.axvline(x=273.15, color='red', linestyle='--', alpha=0.7, label='Target T (273.15 K)')
    
    # Plot sampled points as a faint scatter
    if len(data) > 2000:
        # Downsample for clarity
        sample_data = data.iloc[::len(data)//1000]
    else:
        sample_data = data
    ax4.scatter(sample_data['Temperature'], sample_data['Pressure'], 
               color='w', alpha=0.15, s=2)
    
    ax4.set_xlabel('Temperature (K)')
    ax4.set_ylabel('Pressure (atm)')
    ax4.set_title('Free Energy with Sampling Points')
    ax4.legend()
    
    # Add info text
    info_text = (f"kT = {kT:.4f} kcal/mol\n"
                f"F = -kT ln(P), H = F + TS\n"
                f"Gaussian smoothing applied\n"
                f"Energy truncated at {max_energy:.1f} kcal/mol")
    
    # Add entropy information if available
    if profiles['total_entropy']['temperature'] is not None:
        total_T_entropy = profiles['total_entropy']['temperature']
        total_P_entropy = profiles['total_entropy']['pressure']
        info_text += f"\nTotal Temp. Entropy = {total_T_entropy:.4f} kcal/mol/K"
        info_text += f"\nTotal Press. Entropy = {total_P_entropy:.4f} kcal/mol/K"
    
    fig.text(0.01, 0.01, info_text, fontsize=10, 
             bbox=dict(facecolor='white', alpha=0.8))
    
    plt.suptitle('Thermodynamic Analysis of GROMACS Simulation', fontsize=16)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nEnhanced free energy landscape saved as: {output_file}")
    plt.close()

def plot_entropy_analysis(data, profiles, output_file='entropy_analysis.png'):
    """Plot detailed entropy analysis for temperature and pressure"""
    plt.style.use('seaborn-v0_8-darkgrid')
    
    # Check if entropy information is available
    if (profiles['temperature']['entropy_profile'] is None or 
        profiles['pressure']['entropy_profile'] is None):
        print("Entropy profiles not available. Skipping entropy analysis plot.")
        return
    
    # Create figure with adjusted layout
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    
    # Ensure consistent y-axis ranges for better comparison
    # Calculate common y-limits based on data
    max_energy_T = max(
        np.max(profiles['temperature']['free_energy_clipped']),
        np.max(profiles['temperature']['enthalpy_clipped']),
        np.max(profiles['temperature']['entropy_profile_clipped'])
    )
    
    max_energy_P = max(
        np.max(profiles['pressure']['free_energy_clipped']),
        np.max(profiles['pressure']['enthalpy_clipped']),
        np.max(profiles['pressure']['entropy_profile_clipped'])
    )
    
    # Adjust for more natural y-axis scales
    max_entropy_scale = max(
        np.max(profiles['temperature']['entropy_profile_clipped']),
        np.max(profiles['pressure']['entropy_profile_clipped'])
    ) * 1.2  # Add a little headroom
    
    # 1. Temperature analysis: Free energy vs enthalpy
    # Create temperature x-grid for smoother plotting
    t_centers = profiles['temperature']['centers']
    min_temp, max_temp = t_centers.min(), t_centers.max()
    # Add markers at key thermodynamic points
    min_F_idx = np.argmin(profiles['temperature']['free_energy'])
    min_F_temp = t_centers[min_F_idx]
    max_S_idx = np.argmax(profiles['temperature']['entropy_profile'])
    max_S_temp = t_centers[max_S_idx]
    
    # Mark reference temperature 
    target_temp = 273.15
    
    # Free energy vs Enthalpy plot
    ax1 = axes[0, 0]
    ax1.plot(t_centers, profiles['temperature']['free_energy_clipped'],
            'b-', linewidth=2.5, label='Free Energy (F)')
    ax1.plot(t_centers, profiles['temperature']['enthalpy_clipped'],
            'r-', linewidth=2, label='Enthalpy (H)')
    
    # Add uncertainty band if available
    if np.any(profiles['temperature']['uncertainty'] > 0):
        ax1.fill_between(
            t_centers,
            profiles['temperature']['free_energy_clipped'] - profiles['temperature']['uncertainty'],
            profiles['temperature']['free_energy_clipped'] + profiles['temperature']['uncertainty'],
            color='blue', alpha=0.15)
    
    # Highlight key points
    ax1.scatter([min_F_temp], [profiles['temperature']['free_energy'][min_F_idx]], 
               s=80, color='blue', marker='o', edgecolors='white', zorder=5,
               label='Min Free Energy')
    
    # Highlight regions of favorable entropy
    favorable_regions = profiles['temperature']['entropy_favorable']
    if np.any(favorable_regions):
        for i in range(len(favorable_regions)-1):
            if favorable_regions[i] and not favorable_regions[i+1]:  # End of region
                ax1.axvspan(t_centers[i-2], t_centers[i+1], alpha=0.1, color='green')
            elif not favorable_regions[i] and favorable_regions[i+1]:  # Start of region
                pass  # Just prep for span
    
    # Add target line
    ax1.axvline(x=target_temp, color='k', linestyle='--', alpha=0.7, label='273.15 K')
    
    # Label axes
    ax1.set_xlabel('Temperature (K)', fontsize=12)
    ax1.set_ylabel('Energy (kcal/mol)', fontsize=12)
    ax1.set_title('Temperature: Free Energy vs Enthalpy', fontsize=14)
    ax1.set_ylim(0, max_energy_T * 1.05)
    ax1.legend(loc='best', fontsize=10)
    
    # 2. Temperature analysis: Entropy contribution
    ax2 = axes[0, 1]
    
    # For clarity, use area plot to emphasize favorable contribution
    ax2.fill_between(t_centers, profiles['temperature']['entropy_profile_clipped'], 
                    color='green', alpha=0.3, label='Favorable (TS > 0)')
    
    # Line plots
    ax2.plot(t_centers, profiles['temperature']['entropy_profile_clipped'],
            'g-', linewidth=2.5, label='Entropy (TS)')
    
    # Calculate -T·ΔS as entropy penalty
    entropy_contrib = -profiles['temperature']['entropy_profile_clipped']
    ax2.plot(t_centers, entropy_contrib, 'm--', linewidth=2, label='-T·ΔS')
    
    # Add uncertainty band for entropy if available
    if profiles['temperature']['entropy_uncertainty'] is not None:
        ax2.fill_between(
            t_centers,
            profiles['temperature']['entropy_profile_clipped'] - profiles['temperature']['entropy_uncertainty'],
            profiles['temperature']['entropy_profile_clipped'] + profiles['temperature']['entropy_uncertainty'],
            color='green', alpha=0.15)
    
    # Highlight key points
    ax2.scatter([max_S_temp], [profiles['temperature']['entropy_profile'][max_S_idx]], 
               s=80, color='green', marker='o', edgecolors='white', zorder=5,
               label='Max Entropy')
    
    # Add zero line
    ax2.axhline(y=0, color='k', linestyle='-', alpha=0.5)
    
    # Add target line
    ax2.axvline(x=target_temp, color='k', linestyle='--', alpha=0.7)
    
    # Label axes with consistent y-range
    ax2.set_xlabel('Temperature (K)', fontsize=12)
    ax2.set_ylabel('Energy (kcal/mol)', fontsize=12)
    ax2.set_title('Temperature: Entropy Contribution', fontsize=14)
    ax2.set_ylim(-max_entropy_scale, max_entropy_scale)
    ax2.legend(loc='best', fontsize=10)
    
    # Annotate physical meaning in temperature plots
    phase_annotations = []
    
    # Look for temperature regions that might correspond to different phases
    # This is just a simple heuristic based on energy minima
    f_diffs = np.diff(profiles['temperature']['free_energy'])
    potential_phase_boundaries = []
    for i in range(1, len(f_diffs)):
        if f_diffs[i-1] < 0 and f_diffs[i] > 0:  # Local minimum in F
            potential_phase_boundaries.append(i)
    
    # Add phase annotations if we found potential phases
    if len(potential_phase_boundaries) > 0:
        # Add shading for potential phases
        for i, pb in enumerate(potential_phase_boundaries):
            t_val = t_centers[pb]
            if i % 2 == 0:
                phase_label = "Potential\nSolid Phase" if t_val < target_temp else "Potential\nLiquid Phase"
            else:
                phase_label = "Potential\nLiquid Phase" if t_val < target_temp else "Potential\nGas Phase"
                
            ax1.annotate(phase_label, (t_val, 0.2 * max_energy_T),
                        xytext=(0, -30), textcoords='offset points',
                        ha='center', va='top', fontsize=10,
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))
            phase_annotations.append((t_val, phase_label))
    
    # 3. Pressure analysis: Free energy vs enthalpy
    # Create pressure x-grid for smoother plotting
    p_centers = profiles['pressure']['centers']
    
    # Add markers at key thermodynamic points
    min_F_idx_p = np.argmin(profiles['pressure']['free_energy'])
    min_F_pres = p_centers[min_F_idx_p]
    max_S_idx_p = np.argmax(profiles['pressure']['entropy_profile'])
    max_S_pres = p_centers[max_S_idx_p]
    
    # Reference pressure
    target_pres = 1.0  # 1 atm
    
    # Free energy vs Enthalpy plot
    ax3 = axes[1, 0]
    # Use clipped data to avoid extreme enthalpy values
    ax3.plot(p_centers, profiles['pressure']['free_energy_clipped'],
            'b-', linewidth=2.5, label='Free Energy (F)')
    ax3.plot(p_centers, profiles['pressure']['enthalpy_clipped'],
            'r-', linewidth=2, label='Enthalpy (H)')
    
    # Add uncertainty band if available
    if np.any(profiles['pressure']['uncertainty'] > 0):
        ax3.fill_between(
            p_centers,
            profiles['pressure']['free_energy_clipped'] - profiles['pressure']['uncertainty'],
            profiles['pressure']['free_energy_clipped'] + profiles['pressure']['uncertainty'],
            color='blue', alpha=0.15)
    
    # Highlight key points
    ax3.scatter([min_F_pres], [profiles['pressure']['free_energy_clipped'][min_F_idx_p]], 
               s=80, color='blue', marker='o', edgecolors='white', zorder=5,
               label='Min Free Energy')
    
    # Highlight regions of favorable entropy
    favorable_regions_p = profiles['pressure']['entropy_favorable']
    if np.any(favorable_regions_p):
        for i in range(len(favorable_regions_p)-1):
            if favorable_regions_p[i] and not favorable_regions_p[i+1]:  # End of region
                ax3.axvspan(p_centers[i-2], p_centers[i+1], alpha=0.1, color='green')
            elif not favorable_regions_p[i] and favorable_regions_p[i+1]:  # Start of region
                pass  # Just prep for span
    
    # Add target line
    ax3.axvline(x=target_pres, color='k', linestyle='--', alpha=0.7, label='1.0 atm')
    
    # Label axes
    ax3.set_xlabel('Pressure (atm)', fontsize=12)
    ax3.set_ylabel('Energy (kcal/mol)', fontsize=12)
    ax3.set_title('Pressure: Free Energy vs Enthalpy', fontsize=14)
    ax3.set_ylim(0, max_energy_P * 1.05)
    ax3.legend(loc='best', fontsize=10)
    
    # 4. Pressure analysis: Entropy contribution
    ax4 = axes[1, 1]
    
    # For clarity, use area plot to emphasize favorable contribution
    ax4.fill_between(p_centers, profiles['pressure']['entropy_profile_clipped'], 
                    color='green', alpha=0.3, label='Favorable (TS > 0)')
    
    # Line plots
    ax4.plot(p_centers, profiles['pressure']['entropy_profile_clipped'],
            'g-', linewidth=2.5, label='Entropy (TS)')
    
    # Calculate -T·ΔS as entropy penalty
    entropy_contrib_p = -profiles['pressure']['entropy_profile_clipped']
    ax4.plot(p_centers, entropy_contrib_p, 'm--', linewidth=2, label='-T·ΔS')
    
    # Add uncertainty band for entropy if available
    if profiles['pressure']['entropy_uncertainty'] is not None:
        ax4.fill_between(
            p_centers,
            profiles['pressure']['entropy_profile_clipped'] - profiles['pressure']['entropy_uncertainty'],
            profiles['pressure']['entropy_profile_clipped'] + profiles['pressure']['entropy_uncertainty'],
            color='green', alpha=0.15)
    
    # Highlight key points
    ax4.scatter([max_S_pres], [profiles['pressure']['entropy_profile_clipped'][max_S_idx_p]], 
               s=80, color='green', marker='o', edgecolors='white', zorder=5,
               label='Max Entropy')
    
    # Add zero line
    ax4.axhline(y=0, color='k', linestyle='-', alpha=0.5)
    
    # Add target line
    ax4.axvline(x=target_pres, color='k', linestyle='--', alpha=0.7)
    
    # Label axes with consistent y-range
    ax4.set_xlabel('Pressure (atm)', fontsize=12)
    ax4.set_ylabel('Energy (kcal/mol)', fontsize=12)
    ax4.set_title('Pressure: Entropy Contribution', fontsize=14)
    ax4.set_ylim(-max_entropy_scale, max_entropy_scale)
    ax4.legend(loc='best', fontsize=10)
    
    # Add annotations for pressure effects
    # Regions of extreme pressure often have physical interpretations
    # Identify regions of low/high pressure
    p_mean = np.mean(p_centers)
    if np.min(p_centers) < -50:  # Significant negative pressure
        ax3.annotate("Negative Pressure\n(System Expansion)", 
                    (np.min(p_centers) * 0.9, 0.8 * max_energy_P),
                    xytext=(20, 0), textcoords='offset points', 
                    ha='left', va='center', fontsize=10,
                    arrowprops=dict(arrowstyle="->", color='black'),
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))
    
    if np.max(p_centers) > 50:  # Significant positive pressure
        ax3.annotate("High Pressure\n(System Compression)", 
                    (np.max(p_centers) * 0.9, 0.8 * max_energy_P),
                    xytext=(-20, 0), textcoords='offset points', 
                    ha='right', va='center', fontsize=10,
                    arrowprops=dict(arrowstyle="->", color='black'),
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))
    
    # Add enhanced info text with thermodynamic relationships and physical interpretation
    info_text = "Thermodynamic Relationships:\n"
    info_text += "• F = H - TS  (Helmholtz free energy)\n"
    info_text += "• H = U + PV  (Enthalpy)\n"
    info_text += "• F = -kT ln(P)  (Statistical mechanics)\n\n"
    info_text += "Physical Interpretation:\n"
    info_text += "• Regions with high TS (green shading) are entropically favored\n"
    info_text += "• Minimum in F indicates thermodynamically stable states\n"
    info_text += "• Entropy-enthalpy compensation occurs where F remains flat\n"
    info_text += "  despite changes in TS and H\n"
    info_text += f"• T = {profiles['avg_temp']:.1f}K, Pressure = {np.mean(p_centers):.1f} atm (average values)"
    
    # Create a separate info box with better positioning
    prop = dict(boxstyle='round', facecolor='wheat', alpha=0.4)
    fig.text(0.5, 0.01, info_text, fontsize=10, ha='center', va='bottom',
            bbox=prop)
    
    # Add title with entropy values
    title = 'Entropy Analysis from Simulation Data'
    if profiles['total_entropy']['temperature'] is not None:
        title += f' (Total S = {profiles["total_entropy"]["temperature"]:.4f} kcal/mol/K)'
    
    plt.suptitle(title, fontsize=16, y=0.98)
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])  # Adjust for info text
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nEnhanced entropy analysis saved as: {output_file}")
    plt.close()
    
    # Create an additional plot specifically for entropy-enthalpy compensation
    create_compensation_plot(profiles, output_file.replace('.png', '_compensation.png'))

def create_compensation_plot(profiles, output_file):
    """Create a specialized plot showing entropy-enthalpy compensation effects.
    
    This plot illustrates how entropy and enthalpy changes often compensate 
    each other in thermodynamic systems, resulting in smaller free energy changes.
    """
    plt.figure(figsize=(12, 10))
    
    # Set up grid for two subplots
    gs = plt.GridSpec(2, 1, height_ratios=[2, 1])
    
    # Create two axes
    ax1 = plt.subplot(gs[0])  # Main compensation plot
    ax2 = plt.subplot(gs[1])  # Correlation analysis
    
    # For temperature data
    t_centers = profiles['temperature']['centers']
    
    # Get data
    temp_h = profiles['temperature']['enthalpy_clipped']
    temp_ts = profiles['temperature']['entropy_profile_clipped']
    temp_f = profiles['temperature']['free_energy_clipped']
    
    # Sort by H values for better visualization
    sort_idx = np.argsort(temp_h)
    h_sorted = temp_h[sort_idx]
    ts_sorted = temp_ts[sort_idx]
    f_sorted = temp_f[sort_idx]
    t_sorted = t_centers[sort_idx]
    
    # Plot on main axes with colors based on temperature
    im = ax1.scatter(h_sorted, ts_sorted, 
                    c=t_sorted, cmap='coolwarm', 
                    s=80, alpha=0.7, 
                    edgecolors='black', linewidth=0.5,
                    label='Temperature Points')
    cbar = plt.colorbar(im, ax=ax1)
    cbar.set_label('Temperature (K)', fontsize=12)
    
    # Add linear regression line to show compensation effect
    slope, intercept, r_value, p_value, std_err = linregress(temp_h, temp_ts)
    x_range = np.linspace(min(temp_h), max(temp_h), 100)
    ax1.plot(x_range, slope * x_range + intercept, 'k--', 
             linewidth=2, label=f'Linear Fit (r = {r_value:.2f})')
    
    # Add the diagonal line representing perfect compensation
    max_val = max(max(temp_h), max(temp_ts))
    min_val = min(min(temp_h), min(temp_ts))
    range_val = max_val - min_val
    
    # Get points where H ≈ TS (free energy is minimized)
    close_points = np.where(np.abs(temp_h - temp_ts) < (range_val * 0.1))[0]
    if len(close_points) > 0:
        ax1.scatter(temp_h[close_points], temp_ts[close_points], 
                   s=100, facecolors='none', edgecolors='lime', linewidth=2,
                   label='H ≈ TS Points')
    
    # Add the perfect compensation line
    diagonal = np.linspace(min_val - range_val*0.1, max_val + range_val*0.1, 100)
    ax1.plot(diagonal, diagonal, 'r-', alpha=0.5, linewidth=2, 
             label='Perfect Compensation (H = TS)')
    
    # Add annotations showing degree of compensation
    compensation_ratio = slope
    ax1.annotate(f"Compensation Ratio: {compensation_ratio:.2f}", 
                xy=(0.05, 0.95), xycoords='axes fraction',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
                fontsize=12, ha='left', va='top')
    
    # Add another annotation explaining compensation phenomenon
    if abs(compensation_ratio) > 0.7:
        explanation = "Strong entropy-enthalpy compensation observed"
    elif abs(compensation_ratio) > 0.4:
        explanation = "Moderate entropy-enthalpy compensation observed"
    else:
        explanation = "Weak entropy-enthalpy compensation observed"
        
    ax1.annotate(explanation,
                xy=(0.05, 0.89), xycoords='axes fraction',
                bbox=dict(boxstyle="round,pad=0.3", fc="wheat", alpha=0.8),
                fontsize=12, ha='left', va='top')
    
    # Annotate points of interest
    max_h_idx = np.argmax(temp_h)
    max_ts_idx = np.argmax(temp_ts)
    min_f_idx = np.argmin(temp_f)
    
    ax1.annotate(f"T = {t_centers[max_h_idx]:.1f}K",
                xy=(temp_h[max_h_idx], temp_ts[max_h_idx]),
                xytext=(10, 10), textcoords='offset points',
                arrowprops=dict(arrowstyle="->", color='darkblue'),
                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7),
                fontsize=10)
                
    ax1.annotate(f"T = {t_centers[max_ts_idx]:.1f}K",
                xy=(temp_h[max_ts_idx], temp_ts[max_ts_idx]),
                xytext=(10, -10), textcoords='offset points',
                arrowprops=dict(arrowstyle="->", color='darkgreen'),
                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7),
                fontsize=10)
    
    # Set labels and title for main plot
    ax1.set_xlabel('Enthalpy (H) [kcal/mol]', fontsize=14)
    ax1.set_ylabel('Entropy (TS) [kcal/mol]', fontsize=14)
    ax1.set_title('Entropy-Enthalpy Compensation Analysis', fontsize=16)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper left', fontsize=10)
    
    # Second plot: Show correlation between ΔH and ΔTS
    # Calculate differences between consecutive temperature points
    delta_h = np.diff(temp_h[np.argsort(t_centers)])
    delta_ts = np.diff(temp_ts[np.argsort(t_centers)])
    delta_f = np.diff(temp_f[np.argsort(t_centers)])
    
    # Create bar chart showing the relationship
    bar_width = 0.35
    x = np.arange(len(delta_h))
    
    # Only show a reasonable number of points
    max_bars = 20
    if len(x) > max_bars:
        # Take evenly spaced samples
        indices = np.round(np.linspace(0, len(x)-1, max_bars)).astype(int)
        x = x[indices]
        delta_h = delta_h[indices]
        delta_ts = delta_ts[indices]
        delta_f = delta_f[indices]
    
    # Plot the bars
    ax2.bar(x - bar_width/2, delta_h, bar_width, color='red', alpha=0.6, label='ΔH')
    ax2.bar(x + bar_width/2, delta_ts, bar_width, color='green', alpha=0.6, label='ΔTS')
    ax2.plot(x, delta_f, 'bo-', linewidth=2, label='ΔF')
    
    # Add a horizontal line at y=0
    ax2.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    # Set labels
    ax2.set_xlabel('Temperature Index', fontsize=12)
    ax2.set_ylabel('Energy Changes (kcal/mol)', fontsize=12)
    ax2.set_title('Changes in H, TS, and F Between States', fontsize=14)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{i+1}" for i in range(len(x))], rotation=45)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Add text explaining the phenomenon
    explanation_text = (
        "Entropy-enthalpy compensation is a thermodynamic phenomenon where changes in the "
        "enthalpy (H) of a system are often counterbalanced by changes in the entropy term (TS), "
        "resulting in a smaller change in free energy (F) than would occur without compensation.\n\n"
        "When points fall near the diagonal line (H = TS), the system is in balance between "
        "enthalpy and entropy contributions. Points above the line are entropy-dominated, "
        "while points below are enthalpy-dominated.\n\n"
        f"For this system, the compensation ratio is {compensation_ratio:.2f}, "
        f"with correlation coefficient r = {r_value:.2f}. "
        f"The average free energy is {np.mean(temp_f):.2f} kcal/mol."
    )
    
    # Add text box with explanation
    plt.figtext(0.5, 0.01, explanation_text, fontsize=10, ha='center',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout(rect=[0, 0.07, 1, 0.98])
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Entropy-enthalpy compensation plot saved as: {output_file}")
    plt.close()
    return

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Generate free energy profiles from GROMACS simulation data')
    parser.add_argument('--input_dir', default='data', help='Directory containing input data files')
    parser.add_argument('--output_dir', default='.', help='Directory for output plots')
    parser.add_argument('--temp_bins', type=int, default=30, help='Number of temperature bins')
    parser.add_argument('--pres_bins', type=int, default=30, help='Number of pressure bins')
    parser.add_argument('--bootstrap', type=int, default=100, help='Number of bootstrap samples (0 to disable)')
    parser.add_argument('--smooth', type=float, default=1.0, help='Smoothing factor (0 to disable)')
    parser.add_argument('--clip_percentile', type=float, default=2.0, help='Percentile for clipping extreme values')
    parser.add_argument('--skip_plots', action='store_true', help='Skip generating plots (only calculate data)')
    args = parser.parse_args()
    
    # Set paths based on args
    data_dir = args.input_dir
    output_dir = args.output_dir
    
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Output file paths with proper directory joining
    profiles_plot = os.path.join(output_dir, 'free_energy_profiles.png')
    landscape_plot = os.path.join(output_dir, 'free_energy_2d.png')
    entropy_plot = os.path.join(output_dir, 'entropy_analysis.png')
    
    print("Reading energy data from GROMACS files...")
    data = run_gmx_energy(edr_file, tpr_file)
    
    # Print basic statistics
    print("\nRaw Data Statistics:")
    print(f"Temperature: mean = {np.mean(data['Temperature']):.2f} K, std = {np.std(data['Temperature']):.2f} K")
    print(f"Min: {np.min(data['Temperature']):.2f} K, Max: {np.max(data['Temperature']):.2f} K")
    print(f"Pressure: mean = {np.mean(data['Pressure']):.2f} atm, std = {np.std(data['Pressure']):.2f} atm")
    print(f"Min: {np.min(data['Pressure']):.2f} atm, Max: {np.max(data['Pressure']):.2f} atm")
    print(f"Total number of data points: {len(data['Temperature'])}")
    
    print("\nCalculating free energy profiles with entropy decomposition...")
    
    # Calculate free energy profiles with enhanced entropy analysis
    profiles = calculate_free_energy_profiles(
        data, 
        temp_bins=args.temp_bins, 
        pres_bins=args.pres_bins,
        bootstrap_samples=args.bootstrap,
        smoothing_factor=args.smooth,
        clip_percentile=args.clip_percentile
    )
    
    if args.skip_plots:
        print("\nSkipping plot generation as requested.")
        return
    
    print("\nCreating free energy profile plots...")
    # Create 1D free energy profile plots
    plot_free_energy_profiles(data, profiles, output_file=profiles_plot)
    
    print("\nCreating enhanced free energy landscape with entropy decomposition...")
    # Create 2D free energy landscape
    plot_combined_free_energy(data, profiles, output_file=landscape_plot)
    
    print("\nCreating detailed entropy analysis...")
    # Create entropy analysis plots
    plot_entropy_analysis(data, profiles, output_file=entropy_plot)
    
    print("\nDone!")

if __name__ == "__main__":
    main() 