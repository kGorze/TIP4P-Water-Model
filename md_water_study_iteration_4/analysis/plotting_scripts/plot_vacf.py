#!/usr/bin/python3
import os
import sys
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal, integrate
import warnings
warnings.filterwarnings('ignore')
try:
    import seaborn as sns
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.5)
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    print("Seaborn not available, using matplotlib defaults")

# Reference values for TIP4P water vibrational modes (cm^-1)
REFERENCE_MODES = {
    'libration': [400, 650],        # Librational modes
    'bending': [1600, 1700],        # H-O-H bending
    'stretching': [3400, 3600],     # O-H stretching
    'intermolecular': [50, 250]     # Intermolecular vibrations
}

def read_xvg(filename):
    """Read .xvg files and extract x, y data while skipping comment/label lines."""
    x = []
    y = []
    title = ""
    xlabel = ""
    ylabel = ""
    legend_labels = []
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                elif line.startswith('@'):
                    # Try to extract title and axis labels
                    if 'title' in line:
                        parts = line.split('"')
                        if len(parts) >= 2:
                            title = parts[1]
                    elif 'xaxis label' in line:
                        parts = line.split('"')
                        if len(parts) >= 2:
                            xlabel = parts[1]
                    elif 'yaxis label' in line:
                        parts = line.split('"')
                        if len(parts) >= 2:
                            ylabel = parts[1]
                    elif 's' in line and 'legend' in line:
                        parts = line.split('"')
                        if len(parts) >= 2:
                            legend_labels.append(parts[1])
                else:
                    values = line.strip().split()
                    if len(values) >= 2:
                        x.append(float(values[0]))
                        # If there are multiple y columns, store them all
                        if len(values) > 2:
                            y.append([float(val) for val in values[1:]])
                        else:
                            y.append(float(values[1]))
        
        # Convert to numpy arrays
        x = np.array(x)
        
        # If y is a list of lists, convert to a 2D array
        if y and isinstance(y[0], list):
            y = np.array(y)
        else:
            y = np.array(y)
            
        return x, y, title, xlabel, ylabel, legend_labels
    except Exception as e:
        print(f'Error reading {filename}: {e}')
        return np.array([]), np.array([]), "", "", "", []

def calculate_power_spectrum(t, vacf, zero_padding=True, window=True):
    """
    Calculate the power spectrum (vibrational density of states) from VACF
    
    Parameters:
    -----------
    t : array
        Time values in ps
    vacf : array
        Velocity autocorrelation function values
    zero_padding : bool
        Whether to apply zero padding to improve frequency resolution
    window : bool
        Whether to apply a window function to reduce spectral leakage
        
    Returns:
    --------
    freq : array
        Frequency values in cm^-1
    power : array
        Power spectrum values
    """
    # Ensure the data is evenly spaced
    dt = t[1] - t[0]  # Time step in ps
    
    # Apply window function to reduce spectral leakage
    if window:
        # Hann window is a good general-purpose window
        win = np.hanning(len(vacf))
        vacf_windowed = vacf * win
    else:
        vacf_windowed = vacf
    
    # Zero padding to improve frequency resolution
    if zero_padding:
        n_pad = len(vacf) * 4
        vacf_padded = np.zeros(n_pad)
        vacf_padded[:len(vacf)] = vacf_windowed
    else:
        vacf_padded = vacf_windowed
    
    # Compute the FFT
    fft_result = np.fft.rfft(vacf_padded)
    
    # Get the corresponding frequencies
    # Convert from angular frequency to frequency in THz
    # Then convert from THz to cm^-1 (1 THz = 33.356 cm^-1)
    freq_thz = np.fft.rfftfreq(len(vacf_padded), dt)  # in 1/ps = THz
    freq_cm1 = freq_thz * 33.356  # Convert to cm^-1
    
    # Power spectrum (magnitude squared of the FFT)
    power = np.abs(fft_result)**2
    
    # Normalize the power spectrum
    power = power / np.max(power)
    
    return freq_cm1, power

def plot_vacf(x, y, title, xlabel, ylabel, legend_labels, output_path):
    """Plot VACF data with spectral analysis"""
    # Create a figure with two subplots - VACF and power spectrum
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), dpi=300)
    
    # Check if y has multiple columns
    multi_column = len(y.shape) > 1 and y.shape[1] > 1
    
    # Plot VACF data
    if multi_column:
        # Plot each column
        for i in range(y.shape[1]):
            label = legend_labels[i] if i < len(legend_labels) else f"Series {i+1}"
            if HAS_SEABORN:
                sns.lineplot(x=x, y=y[:, i], label=label, ax=ax1)
            else:
                ax1.plot(x, y[:, i], label=label)
        
        # Use the first column for spectral analysis
        vacf_data = y[:, 0]
    else:
        # Single column plot
        if HAS_SEABORN:
            sns.lineplot(x=x, y=y, color='#1f77b4', linewidth=2, ax=ax1)
        else:
            ax1.plot(x, y, color='#1f77b4', linewidth=2)
        
        vacf_data = y
    
    # Calculate the decay time (time to reach 1/e of initial value)
    try:
        # Find where VACF crosses 1/e
        e_idx = np.where(vacf_data < vacf_data[0]/np.e)[0][0]
        decay_time = x[e_idx]
        
        # Add a marker and annotation for the decay time
        ax1.plot(decay_time, vacf_data[0]/np.e, 'ro', markersize=8)
        ax1.annotate(f'Decay time: {decay_time:.2f} ps',
                    xy=(decay_time, vacf_data[0]/np.e),
                    xytext=(decay_time + 0.1, vacf_data[0]/np.e + 0.1),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                    fontsize=10)
        
        # Add a horizontal line at 1/e
        ax1.axhline(y=vacf_data[0]/np.e, color='#d62728', linestyle='--', alpha=0.7,
                   label='1/e threshold')
    except:
        print("  - Could not determine VACF decay time")
    
    # Add zero line
    ax1.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    # Calculate the integral of VACF (related to diffusion coefficient)
    try:
        # Integrate using Simpson's rule
        integral = integrate.simps(vacf_data, x)
        
        # Add annotation for the integral
        ax1.text(0.05, 0.95, f'Integral: {integral:.3e}',
                transform=ax1.transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    except:
        print("  - Could not calculate VACF integral")
    
    ax1.set_xlabel(xlabel, fontsize=12)
    ax1.set_ylabel(ylabel, fontsize=12)
    ax1.set_title(f'{title}', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10, loc='best')
    
    # Calculate and plot the power spectrum (vibrational density of states)
    freq, power = calculate_power_spectrum(x, vacf_data)
    
    if HAS_SEABORN:
        sns.lineplot(x=freq, y=power, color='#1f77b4', linewidth=1.5, ax=ax2)
    else:
        ax2.plot(freq, y=power, color='#1f77b4', linewidth=1.5)
    
    # Add vertical spans for known vibrational modes
    for mode_name, (freq_min, freq_max) in REFERENCE_MODES.items():
        ax2.axvspan(freq_min, freq_max, alpha=0.2, 
                   label=f'{mode_name} ({freq_min}-{freq_max} cm⁻¹)')
    
    # Find and annotate major peaks
    peak_indices, _ = signal.find_peaks(power, height=0.1, distance=10)
    for idx in peak_indices[:5]:  # Annotate top 5 peaks
        if idx < len(freq):
            ax2.annotate(f'{freq[idx]:.0f} cm⁻¹',
                        xy=(freq[idx], power[idx]),
                        xytext=(freq[idx], power[idx] + 0.05),
                        arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                        fontsize=9, ha='center')
    
    ax2.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax2.set_ylabel('Power (normalized)', fontsize=12)
    ax2.set_title('Vibrational Density of States', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9, loc='best')
    
    # Limit the x-axis to the relevant frequency range
    ax2.set_xlim(0, 4000)
    
    # Add a watermark with simulation details
    fig.text(0.5, 0.01, 'TIP4P Water Model - Vibrational Analysis', 
            ha='center', fontsize=10, style='italic', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')
    
    # Create a detailed spectral analysis plot
    plt.figure(figsize=(12, 8), dpi=300)
    
    # Plot the power spectrum with more detail
    if HAS_SEABORN:
        sns.lineplot(x=freq, y=power, color='#1f77b4', linewidth=2)
    else:
        plt.plot(freq, power, color='#1f77b4', linewidth=2)
    
    # Add vertical spans for known vibrational modes with different colors
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    for i, (mode_name, (freq_min, freq_max)) in enumerate(REFERENCE_MODES.items()):
        plt.axvspan(freq_min, freq_max, alpha=0.2, color=colors[i % len(colors)], 
                   label=f'{mode_name} ({freq_min}-{freq_max} cm⁻¹)')
    
    # Find and annotate major peaks
    for idx in peak_indices:
        if idx < len(freq) and power[idx] > 0.2:  # Only annotate significant peaks
            plt.annotate(f'{freq[idx]:.0f} cm⁻¹',
                        xy=(freq[idx], power[idx]),
                        xytext=(freq[idx], power[idx] + 0.05),
                        arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                        fontsize=10, ha='center')
    
    # Add text box with information about vibrational modes
    info_text = (
        "Vibrational Modes of Water:\n"
        f"Intermolecular: {REFERENCE_MODES['intermolecular'][0]}-{REFERENCE_MODES['intermolecular'][1]} cm⁻¹\n"
        f"Libration: {REFERENCE_MODES['libration'][0]}-{REFERENCE_MODES['libration'][1]} cm⁻¹\n"
        f"Bending: {REFERENCE_MODES['bending'][0]}-{REFERENCE_MODES['bending'][1]} cm⁻¹\n"
        f"Stretching: {REFERENCE_MODES['stretching'][0]}-{REFERENCE_MODES['stretching'][1]} cm⁻¹"
    )
    
    plt.text(0.02, 0.98, info_text, transform=plt.gca().transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    plt.xlabel('Frequency (cm⁻¹)', fontsize=14)
    plt.ylabel('Power (normalized)', fontsize=14)
    plt.title('Detailed Vibrational Spectrum Analysis', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=10, loc='upper right')
    
    # Limit the x-axis to the relevant frequency range
    plt.xlim(0, 4000)
    
    plt.tight_layout()
    plt.savefig(os.path.join(os.path.dirname(output_path), 'vibrational_spectrum_plot.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print('  - vibrational_spectrum_plot.png saved successfully')

def main():
    if len(sys.argv) < 3:
        print("Usage: plot_vacf.py <analysis_dir> <plots_dir>")
        sys.exit(1)
        
    analysis_dir = sys.argv[1]
    plots_dir = sys.argv[2]
    
    # Ensure plots directory exists
    os.makedirs(plots_dir, exist_ok=True)
    
    # Plot VACF data
    vacf_file = os.path.join(analysis_dir, 'vacf.xvg')
    if os.path.exists(vacf_file):
        print('Plotting velocity autocorrelation function...')
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(vacf_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Velocity Autocorrelation Function'
            plot_xlabel = xlabel if xlabel else 'Time (ps)'
            plot_ylabel = ylabel if ylabel else 'VACF'
            
            output_path = os.path.join(plots_dir, 'vacf_plot.png')
            plot_vacf(x, y, plot_title, plot_xlabel, plot_ylabel, legend_labels, output_path)
    else:
        print(f"VACF file not found: {vacf_file}")

if __name__ == "__main__":
    main() 