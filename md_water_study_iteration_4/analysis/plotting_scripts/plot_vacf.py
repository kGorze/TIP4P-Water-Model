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
    
    # Add zero line
    ax1.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    # Calculate the integral of VACF (related to diffusion coefficient)
    try:
        # Integrate using Simpson's rule
        integral = integrate.simps(vacf_data, x)
        
        # Add annotation for the integral with better context
        # Position it to avoid overlap with the plot line
        ax1.text(0.5, 0.95, f'VACF Integral: {integral:.3e}\n(Related to diffusion coefficient)',
                transform=ax1.transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    except:
        print("  - Could not calculate VACF integral")
    
    # Calculate the decay time (time to reach 1/e of initial value)
    try:
        # Find where VACF crosses 1/e
        e_idx = np.where(vacf_data < vacf_data[0]/np.e)[0][0]
        decay_time = x[e_idx]
        
        # Add a marker and annotation for the decay time (standardized style)
        ax1.plot(decay_time, vacf_data[0]/np.e, 'ro', markersize=10)  # Increased marker size
        ax1.annotate(f'Decay time: {decay_time:.2f} ps',
                    xy=(decay_time, vacf_data[0]/np.e),
                    xytext=(decay_time + 0.1, vacf_data[0]/np.e + 0.1),
                    arrowprops=dict(facecolor='#d62728', shrink=0.05, width=1.5, headwidth=8),
                    fontsize=10, bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        
        # Add a horizontal line at 1/e with explanation
        ax1.axhline(y=vacf_data[0]/np.e, color='#d62728', linestyle='--', alpha=0.7,
                   label='1/e threshold')
        
        # Add explanation for 1/e threshold
        ax1.text(0.98, vacf_data[0]/np.e + 0.05, 
                '1/e threshold: time scale of molecular motion decorrelation',
                ha='right', fontsize=9, bbox=dict(boxstyle='round', facecolor='#ffedee', alpha=0.7))
    except:
        print("  - Could not determine VACF decay time")
    
    # Add annotation for the initial spike
    if len(vacf_data) > 0:
        ax1.annotate('Initial correlation\n(t=0 self-correlation)',
                    xy=(0, vacf_data[0]),
                    xytext=(1, vacf_data[0] - 0.3),
                    arrowprops=dict(facecolor='#1f77b4', shrink=0.05, width=1.5, headwidth=8),
                    fontsize=9, bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    ax1.set_xlabel(xlabel, fontsize=12)
    ax1.set_ylabel(ylabel, fontsize=12)
    ax1.set_title(f'{title}', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10, loc='best')
    
    # Calculate and plot the power spectrum (vibrational density of states)
    freq, power = calculate_power_spectrum(x, vacf_data)
    
    # Determine if we're looking at low-frequency data (likely from VACF spectrum file)
    is_low_freq = np.max(freq) < 100  # If max frequency is less than 100 cm⁻¹
    
    # Find the maximum value in the power spectrum to ensure proper y-axis scaling
    max_power = np.max(power) * 1.1  # Add 10% margin
    
    if HAS_SEABORN:
        sns.lineplot(x=freq, y=power, color='#1f77b4', linewidth=1.5, ax=ax2)
    else:
        ax2.plot(freq, y=power, color='#1f77b4', linewidth=1.5)
    
    # Only add vibrational mode spans if we're looking at the full spectrum
    if not is_low_freq:
        # Add vertical spans for known vibrational modes
        for mode_name, (freq_min, freq_max) in REFERENCE_MODES.items():
            ax2.axvspan(freq_min, freq_max, alpha=0.2, 
                       label=f'{mode_name} ({freq_min}-{freq_max} cm⁻¹)')
    
    # Find and annotate all significant peaks
    peak_indices, _ = signal.find_peaks(power, height=0.1, distance=10)
    for idx in peak_indices:  # Annotate all significant peaks
        if idx < len(freq) and power[idx] > 0.1:  # Only annotate peaks above 10% of max
            ax2.annotate(f'{freq[idx]:.1f} cm⁻¹',
                        xy=(freq[idx], power[idx]),
                        xytext=(freq[idx], power[idx] + 0.05),
                        arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                        fontsize=9, ha='center')
    
    # If there's a rising feature at the edge, annotate it
    if is_low_freq and len(freq) > 0:
        right_edge_idx = len(freq) - 1
        if power[right_edge_idx] > 0.5:  # If the edge value is significant
            ax2.annotate('Librational modes\n(begin at ~50 cm⁻¹)',
                        xy=(freq[right_edge_idx], power[right_edge_idx]),
                        xytext=(freq[right_edge_idx] - 3, power[right_edge_idx] - 0.3),
                        arrowprops=dict(facecolor='#1f77b4', shrink=0.05, width=1.5, headwidth=8),
                        fontsize=9, bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    ax2.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
    ax2.set_ylabel('Relative Intensity (normalized)', fontsize=12)  # More specific y-axis label
    
    # Set appropriate title based on frequency range
    if is_low_freq:
        ax2.set_title('Low-Frequency Vibrational Spectrum (Collective Modes)', fontsize=14)
    else:
        ax2.set_title('Vibrational Density of States (Full Spectrum)', fontsize=14)
    
    ax2.grid(True, alpha=0.3)
    
    # Place legend outside the plot area if we have many items
    if not is_low_freq:
        ax2.legend(fontsize=9, loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0)
    else:
        ax2.legend(fontsize=9, loc='best')
    
    # Set appropriate x-axis limits based on the data
    if is_low_freq:
        # Find the maximum x value from the top plot to match axes
        max_x_top = np.max(x)
        # Convert time (ps) to frequency (cm⁻¹) approximately
        # This is a rough conversion to match scales
        max_freq = max_x_top
        ax2.set_xlim(0, max_freq)  # Match the x-axis range with the top plot
    else:
        ax2.set_xlim(0, 4000)  # Full vibrational spectrum range
    
    # Ensure y-axis captures all peaks with margin
    ax2.set_ylim(0, max_power)
    
    # Add a connecting explanation between the plots
    fig.text(0.02, 0.5, 'Fourier\nTransform', fontsize=12, rotation=90, 
             ha='center', va='center', bbox=dict(boxstyle='rarrow', facecolor='lightgray', alpha=0.7))
    
    # Add a watermark with simulation details and explanation
    if is_low_freq:
        explanation = 'VACF (top) → Fourier transform → Vibrational spectrum (bottom)'
    else:
        explanation = 'Full spectrum shows all vibrational modes of water molecules'
    
    fig.text(0.5, 0.01, f'TIP4P Water Model | Vibrational Analysis | {explanation}', 
            ha='center', fontsize=10, style='italic', alpha=0.7)
    
    plt.tight_layout(rect=[0.03, 0.03, 1, 0.97])  # Adjust for the watermark and arrow
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')
    
    # Create a detailed spectral analysis plot with improved visualization
    if not is_low_freq:  # Only create the detailed plot for the full spectrum
        fig = plt.figure(figsize=(12, 10), dpi=300)
        
        # Create a figure with two subplots - one for low frequencies, one for full spectrum
        gs = plt.GridSpec(2, 1, height_ratios=[1, 2], figure=fig)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        
        # Plot the low-frequency part (0-50 cm⁻¹) in the top subplot
        low_freq_mask = freq <= 50
        if np.any(low_freq_mask):
            low_freq = freq[low_freq_mask]
            low_power = power[low_freq_mask]
            
            ax1.plot(low_freq, low_power, color='#1f77b4', linewidth=2.5, label='Low-Frequency Modes')
            
            # Add vertical span for intermolecular vibrations in low-frequency plot
            freq_min, freq_max = REFERENCE_MODES['intermolecular']
            ax1.axvspan(freq_min, freq_max, alpha=0.3, color='#ffb6c1', 
                       label=f'Intermolecular ({freq_min}-{freq_max} cm⁻¹)')
            
            # Annotate any significant peaks in the low-frequency region
            low_peak_indices, _ = signal.find_peaks(low_power, height=0.1, distance=5)
            for idx in low_peak_indices:  # Annotate all significant peaks
                if idx < len(low_freq) and low_power[idx] > 0.1:
                    ax1.annotate(f'{low_freq[idx]:.1f} cm⁻¹',
                                xy=(low_freq[idx], low_power[idx]),
                                xytext=(low_freq[idx], low_power[idx] + 0.05),
                                arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                                fontsize=9, ha='center')
            
            # Add context for the integral
            if 'integral' in locals():
                ax1.text(0.05, 0.95, f'VACF Integral: {integral:.3e}\n(Related to molecular mobility)',
                        transform=ax1.transAxes, fontsize=10, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
            
            ax1.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
            ax1.set_ylabel('Relative Intensity', fontsize=12)
            ax1.set_title('Low-Frequency Collective Modes', fontsize=14)
            ax1.grid(True, linestyle='--', alpha=0.4)
            ax1.legend(fontsize=10, loc='upper right')
            ax1.set_xlim(0, 50)
            
            # Ensure y-axis captures all peaks with margin
            ax1.set_ylim(0, np.max(low_power) * 1.1)
        
        # Plot the full spectrum in the bottom subplot with distinct colored bands
        ax2.plot(freq, power, color='#1f77b4', linewidth=2.5, label='Vibrational Spectrum')
        
        # Add vertical spans for known vibrational modes with different colors
        colors = ['#ffb6c1', '#add8e6', '#ffebcd', '#d8f0d8']  # Distinct colors for each mode
        mode_order = ['intermolecular', 'libration', 'bending', 'stretching']  # Consistent order
        
        # Add the colored bands with appropriate opacity
        for i, mode_name in enumerate(mode_order):
            freq_min, freq_max = REFERENCE_MODES[mode_name]
            ax2.axvspan(freq_min, freq_max, alpha=0.3, color=colors[i], 
                       label=f'{mode_name} ({freq_min}-{freq_max} cm⁻¹)')
        
        # Find and annotate major peaks in the full spectrum
        for idx in peak_indices:
            if idx < len(freq) and power[idx] > 0.2:  # Only annotate significant peaks
                ax2.annotate(f'{freq[idx]:.0f} cm⁻¹',
                            xy=(freq[idx], power[idx]),
                            xytext=(freq[idx], power[idx] + 0.05),
                            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                            fontsize=10, ha='center')
        
        ax2.set_xlabel('Frequency (cm⁻¹)', fontsize=14)
        ax2.set_ylabel('Relative Intensity', fontsize=14)
        ax2.set_title('Full Vibrational Spectrum', fontsize=14)
        ax2.grid(True, linestyle='--', alpha=0.4)
        
        # Place legend outside the plot area to avoid overlap
        ax2.legend(fontsize=10, loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0)
        
        # Set appropriate x-axis limits
        ax2.set_xlim(0, 4000)
        
        # Add a title for the entire figure
        fig.suptitle('Vibrational Spectrum Analysis of TIP4P Water', fontsize=16, fontweight='bold', y=0.98)
        
        # Add simulation conditions and explanation as a footer
        sim_details = (
            f"TIP4P Water Model | T=298K | P=1 bar | 5500 molecules | "
            f"Top: Low-frequency collective modes | Bottom: Full molecular vibrational spectrum"
        )
        fig.text(0.5, 0.01, sim_details, ha='center', fontsize=10, 
               style='italic', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))
        
        plt.tight_layout(rect=[0, 0.05, 0.85, 0.95])  # Adjust for the title, legend, and footer
        plt.savefig(os.path.join(os.path.dirname(output_path), 'vibrational_spectrum_plot.png'), dpi=300, bbox_inches='tight')
        plt.close()
        print('  - vibrational_spectrum_plot.png saved successfully')

def main():
    if len(sys.argv) < 3:
        print("Usage: plot_vacf.py <analysis_dir> <plots_dir>")
        sys.exit(1)
        
    analysis_dir = sys.argv[1]
    plots_dir = sys.argv[2]
    
    # Define data directory
    data_dir = os.path.join(analysis_dir, "data")
    
    # Ensure plots directory exists
    os.makedirs(plots_dir, exist_ok=True)
    
    # VACF plot
    vacf_file = os.path.join(data_dir, 'vacf.xvg')
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
    
    # VACF spectrum plot
    vacf_spectrum_file = os.path.join(data_dir, 'vacf_spectrum.xvg')
    if os.path.exists(vacf_spectrum_file):
        print('Plotting vibrational density of states...')
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(vacf_spectrum_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Vibrational Density of States'
            plot_xlabel = xlabel if xlabel else 'Frequency (cm⁻¹)'
            plot_ylabel = ylabel if ylabel else 'Power (normalized)'
            
            output_path = os.path.join(plots_dir, 'vacf_spectrum_plot.png')
            plot_vacf(x, y, plot_title, plot_xlabel, plot_ylabel, legend_labels, output_path)
    else:
        print(f"VACF spectrum file not found: {vacf_spectrum_file}")

if __name__ == "__main__":
    main() 