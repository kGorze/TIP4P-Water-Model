#!/usr/bin/python3
import os
import sys
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
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

# Reference values for water
REFERENCE_VALUES = {
    'diffusion_coefficient_298K': 2.3e-9,  # m^2/s, self-diffusion coefficient at 298K
    'diffusion_coefficient_273K': 1.2e-9,  # m^2/s, self-diffusion coefficient at 273K (freezing point)
    'viscosity_298K': 0.896e-3,            # Pa·s, viscosity at 298K
    'viscosity_273K': 1.79e-3,             # Pa·s, viscosity at 273K (freezing point)
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

def calculate_diffusion_coefficient(time, msd, fit_start_ps=20, fit_end_ps=200, 
                                   min_fit_points=50, try_multiple_ranges=True,
                                   temperature=273.1, box_length=None):
    """
    Calculate the diffusion coefficient from MSD data using Einstein relation
    
    Parameters:
    -----------
    time : array
        Time values in ps
    msd : array
        MSD values in nm^2
    fit_start_ps : float
        Time (in ps) to start the linear fit (default: 20 ps)
    fit_end_ps : float
        Time (in ps) to end the linear fit (default: 200 ps)
    min_fit_points : int
        Minimum number of points required for fitting
    try_multiple_ranges : bool
        Whether to try multiple fitting ranges to find the best fit
    temperature : float
        Simulation temperature in K
    box_length : float
        Cubic box length in nm for finite-size correction
        
    Returns:
    --------
    D : float
        Diffusion coefficient in m^2/s
    D_corrected : float
        Finite-size corrected diffusion coefficient in m^2/s
    slope : float
        Slope of the linear fit in nm^2/ps
    intercept : float
        Intercept of the linear fit
    r_value : float
        Correlation coefficient
    p_value : float
        Two-sided p-value for a hypothesis test with null hypothesis that the slope is zero
    std_err : float
        Standard error of the slope estimate
    fit_start_idx : int
        Index where the fit starts
    fit_end_idx : int
        Index where the fit ends
    """
    # If we want to try multiple ranges to find the best fit
    if try_multiple_ranges:
        best_r2 = -1
        best_results = None
        
        # Try different fitting ranges
        start_ps_values = [10, 15, 20, 25, 30]
        end_ps_values = [150, 200, 250, 300, 400]
        
        for start_ps in start_ps_values:
            for end_ps in end_ps_values:
                if end_ps <= start_ps:
                    continue
                
                # Calculate with this range
                try:
                    results = calculate_diffusion_coefficient(
                        time, msd, 
                        fit_start_ps=start_ps, 
                        fit_end_ps=end_ps,
                        try_multiple_ranges=False,
                        temperature=temperature,
                        box_length=box_length
                    )
                    
                    # Unpack results
                    D, D_corrected, slope, intercept, r_value, p_value, std_err, fit_start_idx, fit_end_idx = results
                    
                    # Check if this is the best fit so far
                    if r_value**2 > best_r2:
                        best_r2 = r_value**2
                        best_results = results
                except Exception as e:
                    print(f"Error fitting range {start_ps}-{end_ps} ps: {e}")
                    continue
        
        # Return the best fit
        if best_results:
            return best_results
        else:
            # If no good fit was found, try with default values
            print("No good fit found with multiple ranges, using default range")
    
    # Find indices corresponding to the time range
    fit_start_idx = np.searchsorted(time, fit_start_ps)
    fit_end_idx = np.searchsorted(time, fit_end_ps)
    
    # Ensure we have enough points for fitting
    if fit_end_idx - fit_start_idx < min_fit_points:
        print(f"Warning: Not enough points in range {fit_start_ps}-{fit_end_ps} ps, adjusting range")
        fit_start_idx = max(0, len(time) // 5)
        fit_end_idx = min(len(time), len(time) * 4 // 5)
    
    # Extract the data for fitting
    fit_time = time[fit_start_idx:fit_end_idx]
    fit_msd = msd[fit_start_idx:fit_end_idx]
    
    # Perform linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(fit_time, fit_msd)
    
    # Calculate diffusion coefficient using Einstein relation: MSD = 6Dt
    # For 3D diffusion, D = slope / 6
    # Convert from nm^2/ps to m^2/s: 1 nm^2/ps = 1e-18 m^2 / 1e-12 s = 1e-6 m^2/s
    D = slope / 6.0 * 1e-6  # Correct conversion from nm²/ps to m²/s
    
    # Apply Yeh-Hummer finite-size correction if box length is provided
    D_corrected = D
    if box_length is not None:
        # Constants for the Yeh-Hummer correction
        kB = 1.380649e-23  # Boltzmann constant in J/K
        zeta = 2.837297  # Constant for cubic periodic boundary conditions
        
        # Get viscosity based on temperature
        if abs(temperature - 273.1) < 5:  # Close to 273K
            viscosity = REFERENCE_VALUES['viscosity_273K']
        else:  # Default to 298K
            viscosity = REFERENCE_VALUES['viscosity_298K']
        
        # Calculate the correction term
        # D_corrected = D - (kB * T * zeta) / (6 * pi * eta * L)
        correction = (kB * temperature * zeta) / (6 * np.pi * viscosity * box_length * 1e-9)  # Convert box_length from nm to m
        D_corrected = D - correction
    
    return D, D_corrected, slope, intercept, r_value, p_value, std_err, fit_start_idx, fit_end_idx

def plot_msd(x, y, title, xlabel, ylabel, legend_labels, output_path):
    """Plot MSD data with diffusion coefficient calculation"""
    # Create a figure with two subplots - linear and log-log
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7), dpi=300)
    
    # Check if y has multiple columns
    multi_column = len(y.shape) > 1 and y.shape[1] > 1
    
    # Plot data on linear scale
    if multi_column:
        # Plot each column
        for i in range(y.shape[1]):
            label = legend_labels[i] if i < len(legend_labels) else f"Series {i+1}"
            if HAS_SEABORN:
                sns.lineplot(x=x, y=y[:, i], label=label, ax=ax1)
            else:
                ax1.plot(x, y[:, i], label=label)
        
        # Use the first column for diffusion coefficient calculation
        msd_data = y[:, 0]
    else:
        # Single column plot
        if HAS_SEABORN:
            sns.lineplot(x=x, y=y, color='#1f77b4', linewidth=2, ax=ax1)
        else:
            ax1.plot(x, y, color='#1f77b4', linewidth=2)
        
        msd_data = y
    
    # Estimate box length from water density and number of molecules
    # For TIP4P water at 273K, density is ~1000 kg/m³
    # From the summary, we have 5500 water molecules
    num_water_molecules = 5500
    water_density_molecules_per_nm3 = 33.59  # From the summary
    
    # Calculate box volume and length
    box_volume_nm3 = num_water_molecules / water_density_molecules_per_nm3
    box_length_nm = box_volume_nm3 ** (1/3)
    
    # Get the simulation temperature
    temperature = 273.1  # From the summary
    
    # Calculate diffusion coefficient with improved fitting
    # Use a range that excludes the ballistic regime (typically first 10-20 ps)
    D, D_corrected, slope, intercept, r_value, p_value, std_err, fit_start_idx, fit_end_idx = calculate_diffusion_coefficient(
        x, msd_data, fit_start_ps=20, fit_end_ps=200, try_multiple_ranges=True, 
        min_fit_points=50, temperature=temperature, box_length=box_length_nm
    )
    
    # Plot the linear fit
    fit_x = np.array([x[fit_start_idx], x[fit_end_idx]])
    fit_y = slope * fit_x + intercept
    ax1.plot(fit_x, fit_y, 'r--', linewidth=2, 
            label=f'Linear fit (R²={r_value**2:.3f})')
    
    # Get the appropriate reference diffusion coefficient based on temperature
    if abs(temperature - 273.1) < 5:  # Close to 273K
        reference_D = REFERENCE_VALUES['diffusion_coefficient_273K']
        temp_label = "273K"
    else:  # Default to 298K
        reference_D = REFERENCE_VALUES['diffusion_coefficient_298K']
        temp_label = "298K"
    
    # Add annotation for diffusion coefficient
    ax1.text(0.05, 0.95, 
            f'D = {D:.3e} m²/s\n'
            f'D (corrected) = {D_corrected:.3e} m²/s\n'
            f'Slope = {slope:.3f} nm²/ps\n'
            f'Fit range: {x[fit_start_idx]:.0f}-{x[fit_end_idx]:.0f} ps\n'
            f'Reference D ({temp_label}) = {reference_D:.3e} m²/s\n'
            f'Difference: {((D - reference_D) / reference_D * 100):.1f}%',
            transform=ax1.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    ax1.set_xlabel(xlabel, fontsize=12)
    ax1.set_ylabel(ylabel, fontsize=12)
    ax1.set_title(f'{title} (Linear Scale)', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10, loc='lower right')
    
    # Plot data on log-log scale
    if multi_column:
        for i in range(y.shape[1]):
            label = legend_labels[i] if i < len(legend_labels) else f"Series {i+1}"
            ax2.loglog(x, y[:, i], label=label)
    else:
        ax2.loglog(x, y, color='#1f77b4', linewidth=2)
    
    # Add reference lines for different regimes
    # Ballistic regime: MSD ~ t^2
    t_ballistic = np.logspace(np.log10(x[1]), np.log10(x[-1]), 100)
    ballistic_factor = y[5] / (x[5]**2)  # Adjust factor to match data
    msd_ballistic = ballistic_factor * t_ballistic**2
    ax2.loglog(t_ballistic, msd_ballistic, 'g--', linewidth=1.5, alpha=0.7,
              label='Ballistic (t²)')
    
    # Diffusive regime: MSD ~ t
    t_diffusive = np.logspace(np.log10(x[len(x)//4]), np.log10(x[-1]), 100)
    diffusive_factor = slope  # Use the calculated slope
    msd_diffusive = diffusive_factor * t_diffusive
    ax2.loglog(t_diffusive, msd_diffusive, 'r--', linewidth=1.5, alpha=0.7,
              label='Diffusive (t¹)')
    
    # Add annotations for different regimes
    ax2.text(x[3], y[3]*1.5, 'Ballistic\nRegime', fontsize=10, 
            ha='center', va='center', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    ax2.text(x[len(x)//2], y[len(x)//2]*0.5, 'Diffusive\nRegime', fontsize=10, 
            ha='center', va='center', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    ax2.set_xlabel(xlabel, fontsize=12)
    ax2.set_ylabel(ylabel, fontsize=12)
    ax2.set_title(f'{title} (Log-Log Scale)', fontsize=14)
    ax2.grid(True, alpha=0.3, which='both')
    ax2.legend(fontsize=10, loc='lower right')
    
    # Add a watermark with simulation details
    fig.text(0.5, 0.01, 'TIP4P Water Model - Self-Diffusion Analysis', 
            ha='center', fontsize=10, style='italic', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')
    
    # Create a detailed diffusion analysis plot with multiple fitting ranges
    plt.figure(figsize=(12, 10), dpi=300)
    
    # Plot MSD data
    if HAS_SEABORN:
        sns.lineplot(x=x, y=msd_data, color='#1f77b4', linewidth=2, label='MSD')
    else:
        plt.plot(x, msd_data, color='#1f77b4', linewidth=2, label='MSD')
    
    # Plot the linear fit for the entire fitting region
    fit_x_full = x[fit_start_idx:fit_end_idx]
    fit_y_full = slope * fit_x_full + intercept
    plt.plot(fit_x_full, fit_y_full, 'r-', linewidth=2, 
            label=f'Best fit (R²={r_value**2:.3f})')
    
    # Highlight the fitting region
    plt.axvspan(x[fit_start_idx], x[fit_end_idx], alpha=0.2, color='green',
               label=f'Best fitting region ({x[fit_start_idx]:.0f}-{x[fit_end_idx]:.0f} ps)')
    
    # Try alternative fitting ranges for comparison
    alt_ranges = [
        (10, 100, 'Early regime (10-100 ps)'),
        (50, 300, 'Middle regime (50-300 ps)'),
        (100, 500, 'Late regime (100-500 ps)')
    ]
    
    alt_colors = ['#ff7f0e', '#2ca02c', '#d62728']
    alt_results = []
    
    for i, (start_ps, end_ps, label) in enumerate(alt_ranges):
        # Find indices corresponding to the time range
        alt_start_idx = np.searchsorted(x, start_ps)
        alt_end_idx = np.searchsorted(x, end_ps)
        
        # Skip if range is too small or out of bounds
        if alt_end_idx >= len(x) or alt_start_idx >= alt_end_idx or alt_end_idx - alt_start_idx < 20:
            continue
            
        alt_time = x[alt_start_idx:alt_end_idx]
        alt_msd = msd_data[alt_start_idx:alt_end_idx]
        
        # Perform linear regression
        alt_slope, alt_intercept, alt_r_value, alt_p_value, alt_std_err = stats.linregress(alt_time, alt_msd)
        
        # Calculate diffusion coefficient
        alt_D = alt_slope / 6.0 * 1e-6  # Convert to m²/s
        
        # Apply Yeh-Hummer correction
        alt_D_corrected = alt_D
        if box_length_nm is not None:
            # Get viscosity based on temperature
            if abs(temperature - 273.1) < 5:  # Close to 273K
                viscosity = REFERENCE_VALUES['viscosity_273K']
            else:  # Default to 298K
                viscosity = REFERENCE_VALUES['viscosity_298K']
            
            # Calculate the correction term
            kB = 1.380649e-23  # Boltzmann constant in J/K
            zeta = 2.837297  # Constant for cubic periodic boundary conditions
            correction = (kB * temperature * zeta) / (6 * np.pi * viscosity * box_length_nm * 1e-9)
            alt_D_corrected = alt_D - correction
        
        # Store results
        alt_results.append((alt_D, alt_D_corrected, alt_slope, alt_r_value**2, alt_start_idx, alt_end_idx, label))
        
        # Plot alternative fit
        alt_fit_y = alt_slope * alt_time + alt_intercept
        plt.plot(alt_time, alt_fit_y, '--', color=alt_colors[i], linewidth=1.5, alpha=0.7,
                label=f'{label} fit: D={alt_D_corrected:.3e} m²/s (R²={alt_r_value**2:.3f})')
    
    # Calculate the derivative of MSD (numerical)
    if len(x) > 10:
        # Use a window for smoother derivative
        window = min(20, len(x) // 10)
        x_deriv = []
        msd_deriv = []
        
        for i in range(window, len(x) - window):
            x_deriv.append(x[i])
            # Calculate derivative over the window
            dx = x[i + window] - x[i - window]
            dy = msd_data[i + window] - msd_data[i - window]
            msd_deriv.append(dy / dx)
        
        # Plot the derivative divided by 6 (should approach D for diffusive regime)
        deriv_array = np.array(msd_deriv) / 6.0
        plt.plot(x_deriv, deriv_array, 'g--', linewidth=1.5, alpha=0.7,
                label='dMSD/dt ÷ 6')
        
        # Add horizontal line for calculated D (in nm²/ps)
        D_nm2_ps = D * 1e-3  # Convert back to nm²/ps
        plt.axhline(y=D_nm2_ps, color='#d62728', linestyle='--', alpha=0.7,
                   label=f'D = {D_nm2_ps:.3e} nm²/ps')
    
    # Add detailed information
    info_text = (
        f"Diffusion Coefficient Analysis (T = {temperature:.1f} K):\n\n"
        f"Raw D = {D:.3e} m²/s\n"
        f"Finite-size corrected D = {D_corrected:.3e} m²/s\n"
        f"Slope = {slope:.3f} nm²/ps\n"
        f"Fit range: {x[fit_start_idx]:.0f}-{x[fit_end_idx]:.0f} ps\n"
        f"R² = {r_value**2:.3f}\n"
        f"Box length = {box_length_nm:.2f} nm\n\n"
        f"Reference D ({temp_label}) = {reference_D:.3e} m²/s\n"
        f"Difference (raw): {((D - reference_D) / reference_D * 100):.1f}%\n"
        f"Difference (corrected): {((D_corrected - reference_D) / reference_D * 100):.1f}%\n\n"
        f"Note: At 273K, water diffusion is ~50% slower than at 298K\n"
        f"Yeh-Hummer correction applied for finite-size effects\n"
    )
    
    # Add alternative fitting results
    if alt_results:
        info_text += "\nAlternative Fitting Ranges:\n"
        for alt_D, alt_D_corrected, alt_slope, alt_r2, alt_start_idx, alt_end_idx, label in alt_results:
            info_text += f"{label}:\n"
            info_text += f"  Raw D = {alt_D:.3e} m²/s\n"
            info_text += f"  Corrected D = {alt_D_corrected:.3e} m²/s\n"
            info_text += f"  R² = {alt_r2:.3f}\n"
    
    # Add explanation about temperature effects and finite-size correction
    explanation_text = (
        "Temperature Effects on Diffusion:\n"
        "• At 273K (freezing), water diffusion is ~1.0-1.3×10⁻⁹ m²/s\n"
        "• At 298K (room temp), water diffusion is ~2.3×10⁻⁹ m²/s\n"
        "• TIP4P models often overestimate diffusion by 10-40%\n\n"
        "Finite-Size Correction:\n"
        "• Periodic boundary conditions artificially enhance diffusion\n"
        "• Yeh-Hummer correction: D∞ = D(L) - kBTζ/(6πηL)\n"
        "• Correction is larger at lower temperatures (higher viscosity)\n\n"
        "Fitting Considerations:\n"
        "• Exclude ballistic regime (first ~10-20 ps)\n"
        "• Use clear diffusive regime (typically 20-200 ps)\n"
        "• R² value should be >0.99 for reliable fits"
    )
    
    # Position the text boxes
    plt.text(0.02, 0.98, info_text, transform=plt.gca().transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    plt.text(0.02, 0.40, explanation_text, transform=plt.gca().transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
    
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title('Detailed Diffusion Coefficient Analysis', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=10, loc='lower right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(os.path.dirname(output_path), 'diffusion_analysis_plot.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print('  - diffusion_analysis_plot.png saved successfully')

def main():
    if len(sys.argv) < 3:
        print("Usage: plot_msd.py <analysis_dir> <plots_dir>")
        sys.exit(1)
        
    analysis_dir = sys.argv[1]
    plots_dir = sys.argv[2]
    
    # Ensure plots directory exists
    os.makedirs(plots_dir, exist_ok=True)
    
    # Plot MSD data
    msd_file = os.path.join(analysis_dir, 'msd.xvg')
    if os.path.exists(msd_file):
        print('Plotting mean square displacement...')
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(msd_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Mean Square Displacement vs Time'
            plot_xlabel = xlabel if xlabel else 'Time (ps)'
            plot_ylabel = ylabel if ylabel else 'MSD (nm²)'
            
            output_path = os.path.join(plots_dir, 'msd_plot.png')
            plot_msd(x, y, plot_title, plot_xlabel, plot_ylabel, legend_labels, output_path)
    else:
        print(f"MSD file not found: {msd_file}")

if __name__ == "__main__":
    main() 