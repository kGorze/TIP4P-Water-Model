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

def plot_msd_enhanced(x, y, title, xlabel, ylabel, output_path, fit_results=None):
    """Create an enhanced MSD plot with annotations, linear fit, and styling"""
    # Create figure with enhanced styling
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    
    # Set a consistent style
    ax.grid(True, linestyle='--', alpha=0.7)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    
    # Plot MSD data with enhanced styling
    if HAS_SEABORN:
        sns.lineplot(x=x, y=y, color='#1f77b4', linewidth=2.5, label='MSD')
    else:
        ax.plot(x, y, color='#1f77b4', linewidth=2.5, label='MSD')
    
    # Add linear fit if available
    if fit_results:
        start_idx, end_idx, slope, intercept, r_squared, diffusion_coef = fit_results
        
        # Plot the linear fit line
        fit_x = x[start_idx:end_idx+1]
        fit_y = slope * fit_x + intercept
        ax.plot(fit_x, fit_y, color='#d62728', linestyle='--', linewidth=2.5,
               label=f'Linear fit (r² = {r_squared:.3f})')
        
        # Highlight the linear region used for fitting
        ax.fill_between(fit_x, fit_y - 0.1*max(y), fit_y + 0.1*max(y), 
                       color='#d62728', alpha=0.1)
        
        # Add annotation for diffusion coefficient
        ax.text(0.05, 0.95, f"D = {diffusion_coef:.3e} cm²/s", transform=ax.transAxes,
               fontsize=12, verticalalignment='top',
               bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
        
        # Add annotation for the linear region
        ax.annotate('Linear region used for D calculation', 
                   xy=(fit_x[len(fit_x)//2], fit_y[len(fit_x)//2]),
                   xytext=(fit_x[len(fit_x)//2], fit_y[len(fit_x)//2] + 0.2*max(y)),
                   arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
                   fontsize=10, ha='center')
    
    # Add reference line for theoretical behavior (MSD ~ t)
    if len(x) > 1 and len(y) > 1:
        # Create a reference line with slope 1 (in log-log scale)
        ref_x = np.array([x[1], x[-1]])
        ref_y = np.array([y[1], y[1] * (x[-1]/x[1])])
        ax.plot(ref_x, ref_y, color='#2ca02c', linestyle=':', linewidth=2,
               label='Theoretical (MSD ~ t)')
    
    # Add title and labels with enhanced styling
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    
    # Add grid and improve styling
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    # Add explanation text
    explanation = (
        "Mean Square Displacement (MSD) measures how far particles move over time.\n"
        "The slope of the linear region is proportional to the diffusion coefficient (D).\n"
        "For normal diffusion: MSD = 6Dt (3D) or MSD = 4Dt (2D) or MSD = 2Dt (1D)."
    )
    plt.figtext(0.5, 0.01, explanation, ha='center', fontsize=9,
               bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
    
    # Add legend with enhanced styling
    ax.legend(fontsize=11, framealpha=0.8, loc='upper left')
    
    # Save figure with tight layout
    plt.tight_layout(rect=[0, 0.08, 1, 1])  # Adjust for the explanation text
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def create_diffusion_analysis_plot(x, y, fit_results, output_path):
    """Create an enhanced diffusion analysis plot with multiple visualizations"""
    if fit_results is None:
        print("  - No fit results available for diffusion analysis plot")
        return
    
    # Extract fit results
    start_idx, end_idx, slope, intercept, r_squared, diffusion_coef = fit_results
    
    # Create a figure with 2x2 subplots
    fig, axs = plt.subplots(2, 2, figsize=(15, 12), dpi=300)
    
    # Set a consistent style for all subplots
    for ax in axs.flatten():
        ax.grid(True, linestyle='--', alpha=0.7)
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
    
    # 1. MSD vs Time (linear scale)
    axs[0, 0].plot(x, y, color='#1f77b4', linewidth=2.5, label='MSD')
    
    # Add linear fit
    fit_x = x[start_idx:end_idx+1]
    fit_y = slope * fit_x + intercept
    axs[0, 0].plot(fit_x, fit_y, color='#d62728', linestyle='--', linewidth=2.5,
                  label=f'Linear fit (r² = {r_squared:.3f})')
    
    # Highlight the linear region
    axs[0, 0].fill_between(fit_x, fit_y - 0.1*max(y), fit_y + 0.1*max(y), 
                          color='#d62728', alpha=0.1)
    
    axs[0, 0].set_title('MSD vs Time (Linear Scale)', fontsize=16, fontweight='bold')
    axs[0, 0].set_xlabel('Time (ps)', fontsize=14)
    axs[0, 0].set_ylabel('MSD (nm²)', fontsize=14)
    axs[0, 0].legend(fontsize=11, framealpha=0.8, loc='upper left')
    axs[0, 0].tick_params(axis='both', which='major', labelsize=12)
    
    # 2. MSD vs Time (log-log scale)
    axs[0, 1].loglog(x, y, color='#1f77b4', linewidth=2.5, label='MSD')
    
    # Add linear fit in log-log scale
    axs[0, 1].loglog(fit_x, fit_y, color='#d62728', linestyle='--', linewidth=2.5,
                    label=f'Linear fit (r² = {r_squared:.3f})')
    
    # Add reference lines for different diffusion regimes
    if len(x) > 1 and len(y) > 1:
        # Create reference lines with different slopes
        ref_x = np.array([x[1], x[-1]])
        
        # Normal diffusion (slope = 1)
        ref_y1 = np.array([y[1], y[1] * (x[-1]/x[1])])
        axs[0, 1].loglog(ref_x, ref_y1, color='#2ca02c', linestyle=':', linewidth=2,
                        label='Normal (α = 1)')
        
        # Sub-diffusion (slope = 0.5)
        ref_y2 = np.array([y[1], y[1] * np.sqrt(x[-1]/x[1])])
        axs[0, 1].loglog(ref_x, ref_y2, color='#ff7f0e', linestyle=':', linewidth=2,
                        label='Sub-diffusion (α = 0.5)')
        
        # Super-diffusion (slope = 1.5)
        ref_y3 = np.array([y[1], y[1] * (x[-1]/x[1])**1.5])
        axs[0, 1].loglog(ref_x, ref_y3, color='#9467bd', linestyle=':', linewidth=2,
                        label='Super-diffusion (α = 1.5)')
    
    axs[0, 1].set_title('MSD vs Time (Log-Log Scale)', fontsize=16, fontweight='bold')
    axs[0, 1].set_xlabel('Time (ps)', fontsize=14)
    axs[0, 1].set_ylabel('MSD (nm²)', fontsize=14)
    axs[0, 1].legend(fontsize=11, framealpha=0.8, loc='upper left')
    axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
    
    # 3. MSD/t vs Time (to check for normal diffusion)
    msd_over_t = np.zeros_like(y)
    msd_over_t[1:] = y[1:] / x[1:]  # Avoid division by zero
    
    axs[1, 0].plot(x[1:], msd_over_t[1:], color='#1f77b4', linewidth=2.5)
    
    # Add horizontal line for the diffusion coefficient
    d_line = diffusion_coef * 1e7 * 6  # Convert to nm²/ps and multiply by 6 for 3D
    axs[1, 0].axhline(y=d_line, color='#d62728', linestyle='--', linewidth=2.5,
                     label=f'D = {diffusion_coef:.3e} cm²/s')
    
    axs[1, 0].set_title('MSD/t vs Time (Normal Diffusion Check)', fontsize=16, fontweight='bold')
    axs[1, 0].set_xlabel('Time (ps)', fontsize=14)
    axs[1, 0].set_ylabel('MSD/t (nm²/ps)', fontsize=14)
    axs[1, 0].legend(fontsize=11, framealpha=0.8, loc='best')
    axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
    
    # 4. Diffusion coefficient comparison with literature
    # Literature values for water diffusion coefficient at different temperatures
    lit_temps = [273, 298, 310, 323]
    lit_diff = [1.099e-5, 2.299e-5, 3.02e-5, 4.01e-5]  # cm²/s
    
    # Plot literature values
    axs[1, 1].plot(lit_temps, lit_diff, 'o-', color='#2ca02c', linewidth=2.5, 
                  label='Literature values')
    
    # Add our calculated value (assuming 298K)
    axs[1, 1].plot(298, diffusion_coef, 'o', color='#d62728', markersize=10,
                  label=f'Our result: {diffusion_coef:.3e} cm²/s')
    
    # Add error bar (assuming 10% error)
    axs[1, 1].errorbar(298, diffusion_coef, yerr=0.1*diffusion_coef, color='#d62728',
                      capsize=5, capthick=2, elinewidth=2)
    
    # Calculate percent difference from literature
    lit_value_298 = 2.299e-5  # cm²/s at 298K
    percent_diff = (diffusion_coef - lit_value_298) / lit_value_298 * 100
    
    # Add text with comparison
    comparison_text = (
        f"Literature (298K): {lit_value_298:.3e} cm²/s\n"
        f"Our result: {diffusion_coef:.3e} cm²/s\n"
        f"Difference: {percent_diff:.1f}%"
    )
    axs[1, 1].text(0.05, 0.95, comparison_text, transform=axs[1, 1].transAxes,
                  fontsize=11, verticalalignment='top',
                  bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
    
    axs[1, 1].set_title('Diffusion Coefficient Comparison', fontsize=16, fontweight='bold')
    axs[1, 1].set_xlabel('Temperature (K)', fontsize=14)
    axs[1, 1].set_ylabel('Diffusion Coefficient (cm²/s)', fontsize=14)
    axs[1, 1].legend(fontsize=11, framealpha=0.8, loc='best')
    axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
    
    # Add a title for the entire figure
    fig.suptitle('Diffusion Analysis for TIP4P Water', fontsize=18, fontweight='bold', y=0.98)
    
    # Add explanation text
    explanation = (
        "Diffusion coefficient (D) is calculated from the slope of the MSD vs time plot using Einstein relation: MSD = 6Dt.\n"
        "Normal diffusion shows a linear relationship between MSD and time (α = 1).\n"
        "Sub-diffusion (α < 1) and super-diffusion (α > 1) indicate anomalous diffusion behavior."
    )
    fig.text(0.5, 0.01, explanation, ha='center', fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
    
    # Add simulation details
    sim_details = (
        f"System: TIP4P Water\n"
        f"Temperature: 298K\n"
        f"Diffusion Coefficient: {diffusion_coef:.3e} cm²/s"
    )
    fig.text(0.02, 0.02, sim_details, fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
    
    # Save figure with tight layout
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])  # Adjust for the title and footer
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def find_linear_region_and_calculate_diffusion(x, y, fit_start_ps=20, fit_end_ps=200, try_multiple_ranges=True):
    """Find the linear region in MSD data and calculate diffusion coefficient"""
    # Convert to numpy arrays if not already
    x = np.array(x)
    y = np.array(y)
    
    if len(x) < 10 or len(y) < 10:
        print("  - Not enough data points for diffusion coefficient calculation")
        return None
    
    # Find indices corresponding to the time range
    start_idx = np.searchsorted(x, fit_start_ps)
    end_idx = np.searchsorted(x, fit_end_ps)
    
    # Adjust indices if out of bounds
    start_idx = max(1, min(start_idx, len(x) - 2))
    end_idx = max(start_idx + 5, min(end_idx, len(x) - 1))
    
    # If try_multiple_ranges is True, try different ranges to find the best fit
    best_r_squared = 0
    best_start_idx = start_idx
    best_end_idx = end_idx
    best_slope = 0
    best_intercept = 0
    
    if try_multiple_ranges:
        # Try different ranges
        ranges = []
        
        # Generate ranges based on the data
        max_time = x[-1]
        
        # Short range (early)
        early_end = min(100, max_time / 3)
        ranges.append((10, early_end))
        
        # Medium range (middle)
        mid_start = max(20, max_time / 5)
        mid_end = min(200, max_time * 2 / 3)
        ranges.append((mid_start, mid_end))
        
        # Long range (late)
        late_start = max(50, max_time / 3)
        late_end = min(400, max_time * 0.9)
        ranges.append((late_start, late_end))
        
        # Try each range
        for range_start, range_end in ranges:
            try:
                r_start_idx = np.searchsorted(x, range_start)
                r_end_idx = np.searchsorted(x, range_end)
                
                # Adjust indices if out of bounds
                r_start_idx = max(1, min(r_start_idx, len(x) - 2))
                r_end_idx = max(r_start_idx + 5, min(r_end_idx, len(x) - 1))
                
                # Skip if range is too small
                if r_end_idx - r_start_idx < 5:
                    continue
                
                # Perform linear regression
                slope, intercept, r_value, p_value, std_err = stats.linregress(
                    x[r_start_idx:r_end_idx+1], y[r_start_idx:r_end_idx+1]
                )
                
                r_squared = r_value ** 2
                
                # Update best fit if this one is better
                if r_squared > best_r_squared:
                    best_r_squared = r_squared
                    best_start_idx = r_start_idx
                    best_end_idx = r_end_idx
                    best_slope = slope
                    best_intercept = intercept
            except:
                # Skip if regression fails
                continue
    else:
        # Just use the specified range
        try:
            # Perform linear regression
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                x[start_idx:end_idx+1], y[start_idx:end_idx+1]
            )
            
            best_r_squared = r_value ** 2
            best_start_idx = start_idx
            best_end_idx = end_idx
            best_slope = slope
            best_intercept = intercept
        except:
            print("  - Linear regression failed")
            return None
    
    # Calculate diffusion coefficient (D = slope/6 for 3D)
    # Convert from nm²/ps to cm²/s
    diffusion_coef = best_slope / 6.0 * 1e-7  # nm²/ps to cm²/s
    
    # Return the results
    return (best_start_idx, best_end_idx, best_slope, best_intercept, best_r_squared, diffusion_coef)

def main():
    if len(sys.argv) < 3:
        print("Usage: plot_msd.py <analysis_dir> <plots_dir>")
        sys.exit(1)
        
    analysis_dir = sys.argv[1]
    plots_dir = sys.argv[2]
    
    # Define data directory
    data_dir = os.path.join(analysis_dir, "data")
    
    # Ensure plots directory exists
    os.makedirs(plots_dir, exist_ok=True)
    
    print('Plotting mean square displacement...')
    
    # MSD plot
    msd_file = os.path.join(data_dir, 'msd.xvg')
    if os.path.exists(msd_file):
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(msd_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Mean Square Displacement vs Time'
            plot_xlabel = xlabel if xlabel else 'Time (ps)'
            plot_ylabel = ylabel if ylabel else 'MSD (nm²)'
            
            # Find the linear region and calculate diffusion coefficient
            fit_results = find_linear_region_and_calculate_diffusion(x, y, fit_start_ps=20, fit_end_ps=200, try_multiple_ranges=True)
            
            # Create enhanced MSD plot
            output_path = os.path.join(plots_dir, 'msd_plot.png')
            plot_msd_enhanced(x, y, plot_title, plot_xlabel, plot_ylabel, output_path, fit_results)
            
            # Create comprehensive diffusion analysis plot
            if fit_results:
                output_path = os.path.join(plots_dir, 'diffusion_analysis_plot.png')
                create_diffusion_analysis_plot(x, y, fit_results, output_path)
    else:
        print(f"MSD file not found: {msd_file}")

if __name__ == "__main__":
    main() 