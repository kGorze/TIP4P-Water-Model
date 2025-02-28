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

# Reference values for TIP4P water at 298K
REFERENCE_VALUES = {
    'diffusion_coefficient': 2.3e-9,  # m^2/s, self-diffusion coefficient
    'viscosity': 0.896e-3,            # Pa·s, viscosity
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

def calculate_diffusion_coefficient(time, msd, fit_start_fraction=0.4, fit_end_fraction=0.8):
    """
    Calculate the diffusion coefficient from MSD data using Einstein relation
    
    Parameters:
    -----------
    time : array
        Time values in ps
    msd : array
        MSD values in nm^2
    fit_start_fraction : float
        Fraction of the trajectory to start the linear fit (default: 0.4)
    fit_end_fraction : float
        Fraction of the trajectory to end the linear fit (default: 0.8)
        
    Returns:
    --------
    D : float
        Diffusion coefficient in m^2/s
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
    # Determine the indices for fitting
    fit_start_idx = int(len(time) * fit_start_fraction)
    fit_end_idx = int(len(time) * fit_end_fraction)
    
    # Ensure we have enough points for fitting
    if fit_end_idx - fit_start_idx < 10:
        fit_start_idx = max(0, len(time) // 4)
        fit_end_idx = min(len(time), len(time) * 3 // 4)
    
    # Extract the data for fitting
    fit_time = time[fit_start_idx:fit_end_idx]
    fit_msd = msd[fit_start_idx:fit_end_idx]
    
    # Perform linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(fit_time, fit_msd)
    
    # Calculate diffusion coefficient using Einstein relation: MSD = 6Dt
    # For 3D diffusion, D = slope / 6
    # Convert from nm^2/ps to m^2/s: 1 nm^2/ps = 1e-9 m^2 / 1e-12 s = 1e-3 m^2/s
    D = slope / 6.0 * 1e-3
    
    return D, slope, intercept, r_value, p_value, std_err, fit_start_idx, fit_end_idx

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
    
    # Calculate diffusion coefficient
    D, slope, intercept, r_value, p_value, std_err, fit_start_idx, fit_end_idx = calculate_diffusion_coefficient(x, msd_data)
    
    # Plot the linear fit
    fit_x = np.array([x[fit_start_idx], x[fit_end_idx]])
    fit_y = slope * fit_x + intercept
    ax1.plot(fit_x, fit_y, 'r--', linewidth=2, 
            label=f'Linear fit (R²={r_value**2:.3f})')
    
    # Add annotation for diffusion coefficient
    ax1.text(0.05, 0.95, 
            f'D = {D:.3e} m²/s\n'
            f'Slope = {slope:.3f} nm²/ps\n'
            f'Reference D = {REFERENCE_VALUES["diffusion_coefficient"]:.3e} m²/s\n'
            f'Difference: {((D - REFERENCE_VALUES["diffusion_coefficient"]) / REFERENCE_VALUES["diffusion_coefficient"] * 100):.1f}%',
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
    
    # Create a detailed diffusion analysis plot
    plt.figure(figsize=(10, 8), dpi=300)
    
    # Plot MSD data
    if HAS_SEABORN:
        sns.lineplot(x=x, y=msd_data, color='#1f77b4', linewidth=2, label='MSD')
    else:
        plt.plot(x, msd_data, color='#1f77b4', linewidth=2, label='MSD')
    
    # Plot the linear fit for the entire fitting region
    fit_x_full = x[fit_start_idx:fit_end_idx]
    fit_y_full = slope * fit_x_full + intercept
    plt.plot(fit_x_full, fit_y_full, 'r-', linewidth=2, 
            label=f'Linear fit (R²={r_value**2:.3f})')
    
    # Highlight the fitting region
    plt.axvspan(x[fit_start_idx], x[fit_end_idx], alpha=0.2, color='green',
               label=f'Fitting region ({fit_start_idx}-{fit_end_idx})')
    
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
        f"Diffusion Coefficient:\n"
        f"D = {D:.3e} m²/s\n"
        f"Slope = {slope:.3f} nm²/ps\n"
        f"Intercept = {intercept:.3f} nm²\n"
        f"R² = {r_value**2:.3f}\n"
        f"Standard Error = {std_err:.3e}\n\n"
        f"Reference D = {REFERENCE_VALUES['diffusion_coefficient']:.3e} m²/s\n"
        f"Difference: {((D - REFERENCE_VALUES['diffusion_coefficient']) / REFERENCE_VALUES['diffusion_coefficient'] * 100):.1f}%"
    )
    
    plt.text(0.02, 0.98, info_text, transform=plt.gca().transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
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