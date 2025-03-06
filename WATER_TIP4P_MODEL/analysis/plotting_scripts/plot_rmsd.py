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

def calculate_moving_average(data, window_size=20):
    """Calculate moving average with the specified window size."""
    weights = np.ones(window_size) / window_size
    return np.convolve(data, weights, mode='valid')

def detect_equilibration(y, window_size=20, threshold=0.005):
    """
    Detect when the RMSD has equilibrated by looking for a plateau.
    Uses a more robust algorithm that considers both the slope and variance.
    Returns the index where equilibration is detected.
    """
    if len(y) < window_size * 3:
        return len(y) // 3  # Not enough data points, return 1/3 of trajectory
    
    # Calculate rolling average
    rolling_avg = []
    for i in range(len(y) - window_size + 1):
        rolling_avg.append(np.mean(y[i:i+window_size]))
    
    rolling_avg = np.array(rolling_avg)
    
    # Calculate the derivative (slope) of the rolling average
    derivatives = np.abs(np.diff(rolling_avg))
    
    # Calculate rolling variance
    rolling_var = []
    for i in range(len(y) - window_size + 1):
        rolling_var.append(np.var(y[i:i+window_size]))
    
    rolling_var = np.array(rolling_var[:-1])  # Match length with derivatives
    
    # Normalize both metrics for combined scoring
    if np.max(derivatives) > 0:
        norm_derivatives = derivatives / np.max(derivatives)
    else:
        norm_derivatives = derivatives
    
    if np.max(rolling_var) > 0:
        norm_variance = rolling_var / np.max(rolling_var)
    else:
        norm_variance = rolling_var
    
    # Combined score (lower is better)
    combined_score = norm_derivatives + norm_variance
    
    # Find the first window where the score is consistently low
    for i in range(window_size, len(combined_score) - window_size):
        if np.all(combined_score[i:i+window_size] < threshold):
            return i + window_size  # Return the end of the stable window
    
    # If no clear equilibration is detected, use a heuristic
    # Look for the point where the derivative becomes consistently small
    for i in range(len(derivatives) - window_size):
        if np.all(derivatives[i:i+window_size] < np.mean(derivatives) * 0.2):
            return i + window_size
    
    # If still no clear point, return 1/3 of the trajectory
    return len(y) // 3

def calculate_rmsd_statistics(y, equil_idx):
    """Calculate statistics for the equilibrated part of the RMSD."""
    if equil_idx >= len(y):
        equil_idx = 0
    
    equil_data = y[equil_idx:]
    
    stats = {
        'mean': np.mean(equil_data),
        'std': np.std(equil_data),
        'min': np.min(equil_data),
        'max': np.max(equil_data),
        'final': y[-1] if len(y) > 0 else 0
    }
    
    return stats

def plot_rmsd(x, y, title, xlabel, ylabel, output_path):
    """Create an enhanced RMSD plot with improved visualization."""
    # Create figure with subplots - main plot and histogram
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), dpi=300, 
                                  gridspec_kw={'height_ratios': [3, 1]}, sharex=False)
    
    # Detect equilibration
    equil_idx = detect_equilibration(y)
    equil_time = x[equil_idx] if equil_idx < len(x) else x[0]
    
    # Extract equilibrated portion
    x_equil = x[equil_idx:]
    y_equil = y[equil_idx:]
    
    # Plot the raw RMSD data on the main plot
    ax1.plot(x, y, color='#1f77b4', linewidth=2, alpha=0.8, label='Raw RMSD')
    
    # Calculate and plot moving average
    if len(y) > 20:
        window_size = 20
        y_ma = calculate_moving_average(y, window_size)
        x_ma = x[window_size-1:]
        ax1.plot(x_ma, y_ma, color='#ff7f0e', linewidth=2.5, linestyle='-', 
                label=f'Moving Avg (n={window_size})')
    
    # Calculate statistics
    stats = calculate_rmsd_statistics(y, equil_idx)
    
    # Add equilibration region shading
    ax1.axvspan(x[0], equil_time, alpha=0.2, color='#ffb6c1', label='Equilibration')
    
    # Add vertical line at equilibration point
    ax1.axvline(x=equil_time, color='red', linestyle='--', linewidth=2, alpha=0.7, 
                label=f'Equilibration at {equil_time:.1f} {xlabel}')
    
    # Add mean line and standard deviation band for production phase
    ax1.axhline(y=stats['mean'], color='green', linestyle='--', linewidth=2, alpha=0.7,
               label=f'Production Mean: {stats["mean"]:.3f} nm')
    
    # Add standard deviation band
    ax1.fill_between(x_equil, stats['mean'] - stats['std'], stats['mean'] + stats['std'],
                    color='green', alpha=0.2, label=f'Std Dev: ±{stats["std"]:.3f} nm')
    
    # Annotate final RMSD value
    ax1.annotate(f'Final RMSD: {stats["final"]:.3f} nm', 
                xy=(x[-1], y[-1]),
                xytext=(x[-1] - (x[-1] * 0.1), y[-1] + (max(y) * 0.05)),
                arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
                fontsize=12)
    
    # Add statistics text box
    stats_text = (
        f"Production Phase Statistics:\n"
        f"Mean: {stats['mean']:.3f} nm\n"
        f"Std Dev: {stats['std']:.3f} nm\n"
        f"Min: {stats['min']:.3f} nm\n"
        f"Max: {stats['max']:.3f} nm\n"
        f"Final: {stats['final']:.3f} nm"
    )
    ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, fontsize=12,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    # Improve main plot styling
    ax1.set_title(f"{title}\nEquilibration (0-{equil_time:.1f} {xlabel}) vs. Production ({equil_time:.1f}-{x[-1]:.1f} {xlabel})", 
                 fontsize=16, fontweight='bold')
    ax1.set_ylabel(ylabel, fontsize=14)
    ax1.set_xlabel(xlabel, fontsize=14)
    ax1.grid(True, linestyle='--', alpha=0.4)
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax1.legend(fontsize=10, loc='upper left')
    
    # Plot RMSD distribution histogram for production phase only
    if HAS_SEABORN:
        sns.histplot(y_equil, kde=True, ax=ax2, color='#1f77b4', alpha=0.7)
    else:
        # Fallback to matplotlib if seaborn is not available
        n, bins, patches = ax2.hist(y_equil, bins=30, alpha=0.7, color='#1f77b4', density=True)
        # Add KDE curve
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(y_equil)
        x_kde = np.linspace(min(y_equil), max(y_equil), 100)
        ax2.plot(x_kde, kde(x_kde), 'r-', linewidth=2)
    
    # Add vertical lines for mean and std dev
    ax2.axvline(x=stats['mean'], color='green', linestyle='--', linewidth=2, 
               label=f'Mean: {stats["mean"]:.3f} nm')
    ax2.axvline(x=stats['mean'] - stats['std'], color='green', linestyle=':', linewidth=1.5, alpha=0.7)
    ax2.axvline(x=stats['mean'] + stats['std'], color='green', linestyle=':', linewidth=1.5, alpha=0.7)
    
    # Improve histogram styling
    ax2.set_title('RMSD Distribution (Production Phase Only)', fontsize=14)
    ax2.set_xlabel(ylabel, fontsize=14)
    ax2.set_ylabel('Frequency', fontsize=14)
    ax2.grid(True, linestyle='--', alpha=0.4)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax2.legend(fontsize=10, loc='upper right')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def plot_rmsd_equilibrated(x, y, title, xlabel, ylabel, output_path):
    """Create an enhanced plot of the equilibrated RMSD with moving average and distribution."""
    # Detect equilibration
    equil_idx = detect_equilibration(y)
    
    if equil_idx >= len(y) - 5:  # Need at least a few points
        equil_idx = len(y) // 3
    
    # Extract equilibrated portion
    x_equil = x[equil_idx:]
    y_equil = y[equil_idx:]
    
    # Create figure with subplots - main plot and histogram
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), dpi=300, 
                                  gridspec_kw={'height_ratios': [3, 1]}, sharex=False)
    
    # Plot the raw RMSD data on the main plot
    ax1.plot(x_equil, y_equil, color='#1f77b4', linewidth=2, alpha=0.8)
    
    # Calculate and plot moving average
    if len(y_equil) > 20:
        window_size = 20
        y_ma = calculate_moving_average(y_equil, window_size)
        x_ma = x_equil[window_size-1:]
        ax1.plot(x_ma, y_ma, color='#ff7f0e', linewidth=2.5, linestyle='-', 
                label=f'Moving Avg (n={window_size})')
    
    # Calculate statistics
    stats = calculate_rmsd_statistics(y, equil_idx)
    
    # Add equilibration region shading
    ax1.axvspan(x[0], x[equil_idx], alpha=0.2, color='#ffb6c1', label='Equilibration')
    
    # Add mean line and standard deviation band
    ax1.axhline(y=stats['mean'], color='green', linestyle='--', linewidth=2, alpha=0.7,
               label=f'Production Mean: {stats["mean"]:.3f} nm')
    
    # Add standard deviation band
    ax1.fill_between(x_equil, stats['mean'] - stats['std'], stats['mean'] + stats['std'],
                    color='green', alpha=0.2, label=f'Std Dev: ±{stats["std"]:.3f} nm')
    
    # Annotate final RMSD value
    ax1.annotate(f'Final RMSD: {stats["final"]:.3f} nm', 
                xy=(x_equil[-1], y_equil[-1]),
                xytext=(x_equil[-1], y_equil[-1] + 0.1),
                arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
                fontsize=12)
    
    # Improve main plot styling
    ax1.set_title(title, fontsize=16, fontweight='bold')
    ax1.set_ylabel(ylabel, fontsize=14)
    ax1.grid(True, linestyle='--', alpha=0.4)
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax1.legend(fontsize=10, loc='upper left')
    
    # Plot RMSD distribution histogram
    sns.histplot(y_equil, kde=True, ax=ax2, color='#1f77b4', alpha=0.7)
    
    # Add vertical lines for mean and std dev
    ax2.axvline(x=stats['mean'], color='green', linestyle='--', linewidth=2, 
               label=f'Mean: {stats["mean"]:.3f} nm')
    ax2.axvline(x=stats['mean'] - stats['std'], color='green', linestyle=':', linewidth=1.5, alpha=0.7)
    ax2.axvline(x=stats['mean'] + stats['std'], color='green', linestyle=':', linewidth=1.5, alpha=0.7)
    
    # Improve histogram styling
    ax2.set_title('RMSD Distribution (Equilibrated Phase)', fontsize=14)
    ax2.set_xlabel(ylabel, fontsize=14)
    ax2.set_ylabel('Frequency', fontsize=14)
    ax2.grid(True, linestyle='--', alpha=0.4)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax2.legend(fontsize=10, loc='upper right')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def main():
    """Main function to generate RMSD plots."""
    if len(sys.argv) != 3:
        print("Usage: python plot_rmsd.py <analysis_dir> <output_dir>")
        sys.exit(1)
    
    analysis_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Path to RMSD data file
    rmsd_file = os.path.join(analysis_dir, "data", "rmsd.xvg")
    
    if not os.path.exists(rmsd_file):
        print(f"Error: RMSD file not found at {rmsd_file}")
        sys.exit(1)
    
    # Read RMSD data
    x, y, title, xlabel, ylabel, legend_labels = read_xvg(rmsd_file)
    
    if len(x) == 0 or len(y) == 0:
        print("Error: No data found in RMSD file")
        sys.exit(1)
    
    # Set default labels if not provided in the file
    if not title:
        title = "Root Mean Square Deviation (RMSD)"
    if not xlabel:
        xlabel = "Time (ns)"
    if not ylabel:
        ylabel = "RMSD (nm)"
    
    # Generate plots
    plot_rmsd(x, y, title, xlabel, ylabel, os.path.join(output_dir, "rmsd_plot.png"))
    plot_rmsd_equilibrated(x, y, title, xlabel, ylabel, os.path.join(output_dir, "rmsd_equilibrated_plot.png"))

if __name__ == "__main__":
    main()
