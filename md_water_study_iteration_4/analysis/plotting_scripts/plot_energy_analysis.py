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

def calculate_block_averages(x, y, num_blocks=5):
    """Calculate block averages for the data."""
    if len(x) == 0 or len(y) == 0:
        return [], [], []
    
    # Determine block size
    block_size = len(x) // num_blocks
    
    block_means = []
    block_stds = []
    block_centers = []
    
    for i in range(num_blocks):
        start_idx = i * block_size
        end_idx = (i + 1) * block_size if i < num_blocks - 1 else len(x)
        
        block_data = y[start_idx:end_idx]
        block_means.append(np.mean(block_data))
        block_stds.append(np.std(block_data))
        block_centers.append(np.mean(x[start_idx:end_idx]))
    
    return block_centers, block_means, block_stds

def calculate_drift(block_means):
    """Calculate drift metrics from block means."""
    if len(block_means) < 2:
        return 0, 0, False
    
    # Calculate linear regression to detect trend
    x = np.arange(len(block_means))
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, block_means)
    
    # Calculate percent change from first to last block
    percent_change = ((block_means[-1] - block_means[0]) / abs(block_means[0])) * 100
    
    # Determine if drift is significant (more than 1% change and p-value < 0.05)
    significant_drift = (abs(percent_change) > 1.0) and (p_value < 0.05)
    
    return slope, percent_change, significant_drift

def plot_enhanced_energy_analysis(x, y, title, xlabel, ylabel, output_path, num_blocks=5):
    """Create an enhanced energy analysis plot with block averaging and drift detection."""
    # Create a figure with two subplots
    fig = plt.figure(figsize=(15, 12), dpi=300)
    gs = plt.GridSpec(3, 2, figure=fig, height_ratios=[2, 1, 1])
    
    # Main time series plot
    ax1 = fig.add_subplot(gs[0, :])
    
    # Extract data for plotting
    if len(y.shape) > 1 and y.shape[1] > 0:
        data = y[:, 0]  # Use first column if multiple columns
    else:
        data = y
    
    # Calculate basic statistics
    mean_val = np.mean(data)
    std_val = np.std(data)
    min_val = np.min(data)
    max_val = np.max(data)
    final_val = data[-1] if len(data) > 0 else 0
    
    # Plot the time series
    if HAS_SEABORN:
        sns.lineplot(x=x, y=data, color='#1f77b4', linewidth=1.5, ax=ax1)
    else:
        ax1.plot(x, data, color='#1f77b4', linewidth=1.5)
    
    # Add mean line and standard deviation band
    ax1.axhline(y=mean_val, color='#2ca02c', linestyle='--', alpha=0.7, label=f'Mean: {mean_val:.2f}')
    ax1.fill_between(x, mean_val - std_val, mean_val + std_val, color='#2ca02c', alpha=0.2, 
                    label=f'Std Dev: ±{std_val:.2f}')
    
    # Add statistics text box
    stats_text = (
        f"Mean: {mean_val:.2f}\n"
        f"Std Dev: {std_val:.2f}\n"
        f"Min: {min_val:.2f}\n"
        f"Max: {max_val:.2f}\n"
        f"Final: {final_val:.2f}\n"
        f"Relative Fluctuation: ±{(std_val/abs(mean_val))*100:.1f}%"
    )
    ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    ax1.set_xlabel(xlabel, fontsize=12)
    ax1.set_ylabel(ylabel, fontsize=12)
    ax1.set_title(f"{title} - Full Time Series", fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10, loc='lower right')
    
    # Block averaging plot
    ax2 = fig.add_subplot(gs[1, :])
    
    # Calculate block averages
    block_centers, block_means, block_stds = calculate_block_averages(x, data, num_blocks)
    
    # Calculate drift metrics
    slope, percent_change, significant_drift = calculate_drift(block_means)
    
    # Plot block averages with error bars
    ax2.errorbar(block_centers, block_means, yerr=block_stds, fmt='o-', color='#d62728', 
                ecolor='#d62728', elinewidth=1, capsize=4, label='Block Averages')
    
    # Add horizontal line for overall mean
    ax2.axhline(y=mean_val, color='#2ca02c', linestyle='--', alpha=0.7, label=f'Overall Mean: {mean_val:.2f}')
    
    # Add drift information
    drift_text = (
        f"Block Analysis (n={num_blocks}):\n"
        f"First Block Mean: {block_means[0]:.2f}\n"
        f"Last Block Mean: {block_means[-1]:.2f}\n"
        f"Change: {percent_change:+.2f}%\n"
        f"Drift Significant: {'Yes' if significant_drift else 'No'}"
    )
    ax2.text(0.02, 0.98, drift_text, transform=ax2.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    ax2.set_xlabel(xlabel, fontsize=12)
    ax2.set_ylabel(ylabel, fontsize=12)
    ax2.set_title("Block Averages Analysis", fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10, loc='lower right')
    
    # Distribution plot
    ax3 = fig.add_subplot(gs[2, 0])
    
    # Plot histogram with KDE
    if HAS_SEABORN:
        sns.histplot(data, kde=True, color='#1f77b4', ax=ax3)
    else:
        ax3.hist(data, bins=30, alpha=0.7, color='#1f77b4', density=True)
    
    # Add vertical lines for mean and final value
    ax3.axvline(x=mean_val, color='#2ca02c', linestyle='--', alpha=0.7, label=f'Mean: {mean_val:.2f}')
    ax3.axvline(x=final_val, color='#d62728', linestyle=':', alpha=0.7, label=f'Final: {final_val:.2f}')
    
    ax3.set_xlabel(ylabel, fontsize=12)  # Use y-axis label from main plot
    ax3.set_ylabel('Density', fontsize=12)
    ax3.set_title("Energy Distribution", fontsize=14)
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=10, loc='best')
    
    # Running average plot
    ax4 = fig.add_subplot(gs[2, 1])
    
    # Calculate running average with different window sizes
    window_sizes = [len(x)//50, len(x)//20, len(x)//10]
    window_sizes = [max(w, 2) for w in window_sizes]  # Ensure minimum window size of 2
    
    for window in window_sizes:
        running_avg = np.convolve(data, np.ones(window)/window, mode='valid')
        running_x = x[window-1:]
        if len(running_x) == len(running_avg):
            ax4.plot(running_x, running_avg, label=f'{window}-point Moving Avg')
    
    ax4.set_xlabel(xlabel, fontsize=12)
    ax4.set_ylabel(ylabel, fontsize=12)
    ax4.set_title("Running Averages", fontsize=14)
    ax4.grid(True, alpha=0.3)
    ax4.legend(fontsize=10, loc='best')
    
    # Add a title for the entire figure
    plt.suptitle(f"Enhanced Energy Analysis - TIP4P Water Model", fontsize=16, y=0.98)
    
    # Add a watermark with simulation details
    plt.figtext(0.5, 0.01, 'TIP4P Water Model - Detailed Energy Terms', 
               ha='center', fontsize=10, style='italic', alpha=0.7)
    
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def main():
    if len(sys.argv) < 3:
        print("Usage: plot_energy_analysis.py <analysis_dir> <plots_dir>")
        sys.exit(1)
        
    analysis_dir = sys.argv[1]
    plots_dir = sys.argv[2]
    
    # Ensure plots directory exists
    os.makedirs(plots_dir, exist_ok=True)
    
    # Plot energy terms data with enhanced analysis
    energy_terms_file = os.path.join(analysis_dir, 'energy_terms.xvg')
    if os.path.exists(energy_terms_file):
        print('Creating enhanced energy terms analysis...')
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(energy_terms_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'GROMACS Energies'
            plot_xlabel = xlabel if xlabel else 'Time (ns)'
            plot_ylabel = ylabel if ylabel else 'Energy (kJ/mol)'
            
            output_path = os.path.join(plots_dir, 'energy_terms_enhanced_analysis.png')
            plot_enhanced_energy_analysis(x, y, plot_title, plot_xlabel, plot_ylabel, output_path)
    else:
        print(f"Energy terms file not found: {energy_terms_file}")
        
        # Try alternative energy file
        energy_file = os.path.join(analysis_dir, 'energy.xvg')
        if os.path.exists(energy_file):
            print('Using general energy file for enhanced analysis...')
            x, y, title, xlabel, ylabel, legend_labels = read_xvg(energy_file)
            
            if len(x) > 0 and len(y) > 0:
                # Use file metadata if available, otherwise use defaults
                plot_title = title if title else 'GROMACS Energies'
                plot_xlabel = xlabel if xlabel else 'Time (ns)'
                plot_ylabel = ylabel if ylabel else 'Energy (kJ/mol)'
                
                output_path = os.path.join(plots_dir, 'energy_enhanced_analysis.png')
                plot_enhanced_energy_analysis(x, y, plot_title, plot_xlabel, plot_ylabel, output_path)
        else:
            print(f"No energy files found for enhanced analysis")

if __name__ == "__main__":
    main() 