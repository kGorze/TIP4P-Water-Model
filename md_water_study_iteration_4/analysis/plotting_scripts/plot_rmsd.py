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
                else:
                    values = line.strip().split()
                    if len(values) >= 2:
                        x.append(float(values[0]))
                        y.append(float(values[1]))
        return np.array(x), np.array(y), title, xlabel, ylabel
    except Exception as e:
        print(f'Error reading {filename}: {e}')
        return np.array([]), np.array([]), "", "", ""

def detect_equilibration(y, window_size=20, threshold=0.01):
    """
    Detect when the system has equilibrated based on the slope of the RMSD
    
    Parameters:
    -----------
    y : array
        RMSD values
    window_size : int
        Size of the window to use for slope calculation
    threshold : float
        Threshold for the absolute slope to consider the system equilibrated
        
    Returns:
    --------
    equil_idx : int
        Index at which the system is considered equilibrated
    """
    if len(y) < window_size * 2:
        return len(y) // 2  # Default to middle if not enough data
    
    # Calculate slopes in windows
    slopes = []
    for i in range(len(y) - window_size):
        window = y[i:i+window_size]
        x_window = np.arange(window_size)
        slope, _, _, _, _ = stats.linregress(x_window, window)
        slopes.append(abs(slope))
    
    # Find where the slope consistently stays below threshold
    for i in range(len(slopes) - window_size):
        if all(s < threshold for s in slopes[i:i+window_size]):
            return i + window_size
    
    # If no clear equilibration is found, use a default (e.g., 20% of trajectory)
    return int(len(y) * 0.2)

def calculate_rmsd_statistics(y, equil_idx):
    """Calculate statistics for the equilibrated portion of the RMSD"""
    equil_data = y[equil_idx:]
    
    stats_dict = {
        'mean': np.mean(equil_data),
        'std': np.std(equil_data),
        'min': np.min(equil_data),
        'max': np.max(equil_data),
        'final': equil_data[-1] if len(equil_data) > 0 else 0
    }
    
    return stats_dict

def main():
    if len(sys.argv) < 3:
        print("Usage: plot_rmsd.py <analysis_dir> <plots_dir>")
        sys.exit(1)
        
    analysis_dir = sys.argv[1]
    plots_dir = sys.argv[2]
    
    # Ensure plots directory exists
    os.makedirs(plots_dir, exist_ok=True)
    
    # RMSD Plot
    rmsd_file = os.path.join(analysis_dir, 'rmsd.xvg')
    if os.path.exists(rmsd_file):
        print('Plotting RMSD...')
        x, y, file_title, file_xlabel, file_ylabel = read_xvg(rmsd_file)
        
        # Use file metadata if available, otherwise use defaults
        title = file_title if file_title else 'Root Mean Square Deviation vs Time'
        xlabel = file_xlabel if file_xlabel else 'Time (ns)'
        ylabel = file_ylabel if file_ylabel else 'RMSD (nm)'
        
        if len(x) > 0 and len(y) > 0:
            # Create figure with two subplots - main plot and histogram
            fig = plt.figure(figsize=(12, 8), dpi=300)
            
            # Main RMSD plot (larger)
            ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=3, rowspan=2)
            
            # Plot with seaborn if available
            if HAS_SEABORN:
                sns.lineplot(x=x, y=y, color='#1f77b4', linewidth=2, ax=ax1)
            else:
                ax1.plot(x, y, color='#1f77b4', linewidth=2)
            
            # Detect equilibration point
            equil_idx = detect_equilibration(y)
            equil_time = x[equil_idx]
            
            # Add a rolling average to smooth the curve
            window = min(len(x) // 20, 30)  # 5% of data points or 30, whichever is smaller
            if window > 1:
                rolling_avg = np.convolve(y, np.ones(window)/window, mode='valid')
                rolling_x = x[window-1:]
                ax1.plot(rolling_x, rolling_avg, color='#ff7f0e', linewidth=2, 
                        linestyle='--', label=f'Moving Avg (n={window})')
            
            # Highlight equilibration and production phases
            ax1.axvspan(x[0], equil_time, alpha=0.2, color='#d62728', label='Equilibration')
            
            # Calculate statistics for production phase
            rmsd_stats = calculate_rmsd_statistics(y, equil_idx)
            
            # Add horizontal line for production mean
            ax1.axhline(y=rmsd_stats['mean'], color='#2ca02c', linestyle='--', alpha=0.7, 
                       label=f'Production Mean: {rmsd_stats["mean"]:.3f} nm')
            
            # Add shaded area for standard deviation
            ax1.fill_between(x[equil_idx:], rmsd_stats['mean'] - rmsd_stats['std'], 
                           rmsd_stats['mean'] + rmsd_stats['std'], 
                           color='#2ca02c', alpha=0.2, 
                           label=f'Std Dev: ±{rmsd_stats["std"]:.3f} nm')
            
            # Add annotation for final RMSD
            final_rmsd = y[-1]
            ax1.annotate(f'Final RMSD: {final_rmsd:.3f} nm', 
                        xy=(x[-1], final_rmsd),
                        xytext=(x[-1] - (x[-1]-x[0])*0.2, final_rmsd + rmsd_stats['std']),
                        arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                        fontsize=10)
            
            ax1.set_xlabel(xlabel, fontsize=14)
            ax1.set_ylabel(ylabel, fontsize=14)
            ax1.set_title(title, fontsize=16)
            ax1.grid(True, alpha=0.3)
            ax1.legend(fontsize=10, loc='upper left')
            
            # Set y-axis to start from 0 for better visualization
            ax1.set_ylim(bottom=0)
            
            # Add histogram of RMSD values for equilibrated portion
            ax2 = plt.subplot2grid((3, 3), (2, 0), colspan=3, rowspan=1)
            equil_data = y[equil_idx:]
            
            if HAS_SEABORN:
                sns.histplot(equil_data, kde=True, color='#1f77b4', ax=ax2)
            else:
                ax2.hist(equil_data, bins=20, alpha=0.7, color='#1f77b4', density=True)
                
                # Add a simple KDE if seaborn is not available
                from scipy.stats import gaussian_kde
                kde = gaussian_kde(equil_data)
                x_kde = np.linspace(min(equil_data), max(equil_data), 100)
                ax2.plot(x_kde, kde(x_kde), 'r-', linewidth=2)
            
            # Add vertical lines for statistics
            ax2.axvline(x=rmsd_stats['mean'], color='#2ca02c', linestyle='--', 
                       label=f'Mean: {rmsd_stats["mean"]:.3f} nm')
            ax2.axvline(x=rmsd_stats['mean'] - rmsd_stats['std'], color='#2ca02c', 
                       linestyle=':', alpha=0.7)
            ax2.axvline(x=rmsd_stats['mean'] + rmsd_stats['std'], color='#2ca02c', 
                       linestyle=':', alpha=0.7)
            
            ax2.set_xlabel('RMSD (nm)', fontsize=12)
            ax2.set_ylabel('Frequency', fontsize=12)
            ax2.set_title('RMSD Distribution (Equilibrated Phase)', fontsize=14)
            ax2.grid(True, alpha=0.3)
            ax2.legend(fontsize=10)
            
            # Add a text box with statistics
            stats_text = (
                f"Equilibration time: {equil_time:.1f} ns\n"
                f"Mean RMSD: {rmsd_stats['mean']:.3f} nm\n"
                f"Std Dev: {rmsd_stats['std']:.3f} nm\n"
                f"Min: {rmsd_stats['min']:.3f} nm\n"
                f"Max: {rmsd_stats['max']:.3f} nm\n"
                f"Final: {rmsd_stats['final']:.3f} nm"
            )
            
            # Add text box to the main plot
            props = dict(boxstyle='round', facecolor='white', alpha=0.7)
            ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, fontsize=10,
                    verticalalignment='top', bbox=props)
            
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, 'rmsd_plot.png'), dpi=300, bbox_inches='tight')
            plt.close()
            print('  - rmsd_plot.png saved successfully')
            
            # Create a second plot focusing on just the equilibrated portion
            plt.figure(figsize=(10, 6), dpi=300)
            
            if HAS_SEABORN:
                sns.lineplot(x=x[equil_idx:], y=y[equil_idx:], color='#1f77b4', linewidth=2)
            else:
                plt.plot(x[equil_idx:], y[equil_idx:], color='#1f77b4', linewidth=2)
            
            plt.axhline(y=rmsd_stats['mean'], color='#2ca02c', linestyle='--', alpha=0.7,
                       label=f'Mean: {rmsd_stats["mean"]:.3f} nm')
            plt.fill_between(x[equil_idx:], rmsd_stats['mean'] - rmsd_stats['std'], 
                           rmsd_stats['mean'] + rmsd_stats['std'], 
                           color='#2ca02c', alpha=0.2)
            
            plt.xlabel(xlabel, fontsize=14)
            plt.ylabel(ylabel, fontsize=14)
            plt.title('Equilibrated RMSD vs Time', fontsize=16)
            plt.grid(True, alpha=0.3)
            plt.legend(fontsize=12)
            
            # Set y-axis limits to focus on the equilibrated values
            y_min = max(0, rmsd_stats['mean'] - 3*rmsd_stats['std'])
            y_max = rmsd_stats['mean'] + 3*rmsd_stats['std']
            plt.ylim(y_min, y_max)
            
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, 'rmsd_equilibrated_plot.png'), dpi=300, bbox_inches='tight')
            plt.close()
            print('  - rmsd_equilibrated_plot.png saved successfully')
    else:
        print(f"RMSD file not found: {rmsd_file}")

if __name__ == "__main__":
    main()
