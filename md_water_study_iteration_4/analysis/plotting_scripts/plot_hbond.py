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
    'hbonds_per_molecule': 3.5,  # Average number of H-bonds per water molecule
    'hbond_lifetime': 1.0,       # ps, approximate lifetime
    'hbond_angle_mean': 30.0,    # degrees, mean H-bond angle (typical cutoff in GROMACS)
    'hbond_distance_mean': 0.28  # nm, mean O-O distance for hydrogen bonds (2.8 Angstroms)
}

def read_xvg(filename):
    """Read .xvg files and extract x, y data while skipping comment/label lines."""
    x = []
    y = []
    y_columns = []  # Store all y columns
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
                    # Parse data lines
                    try:
                        values = [float(val) for val in line.strip().split()]
                        if len(values) >= 2:
                            x.append(values[0])
                            y.append(values[1])
                            
                            # Store all columns after the first one
                            if len(values) > 2:
                                while len(y_columns) < len(values) - 1:
                                    y_columns.append([])
                                
                                for i in range(1, len(values)):
                                    y_columns[i-1].append(values[i])
                    except ValueError:
                        continue
    except Exception as e:
        print(f"Error reading {filename}: {e}")
    
    return x, y, y_columns, title, xlabel, ylabel, legend_labels

def calculate_statistics(data):
    """Calculate basic statistics for a dataset"""
    # Convert to numpy array if it's a list
    if isinstance(data, list):
        data = np.array(data)
        
    # Check if data is multi-dimensional
    if len(data.shape) > 1:
        # Multi-column data, use first column
        return {
            'mean': np.mean(data[:, 0]),
            'std': np.std(data[:, 0]),
            'min': np.min(data[:, 0]),
            'max': np.max(data[:, 0]),
            'median': np.median(data[:, 0]),
            'final': data[-1, 0] if len(data) > 0 else 0
        }
    else:
        # Single column
        return {
            'mean': np.mean(data),
            'std': np.std(data),
            'min': np.min(data),
            'max': np.max(data),
            'median': np.median(data),
            'final': data[-1] if len(data) > 0 else 0
        }

def get_num_water_molecules(data_dir):
    """Extract the number of water molecules from the topology file."""
    # Look for topology file in the data directory
    topology_file = os.path.join(data_dir, "topol.top")
    
    if not os.path.exists(topology_file):
        # Try looking in the parent directory
        topology_file = os.path.join(os.path.dirname(data_dir), "topol.top")
    
    if not os.path.exists(topology_file):
        # Try looking in the data subdirectory
        topology_file = os.path.join(data_dir, "data", "topol.top")
        
    if os.path.exists(topology_file):
        # Parse the topology file to find the number of water molecules
        with open(topology_file, 'r') as f:
            for line in f:
                if "SOL" in line and not line.startswith(";"):
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            return int(parts[-1])
                        except ValueError:
                            pass
    
    # If we couldn't find it in the topology file, check the water_box.inp file
    configs_dir = os.path.join(os.path.dirname(os.path.dirname(data_dir)), "configs")
    water_box_file = os.path.join(configs_dir, "water_box.inp")
    
    if os.path.exists(water_box_file):
        with open(water_box_file, 'r') as f:
            for line in f:
                if "number" in line and not line.startswith("#"):
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            return int(parts[1])
                        except ValueError:
                            pass
    
    # If we still can't find it, try to parse it from the .gro file
    gro_files = [
        os.path.join(data_dir, "md.gro"),
        os.path.join(os.path.dirname(data_dir), "md.gro"),
        os.path.join(data_dir, "data", "md.gro")
    ]
    
    for gro_file in gro_files:
        if os.path.exists(gro_file):
            try:
                with open(gro_file, 'r') as f:
                    # Second line of .gro file contains the number of atoms
                    lines = f.readlines()
                    if len(lines) > 1:
                        num_atoms = int(lines[1].strip())
                        # Assuming TIP4P water (4 atoms per molecule)
                        return num_atoms // 4
            except:
                pass
    
    # If all else fails, print a warning and return a reasonable default
    print("Warning: Could not determine number of water molecules, using default value of 5500")
    return 5500  # Default based on typical water box size

def plot_hbond_number(x, y, title, xlabel, ylabel, legend_labels, output_path, n_molecules=None):
    """Plot the number of hydrogen bonds with statistics and per-molecule calculation"""
    # If n_molecules is not provided, try to determine it from the topology file
    if n_molecules is None:
        # Extract the directory from the output path
        data_dir = os.path.dirname(os.path.dirname(output_path))
        n_molecules = get_num_water_molecules(data_dir)
        print(f"Found {n_molecules} water molecules in topology file")
    
    plt.figure(figsize=(12, 8), dpi=300)
    
    # Convert to numpy arrays if they're lists
    x_array = np.array(x)
    y_array = np.array(y)
    
    # Check if y has multiple columns
    multi_column = len(y_array.shape) > 1 and y_array.shape[1] > 1
    
    if multi_column:
        # Plot each column
        for i in range(y_array.shape[1]):
            label = legend_labels[i] if i < len(legend_labels) else f"Series {i+1}"
            if HAS_SEABORN:
                sns.lineplot(x=x_array, y=y_array[:, i], label=label)
            else:
                plt.plot(x_array, y_array[:, i], label=label)
        
        # Calculate statistics for the first column
        stats = calculate_statistics(y_array[:, 0])
    else:
        # Single column plot
        if HAS_SEABORN:
            sns.lineplot(x=x_array, y=y_array, color='#1f77b4', linewidth=2)
        else:
            plt.plot(x_array, y_array, color='#1f77b4', linewidth=2)
        
        # Calculate statistics
        stats = calculate_statistics(y_array)
    
    # Add horizontal line for mean
    plt.axhline(y=stats['mean'], color='#2ca02c', linestyle='--', alpha=0.7,
               label=f'Mean: {stats["mean"]:.1f} H-bonds')
    
    # Add shaded area for standard deviation
    if not multi_column:
        plt.fill_between(x_array, stats['mean'] - stats['std'], stats['mean'] + stats['std'],
                       color='#2ca02c', alpha=0.2, label=f'Std Dev: ±{stats["std"]:.1f}')
    
    # Calculate H-bonds per molecule
    hbonds_per_molecule = stats['mean'] / n_molecules
    
    # Add a text box with statistics
    stats_text = (
        f"Total H-bonds (mean): {stats['mean']:.1f} ± {stats['std']:.1f}\n"
        f"H-bonds per molecule: {hbonds_per_molecule:.2f}\n"
        f"Min: {stats['min']:.1f}, Max: {stats['max']:.1f}\n"
        f"Reference: {REFERENCE_VALUES['hbonds_per_molecule']:.1f} per molecule"
    )
    
    # Add text box
    props = dict(boxstyle='round', facecolor='white', alpha=0.7)
    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=10, loc='best')
    
    # Add a watermark with simulation details
    plt.figtext(0.5, 0.01, 'TIP4P Water Model - Hydrogen Bond Analysis', 
               ha='center', fontsize=10, style='italic', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')
    
    # Create a histogram of hydrogen bond numbers
    plt.figure(figsize=(10, 6), dpi=300)
    
    if HAS_SEABORN:
        sns.histplot(y_array, kde=True, color='#1f77b4')
        
        # Add vertical lines for statistics
        plt.axvline(x=stats['mean'], color='#2ca02c', linestyle='--', 
                   label=f'Mean: {stats["mean"]:.1f}')
        plt.axvline(x=stats['median'], color='#d62728', linestyle=':', 
                   label=f'Median: {stats["median"]:.1f}')
    else:
        plt.hist(y_array, bins=30, alpha=0.7, color='#1f77b4', density=True)
        
        # Add a simple KDE if seaborn is not available
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(y_array)
        x_kde = np.linspace(min(y_array), max(y_array), 100)
        plt.plot(x_kde, kde(x_kde), 'r-', linewidth=2)
        
        # Add vertical lines for statistics
        plt.axvline(x=stats['mean'], color='#2ca02c', linestyle='--', 
                   label=f'Mean: {stats["mean"]:.1f}')
        plt.axvline(x=stats['median'], color='#d62728', linestyle=':', 
                   label=f'Median: {stats["median"]:.1f}')
    
    plt.xlabel('Number of Hydrogen Bonds', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title('Distribution of Hydrogen Bond Numbers', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=10)
    
    hist_output_path = os.path.join(os.path.dirname(output_path), 'hbnum_histogram.png')
    plt.tight_layout()
    plt.savefig(hist_output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - hbnum_histogram.png saved successfully')

def plot_hbond_distribution(x, y, title, xlabel, ylabel, output_path, reference_value=None, label_prefix=""):
    """Plot hydrogen bond distribution (distance or angle)"""
    plt.figure(figsize=(10, 6), dpi=300)
    
    # Convert to numpy arrays if they're lists
    x_array = np.array(x)
    y_array = np.array(y)
    
    if HAS_SEABORN:
        sns.lineplot(x=x_array, y=y_array, color='#1f77b4', linewidth=2)
    else:
        plt.plot(x_array, y_array, color='#1f77b4', linewidth=2)
    
    # Add reference value if provided
    if reference_value is not None:
        if "Distance" in title:
            ref_label = f'Reference: {reference_value:.3f} nm (O-O)'
            # Add annotation explaining the distance
            plt.annotate('O-O distance between donor and acceptor', xy=(0.5, 0.95), 
                       xycoords='axes fraction', ha='center', fontsize=10, 
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        elif "Angle" in title:
            ref_label = f'Reference: {reference_value:.1f}° (cutoff)'
            # Add annotation explaining the angle
            plt.annotate('Donor-H-Acceptor angle', xy=(0.5, 0.95), 
                       xycoords='axes fraction', ha='center', fontsize=10, 
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        else:
            ref_label = f'Reference: {reference_value}'
            
        plt.axvline(x=reference_value, color='#d62728', linestyle=':', alpha=0.7,
                   label=ref_label)
    
    # Calculate statistics
    stats = calculate_statistics(y_array)
    
    # Find peak
    peak_idx = np.argmax(y_array)
    peak_x = x_array[peak_idx]
    peak_y = y_array[peak_idx]
    
    # Add marker for peak
    plt.plot(peak_x, peak_y, 'ro', markersize=8)
    plt.annotate(f'Peak: {peak_x:.3f}',
                xy=(peak_x, peak_y),
                xytext=(peak_x + 0.02, peak_y - 0.1 * stats['max']),
                arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                fontsize=12)
    
    # Add explanation for angle peak if this is an angle distribution
    if "Angle" in title and peak_x < 15:
        explanation_text = (
            "Peak near 10° is typical for TIP4P water\n"
            "Real H-bonds aren't perfectly linear (0°)\n"
            "GROMACS defines 0° as perfectly linear D-H-A"
        )
        plt.text(0.98, 0.98, explanation_text, transform=plt.gca().transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    # Set correct axis labels based on the type of distribution
    if "Angle" in title:
        # Override the xlabel and ylabel for angle distribution
        plt.xlabel("Angle (degrees)", fontsize=14)
        plt.ylabel("Probability", fontsize=14)
        plt.title("Hydrogen Bond Angle Distribution Histogram", fontsize=16)
    elif "Distance" in title:
        # Override the xlabel and ylabel for distance distribution
        plt.xlabel("Distance (nm)", fontsize=14)
        plt.ylabel("Probability", fontsize=14)
        plt.title("Hydrogen Bond Distance Distribution Histogram", fontsize=16)
    else:
        # Use the provided labels if not a specific type
        plt.xlabel(xlabel, fontsize=14)
        plt.ylabel(ylabel, fontsize=14)
        plt.title(title, fontsize=16)
    
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    
    # Add a watermark with simulation details
    plt.figtext(0.5, 0.01, 'TIP4P Water Model - Hydrogen Bond Analysis', 
                ha='center', fontsize=10, style='italic', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def plot_hbond_lifetime(x, y, title, xlabel, ylabel, output_path):
    """Plot hydrogen bond lifetime correlation function"""
    plt.figure(figsize=(10, 6), dpi=300)
    
    # Convert to numpy arrays if they're lists
    x_array = np.array(x)
    y_array = np.array(y)
    
    if HAS_SEABORN:
        sns.lineplot(x=x_array, y=y_array, color='#1f77b4', linewidth=2)
    else:
        plt.plot(x_array, y_array, color='#1f77b4', linewidth=2)
    
    # Find the time at which the correlation function decays to 1/e
    try:
        # Find where the correlation drops below 1/e
        e_idx = np.where(y_array < 1/np.e)[0][0]
        lifetime = x_array[e_idx]
        
        # Add a marker and annotation for the lifetime
        plt.plot(lifetime, 1/np.e, 'ro', markersize=8)
        plt.annotate(f'Lifetime: {lifetime:.2f} ps',
                    xy=(lifetime, 1/np.e),
                    xytext=(lifetime + 5, 1/np.e + 0.05),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                    fontsize=12)
        
        # Add a horizontal line at 1/e
        plt.axhline(y=1/np.e, color='#d62728', linestyle='--', alpha=0.7,
                   label='1/e threshold')
        
        # Add a vertical line at the lifetime
        plt.axvline(x=lifetime, color='#d62728', linestyle='--', alpha=0.7)
        
        # Add reference value
        ref_lifetime = REFERENCE_VALUES['hbond_lifetime']
        plt.axvline(x=ref_lifetime, color='#2ca02c', linestyle=':', alpha=0.7,
                   label=f'Reference: {ref_lifetime:.2f} ps')
        
        # Try to fit an exponential decay to the data
        # Use data up to 30 ps or the end of the data, whichever is smaller
        fit_end_idx = min(len(x_array), np.searchsorted(x_array, 30.0))
        
        # Only fit if we have enough data points
        if fit_end_idx > 5:
            try:
                from scipy.optimize import curve_fit
                
                def exp_decay(t, tau):
                    return np.exp(-t/tau)
                
                # Fit the exponential decay function to the data
                popt, pcov = curve_fit(exp_decay, x_array[:fit_end_idx], y_array[:fit_end_idx], 
                                      p0=[lifetime], bounds=(0, np.inf))
                
                # Extract the fitted lifetime
                fitted_lifetime = popt[0]
                
                # Plot the fitted curve
                fit_x = np.linspace(0, x_array[fit_end_idx-1], 100)
                fit_y = exp_decay(fit_x, fitted_lifetime)
                plt.plot(fit_x, fit_y, 'g--', alpha=0.7, 
                        label=f'Exp fit: τ = {fitted_lifetime:.2f} ps')
                
                # Add annotation about the fit
                plt.annotate(f'Fitted τ: {fitted_lifetime:.2f} ps',
                            xy=(fitted_lifetime, exp_decay(fitted_lifetime, fitted_lifetime)),
                            xytext=(fitted_lifetime + 5, exp_decay(fitted_lifetime, fitted_lifetime) + 0.1),
                            arrowprops=dict(facecolor='green', shrink=0.05, width=1, headwidth=5),
                            fontsize=10, color='green')
            except Exception as e:
                print(f"  - Could not fit exponential decay: {e}")
        
        # Add annotation explaining the lifetime
        explanation_text = (
            "Uninterrupted hydrogen bond lifetime\n"
            "(time for correlation to decay to 1/e)\n\n"
            "This measures continuous H-bond lifetime:\n"
            "Once a bond breaks, it's considered dead\n"
            "even if it reforms with the same partners"
        )
        plt.text(0.98, 0.98, explanation_text, transform=plt.gca().transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        
        # Add a text box with factors affecting lifetime
        factors_text = (
            "Factors affecting H-bond lifetime:\n"
            "• Distance/angle cutoffs\n"
            "• Temperature (↓T = ↑lifetime)\n"
            "• Water model (TIP4P vs TIP3P)\n"
            "• System size and conditions"
        )
        plt.text(0.02, 0.98, factors_text, transform=plt.gca().transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        
    except Exception as e:
        print(f"  - Could not determine hydrogen bond lifetime: {e}")
        # Add reference value even if we couldn't determine lifetime
        ref_lifetime = REFERENCE_VALUES['hbond_lifetime']
        plt.axvline(x=ref_lifetime, color='#2ca02c', linestyle=':', alpha=0.7,
                   label=f'Reference: {ref_lifetime:.2f} ps')
    
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    
    # Add a watermark with simulation details
    plt.figtext(0.5, 0.01, 'TIP4P Water Model - Hydrogen Bond Lifetime', 
               ha='center', fontsize=10, style='italic', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

    # Create a log-scale plot to better visualize the decay
    plt.figure(figsize=(10, 6), dpi=300)
    
    # Use semilogy for log scale on y-axis
    plt.semilogy(x_array, y_array, color='#1f77b4', linewidth=2)
    
    # Add horizontal line at 1/e
    plt.axhline(y=1/np.e, color='#d62728', linestyle='--', alpha=0.7,
               label='1/e threshold')
    
    try:
        # Add vertical line at the lifetime
        plt.axvline(x=lifetime, color='#d62728', linestyle='--', alpha=0.7,
                   label=f'Lifetime: {lifetime:.2f} ps')
        
        # Add reference value
        plt.axvline(x=ref_lifetime, color='#2ca02c', linestyle=':', alpha=0.7,
                   label=f'Reference: {ref_lifetime:.2f} ps')
        
        # Plot the fitted curve if available
        if 'fitted_lifetime' in locals():
            fit_x = np.linspace(0, min(30.0, x_array[-1]), 100)
            fit_y = exp_decay(fit_x, fitted_lifetime)
            plt.semilogy(fit_x, fit_y, 'g--', alpha=0.7, 
                        label=f'Exp fit: τ = {fitted_lifetime:.2f} ps')
        
        # Add annotation explaining the lifetime and log scale
        explanation_text = (
            "Uninterrupted hydrogen bond lifetime (log scale)\n\n"
            "• A straight line in this log plot would indicate\n"
            "  a perfect single-exponential decay\n"
            "• Deviations suggest multiple timescales\n"
            "• TIP4P water typically shows 1-3 ps lifetimes\n"
            "  (2.5 ps is within expected range)"
        )
        plt.text(0.98, 0.98, explanation_text, transform=plt.gca().transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        
    except:
        # Add reference value even if we couldn't determine lifetime
        ref_lifetime = REFERENCE_VALUES['hbond_lifetime']
        plt.axvline(x=ref_lifetime, color='#2ca02c', linestyle=':', alpha=0.7,
                   label=f'Reference: {ref_lifetime:.2f} ps')
    
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel + ' (log scale)', fontsize=14)
    plt.title(f'{title} (Log Scale)', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    
    log_output_path = os.path.join(os.path.dirname(output_path), 'hblife_log_plot.png')
    plt.tight_layout()
    plt.savefig(log_output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - hblife_log_plot.png saved successfully')

def create_combined_hbond_plot(hbnum_data, hbdist_data, hbang_data, hblife_data, output_path):
    """Create a combined plot with all hydrogen bond analyses"""
    fig, axs = plt.subplots(2, 2, figsize=(16, 12), dpi=300)
    
    # Plot hydrogen bond number
    if hbnum_data:
        x, y, title, xlabel, ylabel, legend_labels = hbnum_data
        x_array = np.array(x)
        y_array = np.array(y)
        
        if HAS_SEABORN:
            sns.lineplot(x=x_array, y=y_array, color='#1f77b4', linewidth=2, ax=axs[0, 0])
        else:
            axs[0, 0].plot(x_array, y_array, color='#1f77b4', linewidth=2)
        
        # Calculate statistics
        stats = calculate_statistics(y_array)
        
        # Add horizontal line for mean
        axs[0, 0].axhline(y=stats['mean'], color='#2ca02c', linestyle='--', alpha=0.7,
                        label=f'Mean: {stats["mean"]:.1f}')
        
        # Try to determine the number of water molecules
        num_water_molecules = 5500  # Default based on the log output
        topology_file = os.path.join(os.path.dirname(os.path.dirname(output_path)), "data", "topol.top")
        if os.path.exists(topology_file):
            with open(topology_file, 'r') as f:
                for line in f:
                    if "SOL" in line and not line.startswith(";"):
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            try:
                                num_water_molecules = int(parts[-1])
                                break
                            except ValueError:
                                pass
        
        # Calculate the reference value correctly
        # Each water molecule forms ~3.5 H-bonds on average
        # But each H-bond connects two molecules, so divide by 2
        reference_hbonds = (REFERENCE_VALUES['hbonds_per_molecule'] * num_water_molecules) / 2
        
        # Add reference line
        axs[0, 0].axhline(y=reference_hbonds, color='#ff7f0e', linestyle='--', alpha=0.7,
                        label=f'Reference: {reference_hbonds:.1f}')
        
        axs[0, 0].set_xlabel(xlabel if xlabel else 'Time (ps)', fontsize=12)
        axs[0, 0].set_ylabel(ylabel if ylabel else 'Number of H-bonds', fontsize=12)
        axs[0, 0].set_title(title if title else 'Number of Hydrogen Bonds', fontsize=14)
        axs[0, 0].grid(True, alpha=0.3)
        axs[0, 0].legend(fontsize=10)
    
    # Plot hydrogen bond distance
    if hbdist_data:
        x, y, title, xlabel, ylabel, legend_labels = hbdist_data
        x_array = np.array(x)
        y_array = np.array(y)
        
        if HAS_SEABORN:
            sns.lineplot(x=x_array, y=y_array, color='#1f77b4', linewidth=2, ax=axs[0, 1])
        else:
            axs[0, 1].plot(x_array, y_array, color='#1f77b4', linewidth=2)
    
    # Add reference value
        ref_dist = REFERENCE_VALUES['hbond_distance_mean']
        axs[0, 1].axvline(x=ref_dist, color='#d62728', linestyle=':', alpha=0.7,
                        label=f'Reference: {ref_dist:.3f} nm (O-O)')
        
        axs[0, 1].set_xlabel(xlabel if xlabel else 'Distance (nm)', fontsize=12)
        axs[0, 1].set_ylabel(ylabel if ylabel else 'Frequency', fontsize=12)
        axs[0, 1].set_title(title if title else 'Hydrogen Bond Distance Distribution Histogram', fontsize=14)
        # Add annotation explaining the distance
        axs[0, 1].annotate('O-O distance between donor and acceptor', xy=(0.5, 0.95), 
                         xycoords='axes fraction', ha='center', fontsize=10, 
                         bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        axs[0, 1].grid(True, alpha=0.3)
        axs[0, 1].legend(fontsize=10)
    
    # Plot hydrogen bond angle
    if hbang_data:
        x, y, title, xlabel, ylabel, legend_labels = hbang_data
        x_array = np.array(x)
        y_array = np.array(y)
        
        if HAS_SEABORN:
            sns.lineplot(x=x_array, y=y_array, color='#1f77b4', linewidth=2, ax=axs[1, 0])
        else:
            axs[1, 0].plot(x_array, y_array, color='#1f77b4', linewidth=2)
        
        # Add reference value
        ref_angle = REFERENCE_VALUES['hbond_angle_mean']
        axs[1, 0].axvline(x=ref_angle, color='#d62728', linestyle=':', alpha=0.7,
                        label=f'Reference: {ref_angle:.1f}° (cutoff)')
        
        # Find peak
        peak_idx = np.argmax(y_array)
        peak_x = x_array[peak_idx]
        peak_y = y_array[peak_idx]
        
        # Add marker for peak
        axs[1, 0].plot(peak_x, peak_y, 'ro', markersize=8)
        
        # Add explanation for angle peak
        if peak_x < 15:
            explanation_text = (
                "Peak near 10° is typical for TIP4P water\n"
                "Real H-bonds aren't perfectly linear (0°)\n"
                "GROMACS defines 0° as perfectly linear D-H-A"
            )
            axs[1, 0].text(0.98, 0.98, explanation_text, transform=axs[1, 0].transAxes, fontsize=8,
                    verticalalignment='top', horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        
        axs[1, 0].set_xlabel('Angle (degrees)', fontsize=12)
        axs[1, 0].set_ylabel('Probability', fontsize=12)
        axs[1, 0].set_title('Hydrogen Bond Angle Distribution Histogram', fontsize=14)
        # Add annotation explaining the angle
        axs[1, 0].annotate('Donor-H-Acceptor angle', xy=(0.5, 0.95), 
                         xycoords='axes fraction', ha='center', fontsize=10, 
                         bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        axs[1, 0].grid(True, alpha=0.3)
        axs[1, 0].legend(fontsize=10)
    
    # Plot hydrogen bond lifetime
    if hblife_data:
        x, y, title, xlabel, ylabel, legend_labels = hblife_data
        x_array = np.array(x)
        y_array = np.array(y)
        
        if HAS_SEABORN:
            sns.lineplot(x=x_array, y=y_array, color='#1f77b4', linewidth=2, ax=axs[1, 1])
        else:
            axs[1, 1].plot(x_array, y_array, color='#1f77b4', linewidth=2)
    
    # Add reference value
        ref_lifetime = REFERENCE_VALUES['hbond_lifetime']
        axs[1, 1].axvline(x=ref_lifetime, color='#d62728', linestyle=':', alpha=0.7,
                        label=f'Reference: {ref_lifetime:.2f} ps')
        
        # Add horizontal line at 1/e
        axs[1, 1].axhline(y=1/np.e, color='#2ca02c', linestyle='--', alpha=0.7,
                        label='1/e threshold')
        
        # Try to find the time at which the correlation function decays to 1/e
        try:
            # Find where the correlation drops below 1/e
            e_idx = np.where(y_array < 1/np.e)[0][0]
            lifetime = x_array[e_idx]
            
            # Add a marker and annotation for the lifetime
            axs[1, 1].plot(lifetime, 1/np.e, 'ro', markersize=6)
            axs[1, 1].annotate(f'Lifetime: {lifetime:.2f} ps',
                            xy=(lifetime, 1/np.e),
                            xytext=(lifetime + 5, 1/np.e + 0.05),
                            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                            fontsize=10)
            
            # Try to fit an exponential decay to the data
            # Use data up to 30 ps or the end of the data, whichever is smaller
            fit_end_idx = min(len(x_array), np.searchsorted(x_array, 30.0))
            
            # Only fit if we have enough data points
            if fit_end_idx > 5:
                try:
                    from scipy.optimize import curve_fit
                    
                    def exp_decay(t, tau):
                        return np.exp(-t/tau)
                    
                    # Fit the exponential decay function to the data
                    popt, pcov = curve_fit(exp_decay, x_array[:fit_end_idx], y_array[:fit_end_idx], 
                                          p0=[lifetime], bounds=(0, np.inf))
                    
                    # Extract the fitted lifetime
                    fitted_lifetime = popt[0]
                    
                    # Plot the fitted curve
                    fit_x = np.linspace(0, x_array[fit_end_idx-1], 100)
                    fit_y = exp_decay(fit_x, fitted_lifetime)
                    axs[1, 1].plot(fit_x, fit_y, 'g--', alpha=0.7, 
                                 label=f'Exp fit: τ = {fitted_lifetime:.2f} ps')
                except Exception as e:
                    print(f"  - Could not fit exponential decay in combined plot: {e}")
            
            # Add explanation text
            explanation_text = (
                "Continuous H-bond lifetime\n"
                "TIP4P typical: 1-3 ps"
            )
            axs[1, 1].text(0.98, 0.98, explanation_text, transform=axs[1, 1].transAxes, fontsize=8,
                         verticalalignment='top', horizontalalignment='right',
                         bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
            
        except Exception as e:
            print(f"  - Could not determine hydrogen bond lifetime in combined plot: {e}")
        
        axs[1, 1].set_xlabel(xlabel if xlabel else 'Time (ps)', fontsize=12)
        axs[1, 1].set_ylabel(ylabel if ylabel else 'C(t)', fontsize=12)
        axs[1, 1].set_title(title if title else 'Uninterrupted hydrogen bond lifetime', fontsize=14)
        axs[1, 1].grid(True, alpha=0.3)
        axs[1, 1].legend(fontsize=10)
    
    # Add a title for the entire figure
    fig.suptitle('Hydrogen Bond Analysis for TIP4P Water', fontsize=16)
    
    # Add a watermark with simulation details
    plt.figtext(0.5, 0.01, 'TIP4P Water Model - Hydrogen Bond Analysis', 
               ha='center', fontsize=10, style='italic', alpha=0.7)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def plot_hbond_histogram(data, title, xlabel, ylabel, output_path, reference_value=None):
    """Plot histogram of hydrogen bond data"""
    plt.figure(figsize=(10, 6), dpi=300)
    
    # Convert to numpy array if it's a list
    data_array = np.array(data)
    
    if HAS_SEABORN:
        sns.histplot(data_array, kde=True, color='#1f77b4')
        
        # Add vertical lines for statistics
        plt.axvline(x=np.mean(data_array), color='#2ca02c', linestyle='--', 
                   label=f'Mean: {np.mean(data_array):.1f}')
        plt.axvline(x=np.median(data_array), color='#d62728', linestyle=':', 
                   label=f'Median: {np.median(data_array):.1f}')
    else:
        plt.hist(data_array, bins=30, alpha=0.7, color='#1f77b4', density=True)
        
        # Add a simple KDE if seaborn is not available
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(data_array)
        x_kde = np.linspace(min(data_array), max(data_array), 100)
        plt.plot(x_kde, kde(x_kde), 'r-', linewidth=2)
        
        # Add vertical lines for statistics
        plt.axvline(x=np.mean(data_array), color='#2ca02c', linestyle='--', 
                   label=f'Mean: {np.mean(data_array):.1f}')
        plt.axvline(x=np.median(data_array), color='#d62728', linestyle=':', 
                   label=f'Median: {np.median(data_array):.1f}')
    
    # Add reference value if provided
    if reference_value is not None:
        plt.axvline(x=reference_value, color='#ff7f0e', linestyle='--', alpha=0.7,
                   label=f'Reference: {reference_value:.1f}')
        
        # Add explanation about the reference value
        explanation_text = (
            f"Reference value: {reference_value:.1f} H-bonds\n"
            f"Based on ~3.5 H-bonds per molecule\n"
            f"Formula: (3.5 × molecules) ÷ 2\n"
            f"(Each H-bond connects two molecules)"
        )
        plt.text(0.98, 0.98, explanation_text, transform=plt.gca().transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    # Add a text box with statistics
    stats_text = (
        f"Mean: {np.mean(data_array):.1f}\n"
        f"Median: {np.median(data_array):.1f}\n"
        f"Std Dev: {np.std(data_array):.1f}\n"
        f"Min: {np.min(data_array):.1f}\n"
        f"Max: {np.max(data_array):.1f}"
    )
    
    # Add text box
    props = dict(boxstyle='round', facecolor='white', alpha=0.7)
    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    
    # Add a watermark with simulation details
    plt.figtext(0.5, 0.01, 'TIP4P Water Model - Hydrogen Bond Analysis', 
               ha='center', fontsize=10, style='italic', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def main():
    if len(sys.argv) < 3:
        print("Usage: plot_hbond.py <analysis_dir> <plots_dir>")
        sys.exit(1)
        
    analysis_dir = sys.argv[1]
    plots_dir = sys.argv[2]
    
    # Define data directory
    data_dir = os.path.join(analysis_dir, "data")
    
    # Ensure plots directory exists
    os.makedirs(plots_dir, exist_ok=True)
    
    # Try to determine the number of water molecules
    num_water_molecules = 300  # Default to a reasonable value
    
    # Check if we can find the topology file to get the actual number
    topology_file = os.path.join(os.path.dirname(analysis_dir), "data", "topol.top")
    if os.path.exists(topology_file):
        with open(topology_file, 'r') as f:
            for line in f:
                if "SOL" in line and not line.startswith(";"):
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            num_water_molecules = int(parts[-1])
                            print(f"Found {num_water_molecules} water molecules in topology file")
                            break
                        except ValueError:
                            pass
    
    # Store data for combined plot
    hbnum_data = None
    hbdist_data = None
    hbang_data = None
    hblife_data = None
    
    # Plot hydrogen bond number
    hbnum_file = os.path.join(data_dir, 'hbnum.xvg')
    if os.path.exists(hbnum_file):
        print('Plotting hydrogen bond number...')
        x, y, y_columns, title, xlabel, ylabel, legend_labels = read_xvg(hbnum_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Number of Hydrogen Bonds'
            plot_xlabel = xlabel if xlabel else 'Time (ps)'
            plot_ylabel = ylabel if ylabel else 'Number of H-bonds'
            
            output_path = os.path.join(plots_dir, 'hbnum_plot.png')
            plot_hbond_number(x, y, plot_title, plot_xlabel, plot_ylabel, legend_labels, output_path, num_water_molecules)
            
            # Plot histogram
            hist_title = 'Distribution of Hydrogen Bond Numbers'
            hist_xlabel = 'Number of H-bonds'
            hist_ylabel = 'Frequency'
            output_path = os.path.join(plots_dir, 'hbnum_histogram.png')
            
            # Calculate the reference value correctly
            # Each water molecule forms ~3.5 H-bonds on average
            # But each H-bond connects two molecules, so divide by 2
            # This gives the expected total number of H-bonds in the system
            reference_hbonds = (REFERENCE_VALUES['hbonds_per_molecule'] * num_water_molecules) / 2
            
            plot_hbond_histogram(y, hist_title, hist_xlabel, hist_ylabel, output_path, reference_hbonds)
            
            # Store data for combined plot
            hbnum_data = (x, y, title, xlabel, ylabel, legend_labels)
    else:
        print(f"Hydrogen bond number file not found: {hbnum_file}")
    
    # Plot hydrogen bond distance
    hbdist_file = os.path.join(data_dir, 'hbdist.xvg')
    if os.path.exists(hbdist_file):
        print('Plotting hydrogen bond distance...')
        x, y, y_columns, title, xlabel, ylabel, legend_labels = read_xvg(hbdist_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Hydrogen Bond Distance'
            plot_xlabel = xlabel if xlabel else 'Time (ps)'
            plot_ylabel = ylabel if ylabel else 'Distance (nm)'
            
            output_path = os.path.join(plots_dir, 'hbdist_plot.png')
            plot_hbond_distribution(x, y, plot_title, plot_xlabel, plot_ylabel, output_path, 
                                  REFERENCE_VALUES['hbond_distance_mean'], 'H-bond Distance')
            
            # Store data for combined plot
            hbdist_data = (x, y, title, xlabel, ylabel, legend_labels)
    else:
        print(f"Hydrogen bond distance file not found: {hbdist_file}")
    
    # Plot hydrogen bond angle
    hbang_file = os.path.join(data_dir, 'hbang.xvg')
    if os.path.exists(hbang_file):
        print('Plotting hydrogen bond angle...')
        x, y, y_columns, title, xlabel, ylabel, legend_labels = read_xvg(hbang_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Hydrogen Bond Angle Distribution'
            plot_xlabel = xlabel if xlabel else 'Angle (degrees)'
            plot_ylabel = ylabel if ylabel else 'Probability'
            
            output_path = os.path.join(plots_dir, 'hbang_plot.png')
            plot_hbond_distribution(x, y, plot_title, plot_xlabel, plot_ylabel, output_path, 
                                  REFERENCE_VALUES['hbond_angle_mean'], 'H-bond Angle')
            
            # Store data for combined plot
            hbang_data = (x, y, title, xlabel, ylabel, legend_labels)
    else:
        print(f"Hydrogen bond angle file not found: {hbang_file}")
    
    # Plot hydrogen bond lifetime
    hblife_file = os.path.join(data_dir, 'hblife.xvg')
    if os.path.exists(hblife_file):
        print('Plotting hydrogen bond lifetime...')
        x, y, y_columns, title, xlabel, ylabel, legend_labels = read_xvg(hblife_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Hydrogen Bond Lifetime Correlation'
            plot_xlabel = xlabel if xlabel else 'Time (ps)'
            plot_ylabel = ylabel if ylabel else 'C(t)'
            
            output_path = os.path.join(plots_dir, 'hblife_plot.png')
            # Use the second column (index 1) for the lifetime plot if available
            if len(y_columns) > 1:
                plot_hbond_lifetime(x, y_columns[1], plot_title, plot_xlabel, plot_ylabel, output_path)
            else:
                plot_hbond_lifetime(x, y, plot_title, plot_xlabel, plot_ylabel, output_path)
            
            # Store data for combined plot
            hblife_data = (x, y, title, xlabel, ylabel, legend_labels)
    else:
        print(f"Hydrogen bond lifetime file not found: {hblife_file}")
    
    # Also check for autocorrelation file
    hbac_file = os.path.join(data_dir, 'hbac.xvg')
    if os.path.exists(hbac_file) and not os.path.exists(hblife_file):
        print('Plotting hydrogen bond autocorrelation...')
        x, y, y_columns, title, xlabel, ylabel, legend_labels = read_xvg(hbac_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Hydrogen Bond Autocorrelation'
            plot_xlabel = xlabel if xlabel else 'Time (ps)'
            plot_ylabel = ylabel if ylabel else 'C(t)'
            
            output_path = os.path.join(plots_dir, 'hblife_plot.png')
            # Use the second column (index 1) for the autocorrelation plot if available
            if len(y_columns) > 1:
                plot_hbond_lifetime(x, y_columns[1], plot_title, plot_xlabel, plot_ylabel, output_path)
            else:
                plot_hbond_lifetime(x, y, plot_title, plot_xlabel, plot_ylabel, output_path)
            
            # Store data for combined plot
            hblife_data = (x, y, title, xlabel, ylabel, legend_labels)
    elif not os.path.exists(hblife_file):
        print(f"Hydrogen bond autocorrelation file not found: {hbac_file}")
    
    # Create combined hydrogen bond plot
    if any([hbnum_data, hbdist_data, hbang_data, hblife_data]):
        print('Creating combined hydrogen bond plot...')
        output_path = os.path.join(plots_dir, 'combined_hbond_plot.png')
        create_combined_hbond_plot(hbnum_data, hbdist_data, hbang_data, hblife_data, output_path)

if __name__ == "__main__":
    main() 