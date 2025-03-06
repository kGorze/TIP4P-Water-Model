#!/usr/bin/python3
import os
import sys
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import warnings
from scipy import stats
warnings.filterwarnings('ignore')
try:
    import seaborn as sns
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.5)
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    print("Seaborn not available, using matplotlib defaults")

# Reference values for TIP4P water at 298K, 1 bar
REFERENCE_VALUES = {
    'temperature': 298.0,  # K
    'pressure': 1.0,       # bar
    'density': 997.0,      # kg/m^3
    'potential_energy': -41.5,  # kJ/mol per molecule
    'potential_energy_per_molecule': -41.5,  # kJ/mol per molecule
}

# Number of water molecules in the simulation
NUM_WATER_MOLECULES = 5500

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

def calculate_statistics(data):
    """Calculate basic statistics for the data"""
    if len(data.shape) > 1:
        # Multiple columns
        stats = []
        for i in range(data.shape[1]):
            col_data = data[:, i]
            stats.append({
                'mean': np.mean(col_data),
                'std': np.std(col_data),
                'min': np.min(col_data),
                'max': np.max(col_data),
                'final': col_data[-1] if len(col_data) > 0 else 0
            })
        return stats
    else:
        # Single column
        return {
            'mean': np.mean(data),
            'std': np.std(data),
            'min': np.min(data),
            'max': np.max(data),
            'final': data[-1] if len(data) > 0 else 0
        }

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

def plot_time_series(x, y, title, xlabel, ylabel, legend_labels, output_path, reference_value=None, property_name=None):
    """Create a time series plot with statistics and annotations"""
    # Special handling for pressure plots
    if property_name == 'Pressure':
        plot_pressure_time_series(x, y, title, xlabel, ylabel, legend_labels, output_path, reference_value)
        return
    
    # Special handling for potential energy plots
    if property_name == 'Potential Energy':
        plot_potential_energy_time_series(x, y, title, xlabel, ylabel, legend_labels, output_path, reference_value)
        return
    
    plt.figure(figsize=(12, 7), dpi=300)
    
    # Check if y has multiple columns
    multi_column = len(y.shape) > 1 and y.shape[1] > 1
    
    if multi_column:
        # Plot each column
        for i in range(y.shape[1]):
            label = legend_labels[i] if i < len(legend_labels) else f"Series {i+1}"
            if HAS_SEABORN:
                sns.lineplot(x=x, y=y[:, i], label=label)
            else:
                plt.plot(x, y[:, i], label=label)
        
        # Calculate statistics for the first column (usually the main property)
        stats = calculate_statistics(y[:, 0])
    else:
        # Single column plot
        if HAS_SEABORN:
            sns.lineplot(x=x, y=y, color='#1f77b4', linewidth=2)
        else:
            plt.plot(x, y, color='#1f77b4', linewidth=2)
        
        # Calculate statistics
        stats = calculate_statistics(y)
    
    # Add horizontal line for mean
    plt.axhline(y=stats['mean'], color='#2ca02c', linestyle='--', alpha=0.7,
               label=f'Mean: {stats["mean"]:.2f}')
    
    # Add reference value if provided
    if reference_value is not None:
        # Standard reference line for other properties
        plt.axhline(y=reference_value, color='#d62728', linestyle=':', alpha=0.7,
                   label=f'Reference: {reference_value:.2f}')
        
        # Add percent difference annotation
        percent_diff = ((stats['mean'] - reference_value) / reference_value) * 100
        plt.annotate(f'Diff from ref: {percent_diff:.2f}%',
                    xy=(x[-1], reference_value),
                    xytext=(x[-1] - (x[-1]-x[0])*0.2, reference_value + stats['std']),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                    fontsize=10)
    
    # Add shaded area for standard deviation
    if not multi_column:
        plt.fill_between(x, stats['mean'] - stats['std'], stats['mean'] + stats['std'],
                       color='#2ca02c', alpha=0.2, label=f'Std Dev: ±{stats["std"]:.2f}')
    
    # Add a text box with statistics
    stats_text = (
        f"Mean: {stats['mean']:.2f}\n"
        f"Std Dev: {stats['std']:.2f}\n"
        f"Min: {stats['min']:.2f}\n"
        f"Max: {stats['max']:.2f}\n"
        f"Final: {stats['final']:.2f}"
    )
    
    # Add text box
    props = dict(boxstyle='round', facecolor='white', alpha=0.7)
    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=16)
    plt.grid(True, alpha=0.3)
    
    # Add legend if we have multiple series or reference values
    if multi_column or reference_value is not None:
        plt.legend(fontsize=10, loc='best')
    
    # Add a watermark with simulation details
    if property_name:
        plt.figtext(0.5, 0.01, f'TIP4P Water Model - {property_name}', 
                   ha='center', fontsize=10, style='italic', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def plot_pressure_time_series(x, y, title, xlabel, ylabel, legend_labels, output_path, reference_value):
    """Create an enhanced pressure time series plot with block averaging and context about MD pressure fluctuations"""
    # Create a figure with two subplots
    fig = plt.figure(figsize=(15, 12), dpi=300)
    gs = plt.GridSpec(3, 1, figure=fig, height_ratios=[2, 1, 1])
    
    # Main time series plot
    ax1 = fig.add_subplot(gs[0])
    
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
        sns.lineplot(x=x, y=data, color='#1f77b4', linewidth=1.0, alpha=0.7, ax=ax1)
    else:
        ax1.plot(x, data, color='#1f77b4', linewidth=1.0, alpha=0.7)
    
    # Add mean line and standard deviation band
    ax1.axhline(y=mean_val, color='#2ca02c', linestyle='--', alpha=0.7, label=f'Mean: {mean_val:.2f} bar')
    ax1.fill_between(x, mean_val - std_val, mean_val + std_val, color='#2ca02c', alpha=0.2, 
                    label=f'Std Dev: ±{std_val:.2f} bar')
    
    # Add reference line
    ax1.axhline(y=reference_value, color='#d62728', linestyle=':', alpha=0.7, 
               label=f'Reference: {reference_value:.1f} bar')
    
    # Add statistics text box
    stats_text = (
        f"Mean: {mean_val:.2f} bar\n"
        f"Std Dev: {std_val:.2f} bar\n"
        f"Min: {min_val:.2f} bar\n"
        f"Max: {max_val:.2f} bar\n"
        f"Final: {final_val:.2f} bar\n"
        f"Absolute Diff from Ref: {abs(mean_val - reference_value):.2f} bar"
    )
    ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    ax1.set_xlabel(xlabel, fontsize=12)
    ax1.set_ylabel(ylabel, fontsize=12)
    ax1.set_title(f"{title} - Full Time Series", fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10, loc='lower right')
    
    # Block averaging plot
    ax2 = fig.add_subplot(gs[1])
    
    # Calculate block averages
    num_blocks = 10  # Use 10 blocks for better resolution
    block_centers, block_means, block_stds = calculate_block_averages(x, data, num_blocks)
    
    # Calculate drift metrics
    slope, percent_change, significant_drift = calculate_drift(block_means)
    
    # Plot block averages with error bars
    ax2.errorbar(block_centers, block_means, yerr=block_stds, fmt='o-', color='#d62728', 
                ecolor='#d62728', elinewidth=1, capsize=4, label='Block Averages')
    
    # Add horizontal line for overall mean
    ax2.axhline(y=mean_val, color='#2ca02c', linestyle='--', alpha=0.7, label=f'Overall Mean: {mean_val:.2f} bar')
    
    # Add reference line
    ax2.axhline(y=reference_value, color='#d62728', linestyle=':', alpha=0.7, 
               label=f'Reference: {reference_value:.1f} bar')
    
    # Add drift information
    drift_text = (
        f"Block Analysis (n={num_blocks}):\n"
        f"First Block Mean: {block_means[0]:.2f} bar\n"
        f"Last Block Mean: {block_means[-1]:.2f} bar\n"
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
    
    # Context information
    ax3 = fig.add_subplot(gs[2])
    ax3.axis('off')  # Turn off axis
    
    # Add context about pressure fluctuations in MD simulations
    context_text = (
        "INTERPRETING PRESSURE IN MD SIMULATIONS:\n\n"
        "• Large instantaneous pressure fluctuations (±100-200 bar) are NORMAL in molecular dynamics simulations\n"
        "• For water at 273K, the TIP4P model may not reproduce exactly 1 bar pressure\n"
        "• Standard deviations of tens of bars are expected due to the microscopic system size\n"
        "• The absolute difference from the reference (1 bar) is more meaningful than percentage difference\n"
        "• Block averages show if the pressure has converged over the simulation time\n"
        "• If block averages are stable, the system is likely well-equilibrated despite fluctuations"
    )
    
    ax3.text(0.5, 0.5, context_text, transform=ax3.transAxes, fontsize=11,
            horizontalalignment='center', verticalalignment='center', 
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))
    
    # Add a title for the entire figure
    fig.suptitle('Pressure Analysis for TIP4P Water Simulation', fontsize=16)
    
    plt.tight_layout(rect=[0, 0, 1, 0.97])  # Adjust for the suptitle
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def plot_potential_energy_time_series(x, y, title, xlabel, ylabel, legend_labels, output_path, reference_value):
    """Create an enhanced potential energy plot with proper scaling and unit handling"""
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 7), dpi=300)
    
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
    
    # Calculate per-molecule values
    per_molecule_mean = mean_val / NUM_WATER_MOLECULES
    
    # Plot the time series with enhanced styling
    if HAS_SEABORN:
        sns.lineplot(x=x, y=data, color='#1f77b4', linewidth=2.5, ax=ax)
    else:
        ax.plot(x, data, color='#1f77b4', linewidth=2.5)
    
    # Add mean line and standard deviation band
    ax.axhline(y=mean_val, color='#2ca02c', linestyle='--', alpha=0.7, 
              label=f'Mean: {mean_val:.2f} kJ/mol')
    ax.fill_between(x, mean_val - std_val, mean_val + std_val, color='#2ca02c', alpha=0.2, 
                   label=f'Std Dev: ±{std_val:.2f} kJ/mol')
    
    # Determine if the data is likely total system energy or per-molecule energy
    # If mean value is close to reference_value, it's likely per-molecule
    # If mean value is close to reference_value * NUM_WATER_MOLECULES, it's likely total
    is_total_energy = abs(mean_val) > 1000  # Arbitrary threshold, adjust as needed
    
    # Calculate the appropriate reference value and percent difference
    if is_total_energy:
        # Data is total system energy
        total_reference = reference_value * NUM_WATER_MOLECULES
        percent_diff = ((mean_val - total_reference) / abs(total_reference)) * 100
        
        # Add reference line for total energy
        ax.axhline(y=total_reference, color='#d62728', linestyle=':', alpha=0.7,
                  label=f'Reference (total): {total_reference:.2f} kJ/mol')
    else:
        # Data is per-molecule energy
        percent_diff = ((per_molecule_mean - reference_value) / abs(reference_value)) * 100
        
        # Add reference line for per-molecule energy
        ax.axhline(y=reference_value * NUM_WATER_MOLECULES, color='#d62728', linestyle=':', alpha=0.7,
                  label=f'Reference (scaled): {reference_value * NUM_WATER_MOLECULES:.2f} kJ/mol')
    
    # Add explanation text box
    explanation_text = (
        f"Reference value ({reference_value:.1f} kJ/mol) is per water molecule\n"
        f"For {NUM_WATER_MOLECULES} molecules, total reference = {reference_value * NUM_WATER_MOLECULES:.2f} kJ/mol\n"
        f"Observed mean = {mean_val:.2f} kJ/mol (total)\n"
        f"Observed per molecule = {per_molecule_mean:.2f} kJ/mol\n"
        f"Difference: {percent_diff:.2f}%"
    )
    
    # Add text box for explanation
    props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.7)
    ax.text(0.5, 0.25, explanation_text, transform=ax.transAxes, fontsize=10,
           horizontalalignment='center', verticalalignment='center', bbox=props)
    
    # Add statistics text box
    stats_text = (
        f"Mean: {mean_val:.2f} kJ/mol\n"
        f"Std Dev: {std_val:.2f} kJ/mol\n"
        f"Min: {min_val:.2f} kJ/mol\n"
        f"Max: {max_val:.2f} kJ/mol\n"
        f"Final: {final_val:.2f} kJ/mol"
    )
    
    # Add text box
    props = dict(boxstyle='round', facecolor='white', alpha=0.7)
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', bbox=props)
    
    # Set axis labels and title
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=16)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Set y-axis limits to focus on the data range with some padding
    data_range = max_val - min_val
    ax.set_ylim(min_val - 0.1 * data_range, max_val + 0.1 * data_range)
    
    # Add legend
    ax.legend(fontsize=10, loc='best')
    
    # Add a watermark with simulation details
    fig.text(0.5, 0.01, f'TIP4P Water Model - Potential Energy', 
            ha='center', fontsize=10, style='italic', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def plot_energy_components(x, y, legend_labels, output_path):
    """Create a specialized plot for energy components with enhanced visualization"""
    if len(y.shape) <= 1 or y.shape[1] <= 1:
        print("  - Not enough energy components to plot")
        return
    
    # Create a figure with two subplots - individual components and energy drift
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12), dpi=300, gridspec_kw={'height_ratios': [2, 1]})
    
    # Define a distinct color palette with high contrast
    if HAS_SEABORN:
        colors = sns.color_palette("tab10", n_colors=y.shape[1])
    else:
        colors = plt.cm.tab10.colors
    
    # Calculate statistics for each component
    stats_data = []
    for i in range(y.shape[1]):
        component_data = y[:, i]
        mean_val = np.mean(component_data)
        std_val = np.std(component_data)
        min_val = np.min(component_data)
        max_val = np.max(component_data)
        final_val = component_data[-1] if len(component_data) > 0 else 0
        
        # Store statistics for table
        label = legend_labels[i] if i < len(legend_labels) else f"Component {i+1}"
        stats_data.append({
            'Component': label,
            'Mean': mean_val,
            'Std Dev': std_val,
            'Min': min_val,
            'Max': max_val,
            'Final': final_val
        })
    
    # Plot individual components on the first subplot with increased line thickness
    for i in range(y.shape[1]):
        label = legend_labels[i] if i < len(legend_labels) else f"Component {i+1}"
        ax1.plot(x, y[:, i], label=label, linewidth=2.5, color=colors[i % len(colors)])
    
    # Add vertical markers for important time points
    if len(x) > 5:
        # Mark equilibration end (assuming first 20% is equilibration)
        equil_time = x[len(x)//5]
        ax1.axvline(x=equil_time, color='gray', linestyle=':', linewidth=2, alpha=0.7)
        ax1.text(equil_time + 0.02*(x[-1]-x[0]), ax1.get_ylim()[0] + 0.05*(ax1.get_ylim()[1]-ax1.get_ylim()[0]), 
                'Production Phase Start', fontsize=10, rotation=90, va='bottom')
        
        # Mark middle of production
        mid_time = x[len(x)//2]
        ax1.axvline(x=mid_time, color='gray', linestyle=':', linewidth=1.5, alpha=0.5)
        
        # Mark end of simulation
        end_time = x[-1] - 0.05*(x[-1]-x[0])
        ax1.axvline(x=end_time, color='gray', linestyle=':', linewidth=1.5, alpha=0.5)
    
    # Add gridlines to y-axis for easier reading
    ax1.grid(True, axis='y', linestyle='--', alpha=0.7)
    ax1.grid(True, axis='x', linestyle=':', alpha=0.3)
    
    # Add a secondary y-axis for kinetic energy if it exists
    kinetic_idx = None
    for i, label in enumerate(legend_labels):
        if 'kinetic' in label.lower():
            kinetic_idx = i
            break
    
    if kinetic_idx is not None:
        ax1_twin = ax1.twinx()
        ax1_twin.plot(x, y[:, kinetic_idx], label=f"{legend_labels[kinetic_idx]} (right axis)", 
                     linewidth=2.5, linestyle='--', color='darkred')
        ax1_twin.set_ylabel('Kinetic Energy (kJ/mol)', fontsize=12, color='darkred')
        ax1_twin.tick_params(axis='y', labelcolor='darkred')
        ax1_twin.spines['right'].set_color('darkred')
        
        # Add legend for secondary axis
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax1_twin.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=10, loc='best', ncol=2)
    else:
        ax1.legend(fontsize=10, loc='best', ncol=2)
    
    ax1.set_xlabel('Time (ns)', fontsize=14)
    ax1.set_ylabel('Energy (kJ/mol)', fontsize=14)
    ax1.set_title('Individual Energy Components', fontsize=16, fontweight='bold')
    
    # Create a table with statistics for each component
    table_data = []
    table_colors = []
    for i, stat in enumerate(stats_data):
        table_data.append([
            stat['Component'],
            f"{stat['Mean']:.2f}",
            f"{stat['Std Dev']:.2f}",
            f"{stat['Min']:.2f}",
            f"{stat['Max']:.2f}"
        ])
        table_colors.append(colors[i % len(colors)])
    
    # Plot energy drift as percentage over time in the bottom subplot
    if len(x) > 10:  # Only if we have enough data points
        # Calculate drift for total energy if available
        total_idx = None
        for i, label in enumerate(legend_labels):
            if 'total' in label.lower():
                total_idx = i
                break
        
        if total_idx is not None:
            # Calculate percentage drift relative to initial value
            initial_energy = y[0, total_idx]
            energy_drift_percent = [(e - initial_energy) / abs(initial_energy) * 100 for e in y[:, total_idx]]
            
            # Plot drift
            ax2.plot(x, energy_drift_percent, linewidth=2.5, color='darkblue', label='Total Energy Drift')
            ax2.set_xlabel('Time (ns)', fontsize=14)
            ax2.set_ylabel('Energy Drift (%)', fontsize=14)
            ax2.set_title('Energy Conservation (Drift)', fontsize=16, fontweight='bold')
            ax2.grid(True, linestyle='--', alpha=0.7)
            
            # Add horizontal line at zero
            ax2.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.5)
            
            # Add drift statistics
            final_drift = energy_drift_percent[-1]
            drift_per_ns = final_drift / x[-1] if x[-1] > 0 else 0
            
            drift_text = (
                f"Final Drift: {final_drift:.3f}%\n"
                f"Drift Rate: {drift_per_ns:.3f}% per ns\n"
                f"Drift Quality: {'Excellent' if abs(final_drift) < 0.01 else 'Good' if abs(final_drift) < 0.1 else 'Acceptable' if abs(final_drift) < 1 else 'Poor'}"
            )
            
            ax2.text(0.02, 0.95, drift_text, transform=ax2.transAxes, fontsize=11,
                    verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            # Add vertical markers matching the top plot
            if len(x) > 5:
                ax2.axvline(x=equil_time, color='gray', linestyle=':', linewidth=2, alpha=0.7)
                ax2.axvline(x=mid_time, color='gray', linestyle=':', linewidth=1.5, alpha=0.5)
                ax2.axvline(x=end_time, color='gray', linestyle=':', linewidth=1.5, alpha=0.5)
        else:
            # If no total energy, plot potential energy components drift
            for i in range(min(3, y.shape[1])):  # Limit to 3 components for clarity
                component_data = y[:, i]
                initial_value = component_data[0]
                drift_percent = [(e - initial_value) / abs(initial_value) * 100 for e in component_data]
                
                label = legend_labels[i] if i < len(legend_labels) else f"Component {i+1}"
                ax2.plot(x, drift_percent, linewidth=2.5, color=colors[i % len(colors)], 
                        label=f"{label} Drift")
            
            ax2.set_xlabel('Time (ns)', fontsize=14)
            ax2.set_ylabel('Energy Component Drift (%)', fontsize=14)
            ax2.set_title('Energy Components Stability', fontsize=16, fontweight='bold')
            ax2.grid(True, linestyle='--', alpha=0.7)
            ax2.legend(fontsize=10, loc='best')
            
            # Add horizontal line at zero
            ax2.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.5)
    else:
        # If not enough data for drift analysis, show statistics table instead
        ax2.axis('off')  # Turn off axis
        
        # Create a table with statistics
        columns = ['Component', 'Mean (kJ/mol)', 'Std Dev', 'Min', 'Max']
        table = ax2.table(cellText=table_data, colLabels=columns, loc='center', cellLoc='center')
        
        # Style the table
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)
        
        # Color the component names according to their plot colors
        for i, color in enumerate(table_colors):
            table[(i+1, 0)].set_facecolor(matplotlib.colors.to_rgba(color, alpha=0.3))
    
    # Add a title for the entire figure
    plt.suptitle('Energy Components Analysis - TIP4P Water Model', fontsize=18, fontweight='bold', y=0.98)
    
    # Add a caption with simulation conditions
    sim_details = (
        f"TIP4P Water Model | {NUM_WATER_MOLECULES} molecules | "
        f"Temperature: {REFERENCE_VALUES['temperature']} K | "
        f"Pressure: {REFERENCE_VALUES['pressure']} bar | "
        f"Simulation Time: {x[-1]:.1f} ns"
    )
    plt.figtext(0.5, 0.01, sim_details, ha='center', fontsize=11, 
               bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8, edgecolor='gray'))
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust for the suptitle and caption
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def create_thermodynamic_properties_plot(analysis_dir, plots_dir):
    """Create a combined plot of thermodynamic properties with enhanced visualization"""
    # Create a figure with 2x2 subplots
    fig, axs = plt.subplots(2, 2, figsize=(15, 12), dpi=300)
    
    # Define data directory
    data_dir = os.path.join(analysis_dir, "data")
    
    # Set a consistent style for all subplots
    for ax in axs.flatten():
        ax.grid(True, linestyle='--', alpha=0.7)
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
    
    # Temperature subplot
    temp_file = os.path.join(data_dir, 'temperature.xvg')
    if os.path.exists(temp_file):
        x, y, _, _, _, _ = read_xvg(temp_file)
        if len(x) > 0 and len(y) > 0:
            # Plot with enhanced styling
            if HAS_SEABORN:
                sns.lineplot(x=x, y=y, color='#1f77b4', ax=axs[0, 0], linewidth=2.5)
            else:
                axs[0, 0].plot(x, y, color='#1f77b4', linewidth=2.5)
            
            # Add reference line with annotation
            ref_temp = REFERENCE_VALUES['temperature']
            axs[0, 0].axhline(y=ref_temp, color='#d62728', linestyle='--', alpha=0.7, 
                             linewidth=2, label=f'Reference: {ref_temp} K')
            
            # Calculate statistics
            temp_mean = np.mean(y)
            temp_std = np.std(y)
            
            # Add text with statistics and enhanced styling
            stats_text = (
                f"Mean: {temp_mean:.1f} K\n"
                f"Std Dev: {temp_std:.1f} K\n"
                f"Target: {REFERENCE_VALUES['temperature']:.1f} K\n"
                f"Deviation: {abs(temp_mean - REFERENCE_VALUES['temperature']):.1f} K"
            )
            axs[0, 0].text(0.05, 0.95, stats_text, transform=axs[0, 0].transAxes, 
                          fontsize=11, verticalalignment='top', 
                          bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
            
            # Mark equilibration region if simulation is long enough
            if len(x) > 5:
                equil_time = x[len(x)//5]  # Assume first 20% is equilibration
                axs[0, 0].axvline(x=equil_time, color='gray', linestyle=':', linewidth=2, alpha=0.7)
                axs[0, 0].text(equil_time + 0.5, axs[0, 0].get_ylim()[0] + 0.1*(axs[0, 0].get_ylim()[1]-axs[0, 0].get_ylim()[0]), 
                              'Production', fontsize=10, rotation=0, va='bottom')
            
            axs[0, 0].set_title('Temperature', fontsize=16, fontweight='bold')
            axs[0, 0].set_xlabel('Time (ps)', fontsize=14)
            axs[0, 0].set_ylabel('Temperature (K)', fontsize=14)
            axs[0, 0].tick_params(axis='both', which='major', labelsize=12)
            axs[0, 0].legend(fontsize=11, framealpha=0.8)
    
    # Pressure subplot
    pressure_file = os.path.join(data_dir, 'pressure.xvg')
    if os.path.exists(pressure_file):
        x, y, _, _, _, _ = read_xvg(pressure_file)
        if len(x) > 0 and len(y) > 0:
            # Plot with enhanced styling
            if HAS_SEABORN:
                sns.lineplot(x=x, y=y, color='#2ca02c', ax=axs[0, 1], linewidth=2.5)
            else:
                axs[0, 1].plot(x, y, color='#2ca02c', linewidth=2.5)
            
            # Add reference line with annotation
            ref_pressure = REFERENCE_VALUES['pressure']
            axs[0, 1].axhline(y=ref_pressure, color='#d62728', linestyle='--', alpha=0.7, 
                             linewidth=2, label=f'Reference: {ref_pressure} bar')
            
            # Calculate statistics
            pressure_mean = np.mean(y)
            pressure_std = np.std(y)
            
            # Add text with statistics and enhanced styling
            stats_text = (
                f"Mean: {pressure_mean:.1f} bar\n"
                f"Std Dev: {pressure_std:.1f} bar\n"
                f"Target: {REFERENCE_VALUES['pressure']:.1f} bar\n"
                f"Deviation: {abs(pressure_mean - REFERENCE_VALUES['pressure']):.1f} bar"
            )
            axs[0, 1].text(0.05, 0.95, stats_text, transform=axs[0, 1].transAxes, 
                          fontsize=11, verticalalignment='top', 
                          bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
            
            # Add note about pressure fluctuations with enhanced styling
            note_text = "Note: Large pressure fluctuations\nare normal in MD simulations"
            axs[0, 1].text(0.05, 0.75, note_text, transform=axs[0, 1].transAxes,
                          fontsize=10, style='italic',
                          bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.8, edgecolor='gray'))
            
            # Mark equilibration region if simulation is long enough
            if len(x) > 5:
                equil_time = x[len(x)//5]  # Assume first 20% is equilibration
                axs[0, 1].axvline(x=equil_time, color='gray', linestyle=':', linewidth=2, alpha=0.7)
            
            axs[0, 1].set_title('Pressure', fontsize=16, fontweight='bold')
            axs[0, 1].set_xlabel('Time (ps)', fontsize=14)
            axs[0, 1].set_ylabel('Pressure (bar)', fontsize=14)
            axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
            axs[0, 1].legend(fontsize=11, framealpha=0.8)
    
    # Energy subplot - FIXED to avoid recursion
    energy_file = os.path.join(data_dir, 'energy.xvg')
    if os.path.exists(energy_file):
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(energy_file)
        
        if len(x) > 0 and len(y) > 0 and len(y.shape) > 1:
            # Plot with enhanced styling
            for i in range(min(3, y.shape[1])):  # Limit to 3 components for clarity
                label = legend_labels[i] if i < len(legend_labels) else f"Component {i+1}"
                axs[1, 0].plot(x, y[:, i], label=label, linewidth=2.5)
            
            # Calculate statistics for total energy
            if y.shape[1] > 0:
                energy_mean = np.mean(y[:, 0])
                energy_std = np.std(y[:, 0])
                
                # Add text with statistics
                stats_text = (
                    f"Mean Total: {energy_mean:.1f} kJ/mol\n"
                    f"Std Dev: {energy_std:.1f} kJ/mol"
                )
                axs[1, 0].text(0.05, 0.95, stats_text, transform=axs[1, 0].transAxes, 
                              fontsize=11, verticalalignment='top', 
                              bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
            
            # Mark equilibration region if simulation is long enough
            if len(x) > 5:
                equil_time = x[len(x)//5]  # Assume first 20% is equilibration
                axs[1, 0].axvline(x=equil_time, color='gray', linestyle=':', linewidth=2, alpha=0.7)
            
            axs[1, 0].set_title('Energy Components', fontsize=16, fontweight='bold')
            axs[1, 0].set_xlabel('Time (ps)', fontsize=14)
            axs[1, 0].set_ylabel('Energy (kJ/mol)', fontsize=14)
            axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
            axs[1, 0].legend(fontsize=11, framealpha=0.8, loc='best')
    
    # Potential energy subplot - FIXED to avoid recursion
    potential_file = os.path.join(data_dir, 'potential.xvg')
    if os.path.exists(potential_file):
        x, y, _, _, _, _ = read_xvg(potential_file)
        if len(x) > 0 and len(y) > 0:
            # Plot with enhanced styling
            if HAS_SEABORN:
                sns.lineplot(x=x, y=y, color='#d62728', ax=axs[1, 1], linewidth=2.5)
            else:
                axs[1, 1].plot(x, y, color='#d62728', linewidth=2.5)
            
            # Calculate statistics
            potential_mean = np.mean(y)
            potential_std = np.std(y)
            
            # Scale the per-molecule reference to the total system
            total_reference = REFERENCE_VALUES['potential_energy'] * NUM_WATER_MOLECULES
            
            # Add reference line with enhanced styling
            axs[1, 1].axhline(y=total_reference, color='#7f7f7f', 
                             linestyle='--', alpha=0.7, linewidth=2,
                             label=f'Reference: {total_reference:.1f} kJ/mol')
            
            # Add text with statistics and enhanced styling
            stats_text = (
                f"Mean: {potential_mean:.2f} kJ/mol\n"
                f"Std Dev: {potential_std:.2f} kJ/mol\n"
                f"Per molecule: {potential_mean/NUM_WATER_MOLECULES:.2f} kJ/mol\n"
                f"Reference: {REFERENCE_VALUES['potential_energy']:.1f} kJ/mol per molecule"
            )
            axs[1, 1].text(0.05, 0.95, stats_text, transform=axs[1, 1].transAxes, 
                          fontsize=11, verticalalignment='top', 
                          bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
            
            # Mark equilibration region if simulation is long enough
            if len(x) > 5:
                equil_time = x[len(x)//5]  # Assume first 20% is equilibration
                axs[1, 1].axvline(x=equil_time, color='gray', linestyle=':', linewidth=2, alpha=0.7)
            
            axs[1, 1].set_title('Potential Energy', fontsize=16, fontweight='bold')
            axs[1, 1].set_xlabel('Time (ps)', fontsize=14)
            axs[1, 1].set_ylabel('Potential Energy (kJ/mol)', fontsize=14)
            axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
            axs[1, 1].legend(fontsize=11, framealpha=0.8)
    
    # Add a title for the entire figure with enhanced styling
    fig.suptitle('Thermodynamic Properties of TIP4P Water', fontsize=18, fontweight='bold', y=0.98)
    
    # Add a text box with simulation details
    sim_details = (
        f"System: TIP4P Water\n"
        f"Temperature: {REFERENCE_VALUES['temperature']} K\n"
        f"Pressure: {REFERENCE_VALUES['pressure']} bar\n"
        f"Molecules: {NUM_WATER_MOLECULES}"
    )
    fig.text(0.02, 0.02, sim_details, fontsize=10, 
             bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust for the suptitle and footer
    plt.savefig(os.path.join(plots_dir, 'thermodynamic_properties_plot.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print('  - thermodynamic_properties_plot.png saved successfully')

def main():
    if len(sys.argv) < 3:
        print("Usage: plot_temperature.py <analysis_dir> <plots_dir>")
        sys.exit(1)
        
    analysis_dir = sys.argv[1]
    plots_dir = sys.argv[2]
    
    # Define data directory
    data_dir = os.path.join(analysis_dir, "data")
    
    # Ensure plots directory exists
    os.makedirs(plots_dir, exist_ok=True)
    
    # Temperature plot
    temperature_file = os.path.join(data_dir, 'temperature.xvg')
    if os.path.exists(temperature_file):
        print('Plotting temperature data...')
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(temperature_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Temperature vs Time'
            plot_xlabel = xlabel if xlabel else 'Time (ns)'
            plot_ylabel = ylabel if ylabel else 'Temperature (K)'
            
            output_path = os.path.join(plots_dir, 'temperature_plot.png')
            plot_time_series(x, y, plot_title, plot_xlabel, plot_ylabel, legend_labels, 
                           output_path, REFERENCE_VALUES['temperature'], 'Temperature')
    else:
        print(f"Temperature file not found: {temperature_file}")
    
    # Pressure plot
    pressure_file = os.path.join(data_dir, 'pressure.xvg')
    if os.path.exists(pressure_file):
        print('Plotting pressure data...')
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(pressure_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Pressure vs Time'
            plot_xlabel = xlabel if xlabel else 'Time (ns)'
            plot_ylabel = ylabel if ylabel else 'Pressure (bar)'
            
            output_path = os.path.join(plots_dir, 'pressure_plot.png')
            plot_time_series(x, y, plot_title, plot_xlabel, plot_ylabel, legend_labels, 
                           output_path, REFERENCE_VALUES['pressure'], 'Pressure')
    else:
        print(f"Pressure file not found: {pressure_file}")
    
    # Energy plot
    energy_file = os.path.join(data_dir, 'energy.xvg')
    if os.path.exists(energy_file):
        print('Plotting energy data...')
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(energy_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Energy Components vs Time'
            plot_xlabel = xlabel if xlabel else 'Time (ns)'
            plot_ylabel = ylabel if ylabel else 'Energy (kJ/mol)'
            
            output_path = os.path.join(plots_dir, 'energy_plot.png')
            plot_time_series(x, y, plot_title, plot_xlabel, plot_ylabel, legend_labels, output_path, 
                           property_name='Energy Components')
            
            # Create a specialized energy components plot
            components_path = os.path.join(plots_dir, 'energy_components_plot.png')
            plot_energy_components(x, y, legend_labels, components_path)
    else:
        print(f"Energy file not found: {energy_file}")
    
    # Potential energy plot
    potential_file = os.path.join(data_dir, 'potential.xvg')
    if os.path.exists(potential_file):
        print('Plotting potential energy data...')
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(potential_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Potential Energy vs Time'
            plot_xlabel = xlabel if xlabel else 'Time (ns)'
            plot_ylabel = ylabel if ylabel else 'Potential Energy (kJ/mol)'
            
            output_path = os.path.join(plots_dir, 'potential_plot.png')
            plot_time_series(x, y, plot_title, plot_xlabel, plot_ylabel, legend_labels, 
                           output_path, REFERENCE_VALUES['potential_energy'], 'Potential Energy')
    else:
        print(f"Potential energy file not found: {potential_file}")
    
    # Detailed energy terms
    energy_terms_file = os.path.join(data_dir, 'energy_terms.xvg')
    if os.path.exists(energy_terms_file):
        print('Plotting detailed energy terms...')
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(energy_terms_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Detailed Energy Terms vs Time'
            plot_xlabel = xlabel if xlabel else 'Time (ns)'
            plot_ylabel = ylabel if ylabel else 'Energy (kJ/mol)'
            
            output_path = os.path.join(plots_dir, 'energy_terms_plot.png')
            plot_time_series(x, y, plot_title, plot_xlabel, plot_ylabel, legend_labels, 
                           output_path, property_name='Detailed Energy Terms')
            
            # Create a specialized energy terms plot
            terms_path = os.path.join(plots_dir, 'energy_terms_stacked_plot.png')
            plot_energy_components(x, y, legend_labels, terms_path)
    else:
        print(f"Detailed energy terms file not found: {energy_terms_file}")
    
    # Create a combined thermodynamic properties plot
    print('Creating combined thermodynamic properties plot...')
    create_thermodynamic_properties_plot(analysis_dir, plots_dir)

if __name__ == "__main__":
    main()
