#!/usr/bin/python3
import os
import sys
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
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

# Reference values for TIP4P water at 298K, 1 bar
REFERENCE_VALUES = {
    'temperature': 298.0,  # K
    'pressure': 1.0,       # bar
    'density': 997.0,      # kg/m^3
    'potential_energy': -41.5,  # kJ/mol
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

def plot_time_series(x, y, title, xlabel, ylabel, legend_labels, output_path, reference_value=None, property_name=None):
    """Create a time series plot with statistics and annotations"""
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

def plot_energy_components(x, y, legend_labels, output_path):
    """Create a specialized plot for energy components"""
    if len(y.shape) <= 1 or y.shape[1] <= 1:
        print("  - Not enough energy components to plot")
        return
    
    # Create a figure with two subplots - stacked and individual
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), dpi=300)
    
    # Plot individual components on the first subplot
    for i in range(y.shape[1]):
        label = legend_labels[i] if i < len(legend_labels) else f"Component {i+1}"
        ax1.plot(x, y[:, i], label=label, linewidth=1.5)
    
    ax1.set_xlabel('Time (ns)', fontsize=12)
    ax1.set_ylabel('Energy (kJ/mol)', fontsize=12)
    ax1.set_title('Individual Energy Components', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9, loc='best', ncol=2)
    
    # Plot stacked area chart on the second subplot
    # Only include components that make sense in a stacked chart (exclude total energy)
    components_to_stack = []
    stack_labels = []
    
    # Identify which components to include in the stack
    exclude_terms = ['total', 'potential', 'kinetic']
    for i in range(y.shape[1]):
        label = legend_labels[i].lower() if i < len(legend_labels) else ""
        if not any(term in label for term in exclude_terms):
            components_to_stack.append(y[:, i])
            stack_labels.append(legend_labels[i] if i < len(legend_labels) else f"Component {i+1}")
    
    if components_to_stack:
        components_array = np.array(components_to_stack).T  # Transpose for stacking
        
        if HAS_SEABORN:
            # Use a custom colormap for better visualization
            cmap = sns.color_palette("viridis", n_colors=len(components_to_stack))
            ax2.stackplot(x, components_array.T, labels=stack_labels, colors=cmap, alpha=0.7)
        else:
            ax2.stackplot(x, components_array.T, labels=stack_labels, alpha=0.7)
        
        ax2.set_xlabel('Time (ns)', fontsize=12)
        ax2.set_ylabel('Energy (kJ/mol)', fontsize=12)
        ax2.set_title('Stacked Energy Components', fontsize=14)
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=9, loc='best', ncol=2)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def main():
    if len(sys.argv) < 3:
        print("Usage: plot_temperature.py <analysis_dir> <plots_dir>")
        sys.exit(1)
        
    analysis_dir = sys.argv[1]
    plots_dir = sys.argv[2]
    
    # Ensure plots directory exists
    os.makedirs(plots_dir, exist_ok=True)
    
    # Plot temperature data
    temp_file = os.path.join(analysis_dir, 'temperature.xvg')
    if os.path.exists(temp_file):
        print('Plotting temperature data...')
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(temp_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Temperature vs Time'
            plot_xlabel = xlabel if xlabel else 'Time (ns)'
            plot_ylabel = ylabel if ylabel else 'Temperature (K)'
            
            output_path = os.path.join(plots_dir, 'temperature_plot.png')
            plot_time_series(x, y, plot_title, plot_xlabel, plot_ylabel, legend_labels, 
                           output_path, REFERENCE_VALUES['temperature'], 'Temperature')
    else:
        print(f"Temperature file not found: {temp_file}")
    
    # Plot pressure data
    pressure_file = os.path.join(analysis_dir, 'pressure.xvg')
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
    
    # Plot energy data
    energy_file = os.path.join(analysis_dir, 'energy.xvg')
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
    
    # Plot potential energy data
    potential_file = os.path.join(analysis_dir, 'potential.xvg')
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
    
    # Plot detailed energy terms
    energy_terms_file = os.path.join(analysis_dir, 'energy_terms.xvg')
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
    plt.figure(figsize=(15, 10), dpi=300)
    
    # Create a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(15, 10), dpi=300)
    
    # Temperature subplot
    temp_file = os.path.join(analysis_dir, 'temperature.xvg')
    if os.path.exists(temp_file):
        x, y, _, _, _, _ = read_xvg(temp_file)
        if len(x) > 0 and len(y) > 0:
            if HAS_SEABORN:
                sns.lineplot(x=x, y=y, color='#1f77b4', ax=axs[0, 0])
            else:
                axs[0, 0].plot(x, y, color='#1f77b4')
            
            # Add reference line
            axs[0, 0].axhline(y=REFERENCE_VALUES['temperature'], color='#d62728', 
                             linestyle=':', alpha=0.7)
            
            # Calculate statistics
            temp_mean = np.mean(y)
            temp_std = np.std(y)
            
            # Add text with statistics
            stats_text = f"Mean: {temp_mean:.1f} K\nStd Dev: {temp_std:.1f} K"
            axs[0, 0].text(0.05, 0.95, stats_text, transform=axs[0, 0].transAxes, 
                          fontsize=10, verticalalignment='top', 
                          bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
            
            axs[0, 0].set_title('Temperature', fontsize=14)
            axs[0, 0].set_xlabel('Time (ns)', fontsize=12)
            axs[0, 0].set_ylabel('Temperature (K)', fontsize=12)
            axs[0, 0].grid(True, alpha=0.3)
    
    # Pressure subplot
    pressure_file = os.path.join(analysis_dir, 'pressure.xvg')
    if os.path.exists(pressure_file):
        x, y, _, _, _, _ = read_xvg(pressure_file)
        if len(x) > 0 and len(y) > 0:
            if HAS_SEABORN:
                sns.lineplot(x=x, y=y, color='#2ca02c', ax=axs[0, 1])
            else:
                axs[0, 1].plot(x, y, color='#2ca02c')
            
            # Add reference line
            axs[0, 1].axhline(y=REFERENCE_VALUES['pressure'], color='#d62728', 
                             linestyle=':', alpha=0.7)
            
            # Calculate statistics
            pressure_mean = np.mean(y)
            pressure_std = np.std(y)
            
            # Add text with statistics
            stats_text = f"Mean: {pressure_mean:.1f} bar\nStd Dev: {pressure_std:.1f} bar"
            axs[0, 1].text(0.05, 0.95, stats_text, transform=axs[0, 1].transAxes, 
                          fontsize=10, verticalalignment='top', 
                          bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
            
            axs[0, 1].set_title('Pressure', fontsize=14)
            axs[0, 1].set_xlabel('Time (ns)', fontsize=12)
            axs[0, 1].set_ylabel('Pressure (bar)', fontsize=12)
            axs[0, 1].grid(True, alpha=0.3)
    
    # Energy subplot
    energy_file = os.path.join(analysis_dir, 'energy.xvg')
    if os.path.exists(energy_file):
        x, y, _, _, _, legend_labels = read_xvg(energy_file)
        if len(x) > 0 and len(y) > 0 and len(y.shape) > 1:
            # Plot total energy and potential energy if available
            for i in range(min(2, y.shape[1])):
                label = legend_labels[i] if i < len(legend_labels) else f"Series {i+1}"
                axs[1, 0].plot(x, y[:, i], label=label)
            
            axs[1, 0].set_title('Energy Components', fontsize=14)
            axs[1, 0].set_xlabel('Time (ns)', fontsize=12)
            axs[1, 0].set_ylabel('Energy (kJ/mol)', fontsize=12)
            axs[1, 0].grid(True, alpha=0.3)
            axs[1, 0].legend(fontsize=10)
    
    # Potential energy subplot
    potential_file = os.path.join(analysis_dir, 'potential.xvg')
    if os.path.exists(potential_file):
        x, y, _, _, _, _ = read_xvg(potential_file)
        if len(x) > 0 and len(y) > 0:
            if HAS_SEABORN:
                sns.lineplot(x=x, y=y, color='#d62728', ax=axs[1, 1])
            else:
                axs[1, 1].plot(x, y, color='#d62728')
            
            # Add reference line
            axs[1, 1].axhline(y=REFERENCE_VALUES['potential_energy'], color='#7f7f7f', 
                             linestyle=':', alpha=0.7)
            
            # Calculate statistics
            potential_mean = np.mean(y)
            potential_std = np.std(y)
            
            # Add text with statistics
            stats_text = f"Mean: {potential_mean:.1f} kJ/mol\nStd Dev: {potential_std:.1f} kJ/mol"
            axs[1, 1].text(0.05, 0.95, stats_text, transform=axs[1, 1].transAxes, 
                          fontsize=10, verticalalignment='top', 
                          bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
            
            axs[1, 1].set_title('Potential Energy', fontsize=14)
            axs[1, 1].set_xlabel('Time (ns)', fontsize=12)
            axs[1, 1].set_ylabel('Potential Energy (kJ/mol)', fontsize=12)
            axs[1, 1].grid(True, alpha=0.3)
    
    # Add a title for the entire figure
    fig.suptitle('Thermodynamic Properties of TIP4P Water', fontsize=16)
    
    plt.tight_layout(rect=[0, 0, 1, 0.97])  # Adjust for the suptitle
    plt.savefig(os.path.join(plots_dir, 'thermodynamic_properties_plot.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print('  - thermodynamic_properties_plot.png saved successfully')

if __name__ == "__main__":
    main()
