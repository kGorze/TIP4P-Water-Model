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

# Reference values for TIP4P water at 298K
REFERENCE_VALUES = {
    'density': 997.0,  # kg/m^3, bulk density
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

def read_xpm(filename):
    """Read .xpm files and extract matrix data."""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Extract dimensions from the XPM file
        nx, ny = None, None
        for i, line in enumerate(lines):
            if "static char *gromacs_xpm[]" in line:
                # The next line should contain dimensions
                dim_line = lines[i+1].strip()
                # Format is typically: /* nx ny ncolors char_per_pixel */
                parts = dim_line.split()
                if len(parts) >= 4:
                    try:
                        nx = int(parts[0])
                        ny = int(parts[1])
                        break
                    except ValueError:
                        pass
        
        if nx is None or ny is None:
            print(f"Could not determine dimensions in XPM file: {filename}")
            return None, None, None, None
        
        # Find the start of the data
        data_start = None
        for i, line in enumerate(lines):
            if line.strip().startswith('"'):
                data_start = i
                break
        
        if data_start is None:
            print(f"Could not find data start in XPM file: {filename}")
            return None, None, None, None
        
        # Extract the color mapping
        color_map = {}
        color_values = []
        for i in range(data_start, len(lines)):
            line = lines[i].strip()
            if line.startswith('"') and "c " in line:
                parts = line.split('c ')
                if len(parts) >= 2:
                    char = parts[0].strip('" ')
                    color_parts = parts[1].split()
                    if len(color_parts) >= 2:
                        try:
                            value = float(color_parts[1])
                            color_map[char] = value
                            color_values.append(value)
                        except ValueError:
                            pass
            elif line.startswith('"') and len(line) > 2 and not "c " in line:
                # This is likely the start of the data matrix
                break
        
        # Extract the data matrix
        data_matrix = np.zeros((ny, nx))
        row = 0
        for i in range(data_start, len(lines)):
            line = lines[i].strip()
            if line.startswith('"') and not "c " in line:
                if row < ny:
                    data_line = line.strip('" ')
                    for col in range(min(nx, len(data_line))):
                        char = data_line[col]
                        if char in color_map:
                            data_matrix[row, col] = color_map[char]
                    row += 1
        
        # Create x and y axes
        # Assuming the XPM represents a square grid from 0 to 1 in both dimensions
        x = np.linspace(0, 1, nx)
        y = np.linspace(0, 1, ny)
        
        return x, y, data_matrix, color_values
    except Exception as e:
        print(f"Error reading XPM file {filename}: {e}")
        return None, None, None, None

def plot_density_profile(x, y, title, xlabel, ylabel, legend_labels, output_path):
    """Plot density profile with statistics and annotations"""
    plt.figure(figsize=(12, 8), dpi=300)
    
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
        
        # Use the first column for statistics
        density_data = y[:, 0]
    else:
        # Single column plot
        if HAS_SEABORN:
            sns.lineplot(x=x, y=y, color='#1f77b4', linewidth=2)
        else:
            plt.plot(x, y, color='#1f77b4', linewidth=2)
        
        density_data = y
    
    # Calculate statistics
    mean_density = np.mean(density_data)
    std_density = np.std(density_data)
    max_density = np.max(density_data)
    min_density = np.min(density_data)
    
    # Add horizontal line for mean density
    plt.axhline(y=mean_density, color='#2ca02c', linestyle='--', alpha=0.7,
               label=f'Mean: {mean_density:.1f} kg/m³')
    
    # Add reference value
    plt.axhline(y=REFERENCE_VALUES['density'], color='#d62728', linestyle=':', alpha=0.7,
               label=f'Reference: {REFERENCE_VALUES["density"]:.1f} kg/m³')
    
    # Add shaded area for standard deviation
    plt.fill_between(x, mean_density - std_density, mean_density + std_density,
                   color='#2ca02c', alpha=0.2, label=f'Std Dev: ±{std_density:.1f} kg/m³')
    
    # Add annotations for fluctuations
    fluctuation_percent = (std_density / mean_density) * 100
    plt.annotate(f'Fluctuation: {fluctuation_percent:.1f}%',
                xy=(x[len(x)//2], mean_density + std_density),
                xytext=(x[len(x)//2], mean_density + 2*std_density),
                arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                fontsize=10, ha='center')
    
    # Add a text box with statistics
    stats_text = (
        f"Mean Density: {mean_density:.1f} kg/m³\n"
        f"Std Dev: {std_density:.1f} kg/m³\n"
        f"Min: {min_density:.1f} kg/m³\n"
        f"Max: {max_density:.1f} kg/m³\n"
        f"Fluctuation: {fluctuation_percent:.1f}%\n"
        f"Reference: {REFERENCE_VALUES['density']:.1f} kg/m³\n"
        f"Difference: {((mean_density - REFERENCE_VALUES['density']) / REFERENCE_VALUES['density'] * 100):.1f}%"
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
    plt.figtext(0.5, 0.01, 'TIP4P Water Model - Density Analysis', 
               ha='center', fontsize=10, style='italic', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')
    
    # Create a histogram of density values
    plt.figure(figsize=(10, 6), dpi=300)
    
    if HAS_SEABORN:
        sns.histplot(density_data, kde=True, color='#1f77b4')
        
        # Add vertical lines for statistics
        plt.axvline(x=mean_density, color='#2ca02c', linestyle='--', 
                   label=f'Mean: {mean_density:.1f} kg/m³')
        plt.axvline(x=REFERENCE_VALUES['density'], color='#d62728', linestyle=':', 
                   label=f'Reference: {REFERENCE_VALUES["density"]:.1f} kg/m³')
    else:
        plt.hist(density_data, bins=30, alpha=0.7, color='#1f77b4', density=True)
        
        # Add a simple KDE if seaborn is not available
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(density_data)
        x_kde = np.linspace(min(density_data), max(density_data), 100)
        plt.plot(x_kde, kde(x_kde), 'r-', linewidth=2)
        
        # Add vertical lines for statistics
        plt.axvline(x=mean_density, color='#2ca02c', linestyle='--', 
                   label=f'Mean: {mean_density:.1f} kg/m³')
        plt.axvline(x=REFERENCE_VALUES['density'], color='#d62728', linestyle=':', 
                   label=f'Reference: {REFERENCE_VALUES["density"]:.1f} kg/m³')
    
    plt.xlabel('Density (kg/m³)', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title('Distribution of Density Values', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=10)
    
    hist_output_path = os.path.join(os.path.dirname(output_path), 'density_histogram.png')
    plt.tight_layout()
    plt.savefig(hist_output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - density_histogram.png saved successfully')

def plot_radial_density_map(x, y, data, color_values, output_path):
    """Plot radial density map as a heatmap"""
    if data is None or x is None or y is None:
        print("  - Not enough data for radial density map")
        return
    
    plt.figure(figsize=(10, 8), dpi=300)
    
    # Create a heatmap
    if HAS_SEABORN:
        # Use seaborn's heatmap for better visualization
        ax = sns.heatmap(data, cmap='viridis', xticklabels=False, yticklabels=False)
        
        # Add colorbar with proper label
        cbar = ax.collections[0].colorbar
        cbar.set_label('Density (kg/m³)', fontsize=12)
    else:
        # Use matplotlib's imshow
        plt.imshow(data, cmap='viridis', origin='lower', aspect='equal')
        
        # Add colorbar
        cbar = plt.colorbar()
        cbar.set_label('Density (kg/m³)', fontsize=12)
    
    # Add axes labels and title
    plt.xlabel('X (nm)', fontsize=14)
    plt.ylabel('Y (nm)', fontsize=14)
    plt.title('Radial Density Map', fontsize=16)
    
    # Add a circle to indicate the center
    center_x = len(x) // 2
    center_y = len(y) // 2
    circle = plt.Circle((center_x, center_y), 5, color='red', fill=False, linewidth=2)
    plt.gca().add_patch(circle)
    
    # Add annotations for density regions
    if color_values:
        max_density = max(color_values)
        min_density = min(color_values)
        
        # Add text with density range
        plt.text(0.02, 0.98, f"Density Range:\n{min_density:.1f} - {max_density:.1f} kg/m³",
                transform=plt.gca().transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    # Add a watermark with simulation details
    plt.figtext(0.5, 0.01, 'TIP4P Water Model - Radial Density Analysis', 
               ha='center', fontsize=10, style='italic', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def main():
    if len(sys.argv) < 3:
        print("Usage: plot_density.py <analysis_dir> <plots_dir>")
        sys.exit(1)
        
    analysis_dir = sys.argv[1]
    plots_dir = sys.argv[2]
    
    # Ensure plots directory exists
    os.makedirs(plots_dir, exist_ok=True)
    
    # Plot density profile
    density_file = os.path.join(analysis_dir, 'density.xvg')
    if os.path.exists(density_file):
        print('Plotting density profile...')
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(density_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Density Profile'
            plot_xlabel = xlabel if xlabel else 'Position (nm)'
            plot_ylabel = ylabel if ylabel else 'Density (kg/m³)'
            
            output_path = os.path.join(plots_dir, 'density_profile_plot.png')
            plot_density_profile(x, y, plot_title, plot_xlabel, plot_ylabel, legend_labels, output_path)
    else:
        print(f"Density profile file not found: {density_file}")
    
    # Plot radial density map
    radial_density_file = os.path.join(analysis_dir, 'density_radial.xpm')
    if os.path.exists(radial_density_file):
        print('Plotting radial density map...')
        x, y, data, color_values = read_xpm(radial_density_file)
        
        if data is not None:
            output_path = os.path.join(plots_dir, 'radial_density_map.png')
            plot_radial_density_map(x, y, data, color_values, output_path)
        else:
            print("  - Could not read radial density map data")
    else:
        print(f"Radial density map file not found: {radial_density_file}")

if __name__ == "__main__":
    main() 