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

def plot_density_histogram(density_data, output_path):
    """Create an enhanced density histogram plot with annotations and styling"""
    # Extract data
    x, y = density_data
    
    # Create figure with enhanced styling
    plt.figure(figsize=(10, 6), dpi=300)
    
    # Calculate statistics
    mean_density = np.mean(y)
    std_density = np.std(y)
    
    # Determine optimal bin size using Freedman-Diaconis rule
    # This rule is better than Sturges' rule for non-normal distributions
    q75, q25 = np.percentile(y, [75, 25])
    iqr = q75 - q25
    bin_width = 2 * iqr / (len(y) ** (1/3))
    if bin_width > 0:
        num_bins = int(np.ceil((max(y) - min(y)) / bin_width))
        num_bins = min(max(10, num_bins), 50)  # Keep bins between 10 and 50
    else:
        num_bins = 20  # Default if calculation fails
    
    # Plot histogram with enhanced styling
    n, bins, patches = plt.hist(y, bins=num_bins, alpha=0.7, color='#1f77b4', 
                               edgecolor='black', linewidth=1.2)
    
    # Add reference line for mean density
    plt.axvline(x=mean_density, color='#d62728', linestyle='--', linewidth=2, 
               label=f'Mean: {mean_density:.2f} kg/m³')
    
    # Add reference lines for standard deviation
    plt.axvline(x=mean_density + std_density, color='#2ca02c', linestyle=':', linewidth=1.5,
               label=f'±1σ: {std_density:.2f} kg/m³')
    plt.axvline(x=mean_density - std_density, color='#2ca02c', linestyle=':', linewidth=1.5)
    
    # Add reference for experimental density of water at 298K (if applicable)
    exp_density = 997.0  # kg/m³
    plt.axvline(x=exp_density, color='#ff7f0e', linestyle='-', linewidth=1.5,
               label=f'Exp. (298K): {exp_density:.1f} kg/m³')
    
    # Add title and labels with enhanced styling
    plt.title('Density Distribution of TIP4P Water', fontsize=16, fontweight='bold')
    plt.xlabel('Density (kg/m³)', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    
    # Add grid and improve styling
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tick_params(axis='both', which='major', labelsize=12)
    
    # Add statistics text box
    stats_text = (
        f"Mean: {mean_density:.2f} kg/m³\n"
        f"Std Dev: {std_density:.2f} kg/m³\n"
        f"Min: {min(y):.2f} kg/m³\n"
        f"Max: {max(y):.2f} kg/m³\n"
        f"Samples: {len(y)}"
    )
    plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes, fontsize=11,
            verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', 
                                             facecolor='white', alpha=0.8, edgecolor='gray'))
    
    # Add legend with enhanced styling
    plt.legend(fontsize=11, framealpha=0.8, loc='upper right')
    
    # Save figure with tight layout
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def plot_density_profile(x, y, title, xlabel, ylabel, output_path):
    """Create an enhanced density profile plot with annotations and styling"""
    # Create figure with enhanced styling
    plt.figure(figsize=(10, 6), dpi=300)
    
    # Plot density profile with enhanced styling
    if HAS_SEABORN:
        sns.lineplot(x=x, y=y, color='#1f77b4', linewidth=2.5)
    else:
        plt.plot(x, y, color='#1f77b4', linewidth=2.5)
    
    # Calculate statistics
    mean_density = np.mean(y)
    std_density = np.std(y)
    
    # Add reference line for mean density
    plt.axhline(y=mean_density, color='#d62728', linestyle='--', linewidth=2, 
               label=f'Mean: {mean_density:.2f} kg/m³')
    
    # Add reference for experimental density of water at 298K (if applicable)
    exp_density = 997.0  # kg/m³
    plt.axhline(y=exp_density, color='#ff7f0e', linestyle='-', linewidth=1.5,
               label=f'Exp. (298K): {exp_density:.1f} kg/m³')
    
    # Add shaded region for standard deviation
    plt.fill_between(x, mean_density - std_density, mean_density + std_density, 
                    color='#2ca02c', alpha=0.2, label=f'±1σ: {std_density:.2f} kg/m³')
    
    # Identify regions of interest (e.g., bulk regions, interfaces)
    # This is a simple example - for a real system, you might need more sophisticated detection
    try:
        # Find regions where density is significantly different from the mean
        threshold = 0.8 * mean_density
        regions = []
        in_region = False
        start_idx = 0
        
        for i, val in enumerate(y):
            if val < threshold and not in_region:
                in_region = True
                start_idx = i
            elif val >= threshold and in_region:
                in_region = False
                regions.append((start_idx, i))
        
        # Add annotations for identified regions
        for i, (start, end) in enumerate(regions):
            mid = (start + end) // 2
            if mid < len(x):
                plt.axvspan(x[start], x[end], alpha=0.2, color='gray')
                plt.text(x[mid], min(y) + 0.1 * (max(y) - min(y)), 
                        f'Region {i+1}', ha='center', fontsize=10,
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
    except:
        # If region detection fails, continue without annotations
        pass
    
    # Add title and labels with enhanced styling
    plt.title(title, fontsize=16, fontweight='bold')
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    
    # Add grid and improve styling
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tick_params(axis='both', which='major', labelsize=12)
    
    # Add statistics text box
    stats_text = (
        f"Mean: {mean_density:.2f} kg/m³\n"
        f"Std Dev: {std_density:.2f} kg/m³\n"
        f"Min: {min(y):.2f} kg/m³\n"
        f"Max: {max(y):.2f} kg/m³"
    )
    plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes, fontsize=11,
            verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', 
                                             facecolor='white', alpha=0.8, edgecolor='gray'))
    
    # Add explanation text
    explanation = (
        "The density profile shows how density varies across the simulation box.\n"
        "Uniform density indicates a homogeneous system (bulk liquid).\n"
        "Variations may indicate interfaces or structural features."
    )
    plt.figtext(0.5, 0.01, explanation, ha='center', fontsize=9,
               bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
    
    # Add legend with enhanced styling
    plt.legend(fontsize=11, framealpha=0.8, loc='upper right')
    
    # Save figure with tight layout
    plt.tight_layout(rect=[0, 0.08, 1, 1])  # Adjust for the explanation text
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def main():
    if len(sys.argv) < 3:
        print("Usage: plot_density.py <analysis_dir> <plots_dir>")
        sys.exit(1)
        
    analysis_dir = sys.argv[1]
    plots_dir = sys.argv[2]
    
    # Define data directory
    data_dir = os.path.join(analysis_dir, "data")
    
    # Ensure plots directory exists
    os.makedirs(plots_dir, exist_ok=True)
    
    print('Plotting density profile...')
    
    # Density profile
    density_file = os.path.join(data_dir, 'density.xvg')
    if os.path.exists(density_file):
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(density_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Density Profile'
            plot_xlabel = xlabel if xlabel else 'Position (nm)'
            plot_ylabel = ylabel if ylabel else 'Density (kg/m³)'
            
            # Convert units if needed (GROMACS often uses g/L or amu/nm³)
            # Check if units are in the ylabel
            if plot_ylabel and ('g/l' in plot_ylabel.lower() or 'g/dm3' in plot_ylabel.lower()):
                # Convert g/L to kg/m³ (they are the same numerically)
                pass
            elif plot_ylabel and 'amu/nm3' in plot_ylabel.lower():
                # Convert amu/nm³ to kg/m³ (approximate conversion)
                y = y * 1.66053886  # 1 amu/nm³ ≈ 1.66 kg/m³
                plot_ylabel = 'Density (kg/m³)'
            
            output_path = os.path.join(plots_dir, 'density_profile_plot.png')
            plot_density_profile(x, y, plot_title, plot_xlabel, plot_ylabel, output_path)
            
            # Create histogram of density values
            output_path = os.path.join(plots_dir, 'density_histogram.png')
            plot_density_histogram((x, y), output_path)
    else:
        print(f"Density file not found: {density_file}")
    
    # Radial density map (if available)
    density_map_file = os.path.join(data_dir, 'density_radial.xpm')
    if os.path.exists(density_map_file):
        try:
            # Read XPM file
            density_map, extent = read_xpm(density_map_file)
            
            if density_map is not None:
                # Create a figure
                plt.figure(figsize=(10, 8), dpi=300)
                
                # Plot the density map with enhanced styling
                im = plt.imshow(density_map, extent=extent, origin='lower', cmap='viridis')
                
                # Add colorbar
                cbar = plt.colorbar(im)
                cbar.set_label('Density (kg/m³)', fontsize=14)
                
                # Add title and labels with enhanced styling
                plt.title('Radial Density Map', fontsize=16, fontweight='bold')
                plt.xlabel('X (nm)', fontsize=14)
                plt.ylabel('Y (nm)', fontsize=14)
                
                # Add grid and improve styling
                plt.grid(True, linestyle='--', alpha=0.3)
                plt.tick_params(axis='both', which='major', labelsize=12)
                
                # Add explanation text
                explanation = (
                    "The radial density map shows the spatial distribution of density in the system.\n"
                    "Higher values (yellow/white) indicate regions of higher density."
                )
                plt.figtext(0.5, 0.01, explanation, ha='center', fontsize=9,
                           bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
                
                # Save figure with tight layout
                output_path = os.path.join(plots_dir, 'radial_density_map.png')
                plt.tight_layout(rect=[0, 0.08, 1, 1])  # Adjust for the explanation text
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()
                print(f'  - radial_density_map.png saved successfully')
        except Exception as e:
            print(f"Error reading XPM file {density_map_file}: {e}")
    else:
        print(f"Radial density map file not found: {density_map_file}")

if __name__ == "__main__":
    main() 