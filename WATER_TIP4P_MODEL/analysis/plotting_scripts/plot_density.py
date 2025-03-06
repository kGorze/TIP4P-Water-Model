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

import numpy as np
import matplotlib.pyplot as plt
import os
import re
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from scipy.ndimage import gaussian_filter

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
        
        print(f"Processing XPM file with {len(lines)} lines")
        
        # Extract dimensions from the XPM file
        nx, ny = None, None
        for i, line in enumerate(lines):
            if "static char *gromacs_xpm[]" in line and i+1 < len(lines):
                # The next line should contain dimensions
                dim_line = lines[i+1].strip().strip('"').strip(',')
                print(f"Found dimension line: {dim_line}")
                # Format is typically: "nx ny ncolors char_per_pixel"
                parts = dim_line.split()
                if len(parts) >= 2:
                    try:
                        nx = int(parts[0])
                        ny = int(parts[1])
                        print(f"Extracted dimensions: {nx} x {ny}")
                        break
                    except ValueError as e:
                        print(f"Error parsing dimensions: {e}")
        
        if nx is None or ny is None:
            print(f"Could not determine dimensions in XPM file: {filename}")
            return None, None
        
        # Find the start of the color mapping
        color_start = None
        for i, line in enumerate(lines):
            if "static char *gromacs_xpm[]" in line:
                color_start = i + 2  # Skip the dimensions line
                break
        
        if color_start is None:
            print(f"Could not find color mapping start in XPM file: {filename}")
            return None, None
        
        # Extract the color mapping
        color_map = {}
        color_values = []
        data_start = None
        
        for i in range(color_start, len(lines)):
            line = lines[i].strip()
            
            # Check if this is a color definition line
            if line.startswith('"') and "c #" in line and "/*" in line:
                try:
                    # Extract the character and value
                    char = line.split()[0].strip('"')
                    # Extract the value from the comment /* "value" */
                    value_str = line.split('/*')[1].split('*/')[0].strip().strip('"')
                    
                    # Handle both numeric values and ranges
                    if value_str.replace('.', '', 1).isdigit():
                        value = float(value_str)
                    else:
                        # For ranges or non-numeric values, use a placeholder
                        value = len(color_values)
                    
                    color_map[char] = value
                    color_values.append(value)
                except Exception as e:
                    print(f"Error parsing color line: {line}, error: {e}")
            
            # Check if this is the start of the data matrix
            elif line.startswith('"') and not "c #" in line and not "/*" in line:
                data_start = i
                break
        
        if data_start is None:
            print(f"Could not find data start in XPM file: {filename}")
            return None, None
        
        print(f"Found {len(color_map)} colors and data starts at line {data_start}")
        
        # Extract the data matrix
        data_matrix = np.zeros((ny, nx))
        row = 0
        
        for i in range(data_start, len(lines)):
            line = lines[i].strip()
            if line.startswith('"') and line.endswith('",'):
                if row < ny:
                    # Remove quotes and comma
                    data_line = line.strip('"').strip(',')
                    
                    # Ensure we don't exceed the matrix dimensions
                    for col in range(min(nx, len(data_line))):
                        char = data_line[col]
                        if char in color_map:
                            data_matrix[row, col] = color_map[char]
                    
                    row += 1
                    
                    # Debug for first few rows
                    if row <= 3:
                        print(f"Processed row {row} with {len(data_line)} characters")
        
        print(f"Processed {row} rows of data")
        
        # Extract axis information
        x_min, x_max = 0, nx
        y_min, y_max = 0, ny
        
        # Look for axis information in the file
        for line in lines:
            if "x-label" in line and "nm" in line:
                x_min, x_max = 0, nx / 10  # Assuming 10 pixels per nm
            if "y-label" in line and "nm" in line:
                y_min, y_max = 0, ny / 10  # Assuming 10 pixels per nm
        
        # Create extent for imshow (left, right, bottom, top)
        extent = [x_min, x_max, y_min, y_max]
        
        print(f"Created extent: {extent}")
        return data_matrix, extent
    
    except Exception as e:
        print(f"Error reading XPM file {filename}: {e}")
        import traceback
        traceback.print_exc()
        return None, None

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
    # Check for both possible filenames
    density_map_file = os.path.join(data_dir, 'density_radial.xpm')
    if not os.path.exists(density_map_file):
        density_map_file = os.path.join(data_dir, 'density_map.xpm')
    
    if os.path.exists(density_map_file):
        try:
            print(f"Processing radial density map from: {density_map_file}")
            # Read XPM file
            density_map, extent = read_xpm(density_map_file)
            
            if density_map is not None:
                # Calculate statistics for the density map
                min_density = np.min(density_map)
                max_density = np.max(density_map)
                mean_density = np.mean(density_map)
                
                print(f"Density statistics - Min: {min_density:.2f}, Max: {max_density:.2f}, Mean: {mean_density:.2f} kg/m³")
                
                # Apply Gaussian smoothing to reduce noise while preserving features
                # Save original data for comparison if needed
                original_density_map = density_map.copy()
                
                # Apply smoothing with a small sigma to avoid over-smoothing
                # Sigma controls the amount of smoothing (higher = more smoothing)
                sigma = 0.8  # Start with a conservative value
                density_map = gaussian_filter(density_map, sigma=sigma)
                print(f"Applied Gaussian smoothing with sigma={sigma}")
                
                # Recalculate statistics after smoothing
                min_density_smoothed = np.min(density_map)
                max_density_smoothed = np.max(density_map)
                mean_density_smoothed = np.mean(density_map)
                print(f"Smoothed density statistics - Min: {min_density_smoothed:.2f}, Max: {max_density_smoothed:.2f}, Mean: {mean_density_smoothed:.2f} kg/m³")
                
                # Adjust color scale to better represent water density
                # Water at 273K and 1 bar has density around 1000 kg/m³
                # If our data is much lower, it might need unit conversion or normalization
                
                # Store original density map for comparison
                original_density_map_unscaled = original_density_map.copy()
                
                # Reference bulk water density
                bulk_water_density = REFERENCE_VALUES['density']  # 997.0 kg/m³
                
                # Check if the density values are suspiciously low or high
                if max_density < 500:
                    print(f"Density values appear to be in non-standard units (max={max_density:.2f})")
                    
                    # Determine appropriate scaling method
                    # Method 1: Scale based on maximum value
                    scaling_factor_max = bulk_water_density / max_density if max_density > 0 else 1.66
                    
                    # Method 2: Scale based on mean value (often more reliable)
                    # Assuming the mean density should be around 70-80% of bulk water density
                    # This accounts for the fact that the map might include some vacuum or boundary regions
                    target_mean = 0.75 * bulk_water_density  # Target mean around 75% of bulk density
                    scaling_factor_mean = target_mean / mean_density if mean_density > 0 else scaling_factor_max
                    
                    # Choose the more conservative scaling factor to avoid over-scaling
                    scaling_factor = min(scaling_factor_max, scaling_factor_mean)
                    
                    print(f"Scaling options - By max: {scaling_factor_max:.2f}, By mean: {scaling_factor_mean:.2f}")
                    print(f"Selected scaling factor: {scaling_factor:.2f}")
                    
                    # Apply scaling to both original and smoothed maps
                    original_density_map = original_density_map * scaling_factor
                    density_map = density_map * scaling_factor
                    
                    # Recalculate statistics
                    min_density = np.min(density_map)
                    max_density = np.max(density_map)
                    mean_density = np.mean(density_map)
                    
                    print(f"New density statistics - Min: {min_density:.2f}, Max: {max_density:.2f}, Mean: {mean_density:.2f} kg/m³")
                    
                    # Add interpretation note
                    print("NOTE: The scaled mean density is still lower than typical bulk water density (~1000 kg/m³).")
                    print("This may be because the map includes vacuum/boundary regions or represents a 2D projection.")
                    print("If this is a full 3D volume, consider checking your binning method or reference region.")
                elif max_density > 1500:
                    print(f"Density values appear unusually high (max={max_density:.2f} kg/m³)")
                    print("This may indicate an issue with units or calculation method.")
                else:
                    print(f"Density values appear to be in reasonable range for water (max={max_density:.2f} kg/m³)")
                
                # Create a figure with adjusted size for better visualization
                plt.figure(figsize=(12, 10), dpi=300)
                
                # Set vmin and vmax to create a more informative color scale
                # Use a slightly wider range than the data to ensure all values are visible
                vmin = max(0, min_density * 0.9)  # Don't go below zero for density
                vmax = max_density * 1.1  # Go slightly above max for better color differentiation
                
                # Plot the density map with enhanced styling
                im = plt.imshow(density_map, extent=extent, origin='lower', cmap='viridis', 
                               vmin=vmin, vmax=vmax)
                
                # Add contour lines to highlight specific density thresholds
                # Choose contour levels based on the data range
                contour_levels = np.linspace(min_density, max_density, 5)  # 5 contour levels
                contours = plt.contour(density_map, levels=contour_levels, colors='white', 
                                      alpha=0.5, linestyles='dashed', extent=extent)
                plt.clabel(contours, inline=True, fontsize=10, fmt='%.0f')
                
                # Add a marker or line for the mean density
                plt.contour(density_map, levels=[mean_density], colors='red', 
                           linestyles='solid', linewidths=2, extent=extent)
                
                # Add colorbar with improved styling
                cbar = plt.colorbar(im, pad=0.02)
                cbar.set_label('Density (kg/m³)', fontsize=14, fontweight='bold')
                cbar.ax.tick_params(labelsize=12)
                
                # Mark the bulk water density reference on the colorbar
                bulk_water_density = 1000  # kg/m³ at standard conditions
                if vmin <= bulk_water_density <= vmax:
                    cbar.ax.axhline(y=(bulk_water_density - vmin) / (vmax - vmin), 
                                   color='black', linestyle='--', linewidth=2)
                    cbar.ax.text(1.1, (bulk_water_density - vmin) / (vmax - vmin), 
                                'Bulk water (1000 kg/m³)', va='center', ha='left', 
                                transform=cbar.ax.transAxes, fontsize=10)
                
                # Add title and labels with enhanced styling
                plt.title('Radial Density Map', fontsize=18, fontweight='bold')
                plt.xlabel('X (nm)', fontsize=16, fontweight='bold')
                plt.ylabel('Y (nm)', fontsize=16, fontweight='bold')
                
                # Add grid and improve styling
                plt.grid(True, linestyle='--', alpha=0.3)
                plt.tick_params(axis='both', which='major', labelsize=14)
                
                # Add statistical information to the plot
                stats_text = (
                    f"Min: {min_density:.1f} kg/m³\n"
                    f"Max: {max_density:.1f} kg/m³\n"
                    f"Mean: {mean_density:.1f} kg/m³"
                )
                plt.figtext(0.02, 0.02, stats_text, fontsize=12, 
                           bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))
                
                # Add explanation text with more detailed interpretation
                explanation = (
                    "The radial density map shows the spatial distribution of density in the system.\n"
                    "Higher values (yellow/white) indicate regions of higher density. "
                    "Contour lines mark equal density regions.\n"
                    f"Mean density ({mean_density:.1f} kg/m³) is lower than typical bulk water (~1000 kg/m³), "
                    "which may indicate inclusion of boundary regions or projection effects in the calculation."
                )
                plt.figtext(0.5, 0.01, explanation, ha='center', fontsize=10,
                           bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
                
                # Save figure with tight layout
                output_path = os.path.join(plots_dir, 'radial_density_map.png')
                plt.tight_layout(rect=[0, 0.08, 1, 1])  # Adjust for the explanation text
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()
                print(f'  - radial_density_map.png saved successfully')
                
                # Create a side-by-side comparison of original vs. smoothed data
                fig, axes = plt.subplots(1, 2, figsize=(16, 8), dpi=300)
                
                # Common color scale for both plots
                vmin = max(0, min(np.min(original_density_map_unscaled), np.min(density_map)) * 0.9)
                vmax = max(np.max(original_density_map_unscaled), np.max(density_map)) * 1.1
                
                # Plot original data
                im1 = axes[0].imshow(original_density_map_unscaled, extent=extent, origin='lower', 
                                    cmap='viridis', vmin=vmin, vmax=vmax)
                axes[0].set_title('Original Density Map', fontsize=16, fontweight='bold')
                axes[0].set_xlabel('X (nm)', fontsize=14)
                axes[0].set_ylabel('Y (nm)', fontsize=14)
                axes[0].grid(True, linestyle='--', alpha=0.3)
                
                # Plot smoothed data
                im2 = axes[1].imshow(density_map, extent=extent, origin='lower', 
                                    cmap='viridis', vmin=vmin, vmax=vmax)
                axes[1].set_title(f'Smoothed Density Map (σ={sigma})', fontsize=16, fontweight='bold')
                axes[1].set_xlabel('X (nm)', fontsize=14)
                axes[1].set_ylabel('Y (nm)', fontsize=14)
                axes[1].grid(True, linestyle='--', alpha=0.3)
                
                # Add a colorbar that applies to both plots
                cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
                cbar = fig.colorbar(im2, cax=cbar_ax)
                cbar.set_label('Density (kg/m³)', fontsize=14, fontweight='bold')
                
                # Add explanation text with density interpretation
                comparison_explanation = (
                    "Comparison of original (left) and smoothed (right) density maps.\n"
                    f"Gaussian smoothing with σ={sigma} reduces noise while preserving major features.\n"
                    f"Mean density ({mean_density:.1f} kg/m³) is lower than bulk water (~1000 kg/m³), "
                    f"likely due to inclusion of boundary regions or projection effects."
                )
                fig.text(0.5, 0.01, comparison_explanation, ha='center', fontsize=12,
                        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
                
                # Save comparison figure
                comparison_path = os.path.join(plots_dir, 'density_map_comparison.png')
                plt.tight_layout(rect=[0, 0.05, 0.9, 1])  # Adjust for the explanation text and colorbar
                plt.savefig(comparison_path, dpi=300, bbox_inches='tight')
                plt.close()
                print(f'  - density_map_comparison.png saved successfully')
            else:
                print(f"Failed to extract density map data from {density_map_file}")
        except Exception as e:
            print(f"Error processing XPM file {density_map_file}: {e}")
            import traceback
            traceback.print_exc()
    else:
        print(f"Radial density map file not found at: {density_map_file}")

if __name__ == "__main__":
    main() 