#!/usr/bin/python3
import os
import sys
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
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

# Literature values for TIP4P water model at 298K
LITERATURE_VALUES = {
    'OO_first_peak': 0.28,   # nm, first peak in O-O RDF (2.8 Å)
    'OO_second_peak': 0.45,  # nm, second peak in O-O RDF (4.5 Å)
    'OH_peak': 0.18,         # nm, O-H intermolecular peak (1.8 Å)
    'HH_peak': 0.23,         # nm, H-H intermolecular peak (2.3 Å)
    'OO_coordination': 4.5,  # Coordination number for first shell
    'OH_intramolecular': 0.10, # nm, O-H bond length within water molecule (1.0 Å)
    'HH_intramolecular': 0.16, # nm, H-H distance within water molecule (1.6 Å)
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

def calculate_coordination_number(r, g_r, rmin=0, rmax=0.35, density=33.3679):
    """
    Calculate coordination number by integrating the RDF
    
    Parameters:
    -----------
    r : array
        Radial distance values
    g_r : array
        RDF values
    rmin : float
        Lower integration bound (nm)
    rmax : float
        Upper integration bound (nm), typically the first minimum after the first peak
    density : float
        Number density of water molecules (nm^-3), default for TIP4P at 298K
        
    Returns:
    --------
    cn : float
        Coordination number
    """
    # Filter values within the range
    mask = (r >= rmin) & (r <= rmax)
    r_range = r[mask]
    g_r_range = g_r[mask]
    
    if len(r_range) < 2:
        return 0
    
    # Calculate the coordination number
    # CN = 4πρ∫(r² g(r) dr) from rmin to rmax
    integrand = 4 * np.pi * density * r_range**2 * g_r_range
    dr = r_range[1] - r_range[0]
    cn = np.sum(integrand) * dr
    
    return cn

def find_peak_positions(r, g_r, height_threshold=1.5):
    """Find the positions of peaks in the RDF"""
    peaks = []
    for i in range(1, len(g_r)-1):
        if g_r[i] > height_threshold and g_r[i] > g_r[i-1] and g_r[i] > g_r[i+1]:
            peaks.append((r[i], g_r[i]))
    return sorted(peaks, key=lambda x: x[1], reverse=True)  # Sort by height

def find_minima_after_peak(r, g_r, peak_position, search_range=0.2):
    """Find the first minimum after a peak in the RDF"""
    peak_idx = np.abs(r - peak_position).argmin()
    search_end = min(len(r)-1, peak_idx + int(search_range / (r[1] - r[0])))
    
    for i in range(peak_idx, search_end):
        if i+1 < len(g_r) and g_r[i] < g_r[i+1]:
            return r[i], g_r[i]
    
    return None, None

def plot_rdf(x, y, title, xlabel, ylabel, output_path):
    """Create a high-quality RDF plot with annotations and styling"""
    # Create a figure with enhanced styling
    plt.figure(figsize=(10, 6), dpi=300)
    
    # Plot the RDF with enhanced styling
    if HAS_SEABORN:
        sns.lineplot(x=x, y=y, color='#1f77b4', linewidth=2.5)
    else:
        plt.plot(x, y, color='#1f77b4', linewidth=2.5)
    
    # Add reference line at g(r)=1
    plt.axhline(y=1.0, color='gray', linestyle='--', alpha=0.7, linewidth=1.5, label='g(r)=1 (random)')
    
    # Find peaks for annotation
    try:
        peaks, _ = signal.find_peaks(y, height=1.5, distance=10)
        if len(peaks) > 0:
            # Annotate the first peak
            first_peak_x = x[peaks[0]]
            first_peak_y = y[peaks[0]]
            plt.annotate(f'First peak: {first_peak_x:.2f} nm',
                        xy=(first_peak_x, first_peak_y),
                        xytext=(first_peak_x + 0.1, first_peak_y),
                        arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
                        fontsize=10)
    except:
        # If scipy.signal is not available or peak finding fails, continue without annotations
        pass
    
    # Set axis labels and title with enhanced styling
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=16, fontweight='bold')
    
    # Set axis limits to focus on relevant features
    plt.xlim(0, 1.0)  # Focus on 0-1 nm range where most features are
    plt.ylim(bottom=0)  # Start y-axis at 0
    
    # Add grid and improve styling
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tick_params(axis='both', which='major', labelsize=12)
    
    # Add a text box with explanation
    explanation = (
        "The radial distribution function g(r) describes how\n"
        "density varies as a function of distance from a reference particle.\n"
        "g(r) = 1 indicates random distribution (bulk density)."
    )
    plt.figtext(0.5, 0.01, explanation, ha='center', fontsize=9,
               bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
    
    # Save the figure with tight layout
    plt.tight_layout(rect=[0, 0.08, 1, 1])  # Adjust for the explanation text
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  - {os.path.basename(output_path)} saved successfully')

def main():
    if len(sys.argv) < 3:
        print("Usage: plot_rdf.py <analysis_dir> <plots_dir>")
        sys.exit(1)
        
    analysis_dir = sys.argv[1]
    plots_dir = sys.argv[2]
    
    # Define data directory
    data_dir = os.path.join(analysis_dir, "data")
    
    # Ensure plots directory exists
    os.makedirs(plots_dir, exist_ok=True)
    
    print('Plotting RDF data...')
    
    # Read O-O RDF data
    rdf_oo_file = os.path.join(data_dir, 'rdf_OO.xvg')
    rdf_oo_data = None
    if os.path.exists(rdf_oo_file):
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(rdf_oo_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Oxygen-Oxygen Radial Distribution Function'
            plot_xlabel = xlabel if xlabel else 'Distance (nm)'
            plot_ylabel = ylabel if ylabel else 'g(r)'
            
            output_path = os.path.join(plots_dir, 'rdf_OO_plot.png')
            plot_rdf(x, y, plot_title, plot_xlabel, plot_ylabel, output_path)
            rdf_oo_data = (x, y, title, xlabel, ylabel, legend_labels)
    else:
        print(f"O-O RDF file not found: {rdf_oo_file}")
    
    # Read O-H RDF data
    rdf_oh_file = os.path.join(data_dir, 'rdf_OH.xvg')
    rdf_oh_data = None
    if os.path.exists(rdf_oh_file):
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(rdf_oh_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Oxygen-Hydrogen Radial Distribution Function'
            plot_xlabel = xlabel if xlabel else 'Distance (nm)'
            plot_ylabel = ylabel if ylabel else 'g(r)'
            
            output_path = os.path.join(plots_dir, 'rdf_OH_plot.png')
            plot_rdf(x, y, plot_title, plot_xlabel, plot_ylabel, output_path)
            rdf_oh_data = (x, y, title, xlabel, ylabel, legend_labels)
    else:
        print(f"O-H RDF file not found: {rdf_oh_file}")
    
    # Read H-H RDF data
    rdf_hh_file = os.path.join(data_dir, 'rdf_HH.xvg')
    rdf_hh_data = None
    if os.path.exists(rdf_hh_file):
        x, y, title, xlabel, ylabel, legend_labels = read_xvg(rdf_hh_file)
        
        if len(x) > 0 and len(y) > 0:
            # Use file metadata if available, otherwise use defaults
            plot_title = title if title else 'Hydrogen-Hydrogen Radial Distribution Function'
            plot_xlabel = xlabel if xlabel else 'Distance (nm)'
            plot_ylabel = ylabel if ylabel else 'g(r)'
            
            output_path = os.path.join(plots_dir, 'rdf_HH_plot.png')
            plot_rdf(x, y, plot_title, plot_xlabel, plot_ylabel, output_path)
            rdf_hh_data = (x, y, title, xlabel, ylabel, legend_labels)
    else:
        print(f"H-H RDF file not found: {rdf_hh_file}")
    
    # Create combined RDF plot if we have at least one RDF dataset
    if any([rdf_oo_data, rdf_oh_data, rdf_hh_data]):
        output_path = os.path.join(plots_dir, 'combined_rdf_plot.png')
        create_combined_rdf_plot([rdf_oo_data, rdf_oh_data, rdf_hh_data], output_path)
    else:
        print("No RDF data found to create combined plot")

def create_combined_rdf_plot(rdf_data, output_path):
    """Create a combined plot of all RDF data with enhanced visualization"""
    # Extract data
    rdf_oo_data, rdf_oh_data, rdf_hh_data = rdf_data
    
    # Create figure with enhanced styling
    fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
    
    # Set a consistent style
    ax.grid(True, linestyle='--', alpha=0.7)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    
    # Plot each RDF with enhanced styling and increased line thickness
    if rdf_oo_data:
        x_oo, y_oo, _, _, _, _ = rdf_oo_data
        ax.plot(x_oo, y_oo, label='O-O', color='#1f77b4', linewidth=3.0)
        
        # Calculate coordination number for O-O
        try:
            # Find first minimum after first peak to set integration limit
            peak_indices = signal.find_peaks(y_oo, height=1.5, distance=10)[0]
            if len(peak_indices) > 0:
                first_peak_idx = peak_indices[0]
                first_peak_x = x_oo[first_peak_idx]
                
                # Find the first minimum after the peak
                min_after_peak = None
                for i in range(first_peak_idx, len(y_oo)-1):
                    if y_oo[i] < y_oo[i+1]:
                        min_after_peak = x_oo[i]
                        break
                
                if min_after_peak:
                    # Calculate coordination number
                    cn = calculate_coordination_number(x_oo, y_oo, rmin=0, rmax=min_after_peak)
                    
                    # Add coordination number annotation
                    ax.annotate(f'Coordination number: {cn:.1f}',
                              xy=(min_after_peak, 1.0),
                              xytext=(min_after_peak + 0.1, 2.0),
                              arrowprops=dict(facecolor='#1f77b4', shrink=0.05, width=1.5, headwidth=8),
                              fontsize=10, color='#1f77b4')
        except:
            # If peak finding or coordination number calculation fails, continue without it
            pass
    
    if rdf_oh_data:
        x_oh, y_oh, _, _, _, _ = rdf_oh_data
        ax.plot(x_oh, y_oh, label='O-H', color='#ff7f0e', linewidth=3.0)
        
        # Add annotation for intramolecular O-H peak
        try:
            # Find the highest peak (usually the intramolecular O-H bond)
            max_idx = np.argmax(y_oh)
            intra_oh_x = x_oh[max_idx]
            intra_oh_y = y_oh[max_idx]
            
            # Add annotation for intramolecular O-H bond
            ax.annotate(f'O-H bond: {intra_oh_x:.2f} nm',
                       xy=(intra_oh_x, intra_oh_y),
                       xytext=(intra_oh_x + 0.05, intra_oh_y - 20),
                       arrowprops=dict(facecolor='#ff7f0e', shrink=0.05, width=1.5, headwidth=8),
                       fontsize=10, color='#ff7f0e')
        except:
            # If annotation fails, continue without it
            pass
    
    if rdf_hh_data:
        x_hh, y_hh, _, _, _, _ = rdf_hh_data
        ax.plot(x_hh, y_hh, label='H-H', color='#2ca02c', linewidth=3.0)
        
        # Add annotation for intramolecular H-H peak
        try:
            # Find the highest peak (usually the intramolecular H-H distance)
            max_idx = np.argmax(y_hh)
            intra_hh_x = x_hh[max_idx]
            intra_hh_y = y_hh[max_idx]
            
            # Add annotation for intramolecular H-H distance
            ax.annotate(f'H-H distance: {intra_hh_x:.2f} nm',
                       xy=(intra_hh_x, intra_hh_y),
                       xytext=(intra_hh_x + 0.05, intra_hh_y - 5),
                       arrowprops=dict(facecolor='#2ca02c', shrink=0.05, width=1.5, headwidth=8),
                       fontsize=10, color='#2ca02c')
        except:
            # If annotation fails, continue without it
            pass
    
    # Add reference line at g(r)=1
    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.7, linewidth=1.5, label='g(r)=1 (random)')
    
    # Add vertical markers for literature values
    # O-O first peak
    ax.axvline(x=LITERATURE_VALUES['OO_first_peak'], color='#1f77b4', linestyle=':', alpha=0.5)
    ax.text(LITERATURE_VALUES['OO_first_peak'] + 0.01, ax.get_ylim()[1]*0.9, 
           f"O-O (lit): {LITERATURE_VALUES['OO_first_peak']} nm", 
           rotation=90, color='#1f77b4', fontsize=9, alpha=0.7)
    
    # O-O second peak
    ax.axvline(x=LITERATURE_VALUES['OO_second_peak'], color='#1f77b4', linestyle=':', alpha=0.5)
    ax.text(LITERATURE_VALUES['OO_second_peak'] + 0.01, ax.get_ylim()[1]*0.9, 
           f"O-O (lit): {LITERATURE_VALUES['OO_second_peak']} nm", 
           rotation=90, color='#1f77b4', fontsize=9, alpha=0.7)
    
    # O-H intramolecular
    ax.axvline(x=LITERATURE_VALUES['OH_intramolecular'], color='#ff7f0e', linestyle=':', alpha=0.5)
    ax.text(LITERATURE_VALUES['OH_intramolecular'] + 0.01, ax.get_ylim()[1]*0.8, 
           f"O-H bond (lit): {LITERATURE_VALUES['OH_intramolecular']} nm", 
           rotation=90, color='#ff7f0e', fontsize=9, alpha=0.7)
    
    # H-H intramolecular
    ax.axvline(x=LITERATURE_VALUES['HH_intramolecular'], color='#2ca02c', linestyle=':', alpha=0.5)
    ax.text(LITERATURE_VALUES['HH_intramolecular'] + 0.01, ax.get_ylim()[1]*0.7, 
           f"H-H (lit): {LITERATURE_VALUES['HH_intramolecular']} nm", 
           rotation=90, color='#2ca02c', fontsize=9, alpha=0.7)
    
    # Add solvation shell annotations if O-O data is available
    if rdf_oo_data:
        # First solvation shell (typically around 0.28 nm)
        first_shell_x = LITERATURE_VALUES['OO_first_peak']
        ax.axvline(x=first_shell_x, color='gray', linestyle=':', alpha=0.7)
        ax.text(first_shell_x + 0.02, 0.5, '1st shell', rotation=90, alpha=0.7)
        
        # Second solvation shell (typically around 0.45 nm)
        second_shell_x = LITERATURE_VALUES['OO_second_peak']
        ax.axvline(x=second_shell_x, color='gray', linestyle=':', alpha=0.7)
        ax.text(second_shell_x + 0.02, 0.5, '2nd shell', rotation=90, alpha=0.7)
    
    # Add title and labels with enhanced styling
    ax.set_title('Radial Distribution Functions for TIP4P Water', fontsize=16, fontweight='bold')
    ax.set_xlabel('Distance (nm)', fontsize=14)
    ax.set_ylabel('g(r)', fontsize=14)
    
    # Set axis limits to focus on relevant features but extend to show more of second shell
    ax.set_xlim(0, 1.5)  # Extended to 1.5 nm to show more of the second coordination shell
    
    # Add legend with enhanced styling
    ax.legend(fontsize=12, framealpha=0.8, loc='upper right')
    
    # Add text box with explanation and simulation conditions
    explanation = (
        "The radial distribution function g(r) describes how\n"
        "density varies as a function of distance from a reference particle.\n"
        "g(r) = 1 indicates random distribution (bulk density).\n"
        "Peaks indicate preferred coordination distances."
    )
    ax.text(0.5, 0.05, explanation, transform=ax.transAxes, fontsize=10,
           bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'),
           ha='center')
    
    # Add simulation conditions
    sim_conditions = (
        f"TIP4P Water Model | T=298K | P=1 bar | 5500 molecules"
    )
    plt.figtext(0.5, 0.01, sim_conditions, ha='center', fontsize=10, 
               style='italic', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))
    
    # Add an inset for zoomed view of O-O interactions if O-O data is available
    if rdf_oo_data:
        # Create inset axes
        axins = fig.add_axes([0.2, 0.55, 0.3, 0.3])  # [left, bottom, width, height]
        
        # Plot O-O RDF in the inset
        axins.plot(x_oo, y_oo, color='#1f77b4', linewidth=2.0)
        
        # Set limits for the inset (focus on the first and second peaks)
        axins.set_xlim(0.25, 0.6)  # Range covering first and second O-O peaks
        
        # Find appropriate y-limit
        if len(x_oo) > 0:
            mask = (x_oo >= 0.25) & (x_oo <= 0.6)
            if np.any(mask):
                max_y = np.max(y_oo[mask]) * 1.1  # Add 10% margin
                axins.set_ylim(0, max_y)
        
        # Add grid to inset
        axins.grid(True, linestyle='--', alpha=0.5)
        
        # Add title to inset
        axins.set_title('O-O Interactions Detail', fontsize=10)
        
        # Add reference line at g(r)=1 in the inset
        axins.axhline(y=1.0, color='gray', linestyle='--', alpha=0.7, linewidth=1.0)
        
        # Add vertical markers for literature values in the inset
        axins.axvline(x=LITERATURE_VALUES['OO_first_peak'], color='#1f77b4', linestyle=':', alpha=0.5)
        axins.axvline(x=LITERATURE_VALUES['OO_second_peak'], color='#1f77b4', linestyle=':', alpha=0.5)
        
        # Mark the inset area in the main plot
        from matplotlib.patches import ConnectionPatch
        xy = (0.25, 0)
        xy2 = (0.6, max_y if 'max_y' in locals() else 3)
        con = ConnectionPatch(xyA=xy, xyB=xy2, coordsA="data", coordsB="data",
                             axesA=ax, axesB=axins, color="gray", alpha=0.3)
        ax.add_artist(con)
        
        xy = (0.25, max_y if 'max_y' in locals() else 3)
        xy2 = (0.6, 0)
        con = ConnectionPatch(xyA=xy, xyB=xy2, coordsA="data", coordsB="data",
                             axesA=ax, axesB=axins, color="gray", alpha=0.3)
        ax.add_artist(con)
    
    # Add grid and tick parameters
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    # Save figure with tight layout
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])  # Adjust for the simulation conditions text
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print('  - Combined RDF plot saved successfully')

if __name__ == "__main__":
    main()
