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
    # Intermolecular distances
    'OO_first_peak': 0.28,   # nm, first peak in O-O RDF (2.8 Å)
    'OO_second_peak': 0.45,  # nm, second peak in O-O RDF (4.5 Å)
    'OO_third_peak': 0.68,   # nm, third peak in O-O RDF (6.8 Å)
    'OH_intermolecular_peak': 0.18,  # nm, O-H intermolecular peak (1.8 Å)
    'HH_intermolecular_peak': 0.23,  # nm, H-H intermolecular peak (2.3 Å)
    'OO_coordination': 4.5,  # Coordination number for first shell
    
    # Intramolecular distances
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

def find_peak_positions(r, g_r, height_threshold=1.5, distance_threshold=10):
    """Find the positions of peaks in the RDF"""
    try:
        peaks, properties = signal.find_peaks(g_r, height=height_threshold, distance=distance_threshold)
        peak_positions = [(r[idx], g_r[idx]) for idx in peaks]
        return sorted(peak_positions, key=lambda x: x[1], reverse=True)  # Sort by height
    except:
        # Fallback method if scipy.signal fails
        peaks = []
        for i in range(1, len(g_r)-1):
            if g_r[i] > height_threshold and g_r[i] > g_r[i-1] and g_r[i] > g_r[i+1]:
                # Check if this peak is far enough from previously found peaks
                if not peaks or min(abs(r[i] - p[0]) for p in peaks) > (r[1] - r[0]) * distance_threshold:
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

def plot_rdf(x, y, title, xlabel, ylabel, output_path, rdf_type=None):
    """
    Create a high-quality RDF plot with annotations and styling
    
    Parameters:
    -----------
    x : array
        Distance values (nm)
    y : array
        RDF values
    title : str
        Plot title
    xlabel : str
        X-axis label
    ylabel : str
        Y-axis label
    output_path : str
        Output file path
    rdf_type : str
        Type of RDF ('OO', 'OH', or 'HH') for specialized annotations
    """
    # Create a figure with enhanced styling
    plt.figure(figsize=(10, 7), dpi=300)
    
    # Plot the RDF with enhanced styling
    if HAS_SEABORN:
        sns.lineplot(x=x, y=y, color='#1f77b4', linewidth=2.5)
    else:
        plt.plot(x, y, color='#1f77b4', linewidth=2.5)
    
    # Add reference line at g(r)=1
    plt.axhline(y=1.0, color='gray', linestyle='--', alpha=0.7, linewidth=1.5, label='g(r)=1 (random)')
    
    # Find peaks for annotation
    peaks = find_peak_positions(x, y, height_threshold=1.5, distance_threshold=10)
    
    # Specialized annotations based on RDF type
    if rdf_type == 'OO':
        # O-O RDF typically shows intermolecular structure
        if peaks and len(peaks) > 0:
            # First peak (hydrogen bonding distance)
            first_peak_x, first_peak_y = peaks[0]
            plt.annotate(f'First peak: {first_peak_x:.2f} nm\n(H-bonded neighbors)',
                        xy=(first_peak_x, first_peak_y),
                        xytext=(first_peak_x - 0.05, first_peak_y + 0.3),
                        arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
                        fontsize=10, ha='center')
            
            # Try to find second peak (second hydration shell)
            second_peak_candidates = [p for p in peaks if 0.4 < p[0] < 0.5]
            if second_peak_candidates:
                second_peak_x, second_peak_y = second_peak_candidates[0]
                plt.annotate(f'Second peak: {second_peak_x:.2f} nm\n(Second hydration shell)',
                            xy=(second_peak_x, second_peak_y),
                            xytext=(second_peak_x, second_peak_y + 0.2),
                            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=6),
                            fontsize=9, ha='center')
            
            # Calculate coordination number
            first_min_x, _ = find_minima_after_peak(x, y, first_peak_x)
            if first_min_x:
                cn = calculate_coordination_number(x, y, rmin=0, rmax=first_min_x)
                plt.text(0.7, 0.85, f'Coordination number: {cn:.2f}\n(First shell, r < {first_min_x:.2f} nm)',
                        transform=plt.gca().transAxes, fontsize=10,
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    elif rdf_type == 'OH':
        # O-H RDF shows both intramolecular and intermolecular features
        if peaks and len(peaks) > 0:
            # First peak is typically intramolecular O-H bond
            first_peak_x, first_peak_y = peaks[0]
            plt.annotate(f'First peak: {first_peak_x:.2f} nm\n(Intramolecular O-H bond)',
                        xy=(first_peak_x, first_peak_y),
                        xytext=(first_peak_x - 0.05, first_peak_y * 0.8),
                        arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
                        fontsize=10, ha='center')
            
            # Look for intermolecular O-H peak (hydrogen bonding)
            intermolecular_candidates = [p for p in peaks if 0.15 < p[0] < 0.2]
            if intermolecular_candidates:
                inter_peak_x, inter_peak_y = intermolecular_candidates[0]
                plt.annotate(f'Intermolecular peak: {inter_peak_x:.2f} nm\n(H-bond distance)',
                            xy=(inter_peak_x, inter_peak_y),
                            xytext=(inter_peak_x + 0.05, inter_peak_y + 0.2),
                            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=6),
                            fontsize=9, ha='center')
            
            # Add shaded regions to distinguish intra vs intermolecular
            plt.axvspan(0, 0.12, alpha=0.2, color='lightgreen', label='Intramolecular region')
            plt.axvspan(0.12, 0.3, alpha=0.2, color='lightblue', label='Intermolecular H-bonds')
    
    elif rdf_type == 'HH':
        # H-H RDF shows both intramolecular and intermolecular features
        if peaks and len(peaks) > 0:
            # First peak is typically intramolecular H-H distance
            first_peak_x, first_peak_y = peaks[0]
            plt.annotate(f'First peak: {first_peak_x:.2f} nm\n(Intramolecular H-H distance)',
                        xy=(first_peak_x, first_peak_y),
                        xytext=(first_peak_x - 0.05, first_peak_y * 0.8),
                        arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
                        fontsize=10, ha='center')
            
            # Look for intermolecular H-H peak
            intermolecular_candidates = [p for p in peaks if 0.2 < p[0] < 0.3]
            if intermolecular_candidates:
                inter_peak_x, inter_peak_y = intermolecular_candidates[0]
                plt.annotate(f'Intermolecular peak: {inter_peak_x:.2f} nm\n(H-H between molecules)',
                            xy=(inter_peak_x, inter_peak_y),
                            xytext=(inter_peak_x + 0.05, inter_peak_y + 0.2),
                            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=6),
                            fontsize=9, ha='center')
            
            # Add shaded regions to distinguish intra vs intermolecular
            plt.axvspan(0, 0.18, alpha=0.2, color='lightgreen', label='Intramolecular region')
            plt.axvspan(0.18, 0.35, alpha=0.2, color='lightblue', label='Intermolecular region')
    
    else:
        # Generic peak annotation if RDF type not specified
        if peaks and len(peaks) > 0:
            first_peak_x, first_peak_y = peaks[0]
            plt.annotate(f'First peak: {first_peak_x:.2f} nm',
                        xy=(first_peak_x, first_peak_y),
                        xytext=(first_peak_x + 0.1, first_peak_y),
                        arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
                        fontsize=10)
    
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
    
    # Add legend if we have labeled elements
    if rdf_type in ['OH', 'HH'] or plt.gca().get_legend_handles_labels()[0]:
        plt.legend(fontsize=10, loc='upper right')
    
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
            plot_rdf(x, y, plot_title, plot_xlabel, plot_ylabel, output_path, rdf_type='OO')
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
            plot_rdf(x, y, plot_title, plot_xlabel, plot_ylabel, output_path, rdf_type='OH')
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
            plot_rdf(x, y, plot_title, plot_xlabel, plot_ylabel, output_path, rdf_type='HH')
            rdf_hh_data = (x, y, title, xlabel, ylabel, legend_labels)
    else:
        print(f"H-H RDF file not found: {rdf_hh_file}")
    
    # Create a combined RDF plot if all data is available
    if rdf_oo_data and rdf_oh_data and rdf_hh_data:
        plt.figure(figsize=(12, 8), dpi=300)
        
        # Extract data
        x_oo, y_oo = rdf_oo_data[0], rdf_oo_data[1]
        x_oh, y_oh = rdf_oh_data[0], rdf_oh_data[1]
        x_hh, y_hh = rdf_hh_data[0], rdf_hh_data[1]
        
        # Plot all RDFs
        plt.plot(x_oo, y_oo, label='O-O', linewidth=2.5, color='blue')
        plt.plot(x_oh, y_oh, label='O-H', linewidth=2.5, color='green')
        plt.plot(x_hh, y_hh, label='H-H', linewidth=2.5, color='red')
        
        # Add reference line
        plt.axhline(y=1.0, color='gray', linestyle='--', alpha=0.7, linewidth=1.5)
        
        # Add styling
        plt.xlabel('Distance (nm)', fontsize=14)
        plt.ylabel('g(r)', fontsize=14)
        plt.title('Combined Radial Distribution Functions for TIP4P Water', fontsize=16, fontweight='bold')
        plt.xlim(0, 1.0)
        plt.ylim(bottom=0)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.legend(fontsize=12, loc='upper right')
        
        # Add shaded regions to distinguish intra vs intermolecular
        plt.axvspan(0, 0.12, alpha=0.1, color='lightgray', label='_Intramolecular region')
        plt.axvspan(0.12, 0.35, alpha=0.1, color='lightyellow', label='_Intermolecular H-bonds')
        
        # Add annotations for key features
        # O-O first peak (intermolecular)
        try:
            oo_peaks = find_peak_positions(x_oo, y_oo, height_threshold=1.5)
            if oo_peaks:
                oo_peak_x, oo_peak_y = [(x, y) for x, y in oo_peaks if 0.25 < x < 0.3][0]
                plt.annotate(f'O-O: {oo_peak_x:.2f} nm',
                            xy=(oo_peak_x, oo_peak_y),
                            xytext=(oo_peak_x, oo_peak_y + 1),
                            arrowprops=dict(facecolor='blue', shrink=0.05, width=1.5, headwidth=8),
                            fontsize=10, color='blue')
        except:
            pass
        
        # O-H intramolecular peak
        try:
            oh_peaks = find_peak_positions(x_oh, y_oh, height_threshold=1.5)
            if oh_peaks:
                oh_peak_x, oh_peak_y = [(x, y) for x, y in oh_peaks if x < 0.12][0]
                plt.annotate(f'O-H: {oh_peak_x:.2f} nm',
                            xy=(oh_peak_x, oh_peak_y),
                            xytext=(oh_peak_x - 0.05, oh_peak_y - 3),
                            arrowprops=dict(facecolor='green', shrink=0.05, width=1.5, headwidth=8),
                            fontsize=10, color='green')
        except:
            pass
        
        # H-H intramolecular peak
        try:
            hh_peaks = find_peak_positions(x_hh, y_hh, height_threshold=1.5)
            if hh_peaks:
                hh_peak_x, hh_peak_y = [(x, y) for x, y in hh_peaks if 0.14 < x < 0.18][0]
                plt.annotate(f'H-H: {hh_peak_x:.2f} nm',
                            xy=(hh_peak_x, hh_peak_y),
                            xytext=(hh_peak_x + 0.05, hh_peak_y + 1),
                            arrowprops=dict(facecolor='red', shrink=0.05, width=1.5, headwidth=8),
                            fontsize=10, color='red')
        except:
            pass
        
        # Add explanation text
        explanation = (
            "The radial distribution function g(r) shows the probability of finding atoms at distance r.\n"
            "Peaks at short distances (<0.12 nm) represent intramolecular bonds within the same water molecule.\n"
            "Peaks at 0.18-0.35 nm represent intermolecular hydrogen bonding between different water molecules."
        )
        plt.figtext(0.5, 0.01, explanation, ha='center', fontsize=10,
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
        
        # Save the combined plot
        output_path = os.path.join(plots_dir, 'combined_rdf_plot.png')
        plt.tight_layout(rect=[0, 0.08, 1, 1])  # Adjust for the explanation text
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f'  - {os.path.basename(output_path)} saved successfully')
    
    print('RDF plotting completed successfully')

if __name__ == "__main__":
    main()
