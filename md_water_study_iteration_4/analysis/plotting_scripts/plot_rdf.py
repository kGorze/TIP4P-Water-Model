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

def main():
    if len(sys.argv) < 3:
        print("Usage: plot_rdf.py <analysis_dir> <plots_dir>")
        sys.exit(1)
        
    analysis_dir = sys.argv[1]
    plots_dir = sys.argv[2]
    
    # Ensure plots directory exists
    os.makedirs(plots_dir, exist_ok=True)
    
    # RDF Plots (All three: O-O, O-H, H-H)
    print('Plotting RDF data...')
    rdf_files = {
        'rdf_OO.xvg': ('Oxygen-Oxygen Radial Distribution Function', '#1f77b4', 'O-O RDF'),
        'rdf_OH.xvg': ('Oxygen-Hydrogen Radial Distribution Function', '#2ca02c', 'O-H RDF'),
        'rdf_HH.xvg': ('Hydrogen-Hydrogen Radial Distribution Function', '#d62728', 'H-H RDF')
    }

    # Store data for combined plot
    combined_data = {}
    
    # Plot individual RDF plots
    for filename, (title, color, label) in rdf_files.items():
        filepath = os.path.join(analysis_dir, filename)
        if os.path.exists(filepath):
            x, y, file_title, file_xlabel, file_ylabel = read_xvg(filepath)
            
            # Use file metadata if available, otherwise use defaults
            plot_title = file_title if file_title else title
            xlabel = file_xlabel if file_xlabel else 'r (nm)'
            ylabel = file_ylabel if file_ylabel else 'g(r)'
            
            if len(x) > 0 and len(y) > 0:
                output_filename = f'{filename.split(".")[0]}_plot.png'
                
                # Store data for combined plot
                combined_data[filename] = (x, y, label)
                
                # Create a high-quality figure
                plt.figure(figsize=(10, 6), dpi=300)
                if HAS_SEABORN:
                    sns.lineplot(x=x, y=y, color=color, linewidth=2.5, label=label)
                else:
                    plt.plot(x, y, color=color, linewidth=2.5, label=label)
                
                # Add literature value vertical lines
                if 'OO' in filename:
                    lit_peak = LITERATURE_VALUES['OO_first_peak']
                    plt.axvline(x=lit_peak, color='gray', linestyle='--', alpha=0.7,
                               label=f'Literature peak: {lit_peak} nm')
                    
                    # Calculate coordination number
                    # Find first peak
                    peaks = find_peak_positions(x, y)
                    if peaks:
                        first_peak = peaks[0][0]
                        # Find first minimum after peak
                        min_r, _ = find_minima_after_peak(x, y, first_peak)
                        if min_r:
                            cn = calculate_coordination_number(x, y, rmin=0, rmax=min_r)
                            plt.annotate(f'Coordination number: {cn:.1f}',
                                        xy=(min_r, 1.0),
                                        xytext=(min_r+0.05, 1.5),
                                        arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                                        fontsize=10)
                            
                            # Add shaded area for first coordination shell
                            plt.axvspan(0, min_r, alpha=0.2, color='blue', label='First coordination shell')
                            
                            # Plot the running coordination number
                            ax1 = plt.gca()
                            ax2 = ax1.twinx()
                            
                            # Calculate running coordination number
                            cn_values = [calculate_coordination_number(x, y, rmin=0, rmax=r_val) for r_val in x]
                            ax2.plot(x, cn_values, 'r--', alpha=0.7, label='Coordination number')
                            ax2.set_ylabel('Coordination number', color='r', fontsize=12)
                            ax2.tick_params(axis='y', labelcolor='r')
                            ax2.set_ylim(0, max(cn_values)*1.1)
                            
                            # Add horizontal line at literature coordination number
                            ax2.axhline(y=LITERATURE_VALUES['OO_coordination'], color='r', linestyle=':', alpha=0.7,
                                      label=f'Literature CN: {LITERATURE_VALUES["OO_coordination"]}')
                
                elif 'OH' in filename:
                    lit_peak = LITERATURE_VALUES['OH_peak']
                    plt.axvline(x=lit_peak, color='gray', linestyle='--', alpha=0.7,
                               label=f'Literature peak: {lit_peak} nm')
                
                elif 'HH' in filename:
                    lit_peak = LITERATURE_VALUES['HH_peak']
                    plt.axvline(x=lit_peak, color='gray', linestyle='--', alpha=0.7,
                               label=f'Literature peak: {lit_peak} nm')
                
                plt.xlabel(xlabel, fontsize=14)
                plt.ylabel(ylabel, fontsize=14)
                plt.title(plot_title, fontsize=16)
                plt.grid(True, alpha=0.3)
                plt.legend(fontsize=10, loc='best')
                
                # Add a subtle background color if seaborn is not available
                if not HAS_SEABORN:
                    plt.gca().set_facecolor('#f8f8f8')
                
                # Ensure the plot shows the important features
                if 'OO' in filename:
                    plt.ylim(bottom=0)
                    # Zoom in on the interesting region if there's a clear peak
                    if np.max(y) > 2:
                        plt.xlim(0, min(2.0, np.max(x)))
                
                plt.tight_layout()
                plt.savefig(os.path.join(plots_dir, output_filename), dpi=300, bbox_inches='tight')
                plt.close()
                print(f'  - {output_filename} saved successfully')
                
                # Remove any duplicate lowercase files to avoid confusion
                lowercase_file = os.path.join(plots_dir, f'{filename.lower().split(".")[0]}_plot.png')
                if os.path.exists(lowercase_file) and lowercase_file != os.path.join(plots_dir, output_filename):
                    try:
                        os.remove(lowercase_file)
                        print(f'  - Removed duplicate file: {os.path.basename(lowercase_file)}')
                    except:
                        pass

    # Plot combined RDF plot with enhanced styling
    if combined_data:
        plt.figure(figsize=(12, 8), dpi=300)
        colors = ['#1f77b4', '#2ca02c', '#d62728']  # Blue, Green, Red
        
        # Plot all RDFs
        for idx, (filename, (x, y, label)) in enumerate(combined_data.items()):
            if HAS_SEABORN:
                sns.lineplot(x=x, y=y, color=colors[idx % len(colors)], linewidth=2.5, label=label)
            else:
                plt.plot(x, y, color=colors[idx % len(colors)], linewidth=2.5, label=label)
        
        plt.xlabel('r (nm)', fontsize=14)
        plt.ylabel('g(r)', fontsize=14)
        plt.title('Radial Distribution Functions for TIP4P Water', fontsize=16)
        
        # Add vertical lines for intramolecular distances (dashed gray)
        plt.axvline(x=LITERATURE_VALUES['OH_intramolecular'], color='gray', linestyle=':', alpha=0.5)
        plt.axvline(x=LITERATURE_VALUES['HH_intramolecular'], color='gray', linestyle=':', alpha=0.5)
        
        # Add vertical lines for intermolecular peaks (dashed)
        plt.axvline(x=LITERATURE_VALUES['OO_first_peak'], color='gray', linestyle='--', alpha=0.5)
        plt.axvline(x=LITERATURE_VALUES['OH_peak'], color='gray', linestyle='--', alpha=0.5)
        plt.axvline(x=LITERATURE_VALUES['HH_peak'], color='gray', linestyle='--', alpha=0.5)
        
        # Improve the legend
        plt.legend(fontsize=12, framealpha=0.9, loc='upper right')
        
        # Focus on the interesting region (first 2 nm)
        plt.xlim(0, 2.0)
        plt.ylim(bottom=0)
        
        # Add annotations for key features
        if 'rdf_OO.xvg' in combined_data:
            x, y, _ = combined_data['rdf_OO.xvg']
            # Find the first major peak
            peaks = find_peak_positions(x, y)
            if peaks:
                first_peak = peaks[0]
                plt.annotate(f'First peak: {first_peak[0]:.2f} nm',
                            xy=(first_peak[0], first_peak[1]),
                            xytext=(first_peak[0]+0.1, first_peak[1]),
                            arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
                            fontsize=10)
        
        # Add a note about intramolecular vs intermolecular distances
        note_text = (
            "Note: Intramolecular distances (within same water molecule):\n"
            f"O-H bond: ~{LITERATURE_VALUES['OH_intramolecular']:.2f} nm, "
            f"H-H distance: ~{LITERATURE_VALUES['HH_intramolecular']:.2f} nm\n"
            "These should be excluded from intermolecular RDF analysis."
        )
        plt.figtext(0.5, 0.05, note_text, ha='center', fontsize=10, 
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Add a watermark with simulation details
        plt.figtext(0.5, 0.01, 'TIP4P Water Model, 298K, 1 bar', 
                   ha='center', fontsize=10, style='italic', alpha=0.7)
        
        plt.tight_layout(rect=[0, 0.1, 1, 1])  # Adjust layout to make room for the note
        plt.savefig(os.path.join(plots_dir, 'combined_rdf_plot.png'), dpi=300, bbox_inches='tight')
        plt.close()
        print('  - Combined RDF plot saved successfully')

if __name__ == "__main__":
    main()
