#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

def read_xvg(file_path):
    """Read XVG file and return time and data arrays"""
    time = []
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            values = line.strip().split()
            if values:
                time.append(float(values[0]))
                data.append([float(x) for x in values[1:]])
    return np.array(time), np.array(data)

def plot_rdf_article_style(data_dir, model, temp):
    """Plot radial distribution function in the style of the article
    
    Parameters:
    -----------
    data_dir : Path
        Directory containing the RDF data
    model : str
        Water model name (e.g., 'tip4p', 'tip4p-ice')
    temp : int
        Temperature in K
    """
    rdf_file = data_dir / 'rdf_oo.xvg'
    if not rdf_file.exists():
        print(f"RDF file not found: {rdf_file}")
        return

    # Read the RDF data
    distance_nm, data = read_xvg(rdf_file)
    
    # Convert nm to Å (1 nm = 10 Å)
    distance_angstrom = distance_nm * 10.0
    
    # Create the plot with article-like styling
    plt.figure(figsize=(8, 6))
    
    # Set font sizes to match article style
    plt.rcParams.update({
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 14,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
    })
    
    # Plot the data with appropriate line style
    if model.lower() == 'tip4p':
        line_style = 'b-.'  # Blue dash-dot line for TIP4P
        label = 'TIP4P'
    elif model.lower() == 'tip4p-ice':
        line_style = 'b--'  # Blue dashed line for TIP4P/Ice
        label = 'TIP4P/ice'
    else:
        line_style = 'b-'  # Blue solid line for other models
        label = model
    
    plt.plot(distance_angstrom, data, line_style, label=label, linewidth=1.5)
    
    # Set axis labels and title in article style
    plt.xlabel('r/Å')
    plt.ylabel('g')
    plt.title(f'Oxygen-oxygen correlation functions at T={temp} K')
    
    # Set x-axis range to match the article (0-8 Å)
    plt.xlim(0, 8)
    
    # Set y-axis range to match the article (0-4)
    plt.ylim(0, 4)
    
    # Add grid and legend
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    
    # Make the plot more compact like in the article
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(data_dir / 'rdf_plot_article_style.png', dpi=300)
    print(f"Plot saved to {data_dir / 'rdf_plot_article_style.png'}")
    
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Plot RDF in article style")
    parser.add_argument("--model", type=str, choices=["tip4p", "tip4p-ice"], required=True,
                        help="Water model used")
    parser.add_argument("--temp", type=int, choices=[150, 200, 273, 298], required=True,
                        help="Temperature simulated")
    args = parser.parse_args()
    
    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / "data" / args.model / f"{args.temp}K"
    
    if not data_dir.exists():
        print(f"Data directory not found: {data_dir}")
        return 1
    
    print(f"Plotting RDF from {data_dir} in article style")
    plot_rdf_article_style(data_dir, args.model, args.temp)
    
    return 0

if __name__ == "__main__":
    exit(main()) 