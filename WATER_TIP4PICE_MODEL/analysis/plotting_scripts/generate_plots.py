#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import subprocess
import argparse

def main():
    """
    Main function to generate plots for TIP4P water model analysis.
    This script now serves as a wrapper for the more comprehensive run_all_plots.py script.
    """
    parser = argparse.ArgumentParser(description="Generate plots for TIP4P water model analysis")
    parser.add_argument("--analysis-dir", default=None, help="Directory containing analysis files")
    parser.add_argument("--plots-dir", default=None, help="Directory to save plots")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output")
    args = parser.parse_args()
    
    # Get the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # If analysis_dir is not provided, use the directory of this script
    if args.analysis_dir is None:
        analysis_dir = script_dir
    else:
        analysis_dir = args.analysis_dir
    
    # If plots_dir is not provided, create a plots directory in the analysis directory
    if args.plots_dir is None:
        plots_dir = os.path.join(analysis_dir, "plots")
    else:
        plots_dir = args.plots_dir
    
    # Create the plots directory if it doesn't exist
    os.makedirs(plots_dir, exist_ok=True)
    
    print(f"Analysis directory: {analysis_dir}")
    print(f"Plots directory: {plots_dir}")
    
    # Check if the run_all_plots.py script exists
    coordinator_script = os.path.join(script_dir, "plotting_scripts", "run_all_plots.py")
    if os.path.exists(coordinator_script):
        print("Using the comprehensive plotting system...")
        
        # Construct the command
        cmd = [
            sys.executable, 
            coordinator_script, 
            "--analysis-dir", analysis_dir, 
            "--plots-dir", plots_dir
        ]
        
        if args.verbose:
            cmd.append("--verbose")
        
        # Run the coordinator script
        try:
            subprocess.run(cmd, check=True)
            print("All plots generated successfully!")
        except subprocess.CalledProcessError as e:
            print(f"Error running plotting scripts: {e}")
            print("Falling back to basic plotting...")
            generate_basic_plots(analysis_dir, plots_dir)
    else:
        print("Comprehensive plotting system not found, using basic plotting...")
        generate_basic_plots(analysis_dir, plots_dir)

def generate_basic_plots(analysis_dir, plots_dir):
    """Generate basic plots using the original implementation."""
    print("Creating basic plots...")
    
    # Function to read xvg files
    def read_xvg(filename):
        x = []
        y = []
        with open(filename, "r") as f:
            for line in f:
                if line.startswith("#") or line.startswith("@"):
                    continue
                values = line.strip().split()
                if len(values) >= 2:
                    x.append(float(values[0]))
                    y.append(float(values[1]))
        return np.array(x), np.array(y)
    
    # Plot RDF O-O
    print("Creating O-O RDF plot...")
    try:
        x, y = read_xvg(os.path.join(analysis_dir, "rdf_OO.xvg"))
        plt.figure(figsize=(10, 6))
        plt.plot(x, y, "b-", linewidth=2)
        plt.xlabel("r (nm)", fontsize=14)
        plt.ylabel("g(r)", fontsize=14)
        plt.title("Oxygen-Oxygen Radial Distribution Function", fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(plots_dir, "rdf_oo_plot.png"), dpi=300, bbox_inches="tight")
        plt.close()
    except Exception as e:
        print(f"Error plotting O-O RDF: {e}")
    
    # Plot Temperature
    print("Creating temperature plot...")
    try:
        x, y = read_xvg(os.path.join(analysis_dir, "temperature.xvg"))
        plt.figure(figsize=(10, 6))
        plt.plot(x, y, "r-", linewidth=2)
        plt.xlabel("Time (ps)", fontsize=14)
        plt.ylabel("Temperature (K)", fontsize=14)
        plt.title("Temperature vs Time", fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(plots_dir, "temperature_plot.png"), dpi=300, bbox_inches="tight")
        plt.close()
    except Exception as e:
        print(f"Error plotting temperature: {e}")
    
    # Plot RMSD
    print("Creating RMSD plot...")
    try:
        x, y = read_xvg(os.path.join(analysis_dir, "rmsd.xvg"))
        plt.figure(figsize=(10, 6))
        plt.plot(x, y, "b-", linewidth=2)
        plt.xlabel("Time (ns)", fontsize=14)
        plt.ylabel("RMSD (nm)", fontsize=14)
        plt.title("Root Mean Square Deviation vs Time", fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(plots_dir, "rmsd_plot.png"), dpi=300, bbox_inches="tight")
        plt.close()
    except Exception as e:
        print(f"Error plotting RMSD: {e}")
    
    print("Basic plots created successfully!")

if __name__ == "__main__":
    main()
