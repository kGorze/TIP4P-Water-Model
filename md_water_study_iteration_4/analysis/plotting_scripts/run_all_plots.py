#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Run all plotting scripts for the TIP4P water model analysis.
This script serves as a central coordinator for all the plotting scripts.
"""

import os
import sys
import subprocess
import time
import argparse

def run_script(script_name, analysis_dir, plots_dir, verbose=True):
    """Run a plotting script with the given arguments."""
    if verbose:
        print(f"Running {script_name}...")
    
    start_time = time.time()
    
    # Construct the command
    cmd = [sys.executable, script_name, analysis_dir, plots_dir]
    
    # Run the script
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if verbose:
            print(f"  {script_name} completed in {time.time() - start_time:.2f} seconds")
            if result.stdout:
                print(f"  Output: {result.stdout.strip()}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}: {e}")
        if e.stdout:
            print(f"  Output: {e.stdout.strip()}")
        if e.stderr:
            print(f"  Error: {e.stderr.strip()}")
        return False

def main():
    """Main function to run all plotting scripts."""
    parser = argparse.ArgumentParser(description="Run all plotting scripts for TIP4P water model analysis")
    parser.add_argument("--analysis-dir", default=None, help="Directory containing analysis files")
    parser.add_argument("--plots-dir", default=None, help="Directory to save plots")
    parser.add_argument("--data-dir", default=None, help="Directory containing data files")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output")
    args = parser.parse_args()
    
    # Get the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # If analysis_dir is not provided, use the parent directory of the script directory
    if args.analysis_dir is None:
        analysis_dir = os.path.dirname(script_dir)
    else:
        analysis_dir = args.analysis_dir
    
    # If data_dir is not provided, use the data directory inside the analysis directory
    if args.data_dir is None:
        data_dir = os.path.join(analysis_dir, "data")
    else:
        data_dir = args.data_dir
    
    # If plots_dir is not provided, create a plots directory in the analysis directory
    if args.plots_dir is None:
        plots_dir = os.path.join(analysis_dir, "plots")
    else:
        plots_dir = args.plots_dir
    
    # Create the plots directory if it doesn't exist
    os.makedirs(plots_dir, exist_ok=True)
    
    print(f"Analysis directory: {analysis_dir}")
    print(f"Data directory: {data_dir}")
    print(f"Plots directory: {plots_dir}")
    
    # Check if the data directory exists
    if not os.path.exists(data_dir):
        print(f"Error: Data directory not found at {data_dir}")
        print("Please make sure the data directory exists and contains the necessary data files.")
        return False
    
    # Create symbolic links to data files in the analysis directory for backward compatibility
    print("Creating symbolic links to data files for backward compatibility...")
    for file in os.listdir(data_dir):
        if file.endswith(('.xvg', '.xpm')) and not file.startswith('#'):
            source = os.path.join(data_dir, file)
            target = os.path.join(analysis_dir, file)
            if not os.path.exists(target):
                try:
                    os.symlink(source, target)
                    if args.verbose:
                        print(f"  Created symlink: {target} -> {source}")
                except Exception as e:
                    print(f"  Warning: Could not create symlink for {file}: {e}")
    
    # List of plotting scripts to run
    plotting_scripts = [
        "plot_rdf.py",
        "plot_density.py",
        "plot_msd.py",
        "plot_hbond.py",
        "plot_temperature.py",
        "plot_energy_analysis.py",
        "plot_rmsd.py",
    ]
    
    # Check if VACF data exists
    if os.path.exists(os.path.join(data_dir, "vacf.xvg")):
        plotting_scripts.append("plot_vacf.py")
    else:
        print("VACF data not found, skipping VACF plots")
    
    # Run each plotting script
    success_count = 0
    for script in plotting_scripts:
        script_path = os.path.join(script_dir, script)
        if os.path.exists(script_path):
            if run_script(script_path, analysis_dir, plots_dir, args.verbose):
                success_count += 1
        else:
            print(f"Warning: Script {script} not found at {script_path}")
    
    print(f"Successfully ran {success_count} out of {len(plotting_scripts)} plotting scripts")
    
    # Generate summary report
    summary_script = os.path.join(script_dir, "generate_summary_report.py")
    if os.path.exists(summary_script):
        print("Generating summary report...")
        if run_script(summary_script, analysis_dir, plots_dir, args.verbose):
            print("Summary report generated successfully")
            print(f"Summary report saved to {os.path.join(plots_dir, 'tip4p_water_analysis_summary.png')}")
            print(f"Text summary saved to {os.path.join(plots_dir, 'tip4p_water_analysis_summary.txt')}")
        else:
            print("Failed to generate summary report")
    else:
        print(f"Warning: Summary report script not found at {summary_script}")
    
    print("All plotting tasks completed")
    return True

if __name__ == "__main__":
    main() 