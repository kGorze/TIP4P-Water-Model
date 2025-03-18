#!/usr/bin/env python3
"""
TIP4P Water Analysis - Main Script

This script runs all analysis scripts for TIP4P water simulation.
It creates a comprehensive summary of all analyses.

Usage:
    python run_all_analysis.py [tpr_file] [trajectory_file] [edr_file]

Default:
    Uses md.tpr, md.xtc, and md.edr in the current directory
"""

import sys
import os
import subprocess
import time
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

def check_package_versions():
    """Check versions of required packages and print warnings if needed"""
    print("Checking package compatibility...")
    
    warnings = []
    
    # Check NumPy version
    np_version = np.__version__
    print(f"NumPy version: {np_version}")
    
    # Check if NumPy version is compatible with SciPy
    major, minor, *rest = np_version.split('.')
    np_major, np_minor = int(major), int(minor)
    
    if np_major == 1 and np_minor >= 25:
        warnings.append(
            "NumPy version is too new for some SciPy versions. If you get SciPy errors, try:\n"
            "  pip install numpy<1.25.0"
        )
    
    # Check Matplotlib version
    try:
        mpl_version = plt.matplotlib.__version__
    except:
        mpl_version = plt.matplotlib._get_version()
    print(f"Matplotlib version: {mpl_version}")
    
    # Check MDAnalysis version
    try:
        import MDAnalysis as mda
        mda_version = mda.__version__
        print(f"MDAnalysis version: {mda_version}")
        
        # Check if specific modules exist
        try:
            from MDAnalysis.analysis import diffusion
            print("MDAnalysis diffusion module is available")
        except ImportError:
            try:
                from MDAnalysis.analysis import msd
                print("MDAnalysis msd module is available")
                warnings.append(
                    "Using modern MDAnalysis with msd module instead of deprecated diffusion module."
                )
            except ImportError:
                warnings.append(
                    "Neither diffusion nor msd module is available in your MDAnalysis installation.\n"
                    "  Try: pip install MDAnalysis==2.0.0 (for msd) or MDAnalysis==0.20.0 (for diffusion)"
                )
    except ImportError:
        warnings.append("MDAnalysis is not installed. Try: pip install MDAnalysis")
    
    # Print warnings
    if warnings:
        print("\nWARNINGS:")
        for i, warning in enumerate(warnings, 1):
            print(f"{i}. {warning}")
        print("\nYou may want to create a consistent environment with these versions:")
        print("  pip install numpy==1.24.3 scipy==1.10.1 matplotlib==3.7.2 MDAnalysis==2.0.0 pandas==2.0.3 seaborn==0.12.2")
        
        response = input("\nContinue anyway? (y/n): ")
        if response.lower() != 'y':
            print("Exiting. Please fix the package compatibility issues and try again.")
            sys.exit(1)
    else:
        print("Package versions appear compatible!")
    
    return warnings

def check_files_exist(tpr_file, trajectory_file, edr_file):
    """Check if required files exist in current directory or parent directory"""
    missing_files = []
    
    # Get parent directory path
    parent_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
    
    # Check tpr file
    if not os.path.exists(tpr_file):
        # Try in parent directory
        parent_tpr = os.path.join(parent_dir, os.path.basename(tpr_file))
        if os.path.exists(parent_tpr):
            print(f"File {tpr_file} found in parent directory")
            tpr_file = parent_tpr
        else:
            missing_files.append(tpr_file)
    
    # Check trajectory file
    if not os.path.exists(trajectory_file):
        # Try in parent directory
        parent_trajectory = os.path.join(parent_dir, os.path.basename(trajectory_file))
        if os.path.exists(parent_trajectory):
            print(f"File {trajectory_file} found in parent directory")
            trajectory_file = parent_trajectory
        else:
            missing_files.append(trajectory_file)
    
    # Check edr file
    if not os.path.exists(edr_file):
        # Try in parent directory
        parent_edr = os.path.join(parent_dir, os.path.basename(edr_file))
        if os.path.exists(parent_edr):
            print(f"File {edr_file} found in parent directory")
            edr_file = parent_edr
        else:
            missing_files.append(edr_file)
    
    if missing_files:
        print("Error: The following required files are missing:")
        for file in missing_files:
            print(f"  - {file}")
        return False, None, None, None
    
    return True, tpr_file, trajectory_file, edr_file

def ensure_dirs_exist():
    """Ensure that analysis directory exists"""
    if not os.path.exists('../analysis'):
        os.makedirs('../analysis')
        print("Created analysis directory")
    
    return True

def run_script(script_name, *args):
    """Run a Python script with arguments"""
    script_path = os.path.join(os.path.dirname(__file__), script_name)
    
    cmd = [sys.executable, script_path]
    cmd.extend(args)
    
    print(f"Running {script_name}...")
    start_time = time.time()
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"Output: {result.stdout}")
        if result.stderr:
            print(f"Errors: {result.stderr}")
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}: {e}")
        print(f"Stderr: {e.stderr}")
        return False
    
    end_time = time.time()
    print(f"Completed {script_name} in {end_time - start_time:.2f} seconds")
    
    return True

def create_summary_page():
    """Create a summary PDF with key plots from all analyses"""
    print("Creating summary page...")
    
    # List of plots to include in the summary
    plot_files = [
        # RDF plots
        '../analysis/combined_rdf_plot.png',
        
        # Density plots
        '../analysis/density_histogram.png',
        '../analysis/density_profile_z_plot.png',
        
        # Diffusion plots
        '../analysis/diffusion_analysis_plot.png',
        
        # H-bond plots
        '../analysis/combined_hbond_plot.png',
        '../analysis/hblife_plot.png',
        
        # RMSD plots
        '../analysis/rmsd_equilibrated_plot.png',
        
        # Thermodynamic plots
        '../analysis/thermodynamic_properties_plot.png',
        '../analysis/energy_components_plot.png'
    ]
    
    # Check which plots exist
    existing_plots = [f for f in plot_files if os.path.exists(f)]
    
    if not existing_plots:
        print("No plots found for summary page")
        return False
    
    # Create a multi-page PDF with all plots
    with PdfPages('../analysis/tip4p_water_analysis_summary.pdf') as pdf:
        # Title page
        plt.figure(figsize=(8.5, 11))
        plt.text(0.5, 0.5, 'TIP4P Water Analysis Summary', 
                 horizontalalignment='center', verticalalignment='center',
                 fontsize=24, transform=plt.gca().transAxes)
        plt.text(0.5, 0.45, f'Generated on {time.strftime("%Y-%m-%d %H:%M:%S")}',
                 horizontalalignment='center', verticalalignment='center',
                 fontsize=14, transform=plt.gca().transAxes)
        plt.axis('off')
        pdf.savefig()
        plt.close()
        
        # Add each plot to the PDF
        for plot_file in existing_plots:
            try:
                fig = plt.figure(figsize=(8.5, 11))
                img = plt.imread(plot_file)
                plt.imshow(img)
                plt.axis('off')
                plt.title(os.path.basename(plot_file), fontsize=12)
                pdf.savefig(fig)
                plt.close()
            except Exception as e:
                print(f"Error adding {plot_file} to summary: {e}")
    
    # Also create a single summary image
    if existing_plots:
        try:
            n_plots = len(existing_plots)
            n_cols = 2
            n_rows = (n_plots + n_cols - 1) // n_cols
            
            fig = plt.figure(figsize=(12, 6 * n_rows))
            gs = gridspec.GridSpec(n_rows, n_cols)
            
            for i, plot_file in enumerate(existing_plots):
                ax = plt.subplot(gs[i // n_cols, i % n_cols])
                img = plt.imread(plot_file)
                ax.imshow(img)
                ax.axis('off')
                ax.set_title(os.path.basename(plot_file), fontsize=10)
            
            plt.tight_layout()
            plt.savefig('../analysis/tip4p_water_analysis_summary.png', dpi=300)
            plt.close()
            
            print("Created summary image")
        except Exception as e:
            print(f"Error creating summary image: {e}")
    
    return True

def main():
    # Get command line arguments if provided
    if len(sys.argv) > 3:
        tpr_file = sys.argv[1]
        trajectory_file = sys.argv[2]
        edr_file = sys.argv[3]
    else:
        tpr_file = 'md.tpr'
        trajectory_file = 'md.xtc'
        edr_file = 'md.edr'
    
    # Check for package compatibility
    check_package_versions()
    
    # Check if required files exist
    files_exist, tpr_path, trajectory_path, edr_path = check_files_exist(tpr_file, trajectory_file, edr_file)
    if not files_exist:
        return
    
    # Ensure analysis directory exists
    if not ensure_dirs_exist():
        return
    
    # Run all analysis scripts
    print("\n" + "="*50)
    print("Starting TIP4P Water Analysis")
    print("="*50 + "\n")
    
    overall_start_time = time.time()
    
    # 1. Radial Distribution Functions
    print("\n" + "-"*50)
    print("RDF Analysis")
    print("-"*50)
    run_script('rdf_analysis.py', tpr_path, trajectory_path)
    
    # 2. Density Analysis
    print("\n" + "-"*50)
    print("Density Analysis")
    print("-"*50)
    run_script('density_analysis.py', tpr_path, trajectory_path)
    
    # 3. Diffusion Analysis
    print("\n" + "-"*50)
    print("Diffusion Analysis")
    print("-"*50)
    run_script('diffusion_analysis.py', tpr_path, trajectory_path)
    
    # 4. Hydrogen Bond Analysis
    print("\n" + "-"*50)
    print("Hydrogen Bond Analysis")
    print("-"*50)
    run_script('hbond_analysis.py', tpr_path, trajectory_path)
    
    # 5. RMSD Analysis
    print("\n" + "-"*50)
    print("RMSD Analysis")
    print("-"*50)
    run_script('rmsd_analysis.py', tpr_path, trajectory_path)
    
    # 6. Thermodynamic Analysis
    print("\n" + "-"*50)
    print("Thermodynamic Analysis")
    print("-"*50)
    run_script('thermodynamics_analysis.py', edr_path)
    
    # Create summary
    print("\n" + "-"*50)
    print("Creating Analysis Summary")
    print("-"*50)
    create_summary_page()
    
    overall_end_time = time.time()
    elapsed_time = overall_end_time - overall_start_time
    
    print("\n" + "="*50)
    print(f"Analysis complete in {elapsed_time:.2f} seconds")
    print("="*50 + "\n")
    
    print(f"Analysis results saved to the '../analysis' directory")
    print(f"Summary: ../analysis/tip4p_water_analysis_summary.pdf")

if __name__ == '__main__':
    main() 