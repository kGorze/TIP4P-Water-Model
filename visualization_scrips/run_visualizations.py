#!/usr/bin/env python3

import os
import subprocess
import sys

def run_script(script_name):
    """Run a Python script with proper environment settings."""
    print(f"Running {script_name}...")
    
    # Create the command with unset PYTHONHOME and PYTHONPATH
    cmd = f"unset PYTHONHOME && unset PYTHONPATH && python3 {script_name}"
    
    # Run the command
    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"Successfully ran {script_name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}: {e}")
        return False

def main():
    """Run all visualization scripts."""
    # Create output directory
    output_dir = "visualization_outputs"
    os.makedirs(output_dir, exist_ok=True)
    
    # List of scripts to run
    scripts = [
        "water_molecule_constraint.py",
        "minimal_nonbonded_cutoffs.py",  # Using the minimal version that works
        "timestep_visualization.py",
        "sampling_methods_visualization.py",  # Added new sampling methods visualization
        "energy_minimization_visualization.py",
        "nvt_equilibration_visualization.py",
        "npt_equilibration_visualization.py",
        "regular_md_visualization.py"  # Add the new regular MD visualization script
    ]
    
    # Run each script
    success_count = 0
    for script in scripts:
        if run_script(script):
            success_count += 1
    
    # Print summary
    print("\nVisualization Summary:")
    print(f"Successfully ran {success_count} out of {len(scripts)} scripts")
    
    if success_count == len(scripts):
        print(f"\nAll visualizations have been saved to the '{output_dir}' directory.")
        print("The following files were created:")
        
        # List files in the output directory
        files = os.listdir(output_dir)
        for file in sorted(files):
            file_path = os.path.join(output_dir, file)
            file_size = os.path.getsize(file_path) / 1024  # Size in KB
            print(f"  - {file} ({file_size:.1f} KB)")
    else:
        print("\nSome visualizations failed. Please check the error messages above.")

if __name__ == "__main__":
    main() 