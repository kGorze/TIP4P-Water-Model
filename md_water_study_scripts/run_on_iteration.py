#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from pathlib import Path

def list_available_scripts():
    """List all available Python scripts in the current directory."""
    script_dir = Path(__file__).resolve().parent
    scripts = [f for f in os.listdir(script_dir) if f.endswith('.py') and f != os.path.basename(__file__)]
    return scripts

def list_available_iterations():
    """List all available iterations in the project root."""
    current_dir = Path(__file__).resolve().parent
    project_root = current_dir.parent
    iterations = [d for d in os.listdir(project_root) if d.startswith('md_water_study_iteration_')]
    return iterations

def run_script(script_name, iteration_name, args=None):
    """Run a script on a specific iteration with optional arguments."""
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    
    # Validate script exists
    script_path = script_dir / script_name
    if not script_path.exists():
        print(f"Error: Script '{script_name}' not found in {script_dir}")
        return False
    
    # Validate iteration exists
    iteration_path = project_root / iteration_name
    if not iteration_path.exists():
        print(f"Error: Iteration directory '{iteration_name}' not found in {project_root}")
        return False
    
    # Construct command
    cmd = [sys.executable, str(script_path)]
    
    # Add the iteration directory argument
    cmd.extend(["--iteration-dir", str(iteration_path)])
    
    # Add any additional arguments
    if args:
        cmd.extend(args)
    
    # Print the command
    print(f"Running: {' '.join(cmd)}")
    
    # Execute the command
    try:
        subprocess.run(cmd, check=True)
        print(f"Successfully ran {script_name} on {iteration_name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name} on {iteration_name}: {e}")
        return False

def main():
    # Get available scripts and iterations
    available_scripts = list_available_scripts()
    available_iterations = list_available_iterations()
    
    # Parse arguments
    parser = argparse.ArgumentParser(description="Run a script on a specific iteration")
    parser.add_argument("script", choices=available_scripts + ["list"], 
                        help="Script to run or 'list' to see available scripts")
    parser.add_argument("iteration", choices=available_iterations + ["list", "all"], 
                        help="Iteration to run on, 'list' to see available iterations, or 'all' to run on all iterations")
    parser.add_argument("script_args", nargs=argparse.REMAINDER, 
                        help="Additional arguments to pass to the script")
    
    args = parser.parse_args()
    
    # Handle list commands
    if args.script == "list":
        print("Available scripts:")
        for script in sorted(available_scripts):
            print(f"  {script}")
        return 0
    
    if args.iteration == "list":
        print("Available iterations:")
        for iteration in sorted(available_iterations):
            print(f"  {iteration}")
        return 0
    
    # Run on all iterations
    if args.iteration == "all":
        success = True
        for iteration in sorted(available_iterations):
            print(f"\n=== Running {args.script} on {iteration} ===\n")
            if not run_script(args.script, iteration, args.script_args):
                success = False
        return 0 if success else 1
    
    # Run on a specific iteration
    return 0 if run_script(args.script, args.iteration, args.script_args) else 1

if __name__ == "__main__":
    sys.exit(main()) 