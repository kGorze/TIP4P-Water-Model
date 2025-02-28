#!/usr/bin/env python3

import os
import subprocess
import argparse
from pathlib import Path
import time

def run_command(cmd, cwd=None):
    """Run a shell command and handle errors"""
    print(f"Running: {cmd}")
    try:
        subprocess.run(cmd, check=True, shell=True, cwd=cwd)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {cmd}")
        print(f"Error: {e}")
        return False

def run_workflow(model="tip4p", temp=273, iteration_dir=None):
    """Run the complete TIP4P water simulation workflow"""
    
    start_time = time.time()
    print(f"Starting TIP4P water simulation workflow at {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Model: {model}, Temperature: {temp}K")
    
    # Get the base directory
    if iteration_dir:
        base_dir = Path(iteration_dir).resolve()
    else:
        # Default to the parent directory of this script
        base_dir = Path(__file__).resolve().parent.parent
    
    # Get the scripts directory (always use the central scripts directory)
    scripts_dir = Path(__file__).resolve().parent
    
    print(f"Working with iteration directory: {base_dir}")
    print(f"Using scripts from: {scripts_dir}")
    
    # Step 1: Run the simulation (includes energy minimization, NVT, NPT, and production MD)
    print("\n===== STEP 1: RUNNING SIMULATION =====")
    if not run_command(f"python3 {scripts_dir}/run_simulation.py --model {model} --temp {temp} --iteration-dir {base_dir}", cwd=base_dir):
        print("Simulation failed. Exiting workflow.")
        return False
    
    # Step 2: Analyze the results
    print("\n===== STEP 2: ANALYZING RESULTS =====")
    if not run_command(f"python3 {scripts_dir}/analyze_results.py --model {model} --temp {temp} --iteration-dir {base_dir}", cwd=base_dir):
        print("Analysis failed. Exiting workflow.")
        return False
    
    # Step 3: Perform additional hydrogen bond analysis
    print("\n===== STEP 3: HYDROGEN BOND ANALYSIS =====")
    data_dir = base_dir / "data" / model / f"{temp}K"
    analysis_dir = base_dir / "analysis" / model / f"{temp}K"
    
    print("Analyzing hydrogen bonds...")
    if not run_command(f"echo 'a a' | gmx hbond -f {data_dir}/md.xtc -s {data_dir}/md.tpr -num {analysis_dir}/hbnum.xvg -life {analysis_dir}/hblife.xvg -dist {analysis_dir}/hbdist.xvg -ang {analysis_dir}/hbang.xvg", cwd=base_dir):
        print("Warning: Hydrogen bond analysis failed. Continuing workflow...")
    
    # Step 4: Perform velocity autocorrelation function analysis if trajectory with velocities exists
    print("\n===== STEP 4: VELOCITY AUTOCORRELATION ANALYSIS =====")
    if os.path.exists(data_dir / "md.trr"):
        print("Analyzing velocity autocorrelation function...")
        if not run_command(f"echo 'SOL' | gmx velacc -f {data_dir}/md.trr -s {data_dir}/md.tpr -o {analysis_dir}/vacf.xvg -os {analysis_dir}/vacf_spectrum.xvg -acflen 1000 -nonormalize", cwd=base_dir):
            print("Warning: VACF analysis failed. Continuing workflow...")
    else:
        print("md.trr file not found. Skipping VACF analysis.")
        print("Note: For VACF analysis, ensure velocities are saved in trajectory (nstvout in md.mdp).")
    
    # Calculate total runtime
    end_time = time.time()
    runtime = end_time - start_time
    hours, remainder = divmod(runtime, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    print("\n===== WORKFLOW COMPLETED SUCCESSFULLY =====")
    print(f"Total runtime: {int(hours)}h {int(minutes)}m {int(seconds)}s")
    print(f"Model: {model}, Temperature: {temp}K")
    
    # Print location of results
    output_dir = base_dir / "data" / model / f"{temp}K"
    analysis_dir = base_dir / "analysis" / model / f"{temp}K"
    
    print(f"\nSimulation results located at: {output_dir}")
    print(f"Analysis results located at: {analysis_dir}")
    print("\nTo visualize the results, you can use the visualization scripts in the scripts directory.")
    
    return True

def main():
    parser = argparse.ArgumentParser(description="Run complete TIP4P water simulation workflow")
    parser.add_argument("--model", type=str, default="tip4p", help="Water model (default: tip4p)")
    parser.add_argument("--temp", type=int, default=273, help="Temperature in K (default: 273)")
    parser.add_argument("--iteration-dir", type=str, help="Path to the iteration directory (default: parent directory of this script)")
    args = parser.parse_args()
    
    # Run the workflow
    success = run_workflow(args.model, args.temp, args.iteration_dir)
    
    if not success:
        print("\nWorkflow completed with errors. Please check the logs.")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main()) 