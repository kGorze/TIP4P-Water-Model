#!/usr/bin/env python3

import os
import subprocess
import argparse
from pathlib import Path

def run_command(cmd, cwd=None, input_text=None):
    """Run a shell command and handle errors"""
    try:
        if input_text:
            result = subprocess.run(cmd, check=True, shell=True, cwd=cwd, 
                                   input=input_text.encode(), 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
        else:
            result = subprocess.run(cmd, check=True, shell=True, cwd=cwd,
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
        return True, result.stdout.decode()
    except subprocess.CalledProcessError as e:
        return False, f"Error running command: {cmd}\nError: {e}\nOutput: {e.stdout.decode()}\nError: {e.stderr.decode()}"

def create_md_mdp(output_dir, config_file, temperature):
    """Create production MD configuration file"""
    print(f"Creating production MD configuration for {temperature}K...")
    
    # Check if config file exists
    if not os.path.exists(config_file):
        return False, f"Error: MD config file {config_file} not found."
    
    try:
        # Read the config file
        with open(config_file, 'r') as f:
            content = f.read()
        
        # Replace temperature values
        content = content.replace("ref_t                   = 273", f"ref_t                   = {temperature}")
        
        # Modify time step for stability
        content = content.replace("dt                       = 0.002", "dt                       = 0.001")
        
        # Write the modified configuration
        with open(output_dir / "md.mdp", "w") as f:
            f.write(content)
        
        return True, f"Created production MD configuration for {temperature}K"
    
    except Exception as e:
        return False, f"Error creating MD configuration: {e}"

def run_production_md(output_dir):
    """Run production MD"""
    print("Running production MD...")
    
    # Check if nvt.gro and md.mdp exist
    if not (output_dir / "nvt.gro").exists():
        return False, "Error: nvt.gro not found. Run nvt_equilibration.py first."
    
    if not (output_dir / "nvt.cpt").exists():
        return False, "Error: nvt.cpt not found. Run nvt_equilibration.py first."
    
    if not (output_dir / "md.mdp").exists():
        return False, "Error: md.mdp not found. Run create_md_mdp first."
    
    # Prepare for production MD
    success, output = run_command(
        "gmx grompp -f md.mdp -c nvt.gro -t nvt.cpt -p topol.top -o md.tpr", 
        cwd=output_dir
    )
    
    if not success:
        print(output)
        return False, "Failed to prepare for production MD"
    
    # Run production MD
    success, output = run_command(
        "gmx mdrun -deffnm md", 
        cwd=output_dir
    )
    
    if not success:
        print(output)
        return False, "Failed to run production MD"
    
    # Check if md.gro exists
    if not (output_dir / "md.gro").exists():
        return False, "Error: md.gro not found. Production MD failed."
    
    return True, "Successfully completed production MD"

def main():
    parser = argparse.ArgumentParser(description="Run production MD on water box")
    parser.add_argument("--model", type=str, choices=["tip4p", "tip4p-ice"], required=True,
                        help="Water model to use")
    parser.add_argument("--temp", type=int, choices=[150, 200, 273, 298], required=True,
                        help="Temperature for simulation")
    parser.add_argument("--config", type=str, default=None,
                        help="Path to MD config file (default: configs/md.mdp)")
    args = parser.parse_args()
    
    base_dir = Path(__file__).parent.parent
    output_dir = base_dir / "data" / args.model / f"{args.temp}K"
    
    # Use default config file if not specified
    if args.config is None:
        config_file = base_dir / "configs" / "md.mdp"
    else:
        config_file = Path(args.config)
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Create MD configuration
    success, message = create_md_mdp(output_dir, config_file, args.temp)
    if not success:
        print(f"Error: {message}")
        return 1
    
    # Run production MD
    success, message = run_production_md(output_dir)
    if not success:
        print(f"Error: {message}")
        return 1
    
    print(f"Success: {message}")
    return 0

if __name__ == "__main__":
    exit(main()) 