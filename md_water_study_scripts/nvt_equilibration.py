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

def create_nvt_mdp(output_dir, config_file, temperature):
    """Create NVT equilibration configuration file"""
    print(f"Creating NVT equilibration configuration for {temperature}K...")
    
    # Check if config file exists
    if not os.path.exists(config_file):
        return False, f"Error: NVT config file {config_file} not found."
    
    try:
        # Read the config file
        with open(config_file, 'r') as f:
            content = f.read()
        
        # Replace temperature values
        content = content.replace("ref_t                   = 273", f"ref_t                   = {temperature}")
        content = content.replace("gen_temp                = 273", f"gen_temp                = {temperature}")
        
        # Modify time step for stability
        content = content.replace("dt                       = 0.002", "dt                       = 0.001")
        
        # Write the modified configuration
        with open(output_dir / "nvt.mdp", "w") as f:
            f.write(content)
        
        return True, f"Created NVT equilibration configuration for {temperature}K"
    
    except Exception as e:
        return False, f"Error creating NVT configuration: {e}"

def run_nvt_equilibration(output_dir):
    """Run NVT equilibration"""
    print("Running NVT equilibration...")
    
    # Check if em.gro and nvt.mdp exist
    if not (output_dir / "em.gro").exists():
        return False, "Error: em.gro not found. Run energy_minimization.py first."
    
    if not (output_dir / "nvt.mdp").exists():
        return False, "Error: nvt.mdp not found. Run create_nvt_mdp first."
    
    # Prepare for NVT equilibration
    success, output = run_command(
        "gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr", 
        cwd=output_dir
    )
    
    if not success:
        print(output)
        return False, "Failed to prepare for NVT equilibration"
    
    # Run NVT equilibration
    success, output = run_command(
        "gmx mdrun -deffnm nvt", 
        cwd=output_dir
    )
    
    if not success:
        print(output)
        return False, "Failed to run NVT equilibration"
    
    # Check if nvt.gro exists
    if not (output_dir / "nvt.gro").exists():
        return False, "Error: nvt.gro not found. NVT equilibration failed."
    
    return True, "Successfully completed NVT equilibration"

def main():
    parser = argparse.ArgumentParser(description="Run NVT equilibration on water box")
    parser.add_argument("--model", type=str, choices=["tip4p", "tip4p-ice"], required=True,
                        help="Water model to use")
    parser.add_argument("--temp", type=int, choices=[150, 200, 273, 298], required=True,
                        help="Temperature for simulation")
    parser.add_argument("--config", type=str, default=None,
                        help="Path to NVT config file (default: configs/nvt.mdp)")
    args = parser.parse_args()
    
    base_dir = Path(__file__).parent.parent
    output_dir = base_dir / "data" / args.model / f"{args.temp}K"
    
    # Use default config file if not specified
    if args.config is None:
        config_file = base_dir / "configs" / "nvt.mdp"
    else:
        config_file = Path(args.config)
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Create NVT configuration
    success, message = create_nvt_mdp(output_dir, config_file, args.temp)
    if not success:
        print(f"Error: {message}")
        return 1
    
    # Run NVT equilibration
    success, message = run_nvt_equilibration(output_dir)
    if not success:
        print(f"Error: {message}")
        return 1
    
    print(f"Success: {message}")
    return 0

if __name__ == "__main__":
    exit(main()) 