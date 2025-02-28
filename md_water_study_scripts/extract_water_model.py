#!/usr/bin/env python3

import os
import subprocess
import argparse
import shutil
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

def extract_water_model(output_dir, model):
    """Extract the water model directly from GROMACS"""
    print(f"Extracting {model} water model from GROMACS...")
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Create a temporary directory
    temp_dir = Path(output_dir) / "temp_water"
    temp_dir.mkdir(exist_ok=True)
    
    # Create a small box with a single water molecule
    success, output = run_command(f"gmx solvate -cs {model} -box 0.5 0.5 0.5 -o water_box.gro", 
                                 cwd=temp_dir, input_text="1\n")
    if not success:
        print(output)
        return False, "Failed to create water box"
    
    # Convert to PDB format
    success, output = run_command("gmx editconf -f water_box.gro -o water_box.pdb", cwd=temp_dir)
    if not success:
        print(output)
        return False, "Failed to convert water box to PDB format"
    
    # Extract the first water molecule
    try:
        with open(temp_dir / "water_box.pdb", "r") as f:
            lines = f.readlines()
        
        # Find all atoms of the first water molecule
        water_lines = []
        current_res_id = None
        
        for line in lines:
            if line.startswith("ATOM"):
                res_id = line[22:26].strip()
                if current_res_id is None:
                    current_res_id = res_id
                    water_lines.append(line)
                elif res_id == current_res_id:
                    water_lines.append(line)
                else:
                    break
        
        # Write the single water molecule to a PDB file
        with open(Path(output_dir) / "water.pdb", "w") as f:
            for line in water_lines:
                f.write(line)
            f.write("END\n")
        
        print(f"Extracted {model} water model with {len(water_lines)} atoms")
        
        # Clean up
        shutil.rmtree(temp_dir)
        
        return True, f"Successfully extracted {model} water model with {len(water_lines)} atoms"
    
    except Exception as e:
        print(f"Error extracting water model: {e}")
        return False, f"Error extracting water model: {e}"

def main():
    parser = argparse.ArgumentParser(description="Extract water model from GROMACS")
    parser.add_argument("--model", type=str, choices=["tip4p", "tip4p-ice"], required=True,
                        help="Water model to extract")
    parser.add_argument("--output", type=str, default=None,
                        help="Output directory (default: data/<model>)")
    args = parser.parse_args()
    
    base_dir = Path(__file__).parent.parent
    
    if args.output:
        output_dir = Path(args.output)
    else:
        output_dir = base_dir / "data" / args.model
    
    success, message = extract_water_model(output_dir, args.model)
    
    if success:
        print(f"Success: {message}")
        return 0
    else:
        print(f"Error: {message}")
        return 1

if __name__ == "__main__":
    exit(main()) 