#!/usr/bin/env python3

import os
import subprocess
import argparse
from pathlib import Path
import shutil

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

def calculate_box_size(num_molecules, density=0.93):
    """Calculate box size based on number of molecules and density
    density in g/cm³, returns box size in Angstroms"""
    # Constants
    water_molar_mass = 18.015  # g/mol
    avogadro = 6.022e23  # molecules/mol
    
    # Calculate volume in cm³
    volume = (num_molecules * water_molar_mass) / (density * avogadro)
    # Convert to nm³
    volume_nm = volume * 1e24
    # Calculate box size in nm
    box_size_nm = volume_nm ** (1/3)
    # Convert to Angstroms
    box_size_ang = box_size_nm * 10
    
    return box_size_ang

def read_packmol_config(config_file):
    """Read parameters from PackMol config file"""
    num_molecules = 5500  # default
    box_size = 54.0  # default
    tolerance = 2.0  # default
    
    try:
        with open(config_file, 'r') as f:
            content = f.read()
            
        for line in content.split('\n'):
            if 'number' in line and 'structure' not in line:
                try:
                    num_molecules = int(line.strip().split()[1])
                except:
                    pass
            elif 'inside box' in line:
                try:
                    parts = line.strip().split()
                    box_size = float(parts[-1])
                except:
                    pass
            elif 'tolerance' in line and 'filetype' not in line:
                try:
                    tolerance = float(line.strip().split()[1])
                except:
                    pass
        
        return num_molecules, box_size, tolerance
    except Exception as e:
        print(f"Error reading config file: {e}")
        return num_molecules, box_size, tolerance

def write_packmol_input(output_dir, config_file=None):
    """Write PackMol input file using parameters from config"""
    if not config_file or not os.path.exists(config_file):
        print("Warning: Config file not found. Using default parameters.")
        num_molecules = 5500
        box_size = calculate_box_size(num_molecules)
        tolerance = 2.0
    else:
        num_molecules, box_size, tolerance = read_packmol_config(config_file)
        print(f"Using parameters from config file: {config_file}")
    
    print(f"Writing PackMol input file with {num_molecules} molecules in a {box_size}x{box_size}x{box_size} box...")
    
    # Create the input file
    with open(output_dir / "pack.inp", "w") as f:
        f.write(f"""# PackMol input file for generating water box
# Box with {num_molecules} water molecules at density 0.93 g/cm³

tolerance {tolerance}
filetype pdb
output water_box.pdb

structure water.pdb
    number {num_molecules}
    inside box 0. 0. 0. {box_size} {box_size} {box_size}
end structure
""")
    
    return True, f"Created PackMol input file with {num_molecules} molecules in a {box_size}x{box_size}x{box_size} box"

def generate_water_box(output_dir, config_file=None):
    """Generate water box using PackMol"""
    print(f"Generating water box in {output_dir}...")
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Copy water.pdb from model directory to temperature directory
    model_dir = output_dir.parent
    if not (model_dir / "water.pdb").exists():
        return False, "Error: water.pdb not found in model directory. Run extract_water_model.py first."
    
    shutil.copy2(model_dir / "water.pdb", output_dir / "water.pdb")
    
    # Write PackMol input file
    success, message = write_packmol_input(output_dir, config_file)
    if not success:
        return False, message
    
    # Run PackMol
    print("Running PackMol...")
    success, output = run_command("julia -e 'using Packmol; run_packmol(\"pack.inp\")'", cwd=output_dir)
    if not success:
        print(output)
        return False, "Failed to run PackMol"
    
    # Check if water_box.pdb was created
    if not (output_dir / "water_box.pdb").exists():
        return False, "Error: PackMol did not generate water_box.pdb"
    
    return True, "Successfully generated water box"

def main():
    parser = argparse.ArgumentParser(description="Generate water box using PackMol")
    parser.add_argument("--model", type=str, choices=["tip4p", "tip4p-ice"], required=True,
                        help="Water model to use")
    parser.add_argument("--temp", type=int, choices=[150, 200, 273, 298], required=True,
                        help="Temperature for simulation")
    parser.add_argument("--config", type=str, default=None,
                        help="Path to PackMol config file (default: configs/water_box.inp)")
    args = parser.parse_args()
    
    base_dir = Path(__file__).parent.parent
    output_dir = base_dir / "data" / args.model / f"{args.temp}K"
    
    # Use default config file if not specified
    if args.config is None:
        config_file = base_dir / "configs" / "water_box.inp"
    else:
        config_file = Path(args.config)
    
    success, message = generate_water_box(output_dir, config_file)
    
    if success:
        print(f"Success: {message}")
        return 0
    else:
        print(f"Error: {message}")
        return 1

if __name__ == "__main__":
    exit(main()) 