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

def create_em_mdp(output_dir, config_file=None):
    """Create energy minimization configuration file"""
    print("Creating energy minimization configuration...")
    
    # Default energy minimization parameters
    em_content = """; Energy minimization with very gentle parameters
integrator               = steep     ; Steepest descent energy minimization
emtol                    = 1000.0    ; Stop when max force < 1000.0 kJ/mol/nm
emstep                   = 0.0001    ; Very small initial step-size for stability
nsteps                   = 50000     ; Maximum number of steps
nstlist                  = 1         ; Frequency to update neighbor list
cutoff-scheme           = Verlet    ; Neighbor search method
ns_type                 = grid      ; Method to determine neighbor list
coulombtype             = PME       ; Treatment of long range electrostatic interactions
rcoulomb                = 1.0       ; Short-range electrostatic cut-off
rvdw                    = 1.0       ; Short-range van der Waals cut-off
pbc                     = xyz       ; Periodic Boundary Condition
constraints             = none      ; No constraints during minimization
"""
    
    # If config file is provided, use it instead
    if config_file and os.path.exists(config_file):
        try:
            with open(config_file, 'r') as f:
                em_content = f.read()
            print(f"Using energy minimization parameters from {config_file}")
        except Exception as e:
            print(f"Error reading config file: {e}")
            print("Using default parameters instead.")
    
    # Write the configuration file
    with open(output_dir / "em.mdp", "w") as f:
        f.write(em_content)
    
    return True, "Created energy minimization configuration"

def convert_to_gromacs(output_dir, model):
    """Convert PDB to GROMACS format"""
    print("Converting to GROMACS format...")
    
    # Check if water_box.pdb exists
    if not (output_dir / "water_box.pdb").exists():
        return False, "Error: water_box.pdb not found. Run generate_water_box.py first."
    
    # Convert PDB to GROMACS format
    success, output = run_command(
        f"gmx pdb2gmx -f water_box.pdb -o processed.gro -water {model}", 
        cwd=output_dir, 
        input_text="15\n"  # Select OPLS-AA/L force field
    )
    
    if not success:
        print(output)
        return False, "Failed to convert PDB to GROMACS format"
    
    return True, "Successfully converted PDB to GROMACS format"

def run_energy_minimization(output_dir):
    """Run energy minimization"""
    print("Running energy minimization...")
    
    # Check if processed.gro and em.mdp exist
    if not (output_dir / "processed.gro").exists():
        return False, "Error: processed.gro not found. Run convert_to_gromacs first."
    
    if not (output_dir / "em.mdp").exists():
        return False, "Error: em.mdp not found. Run create_em_mdp first."
    
    # Prepare for energy minimization
    success, output = run_command(
        "gmx grompp -f em.mdp -c processed.gro -p topol.top -o em.tpr", 
        cwd=output_dir
    )
    
    if not success:
        print(output)
        return False, "Failed to prepare for energy minimization"
    
    # Run energy minimization
    success, output = run_command(
        "gmx mdrun -v -deffnm em", 
        cwd=output_dir
    )
    
    if not success:
        print(output)
        return False, "Failed to run energy minimization"
    
    # Check if em.gro exists
    if not (output_dir / "em.gro").exists():
        return False, "Error: em.gro not found. Energy minimization failed."
    
    return True, "Successfully completed energy minimization"

def main():
    parser = argparse.ArgumentParser(description="Run energy minimization on water box")
    parser.add_argument("--model", type=str, choices=["tip4p", "tip4p-ice"], required=True,
                        help="Water model to use")
    parser.add_argument("--temp", type=int, choices=[150, 200, 273, 298], required=True,
                        help="Temperature for simulation")
    parser.add_argument("--config", type=str, default=None,
                        help="Path to energy minimization config file (default: configs/em.mdp)")
    args = parser.parse_args()
    
    base_dir = Path(__file__).parent.parent
    output_dir = base_dir / "data" / args.model / f"{args.temp}K"
    
    # Use default config file if not specified
    if args.config is None:
        config_file = base_dir / "configs" / "em.mdp"
    else:
        config_file = Path(args.config)
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Create energy minimization configuration
    success, message = create_em_mdp(output_dir, config_file)
    if not success:
        print(f"Error: {message}")
        return 1
    
    # Convert PDB to GROMACS format
    success, message = convert_to_gromacs(output_dir, args.model)
    if not success:
        print(f"Error: {message}")
        return 1
    
    # Run energy minimization
    success, message = run_energy_minimization(output_dir)
    if not success:
        print(f"Error: {message}")
        return 1
    
    print(f"Success: {message}")
    return 0

if __name__ == "__main__":
    exit(main()) 