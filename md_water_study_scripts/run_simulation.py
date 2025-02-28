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
            subprocess.run(cmd, check=True, shell=True, cwd=cwd, input=input_text.encode())
        else:
            subprocess.run(cmd, check=True, shell=True, cwd=cwd)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {cmd}")
        print(f"Error: {e}")
        exit(1)

def extract_tip4p_from_gromacs(output_dir, model):
    """Extract the TIP4P water model directly from GROMACS"""
    # Create a temporary directory
    temp_dir = output_dir / "temp_water"
    temp_dir.mkdir(exist_ok=True)
    
    # Create a small box with a single water molecule
    print(f"Extracting {model} water model from GROMACS...")
    run_command(f"gmx solvate -cs {model} -box 0.5 0.5 0.5 -o water_box.gro", cwd=temp_dir, input_text="1\n")
    
    # Convert to PDB format
    run_command("gmx editconf -f water_box.gro -o water_box.pdb", cwd=temp_dir)
    
    # Extract the first water molecule
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
    with open(output_dir / "water.pdb", "w") as f:
        for line in water_lines:
            f.write(line)
        f.write("END\n")
    
    # Clean up
    shutil.rmtree(temp_dir)
    
    print(f"Extracted {model} water model from GROMACS with {len(water_lines)} atoms")
    return len(water_lines)

def write_packmol_input(output_dir, config_file=None):
    """Write PackMol input file with parameters from config file"""
    if config_file and os.path.exists(config_file):
        # Copy the config file to the output directory
        shutil.copy(config_file, output_dir / "pack.inp")
        print(f"Using PackMol configuration from {config_file}")
    else:
        # Use default configuration with 5500 water molecules
        with open(output_dir / "pack.inp", "w") as f:
            f.write("""# PackMol input file for generating water box
# Box with 5500 water molecules at density 0.93 g/cm³

tolerance 3.0
filetype pdb
output water_box.pdb

structure water.pdb
    number 5500
    inside box 0. 0. 0. 78. 78. 78.
end structure
""")
        print("Using default PackMol configuration with 5500 water molecules")

def create_em_mdp(output_dir):
    """Create a more thorough energy minimization configuration"""
    with open(output_dir / "em.mdp", "w") as f:
        f.write("""; Energy minimization with very gentle parameters
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
""")

def create_posre_mdp(output_dir):
    """Create a position restraint equilibration step"""
    with open(output_dir / "posre.mdp", "w") as f:
        f.write("""; Position restrained molecular dynamics
integrator               = md        ; leap-frog integrator
dt                       = 0.0005    ; Very small time step (0.5 fs)
nsteps                   = 10000     ; 5 ps
nstxout                  = 500       ; save coordinates every 0.25 ps
nstvout                  = 500       ; save velocities every 0.25 ps
nstenergy                = 500       ; save energies every 0.25 ps
nstlog                   = 500       ; update log file every 0.25 ps
continuation             = no        ; first dynamics run
constraint_algorithm     = lincs     ; holonomic constraints 
constraints              = h-bonds   ; bonds involving H are constrained
lincs_iter               = 1         ; accuracy of LINCS
lincs_order              = 4         ; also related to accuracy
; Nonbonded settings 
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale ; modified Berendsen thermostat
tc-grps                 = System    ; one coupling group
tau_t                   = 0.1       ; time constant, in ps
ref_t                   = 273       ; reference temperature, one for each group, in K
; Pressure coupling
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 273       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
""")

def setup_simulation(temp, model, output_dir, base_dir, config_file=None):
    """Set up and run a single simulation"""
    print(f"Setting up simulation for {model} at {temp}K")
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Extract water model from GROMACS
    extract_tip4p_from_gromacs(output_dir, model)
    
    # Generate initial structure using PackMol
    print("Generating initial structure with PackMol...")
    write_packmol_input(output_dir, config_file)
    run_command("julia -e 'using Packmol; run_packmol(\"pack.inp\")'", cwd=output_dir)
    
    # Convert PDB to GROMACS format
    print("Converting to GROMACS format...")
    run_command(f"gmx pdb2gmx -f water_box.pdb -o processed.gro -water {model}", cwd=output_dir, input_text="15\n")
    
    # Create a more thorough energy minimization configuration
    create_em_mdp(output_dir)
    
    # Energy minimization
    print("Running energy minimization...")
    run_command("gmx grompp -f em.mdp -c processed.gro -p topol.top -o em.tpr", cwd=output_dir)
    run_command("gmx mdrun -v -deffnm em", cwd=output_dir)
    
    # Create position restraint file
    print("Creating position restraint file...")
    create_posre_mdp(output_dir)
    
    # Position restrained equilibration
    print("Running position restrained equilibration...")
    run_command("gmx grompp -f posre.mdp -c em.gro -r em.gro -p topol.top -o posre.tpr -maxwarn 1", cwd=output_dir)
    run_command("gmx mdrun -v -deffnm posre", cwd=output_dir)
    
    # Update NVT configuration for specific temperature
    print("Preparing NVT equilibration...")
    nvt_mdp = Path(output_dir) / "nvt.mdp"
    with open(base_dir / "configs" / "nvt.mdp") as f:
        content = f.read()
    content = content.replace("ref_t                   = 273", f"ref_t                   = {temp}")
    content = content.replace("gen_temp                = 273", f"gen_temp                = {temp}")
    # Modify NVT parameters for stability
    content = content.replace("dt                       = 0.002", "dt                       = 0.001")
    content = content.replace("nsteps                   = 50000", "nsteps                   = 100000")
    with open(nvt_mdp, "w") as f:
        f.write(content)
    
    # Temperature equilibration (NVT)
    print("Running NVT equilibration...")
    run_command("gmx grompp -f nvt.mdp -c posre.gro -r posre.gro -p topol.top -o nvt.tpr", cwd=output_dir)
    run_command("gmx mdrun -deffnm nvt", cwd=output_dir)
    
    # NPT equilibration
    print("Preparing NPT equilibration...")
    npt_mdp = Path(output_dir) / "npt.mdp"
    with open(base_dir / "configs" / "npt.mdp") as f:
        content = f.read()
    content = content.replace("ref_t                   = 273", f"ref_t                   = {temp}")
    with open(npt_mdp, "w") as f:
        f.write(content)
    
    # Pressure equilibration (NPT)
    print("Running NPT equilibration...")
    run_command("gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 1", cwd=output_dir)
    run_command("gmx mdrun -deffnm npt", cwd=output_dir)
    
    # Update production MD configuration for specific temperature
    print("Preparing production run...")
    md_mdp = Path(output_dir) / "md.mdp"
    with open(base_dir / "configs" / "md.mdp") as f:
        content = f.read()
    content = content.replace("ref_t                   = 273", f"ref_t                   = {temp}")
    with open(md_mdp, "w") as f:
        f.write(content)
    
    # Production MD
    print("Running production MD...")
    run_command("gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr", cwd=output_dir)
    run_command("gmx mdrun -deffnm md", cwd=output_dir)
    
    print("Simulation complete!")
    return True

def main():
    parser = argparse.ArgumentParser(description="Run water model simulations")
    parser.add_argument("--temp", type=int, default=273, help="Temperature in K")
    parser.add_argument("--model", type=str, default="tip4p", help="Water model (tip4p)")
    parser.add_argument("--config", type=str, help="Path to PackMol config file")
    parser.add_argument("--iteration-dir", type=str, help="Path to the iteration directory")
    args = parser.parse_args()
    
    if args.iteration_dir:
        base_dir = Path(args.iteration_dir).resolve()
    else:
        base_dir = Path(__file__).resolve().parent.parent
    
    output_dir = base_dir / "data" / args.model / f"{args.temp}K"
    
    config_file = args.config
    if not config_file and os.path.exists(base_dir / "configs" / "water_box.inp"):
        config_file = base_dir / "configs" / "water_box.inp"
    
    print(f"Working with iteration directory: {base_dir}")
    setup_simulation(args.temp, args.model, output_dir, base_dir, config_file)
    
if __name__ == "__main__":
    main() 