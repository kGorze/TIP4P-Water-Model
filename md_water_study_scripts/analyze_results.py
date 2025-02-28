#!/usr/bin/env python3

import os
import subprocess
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def run_command(cmd, cwd=None, input_text=None):
    """Run a shell command and handle errors"""
    print(f"Running: {cmd}")
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
        print(f"Error running command: {cmd}")
        print(f"Error: {e}")
        return False, f"Error running command: {cmd}\nError: {e}\nOutput: {e.stdout.decode()}\nError: {e.stderr.decode()}"

def analyze_energy(output_dir):
    """Analyze energy from simulation"""
    print("Analyzing energy...")
    
    # Check if md.edr exists
    if not (output_dir / "md.edr").exists():
        return False, "Error: md.edr not found. Run production_md.py first."
    
    # Run energy analysis for multiple terms
    success, output = run_command(
        "echo '10 11 12 13 0' | gmx energy -f md.edr -o energy.xvg", 
        cwd=output_dir
    )
    
    if not success:
        print(output)
        return False, "Failed to analyze energy"
    
    # Also extract potential energy separately
    success, output = run_command(
        "echo '10 0' | gmx energy -f md.edr -o potential.xvg", 
        cwd=output_dir
    )
    
    # Extract temperature
    success, output = run_command(
        "echo '14 0' | gmx energy -f md.edr -o temperature.xvg", 
        cwd=output_dir
    )
    
    # Extract pressure
    success, output = run_command(
        "echo '15 0' | gmx energy -f md.edr -o pressure.xvg", 
        cwd=output_dir
    )
    
    # Extract density
    success, output = run_command(
        "echo '22 0' | gmx energy -f md.edr -o density.xvg", 
        cwd=output_dir
    )
    
    # Check if energy.xvg exists
    if not (output_dir / "energy.xvg").exists():
        return False, "Error: energy.xvg not found. Energy analysis failed."
    
    return True, "Successfully analyzed energy"

def analyze_rmsd(output_dir):
    """Analyze RMSD of the system"""
    print("Analyzing RMSD...")
    
    # Check if md.xtc exists, if not check for md.trr and convert
    if not (output_dir / "md.xtc").exists():
        if (output_dir / "md.trr").exists():
            print("Found md.trr, converting to md.xtc...")
            run_command(
                "echo '0' | gmx trjconv -f md.trr -s md.tpr -o md.xtc",
                cwd=output_dir
            )
        else:
            return False, "Error: Neither md.xtc nor md.trr found. Run production_md.py first."
    
    # Run RMSD analysis
    run_command(
        "echo '4 4' | gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns",
        cwd=output_dir
    )
    
    return True, "RMSD analysis completed successfully"

def analyze_rdf(output_dir):
    """Analyze radial distribution functions"""
    print("Analyzing RDF...")
    
    # Check if md.xtc exists, if not check for md.trr and convert
    if not (output_dir / "md.xtc").exists():
        if (output_dir / "md.trr").exists():
            print("Found md.trr, converting to md.xtc...")
            run_command(
                "echo '0' | gmx trjconv -f md.trr -s md.tpr -o md.xtc",
                cwd=output_dir
            )
        else:
            return False, "Error: Neither md.xtc nor md.trr found. Run production_md.py first."
    
    # O-O RDF
    run_command(
        "echo '2 2' | gmx rdf -f md.xtc -s md.tpr -o rdf_oo.xvg -ref 'name OW' -sel 'name OW'",
        cwd=output_dir
    )
    
    # O-H RDF
    run_command(
        "echo '2 3' | gmx rdf -f md.xtc -s md.tpr -o rdf_oh.xvg -ref 'name OW' -sel 'name HW'",
        cwd=output_dir
    )
    
    # H-H RDF
    run_command(
        "echo '3 3' | gmx rdf -f md.xtc -s md.tpr -o rdf_hh.xvg -ref 'name HW' -sel 'name HW'",
        cwd=output_dir
    )
    
    return True, "RDF analysis completed successfully"

def analyze_diffusion(output_dir):
    """Analyze diffusion coefficient"""
    print("Analyzing diffusion coefficient...")
    
    # Check if md.xtc exists, if not check for md.trr and convert
    if not (output_dir / "md.xtc").exists():
        if (output_dir / "md.trr").exists():
            print("Found md.trr, converting to md.xtc...")
            run_command(
                "echo '0' | gmx trjconv -f md.trr -s md.tpr -o md.xtc",
                cwd=output_dir
            )
        else:
            return False, "Error: Neither md.xtc nor md.trr found. Run production_md.py first."
    
    # Run MSD analysis
    run_command(
        "echo '2' | gmx msd -f md.xtc -s md.tpr -o msd.xvg -mol -rmcomm",
        cwd=output_dir
    )
    
    return True, "Diffusion coefficient analysis completed successfully"

def analyze_structure(output_dir):
    """Analyze structural properties"""
    print("Analyzing structural properties...")
    
    # Check if md.xtc exists, if not check for md.trr and convert
    if not (output_dir / "md.xtc").exists():
        if (output_dir / "md.trr").exists():
            print("Found md.trr, converting to md.xtc...")
            run_command(
                "echo '0' | gmx trjconv -f md.trr -s md.tpr -o md.xtc",
                cwd=output_dir
            )
        else:
            return False, "Error: Neither md.xtc nor md.trr found. Run production_md.py first."
    
    # Clean trajectory (remove PBC effects)
    run_command(
        "echo '1' | gmx trjconv -f md.xtc -s md.tpr -o md_clean.xtc -pbc nojump -center",
        cwd=output_dir
    )
    
    # Principal component analysis
    run_command(
        "echo '1' | gmx covar -f md_clean.xtc -s md.tpr -o eigenvalues.xvg -av average.pdb -xpma covmat.xpm",
        cwd=output_dir
    )
    
    return True, "Structural analysis completed successfully"

def analyze_free_energy(output_dir):
    """Estimate free energy from simulation"""
    print("Estimating free energy...")
    
    # Check if md.edr exists
    if not (output_dir / "md.edr").exists():
        return False, "Error: md.edr not found. Run production_md.py first."
    
    # Extract different energy components for free energy estimation
    energy_terms = ['Potential', 'Kinetic', 'Total-Energy', 'Temperature', 'Pressure', 'Volume']
    term_numbers = ['10', '11', '12', '14', '15', '20']
    
    # Get all energy terms
    success, output = run_command(
        f"echo '{' '.join(term_numbers)} 0' | gmx energy -f md.edr -o energy_terms.xvg", 
        cwd=output_dir
    )
    
    if not success:
        print(output)
        return False, "Failed to extract energy terms for free energy estimation"
    
    # Check if file exists
    if not (output_dir / "energy_terms.xvg").exists():
        return False, "Error: energy_terms.xvg not found. Free energy analysis failed."
    
    # Create a simple file with average values
    try:
        # Extract averages from GMX energy output
        with open(output_dir / "free_energy.txt", "w") as f:
            f.write("Free energy estimation for TIP4P water at 273K\n")
            f.write("================================================\n\n")
            f.write("Note: This is a simplified estimation based on average energy terms\n\n")
            
            # Try to parse the GMX energy output to extract averages
            lines = output.split("\n")
            for line in lines:
                if "Average" in line and "=" in line:
                    f.write(f"{line.strip()}\n")
            
            f.write("\nFor more accurate free energy calculations, consider methods like:\n")
            f.write("1. Thermodynamic Integration\n")
            f.write("2. Free Energy Perturbation\n")
            f.write("3. Umbrella Sampling\n")
    except Exception as e:
        print(f"Warning: Could not create free_energy.txt: {e}")
    
    return True, "Successfully estimated free energy components"

def analyze_hydrogen_bonds(output_dir):
    """Analyze hydrogen bonds in the water system"""
    print("Analyzing hydrogen bonds...")
    
    # Check if md.xtc exists, if not check for md.trr and convert
    if not (output_dir / "md.xtc").exists():
        if (output_dir / "md.trr").exists():
            print("Found md.trr, converting to md.xtc...")
            run_command(
                "echo '0' | gmx trjconv -f md.trr -s md.tpr -o md.xtc",
                cwd=output_dir
            )
        else:
            return False, "Error: Neither md.xtc nor md.trr found. Run production_md.py first."
    
    # Run hydrogen bond analysis - for water-water h-bonds, select all atoms twice
    success, output = run_command(
        "echo 'a a' | gmx hbond -f md.xtc -s md.tpr -num hbnum.xvg -life hblife.xvg -dist hbdist.xvg -ang hbang.xvg",
        cwd=output_dir
    )
    
    if not success:
        print(output)
        return False, "Failed to analyze hydrogen bonds"
    
    # Check if hbnum.xvg exists
    if not (output_dir / "hbnum.xvg").exists():
        return False, "Error: hbnum.xvg not found. Hydrogen bond analysis failed."
    
    # Try to extract average number of hydrogen bonds
    try:
        with open(output_dir / "hbond_summary.txt", "w") as f:
            f.write("Hydrogen Bond Analysis for TIP4P water\n")
            f.write("=====================================\n\n")
            
            # Parse the output to find the average number of hydrogen bonds
            lines = output.split("\n")
            for line in lines:
                if "Average" in line and "hydrogen bonds" in line:
                    f.write(f"{line.strip()}\n")
                if "Hydrogen bond lifetime" in line:
                    f.write(f"{line.strip()}\n")
    except Exception as e:
        print(f"Warning: Could not create hbond_summary.txt: {e}")
    
    return True, "Hydrogen bond analysis completed successfully"

def analyze_velocity_autocorrelation(output_dir):
    """Analyze velocity autocorrelation function (VACF)"""
    print("Analyzing velocity autocorrelation function...")
    
    # Check if md.trr exists (we need velocities which are in .trr files)
    if not (output_dir / "md.trr").exists():
        return False, "Error: md.trr not found. VACF analysis requires velocities saved in trajectory (nstvout in md.mdp)."
    
    # Run VACF analysis
    success, output = run_command(
        "echo 'SOL' | gmx velacc -f md.trr -s md.tpr -o vacf.xvg -os vacf_spectrum.xvg -acflen 1000 -nonormalize",
        cwd=output_dir
    )
    
    if not success:
        print(output)
        return False, "Failed to analyze velocity autocorrelation"
    
    # Check if vacf.xvg exists
    if not (output_dir / "vacf.xvg").exists():
        return False, "Error: vacf.xvg not found. VACF analysis failed."
    
    return True, "VACF analysis completed successfully"

def copy_analysis_files(output_dir, analysis_dir, model, temp):
    """Copy analysis files to analysis directory"""
    print("Copying analysis files...")
    
    # Create model/temp subdirectory in analysis
    model_dir = analysis_dir / model
    temp_dir = model_dir / f"{temp}K"
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy XVG files
    xvg_files = list(output_dir.glob("*.xvg"))
    for file in xvg_files:
        print(f"Copying {file.name} to {temp_dir}")
        try:
            import shutil
            shutil.copy2(file, temp_dir / file.name)
        except Exception as e:
            print(f"Error copying {file}: {e}")
    
    # Copy PDB files (average structure)
    pdb_files = list(output_dir.glob("*.pdb"))
    for file in pdb_files:
        if file.name != "water.pdb" and file.name != "water_box.pdb":
            print(f"Copying {file.name} to {temp_dir}")
            try:
                import shutil
                shutil.copy2(file, temp_dir / file.name)
            except Exception as e:
                print(f"Error copying {file}: {e}")
    
    # Copy free energy text file
    txt_files = list(output_dir.glob("*.txt"))
    for file in txt_files:
        print(f"Copying {file.name} to {temp_dir}")
        try:
            import shutil
            shutil.copy2(file, temp_dir / file.name)
        except Exception as e:
            print(f"Error copying {file}: {e}")
    
    return True, f"Successfully copied analysis files to {temp_dir}"

def analyze_results(model, temp, iteration_dir=None):
    # Set base directory
    if iteration_dir:
        base_dir = Path(iteration_dir).resolve()
    else:
        base_dir = Path(__file__).resolve().parent.parent
    
    # Get the scripts directory (always use the central scripts directory)
    scripts_dir = Path(__file__).resolve().parent
    
    print(f"Analyzing results for {model} at {temp}K")
    print(f"Working with iteration directory: {base_dir}")
    print(f"Using scripts from: {scripts_dir}")
    
    # Define directories
    data_dir = base_dir / "data" / model / f"{temp}K"
    sim_dir = data_dir / "outputs"
    analysis_dir = base_dir / "analysis" / model / f"{temp}K"
    plots_dir = analysis_dir / "plots"
    
    # Create analysis directories
    os.makedirs(analysis_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    
    # Check if simulation files exist
    if not os.path.exists(sim_dir / "md.gro") or not os.path.exists(sim_dir / "md.xtc"):
        print(f"Simulation files not found in {sim_dir}")
        return False
    
    # Calculate radial distribution function (RDF)
    print("\nCalculating radial distribution function...")
    if not run_command(f"gmx rdf -f {sim_dir}/md.xtc -s {sim_dir}/md.tpr -o {analysis_dir}/rdf_OO.xvg -ref \"name OW\" -sel \"name OW\" -bin 0.002 -norm rdf", cwd=base_dir):
        print("Error calculating RDF")
    
    # Calculate density
    print("\nCalculating density...")
    if not run_command(f"gmx density -f {sim_dir}/md.xtc -s {sim_dir}/md.tpr -o {analysis_dir}/density.xvg -d Z", cwd=base_dir):
        print("Error calculating density")
    
    # Calculate mean square displacement (MSD)
    print("\nCalculating mean square displacement...")
    if not run_command(f"gmx msd -f {sim_dir}/md.xtc -s {sim_dir}/md.tpr -o {analysis_dir}/msd.xvg -mol -rmcomm", cwd=base_dir):
        print("Error calculating MSD")
    
    # Energy analysis
    print("\nAnalyzing energy components...")
    if not run_command(f"echo \"Potential Kinetic Total\" | gmx energy -f {sim_dir}/md.edr -o {analysis_dir}/energy.xvg", cwd=base_dir):
        print("Error analyzing energy")
    
    # Temperature analysis
    print("\nAnalyzing temperature...")
    if not run_command(f"echo \"Temperature\" | gmx energy -f {sim_dir}/md.edr -o {analysis_dir}/temperature.xvg", cwd=base_dir):
        print("Error analyzing temperature")
    
    # Pressure analysis
    print("\nAnalyzing pressure...")
    if not run_command(f"echo \"Pressure\" | gmx energy -f {sim_dir}/md.edr -o {analysis_dir}/pressure.xvg", cwd=base_dir):
        print("Error analyzing pressure")
    
    # Hydrogen bond analysis
    print("\nAnalyzing hydrogen bonds...")
    if not run_command(f"echo \"a a\" | gmx hbond -f {sim_dir}/md.xtc -s {sim_dir}/md.tpr -num {analysis_dir}/hbnum.xvg -life {analysis_dir}/hblife.xvg -dist {analysis_dir}/hbdist.xvg -ang {analysis_dir}/hbang.xvg", cwd=base_dir):
        print("Error analyzing hydrogen bonds")
    
    # Velocity autocorrelation function analysis
    print("\nAnalyzing velocity autocorrelation function...")
    if os.path.exists(sim_dir / "md.trr"):
        if not run_command(f"echo \"SOL\" | gmx velacc -f {sim_dir}/md.trr -s {sim_dir}/md.tpr -o {analysis_dir}/vacf.xvg -os {analysis_dir}/vacf_spectrum.xvg -acflen 1000 -nonormalize", cwd=base_dir):
            print("Error analyzing velocity autocorrelation function")
    else:
        print(f"Warning: md.trr file not found in {sim_dir}, skipping VACF analysis")
        print("For VACF analysis, ensure velocities are saved in trajectory (nstvout in md.mdp)")
    
    # Generate plots
    print("\nGenerating plots...")
    if not run_command(f"python3 {scripts_dir}/plot_results.py --model {model} --temp {temp} --iteration-dir {base_dir}", cwd=base_dir):
        print("Error generating plots")
    
    print("\nAnalysis completed successfully!")
    print(f"Results are available in {analysis_dir}")
    
    return True

def main():
    parser = argparse.ArgumentParser(description="Analyze water simulation results")
    parser.add_argument("--model", type=str, default="tip4p", help="Water model (default: tip4p)")
    parser.add_argument("--temp", type=int, default=273, help="Temperature in K (default: 273)")
    parser.add_argument("--iteration-dir", type=str, help="Path to the iteration directory (default: parent directory of this script)")
    args = parser.parse_args()
    
    analyze_results(args.model, args.temp, args.iteration_dir)

if __name__ == "__main__":
    main() 