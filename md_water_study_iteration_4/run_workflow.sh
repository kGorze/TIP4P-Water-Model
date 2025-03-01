#!/bin/bash
# Script to run the complete workflow (simulation and analysis) for md_water_study_iteration_4
# This script uses the centralized scripts from md_water_study_scripts

set -e  # Exit on error

# Define directories
CENTRAL_SCRIPTS_DIR="/home/konrad_guest/Documents/research/cursor/md_water_study_scripts"
ITERATION_DIR="/home/konrad_guest/Documents/research/cursor/md_water_study_iteration_4"
SCRIPTS_DIR="${ITERATION_DIR}/scripts"
DATA_DIR="${ITERATION_DIR}/data"
CONFIGS_DIR="${ITERATION_DIR}/configs"
ANALYSIS_DIR="${ITERATION_DIR}/analysis"

# Create necessary directories
mkdir -p "${SCRIPTS_DIR}"
mkdir -p "${DATA_DIR}"
mkdir -p "${ANALYSIS_DIR}"

# Prepare configs for GROMACS
echo "Setting up simulation..."

# Step 1: Generate water box using Julia and Packmol - using the existing config file
echo "Step 1: Generating water box with Packmol via Julia..."
cd "${DATA_DIR}"

# Create water molecule template first
cat > water.pdb << EOF
ATOM      1  OW  SOL     1       0.000   0.000   0.000  1.00  0.00            
ATOM      2  HW1 SOL     1       0.957   0.000   0.000  1.00  0.00            
ATOM      3  HW2 SOL     1      -0.240   0.927   0.000  1.00  0.00            
END
EOF

# Copy the existing water_box.inp file from the configs directory
cp "${CONFIGS_DIR}/water_box.inp" ./

# Run Packmol directly using Julia with the existing config file
julia -e 'using Packmol; run_packmol("water_box.inp")'
echo "Water box generated."

# Step 2: Generate topology with GROMACS using TIP4P water model and OPLS-AA force field
echo "Step 2: Generating topology with TIP4P water model..."
# Explicitly select force field 15 (OPLS-AA/L) with an echo command
echo "15" | gmx pdb2gmx -f water_box.pdb -o water_box.gro -water tip4p -p topol.top -ff oplsaa || { echo "pdb2gmx failed"; exit 1; }
echo "Topology generated."

# Step 3: Perform energy minimization
echo "Step 3: Running energy minimization..."
gmx grompp -f "${CONFIGS_DIR}/em.mdp" -c water_box.gro -p topol.top -o em.tpr -maxwarn 1 || { echo "grompp for em failed"; exit 1; }
gmx mdrun -deffnm em -ntmpi 1 -ntomp 6 || { echo "mdrun for em failed"; exit 1; }
echo "Energy minimization completed."

# Step 4: Perform NVT equilibration
echo "Step 4: Running NVT equilibration..."
gmx grompp -f "${CONFIGS_DIR}/nvt.mdp" -c em.gro -p topol.top -o nvt.tpr || { echo "grompp for nvt failed"; exit 1; }
gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 6 || { echo "mdrun for nvt failed"; exit 1; }
echo "NVT equilibration completed."

# Step 5: Perform NPT equilibration
echo "Step 5: Running NPT equilibration..."
gmx grompp -f "${CONFIGS_DIR}/npt.mdp" -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1 || { echo "grompp for npt failed"; exit 1; }
gmx mdrun -deffnm npt -ntmpi 1 -ntomp 6 || { echo "mdrun for npt failed"; exit 1; }
echo "NPT equilibration completed."

# Step 6: Perform production MD
echo "Step 6: Running production MD (2 ns)..."
gmx grompp -f "${CONFIGS_DIR}/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o md.tpr || { echo "grompp for md failed"; exit 1; }
gmx mdrun -deffnm md -ntmpi 1 -ntomp 6 || { echo "mdrun for md failed"; exit 1; }
echo "Production MD completed."

# Step 7: Run analysis
echo "Step 7: Running analysis..."
cd "${DATA_DIR}"

# Calculate Radial Distribution Functions (RDFs)
echo "1. Calculating RDFs..."
# O-O RDF
echo "   - Oxygen-Oxygen RDF..."
mkdir -p "${ANALYSIS_DIR}"
gmx rdf -f md.xtc -s md.tpr -o "${ANALYSIS_DIR}/rdf_OO.xvg" -ref "name OW" -sel "name OW" || { echo "RDF analysis failed"; exit 1; }

# Calculate Mean Square Displacement (MSD) for diffusion coefficient
echo "2. Calculating Mean Square Displacement for diffusion coefficient..."
# Add selection parameter to track oxygen atoms of water molecules
gmx msd -f md.xtc -s md.tpr -o "${ANALYSIS_DIR}/msd.xvg" -beginfit 1000 -endfit 2000 -sel "name OW" || { echo "MSD analysis failed"; exit 1; }

# Calculate density profile
echo "3. Calculating density profile..."
gmx density -f md.xtc -s md.tpr -o "${ANALYSIS_DIR}/density.xvg" -dens mass -d Z || { echo "Density analysis failed"; exit 1; }

# Calculate thermodynamic properties
echo "4. Extracting thermodynamic properties..."
# Temperature
echo "Temperature" | gmx energy -f md.edr -o "${ANALYSIS_DIR}/temperature.xvg" || { echo "Temperature analysis failed"; exit 1; }

# Pressure
echo "Pressure" | gmx energy -f md.edr -o "${ANALYSIS_DIR}/pressure.xvg" || { echo "Pressure analysis failed"; exit 1; }

# Energy components
echo "Potential Kinetic-En. Total-Energy" | gmx energy -f md.edr -o "${ANALYSIS_DIR}/energy.xvg" || { echo "Energy analysis failed"; exit 1; }

# Hydrogen Bond Analysis
echo "5. Performing hydrogen bond analysis..."
echo "a a" | gmx hbond -f md.xtc -s md.tpr -num "${ANALYSIS_DIR}/hbnum.xvg" -life "${ANALYSIS_DIR}/hblife.xvg" -dist "${ANALYSIS_DIR}/hbdist.xvg" -ang "${ANALYSIS_DIR}/hbang.xvg" || { echo "Hydrogen bond analysis failed"; exit 1; }
echo "Hydrogen bond analysis completed."

# Velocity Autocorrelation Function (VACF) Analysis
echo "6. Performing velocity autocorrelation function analysis..."
# Check if md.trr exists (we need velocities which are in .trr files, not in .xtc)
if [ -f "md.trr" ]; then
    echo "SOL" | gmx velacc -f md.trr -s md.tpr -o "${ANALYSIS_DIR}/vacf.xvg" -os "${ANALYSIS_DIR}/vacf_spectrum.xvg" -acflen 1000 -nonormalize || { echo "VACF analysis failed"; exit 1; }
    echo "VACF analysis completed."
else
    echo "Warning: md.trr file not found, cannot perform VACF analysis."
    echo "For VACF analysis, ensure velocities are saved in trajectory (nstvout in md.mdp)."
fi

echo "Workflow completed! Results are in ${DATA_DIR} and ${ANALYSIS_DIR}" 