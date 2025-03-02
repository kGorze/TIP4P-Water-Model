#!/bin/bash
# Script to run the complete workflow (simulation and analysis) for md_water_study_iteration_4
# This script uses the centralized scripts from md_water_study_scripts

# Define directories
CENTRAL_SCRIPTS_DIR="/home/konrad_guest/Documents/research/cursor/md_water_study_scripts"
ITERATION_DIR="/home/konrad_guest/Documents/research/cursor/md_water_study_iteration_4"
SCRIPTS_DIR="${ITERATION_DIR}/scripts"
DATA_DIR="${ITERATION_DIR}/data"
CONFIGS_DIR="${ITERATION_DIR}/configs"
ANALYSIS_DIR="${ITERATION_DIR}/analysis"
LOGS_DIR="${ITERATION_DIR}/logs"

# Create necessary directories
mkdir -p "${SCRIPTS_DIR}"
mkdir -p "${DATA_DIR}"
mkdir -p "${ANALYSIS_DIR}"
mkdir -p "${LOGS_DIR}"

# Create log file with timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOGS_DIR}/workflow_${TIMESTAMP}.log"

# Function to log messages to both console and log file
log() {
    echo "$@" | tee -a "${LOG_FILE}"
}

# Redirect all output to log file while also displaying on console
exec > >(tee -a "${LOG_FILE}") 2>&1

log "Starting workflow at $(date)"
log "======================================"

# Prepare configs for GROMACS
log "Setting up simulation..."

# Step 1: Generate water box using Julia and Packmol - using the existing config file
log "Step 1: Generating water box with Packmol via Julia..."
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
julia -e 'using Packmol; run_packmol("water_box.inp")' || { log "Packmol failed"; exit 1; }
log "Water box generated."

# Step 2: Generate topology with GROMACS using TIP4P water model and OPLS-AA force field
log "Step 2: Generating topology with TIP4P water model..."
# Explicitly select force field 15 (OPLS-AA/L) with an echo command
echo "15" | gmx pdb2gmx -f water_box.pdb -o water_box.gro -water tip4p -p topol.top -ff oplsaa || { log "pdb2gmx failed"; exit 1; }
log "Topology generated."

# Step 3: Perform energy minimization
log "Step 3: Running energy minimization..."
gmx grompp -f "${CONFIGS_DIR}/em.mdp" -c water_box.gro -p topol.top -o em.tpr -maxwarn 1 || { log "grompp for em failed"; exit 1; }
gmx mdrun -deffnm em -ntmpi 1 -ntomp 6 || { log "mdrun for em failed"; exit 1; }
log "Energy minimization completed."

# Step 4: Perform NVT equilibration
log "Step 4: Running NVT equilibration..."
gmx grompp -f "${CONFIGS_DIR}/nvt.mdp" -c em.gro -p topol.top -o nvt.tpr || { log "grompp for nvt failed"; exit 1; }
gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 6 || { log "mdrun for nvt failed"; exit 1; }
log "NVT equilibration completed."

# Step 5: Perform NPT equilibration
log "Step 5: Running NPT equilibration..."
gmx grompp -f "${CONFIGS_DIR}/npt.mdp" -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1 || { log "grompp for npt failed"; exit 1; }
gmx mdrun -deffnm npt -ntmpi 1 -ntomp 6 || { log "mdrun for npt failed"; exit 1; }
log "NPT equilibration completed."

# Step 6: Perform production MD
log "Step 6: Running production MD..."
gmx grompp -f "${CONFIGS_DIR}/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o md.tpr || { log "grompp for md failed"; exit 1; }
gmx mdrun -deffnm md -ntmpi 1 -ntomp 6 || { log "mdrun for md failed"; exit 1; }
log "Production MD completed."

# Step 7: Run analysis
log "Step 7: Running analysis..."
cd "${DATA_DIR}"

# Calculate Radial Distribution Functions (RDFs)
log "1. Calculating RDFs..."
# O-O RDF
log "   - Oxygen-Oxygen RDF..."
mkdir -p "${ANALYSIS_DIR}"
gmx rdf -f md.xtc -s md.tpr -o "${ANALYSIS_DIR}/rdf_OO.xvg" -ref "name OW" -sel "name OW" || { log "RDF analysis failed"; exit 1; }

# Calculate Mean Square Displacement (MSD) for diffusion coefficient
log "2. Calculating Mean Square Displacement for diffusion coefficient..."
# Add selection parameter to track oxygen atoms of water molecules
gmx msd -f md.xtc -s md.tpr -o "${ANALYSIS_DIR}/msd.xvg" -beginfit 1000 -endfit 2000 -sel "name OW" || { log "MSD analysis failed"; exit 1; }

# Calculate density profile
log "3. Calculating density profile..."
# Automatically select the Water group (group 1)
echo "1" | gmx density -f md.xtc -s md.tpr -o "${ANALYSIS_DIR}/density.xvg" -dens mass -d Z || { log "Density analysis failed"; exit 1; }

# Calculate thermodynamic properties
log "4. Extracting thermodynamic properties..."
# Temperature
echo "9" | gmx energy -f md.edr -o "${ANALYSIS_DIR}/temperature.xvg" || { log "Temperature analysis failed"; exit 1; }

# Pressure
echo "11" | gmx energy -f md.edr -o "${ANALYSIS_DIR}/pressure.xvg" || { log "Pressure analysis failed"; exit 1; }

# Energy components - select by numbers instead of names
echo "5 6 7" | gmx energy -f md.edr -o "${ANALYSIS_DIR}/energy.xvg" || { log "Energy analysis failed"; exit 1; }

# Hydrogen Bond Analysis
log "5. Performing hydrogen bond analysis..."
# First, let's check what groups are available in our system
log "   - Checking available groups for hydrogen bond analysis..."
# Create a temporary index file to see available groups
gmx make_ndx -f md.tpr -o hbond_groups.ndx << EOF
q
EOF

# For a water-only system, we want to analyze hydrogen bonds between water molecules
# From the output of make_ndx, we can see that SOL is group 2
log "   - Performing hydrogen bond analysis between water molecules..."
echo "2 2" | gmx hbond -f md.xtc -s md.tpr -num "${ANALYSIS_DIR}/hbnum.xvg" -dist "${ANALYSIS_DIR}/hbdist.xvg" -ang "${ANALYSIS_DIR}/hbang.xvg" || { log "Hydrogen bond analysis failed"; exit 1; }
log "Hydrogen bond analysis completed."

# Velocity Autocorrelation Function (VACF) Analysis
log "6. Performing velocity autocorrelation function analysis..."
# Check if md.trr exists (we need velocities which are in .trr files, not in .xtc)
if [ -f "md.trr" ]; then
    echo "SOL" | gmx velacc -f md.trr -s md.tpr -o "${ANALYSIS_DIR}/vacf.xvg" -os "${ANALYSIS_DIR}/vacf_spectrum.xvg" -acflen 1000 -nonormalize || { log "VACF analysis failed"; exit 1; }
    log "VACF analysis completed."
else
    log "Warning: md.trr file not found, cannot perform VACF analysis."
    log "For VACF analysis, ensure velocities are saved in trajectory (nstvout in md.mdp)."
fi

log "======================================"
log "Workflow completed at $(date)!"
log "Results are in ${DATA_DIR} and ${ANALYSIS_DIR}"
log "Log file saved to ${LOG_FILE}" 