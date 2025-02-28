#!/bin/bash
# Script to run the workflow with checkpoints for md_water_study_iteration_4

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

# Source the checkpoint system
source "${SCRIPTS_DIR}/checkpoint_system.sh"

echo "Setting up simulation with checkpoint system..."

# Step 1: Generate water box
run_step "water_box" "Generate water box with Packmol" "
  cd \"${DATA_DIR}\"
  # Create water molecule template
  cat > water.pdb << EOF
ATOM      1  OW  SOL     1       0.000   0.000   0.000  1.00  0.00            
ATOM      2  HW1 SOL     1       0.957   0.000   0.000  1.00  0.00            
ATOM      3  HW2 SOL     1      -0.240   0.927   0.000  1.00  0.00            
END
EOF
  # Copy the water_box.inp file from configs
  cp \"${CONFIGS_DIR}/water_box.inp\" ./
  # Run Packmol
  julia -e 'using Packmol; run_packmol(\"water_box.inp\")'
"

# Step 2: Generate topology
run_step "topology" "Generate topology with TIP4P water model" "
  cd \"${DATA_DIR}\"
  echo \"15\" | gmx pdb2gmx -f water_box.pdb -o water_box.gro -water tip4p -p topol.top -ff oplsaa
"

# Step 3: Energy minimization
run_step "energy_minimization" "Run energy minimization" "
  cd \"${DATA_DIR}\"
  gmx grompp -f \"${CONFIGS_DIR}/em.mdp\" -c water_box.gro -p topol.top -o em.tpr -maxwarn 1
  gmx mdrun -deffnm em
"

# Step 4: NVT equilibration
run_step "nvt_equilibration" "Run NVT equilibration" "
  cd \"${DATA_DIR}\"
  gmx grompp -f \"${CONFIGS_DIR}/nvt.mdp\" -c em.gro -p topol.top -o nvt.tpr
  gmx mdrun -deffnm nvt
"

# Step 5: NPT equilibration
run_step "npt_equilibration" "Run NPT equilibration" "
  cd \"${DATA_DIR}\"
  gmx grompp -f \"${CONFIGS_DIR}/npt.mdp\" -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1
  gmx mdrun -deffnm npt
"

# Step 6: Production MD
run_step "production_md" "Run production MD (2 ns)" "
  cd \"${DATA_DIR}\"
  gmx grompp -f \"${CONFIGS_DIR}/md.mdp\" -c npt.gro -t npt.cpt -p topol.top -o md.tpr
  gmx mdrun -deffnm md
"

# Step 7: Analysis
echo "Step 7: Running analysis..."

# Calculate Radial Distribution Functions (RDFs)
run_step "rdf_analysis" "Calculate Oxygen-Oxygen RDF" "
  cd \"${DATA_DIR}\"
  mkdir -p \"${ANALYSIS_DIR}\"
  gmx rdf -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/rdf_OO.xvg\" -ref \"name OW\" -sel \"name OW\"
"

# Calculate Mean Square Displacement (MSD)
run_step "msd_analysis" "Calculate Mean Square Displacement for diffusion coefficient" "
  cd \"${DATA_DIR}\"
  gmx msd -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/msd.xvg\" -beginfit 1000 -endfit 2000 -sel \"name OW\"
"

# Calculate density profile
run_step "density_analysis" "Calculate density profile" "
  cd \"${DATA_DIR}\"
  gmx density -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/density.xvg\" -dens mass -d Z
"

# Calculate thermodynamic properties
run_step "temperature_analysis" "Extract temperature" "
  cd \"${DATA_DIR}\"
  echo \"Temperature\" | gmx energy -f md.edr -o \"${ANALYSIS_DIR}/temperature.xvg\"
"

run_step "pressure_analysis" "Extract pressure" "
  cd \"${DATA_DIR}\"
  echo \"Pressure\" | gmx energy -f md.edr -o \"${ANALYSIS_DIR}/pressure.xvg\"
"

run_step "energy_analysis" "Extract energy components" "
  cd \"${DATA_DIR}\"
  echo \"Potential Kinetic-En. Total-Energy\" | gmx energy -f md.edr -o \"${ANALYSIS_DIR}/energy.xvg\"
"

# Hydrogen Bond Analysis
run_step "hbond_number" "Calculate number of hydrogen bonds" "
  cd \"${DATA_DIR}\"
  echo \"1 1\" | gmx hbond -f md.xtc -s md.tpr -num \"${ANALYSIS_DIR}/hbnum.xvg\"
"

run_step "hbond_dist" "Calculate hydrogen bond distance distribution" "
  cd \"${DATA_DIR}\"
  echo \"1 1\" | gmx hbond -f md.xtc -s md.tpr -dist \"${ANALYSIS_DIR}/hbdist.xvg\"
"

run_step "hbond_ang" "Calculate hydrogen bond angle distribution" "
  cd \"${DATA_DIR}\"
  echo \"1 1\" | gmx hbond -f md.xtc -s md.tpr -ang \"${ANALYSIS_DIR}/hbang.xvg\"
"

run_step "hbond_lifetime" "Calculate hydrogen bond lifetime correlation" "
  cd \"${DATA_DIR}\"
  echo \"1 1\" | gmx hbond-legacy -f md.xtc -s md.tpr -nopbc -life \"${ANALYSIS_DIR}/hblife.xvg\" -ac \"${ANALYSIS_DIR}/hbac.xvg\"
"

# VACF Analysis if md.trr exists
if [ -f "${DATA_DIR}/md.trr" ]; then
  run_step "vacf_analysis" "Perform velocity autocorrelation function analysis" "
    cd \"${DATA_DIR}\"
    echo \"SOL\" | gmx velacc -f md.trr -s md.tpr -o \"${ANALYSIS_DIR}/vacf.xvg\" -os \"${ANALYSIS_DIR}/vacf_spectrum.xvg\" -acflen 1000 -nonormalize
  "
else
  echo "Warning: md.trr file not found, cannot perform VACF analysis."
  echo "For VACF analysis, ensure velocities are saved in trajectory (nstvout in md.mdp)."
fi

# List all completed steps
list_checkpoints

echo "Workflow completed! Results are in ${DATA_DIR} and ${ANALYSIS_DIR}" 