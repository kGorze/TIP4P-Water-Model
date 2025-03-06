#!/bin/bash
# Script to run the complete workflow (simulation and analysis) for md_water_study_iteration_4
# This script includes a checkpoint system to allow resuming from interruptions

set -e  # Exit on error

# Parse command line arguments
ONLY_STEP=""
RERUN_STEP=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --only)
            ONLY_STEP="$2"
            shift 2
            ;;
        --rerun)
            RERUN_STEP="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --only STEP_NAME    Run only the specified step"
            echo "  --rerun STEP_NAME   Rerun the specified step (and continue from there)"
            echo "  --help, -h          Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0                  Run the full workflow"
            echo "  $0 --only generate_plots  Run only the plotting step"
            echo "  $0 --rerun msd_analysis   Rerun the MSD analysis step and continue"
            echo ""
            echo "Available steps:"
            echo "  water_box, topology, energy_minimization, nvt_equilibration,"
            echo "  npt_equilibration, production_md, rdf_analysis_OO, rdf_analysis_OH,"
            echo "  rdf_analysis_HH, msd_analysis, density_analysis, temperature_analysis,"
            echo "  pressure_analysis, energy_analysis, hbond_index, hbond_analysis,"
            echo "  hbond_lifetime, vacf_analysis, additional_analysis, generate_plots"
            echo ""
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Define directories
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CENTRAL_SCRIPTS_DIR="${SCRIPT_DIR}/../md_water_study_scripts"
ITERATION_DIR="${SCRIPT_DIR}"
SCRIPTS_DIR="${ITERATION_DIR}/scripts"
DATA_DIR="${ITERATION_DIR}/data"
CONFIGS_DIR="${ITERATION_DIR}/configs"
ANALYSIS_DIR="${ITERATION_DIR}/analysis"
LOGS_DIR="${ITERATION_DIR}/logs"
CHECKPOINT_DIR="${ITERATION_DIR}/checkpoints"

# Create necessary directories
mkdir -p "${SCRIPTS_DIR}"
mkdir -p "${DATA_DIR}"
mkdir -p "${ANALYSIS_DIR}"
mkdir -p "${ANALYSIS_DIR}/data"
mkdir -p "${ANALYSIS_DIR}/plots"
mkdir -p "${LOGS_DIR}"
mkdir -p "${CHECKPOINT_DIR}"

# Create log file with timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOGS_DIR}/workflow_${TIMESTAMP}.log"

# Function to log messages to both console and log file
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "${LOG_FILE}"
}

# Redirect all output to log file while also displaying on console
exec > >(tee -a "${LOG_FILE}") 2>&1

# Checkpoint system functions
# Function to create a checkpoint
create_checkpoint() {
    local step_name="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S')" > "${CHECKPOINT_DIR}/${step_name}.done"
    log "Checkpoint created: ${step_name}"
}

# Function to check if a checkpoint exists
check_checkpoint() {
    local step_name="$1"
    if [ -f "${CHECKPOINT_DIR}/${step_name}.done" ]; then
        return 0  # Checkpoint exists
    else
        return 1  # Checkpoint doesn't exist
    fi
}

# Function to remove a checkpoint
remove_checkpoint() {
    local step_name="$1"
    if [ -f "${CHECKPOINT_DIR}/${step_name}.done" ]; then
        rm "${CHECKPOINT_DIR}/${step_name}.done"
        log "Checkpoint removed: ${step_name}"
    else
        log "No checkpoint found for: ${step_name}"
    fi
}

# Function to list all checkpoints
list_checkpoints() {
    log "Completed steps:"
    for checkpoint in "${CHECKPOINT_DIR}"/*.done; do
        if [ -f "$checkpoint" ]; then
            step_name=$(basename "$checkpoint" .done)
            completion_time=$(cat "$checkpoint")
            log "  - ${step_name} (completed at ${completion_time})"
        fi
    done
}

# Function to run a step with checkpoint checking
run_step() {
    local step_name="$1"
    local step_description="$2"
    local command="$3"
    
    # If we're only running a specific step, skip others
    if [ -n "${ONLY_STEP}" ] && [ "${ONLY_STEP}" != "${step_name}" ]; then
        return 0
    fi
    
    # If we're rerunning a specific step, remove its checkpoint
    if [ "${RERUN_STEP}" = "${step_name}" ]; then
        remove_checkpoint "${step_name}"
        # Clear the rerun flag so subsequent steps run normally
        RERUN_STEP=""
    fi
    
    if check_checkpoint "${step_name}" && [ -z "${ONLY_STEP}" ]; then
        log "Step '${step_description}' already completed. Skipping..."
    else
        log "Running step: ${step_description}..."
        eval "${command}"
        if [ $? -eq 0 ]; then
            create_checkpoint "${step_name}"
            log "Step completed: ${step_description}"
        else
            log "Step failed: ${step_description}"
            return 1
        fi
    fi
    return 0
}

log "Starting workflow at $(date)"
log "======================================"

# Step 1: Generate water box using Julia and Packmol
run_step "water_box" "Generate water box with Packmol" "
  cd \"${DATA_DIR}\"
  
  # Create water molecule template
  cat > water.pdb << 'EOF_WATER'
ATOM      1  OW  SOL     1       0.000   0.000   0.000  1.00  0.00            
ATOM      2  HW1 SOL     1       0.957   0.000   0.000  1.00  0.00            
ATOM      3  HW2 SOL     1      -0.240   0.927   0.000  1.00  0.00            
END
EOF_WATER
  
  # Copy the water_box.inp file from configs
  cp \"${CONFIGS_DIR}/water_box.inp\" ./
  
  # Run Packmol
  julia -e 'using Packmol; run_packmol(\"water_box.inp\")' || { log \"Packmol failed\"; exit 1; }
  log \"Water box generated.\"
"

# Step 2: Generate topology with GROMACS using TIP4P water model and OPLS-AA force field
run_step "topology" "Generate topology with TIP4P water model" "
  cd \"${DATA_DIR}\"
  # Explicitly select force field 15 (OPLS-AA/L) with an echo command
  echo \"15\" | gmx pdb2gmx -f water_box.pdb -o water_box.gro -water tip4p -p topol.top -ff oplsaa || { log \"pdb2gmx failed\"; exit 1; }
  log \"Topology generated.\"
"

# Step 3: Perform energy minimization
run_step "energy_minimization" "Run energy minimization" "
  cd \"${DATA_DIR}\"
  gmx grompp -f \"${CONFIGS_DIR}/em.mdp\" -c water_box.gro -p topol.top -o em.tpr -maxwarn 1 || { log \"grompp for em failed\"; exit 1; }
  gmx mdrun -deffnm em -ntmpi 1 -ntomp 6 || { log \"mdrun for em failed\"; exit 1; }
  log \"Energy minimization completed.\"
"

# Step 4: Perform NVT equilibration
run_step "nvt_equilibration" "Run NVT equilibration" "
  cd \"${DATA_DIR}\"
  gmx grompp -f \"${CONFIGS_DIR}/nvt.mdp\" -c em.gro -p topol.top -o nvt.tpr || { log \"grompp for nvt failed\"; exit 1; }
  gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 6 || { log \"mdrun for nvt failed\"; exit 1; }
  log \"NVT equilibration completed.\"
"

# Step 5: Perform NPT equilibration
run_step "npt_equilibration" "Run NPT equilibration" "
  cd \"${DATA_DIR}\"
  gmx grompp -f \"${CONFIGS_DIR}/npt.mdp\" -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1 || { log \"grompp for npt failed\"; exit 1; }
  gmx mdrun -deffnm npt -ntmpi 1 -ntomp 6 || { log \"mdrun for npt failed\"; exit 1; }
  log \"NPT equilibration completed.\"
"

# Step 6: Perform production MD
run_step "production_md" "Run production MD" "
  cd \"${DATA_DIR}\"
  gmx grompp -f \"${CONFIGS_DIR}/md.mdp\" -c npt.gro -t npt.cpt -p topol.top -o md.tpr || { log \"grompp for md failed\"; exit 1; }
  gmx mdrun -deffnm md -ntmpi 1 -ntomp 6 || { log \"mdrun for md failed\"; exit 1; }
  log \"Production MD completed.\"
"

# Step 7: Run analysis
log "Step 7: Running analysis..."

# Calculate Radial Distribution Functions (RDFs)
run_step "rdf_analysis_OO" "Calculate Oxygen-Oxygen RDF" "
  cd \"${DATA_DIR}\"
  gmx rdf -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/data/rdf_OO.xvg\" -ref \"name OW\" -sel \"name OW\" || { log \"RDF O-O analysis failed\"; exit 1; }
"

run_step "rdf_analysis_OH" "Calculate Oxygen-Hydrogen RDF" "
  cd \"${DATA_DIR}\"
  gmx rdf -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/data/rdf_OH.xvg\" -ref \"name OW\" -sel \"name HW1 or name HW2\" || { log \"RDF O-H analysis failed\"; exit 1; }
"

run_step "rdf_analysis_HH" "Calculate Hydrogen-Hydrogen RDF" "
  cd \"${DATA_DIR}\"
  gmx rdf -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/data/rdf_HH.xvg\" -ref \"name HW1 or name HW2\" -sel \"name HW1 or name HW2\" || { log \"RDF H-H analysis failed\"; exit 1; }
"

# Calculate Mean Square Displacement (MSD) for diffusion coefficient
run_step "msd_analysis" "Calculate Mean Square Displacement for diffusion coefficient" "
  cd \"${DATA_DIR}\"
  # Add selection parameter to track oxygen atoms of water molecules
  # Add -dt parameter to match the trajectory output frequency (50 ps)
  # And adjust -trestart to be larger than dt
  gmx msd -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/data/msd.xvg\" -beginfit 1000 -endfit 2000 -sel \"name OW\" -dt 50 -trestart 100 || { log \"MSD analysis failed\"; exit 1; }
"

# Calculate density profile
run_step "density_analysis" "Calculate density profile" "
  cd \"${DATA_DIR}\"
  # Automatically select the Water group (group 1)
  echo \"1\" | gmx density -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/data/density.xvg\" -dens mass -d Z || { log \"Density analysis failed\"; exit 1; }
"

# Calculate thermodynamic properties
run_step "temperature_analysis" "Extract temperature" "
  cd \"${DATA_DIR}\"
  # Temperature
  echo \"9\" | gmx energy -f md.edr -o \"${ANALYSIS_DIR}/data/temperature.xvg\" || { log \"Temperature analysis failed\"; exit 1; }
"

run_step "pressure_analysis" "Extract pressure" "
  cd \"${DATA_DIR}\"
  # Pressure
  echo \"11\" | gmx energy -f md.edr -o \"${ANALYSIS_DIR}/data/pressure.xvg\" || { log \"Pressure analysis failed\"; exit 1; }
"

run_step "energy_analysis" "Extract energy components" "
  cd \"${DATA_DIR}\"
  # Energy components - select by numbers instead of names
  echo \"5 6 7\" | gmx energy -f md.edr -o \"${ANALYSIS_DIR}/data/energy.xvg\" || { log \"Energy analysis failed\"; exit 1; }
"

# Hydrogen Bond Analysis
run_step "hbond_index" "Create index file for hydrogen bond analysis" "
  cd \"${DATA_DIR}\"
  # Create a proper index file for hydrogen bond analysis
  log \"   - Creating index file for hydrogen bond analysis...\"
  # Create a custom index file with water groups
  cat > hbond.ndx.input << 'EOF_NDX'
q
EOF_NDX
  gmx make_ndx -f md.tpr -o hbond.ndx < hbond.ndx.input || { log \"Failed to create index file\"; exit 1; }
"

run_step "hbond_analysis" "Perform hydrogen bond analysis" "
  cd \"${DATA_DIR}\"
  # For a water-only system, we want to analyze hydrogen bonds between water molecules
  # Use the index file with proper group selection
  log \"   - Performing hydrogen bond analysis between water molecules...\"
  # Use the Water group (group 1) for donor and acceptor
  cat > hbond.input << 'EOF_HBOND'
1
1
EOF_HBOND
  gmx hbond -f md.xtc -s md.tpr -n hbond.ndx -num \"${ANALYSIS_DIR}/data/hbnum.xvg\" -dist \"${ANALYSIS_DIR}/data/hbdist.xvg\" -ang \"${ANALYSIS_DIR}/data/hbang.xvg\" < hbond.input || { log \"Hydrogen bond analysis failed\"; exit 1; }
  log \"Hydrogen bond analysis completed.\"
"

# Hydrogen Bond Lifetime Correlation
run_step "hbond_lifetime" "Calculate hydrogen bond lifetime correlation" "
  cd \"${DATA_DIR}\"
  # Use the Water group (group 1) for donor and acceptor
  cat > hbond.input << 'EOF_HBLIFE'
1
1
EOF_HBLIFE
  
  # Run gmx hbond command for hydrogen bond lifetime analysis
  echo \"1\n1\" | gmx hbond -f md.xtc -s md.tpr -n hbond.ndx -life \"${ANALYSIS_DIR}/data/hblife.xvg\" -ac \"${ANALYSIS_DIR}/data/hbac.xvg\" || {
    # If that fails, try the legacy version
    gmx hbond-legacy -f md.xtc -s md.tpr -life \"${ANALYSIS_DIR}/data/hblife.xvg\" -ac \"${ANALYSIS_DIR}/data/hbac.xvg\" < hbond.input || {
      log \"Error: Both hbond and hbond-legacy commands failed. Cannot perform hydrogen bond lifetime analysis.\"
      exit 1
    }
  }
  
  log \"Hydrogen bond lifetime analysis completed.\"
"

# Velocity Autocorrelation Function (VACF) Analysis
run_step "vacf_analysis" "Perform velocity autocorrelation function analysis" "
  cd \"${DATA_DIR}\"
  # Check if md.trr exists (we need velocities which are in .trr files, not in .xtc)
  if [ -f \"md.trr\" ]; then
      # Create index file for VACF analysis if it doesn't exist
      if [ ! -f \"index.ndx\" ]; then
          echo \"q\" | gmx make_ndx -f md.tpr -o index.ndx || { log \"Failed to create index file for VACF analysis\"; exit 1; }
      fi
      
      # Run VACF analysis with the index file
      echo \"SOL\" | gmx velacc -f md.trr -s md.tpr -n index.ndx -o \"${ANALYSIS_DIR}/data/vacf.xvg\" -os \"${ANALYSIS_DIR}/data/vacf_spectrum.xvg\" -acflen 1000 -nonormalize || { log \"VACF analysis failed\"; exit 1; }
      log \"VACF analysis completed.\"
  else
      log \"Warning: md.trr file not found, generating velocities from positions...\"
      
      # Generate velocities using a more robust approach
      # First, create a new TPR file with velocity generation enabled
      gmx grompp -f \"${CONFIGS_DIR}/md.mdp\" -c md.gro -p topol.top -o md_with_vel.tpr -maxwarn 2 || { log \"TPR generation failed\"; exit 1; }
      
      # Then run a short simulation to generate velocities
      # We'll use -rerun to avoid changing the original trajectory
      gmx mdrun -s md_with_vel.tpr -rerun md.xtc -o md_with_vel.trr -v || { log \"TRR generation failed\"; exit 1; }
      
      # Now try VACF analysis with the generated TRR file
      if [ -f \"md_with_vel.trr\" ]; then
          # Create index file for VACF analysis
          echo \"q\" | gmx make_ndx -f md_with_vel.tpr -o index.ndx || { log \"Failed to create index file for VACF analysis\"; exit 1; }
          
          # Run VACF analysis with the index file
          echo \"SOL\" | gmx velacc -f md_with_vel.trr -s md_with_vel.tpr -n index.ndx -o \"${ANALYSIS_DIR}/data/vacf.xvg\" -os \"${ANALYSIS_DIR}/data/vacf_spectrum.xvg\" -acflen 1000 -nonormalize || { log \"VACF analysis failed\"; exit 1; }
          log \"VACF analysis completed after generating velocities.\"
      else
          log \"Failed to generate velocity data. VACF analysis cannot be performed.\"
          exit 1
      fi
  fi
"

# Additional Analysis Steps for Missing Data
run_step "additional_analysis" "Generate additional data files for plotting" "
  cd \"${DATA_DIR}\"
  
  # 1. RMSD Analysis
  log \"Generating RMSD data...\"
  echo -e \"System\nSystem\" | gmx rms -s md.tpr -f md.xtc -o \"${ANALYSIS_DIR}/data/rmsd.xvg\" -tu ns || { log \"RMSD analysis failed\"; exit 1; }
  
  # 2. Potential Energy
  log \"Extracting potential energy...\"
  echo \"10\" | gmx energy -f md.edr -o \"${ANALYSIS_DIR}/data/potential.xvg\" || { log \"Potential energy extraction failed\"; exit 1; }
  
  # 3. Energy Terms
  log \"Extracting energy terms...\"
  echo \"9 10 11 12 13 14\" | gmx energy -f md.edr -o \"${ANALYSIS_DIR}/data/energy_terms.xvg\" || { log \"Energy terms extraction failed\"; exit 1; }
  
  # 4. Radial Density Map
  log \"Generating radial density map...\"
  
  # Try with gmx spatial (which is more accurate for radial density)
  echo \"1\" | gmx spatial -s md.tpr -f md.xtc -nab 50 -b 0 -e 2000 -bin 0.05 -xvg none -od \"${ANALYSIS_DIR}/data/density_radial.dat\" || {
    log \"Standard spatial density map generation failed, trying densmap approach...\"
    
    # Try with gmx densmap which is more reliable
    echo \"1\" | gmx densmap -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/data/density_map.xpm\" -bin 0.05 || {
      log \"Both spatial and densmap failed, trying alternative densmap approach...\"
      
      # Try with a different selection and parameters
      echo \"SOL\" | gmx densmap -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/data/density_map.xpm\" -bin 0.1 -unit nm-3 || {
        log \"All density map generation approaches failed.\"
        exit 1
      }
    }
  }
  
  log \"Additional data files generated successfully.\"
"

# Step 8: Generate plots and summary report
run_step "generate_plots" "Generate plots and summary report" "
  # Set up environment variables
  export PYTHONPATH=\"\"
  unset PYTHONHOME
  
  # Create plots directory if it doesn't exist
  mkdir -p \"${ANALYSIS_DIR}/plots\"
  
  # Run each plotting script individually
  echo \"Running RDF plotting script...\"
  python3 \"${ANALYSIS_DIR}/plotting_scripts/plot_rdf.py\" \"${ANALYSIS_DIR}\" \"${ANALYSIS_DIR}/plots\"
  
  echo \"Running density plotting script...\"
  python3 \"${ANALYSIS_DIR}/plotting_scripts/plot_density.py\" \"${ANALYSIS_DIR}\" \"${ANALYSIS_DIR}/plots\"
  
  echo \"Running MSD plotting script...\"
  python3 \"${ANALYSIS_DIR}/plotting_scripts/plot_msd.py\" \"${ANALYSIS_DIR}\" \"${ANALYSIS_DIR}/plots\"
  
  echo \"Running hydrogen bond plotting script...\"
  python3 \"${ANALYSIS_DIR}/plotting_scripts/plot_hbond.py\" \"${ANALYSIS_DIR}\" \"${ANALYSIS_DIR}/plots\"
  
  echo \"Running temperature plotting script...\"
  python3 \"${ANALYSIS_DIR}/plotting_scripts/plot_temperature.py\" \"${ANALYSIS_DIR}\" \"${ANALYSIS_DIR}/plots\"
  
  echo \"Running energy analysis plotting script...\"
  python3 \"${ANALYSIS_DIR}/plotting_scripts/plot_energy_analysis.py\" \"${ANALYSIS_DIR}\" \"${ANALYSIS_DIR}/plots\"
  
  echo \"Running RMSD plotting script...\"
  python3 \"${ANALYSIS_DIR}/plotting_scripts/plot_rmsd.py\" \"${ANALYSIS_DIR}\" \"${ANALYSIS_DIR}/plots\"
  
  # Check if VACF data exists
  if [ -f \"${ANALYSIS_DIR}/data/vacf.xvg\" ]; then
      echo \"Running VACF plotting script...\"
      python3 \"${ANALYSIS_DIR}/plotting_scripts/plot_vacf.py\" \"${ANALYSIS_DIR}\" \"${ANALYSIS_DIR}/plots\"
  else
      echo \"VACF data not found, skipping VACF plots\"
  fi
  
  # Generate summary report
  echo \"Generating summary report...\"
  python3 \"${ANALYSIS_DIR}/plotting_scripts/generate_summary_report.py\" \"${ANALYSIS_DIR}\" \"${ANALYSIS_DIR}/plots\"
"

# List all completed steps
list_checkpoints

log "======================================"
log "Workflow completed at $(date)!"
log "Results are in ${DATA_DIR} and ${ANALYSIS_DIR}"
log "Plots and summary report are in ${ANALYSIS_DIR}/plots"
log "Summary report is available at ${ANALYSIS_DIR}/plots/tip4p_water_analysis_summary.png"
log "Text summary is available at ${ANALYSIS_DIR}/plots/tip4p_water_analysis_summary.txt"
log "Log file saved to ${LOG_FILE}"

# Print a helpful message about how to rerun specific steps
log ""
log "To rerun a specific step, use: ./run_workflow.sh --rerun STEP_NAME"
log "To run only a specific step, use: ./run_workflow.sh --only STEP_NAME"
log ""
log "For example, to regenerate only the plots:"
log "./run_workflow.sh --only generate_plots" 