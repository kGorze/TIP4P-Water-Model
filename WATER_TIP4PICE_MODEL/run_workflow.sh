#!/bin/bash
# Script to run the complete workflow (simulation and analysis) for md_water_study_iteration_4
# This script includes a checkpoint system to allow resuming from interruptions
# Using TIP4P/Ice water model with OPLS-AA force field

set -e  # Exit on error

# Parse command line arguments
ONLY_STEP=""
RERUN_STEP=""
CLEANUP="no"

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
        --cleanup)
            CLEANUP="yes"
            shift
            ;;
        --help|-h)
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --cleanup          Remove temporary files after completion"
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

# Set GMXLIB to include the current directory for force field searching
export GMXLIB="${ITERATION_DIR}"

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

# If you are sure you're running this under bash and want combined output,
# uncomment the next line. Otherwise, if you see duplicate log messages,
# comment it out.
# exec > >(tee -a "${LOG_FILE}") 2>&1

log "Starting workflow at $(date)"
log "======================================"
log "Set GMXLIB=${GMXLIB} for force field searching"

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
            log "  - ${step_name} \(completed at ${completion_time}\)"
        fi
    done
}

# Function to run a step with checkpoint checking
run_step() {
    local step_name="$1"
    local step_description="$2"
    local step_commands="$3"
    
    # Debug: Print step information
    echo "DEBUG: Running step $step_name, RERUN_STEP=$RERUN_STEP, checkpoint=$CHECKPOINT_DIR/$step_name.done"
    
    # Check if we should skip this step
    if [ -n "${ONLY_STEP}" ] && [ "${ONLY_STEP}" != "${step_name}" ]; then
        log "Step '${step_description}' skipped (--only specified)"
        return
    fi
    
    # Check if we should force rerun this step
    if [ -n "${RERUN_STEP}" ] && [ "${RERUN_STEP}" == "${step_name}" ]; then
        log "Forcing rerun of step '${step_description}' (--rerun specified)"
        remove_checkpoint "${step_name}"
    fi
    
    # Check if step is already completed and required files exist
    if [ -f "${CHECKPOINT_DIR}/${step_name}.done" ]; then
        # For additional_analysis step, check if required files exist
        if [ "${step_name}" == "additional_analysis" ]; then
            if [ ! -f "${ANALYSIS_DIR}/data/density_map.xpm" ]; then
                log "Checkpoint exists but required files are missing. Rerunning step '${step_description}'"
                remove_checkpoint "${step_name}"
            else
        log "Step '${step_description}' already completed. Skipping..."
                return
            fi
        else
            log "Step '${step_description}' already completed. Skipping..."
            return
        fi
    fi
    
    # Run the step
        log "Running step: ${step_description}..."
    eval "${step_commands}"
    
    # Create checkpoint file
            create_checkpoint "${step_name}"
}

# Function to clean up temporary files
cleanup() {
    if [ "${CLEANUP}" = "yes" ]; then
        log "Performing selective cleanup..."
        
        # Only remove truly temporary files - properly handle wildcards
        find "${DATA_DIR}" -name "*.log" -type f -delete
        find "${DATA_DIR}" -name "*.mdp" -type f -delete
        find "${DATA_DIR}" -name "*.edr" -type f -delete
        
        # Preserve important files
        log "Preserving simulation files:"
        log "  - Trajectories (.xtc, .trr)"
        log "  - Topology files (.top, .itp)"
        log "  - Checkpoint files (.cpt)"
        log "  - Coordinate files (.gro, .pdb)"
        log "  - Analysis data files (.xvg, .dat)"
        
        log "Selective cleanup completed. Use --cleanup to remove temporary files."
    else
        log "Skipping cleanup. All simulation files preserved."
        log "To remove temporary files, run with --cleanup flag."
    fi
}

# Step 1: Generate water box using Julia and Packmol
run_step "water_box" "Generate water box with Packmol" "
  cd "${DATA_DIR}"
  
  # Check if water.pdb already exists and use it
  if [ ! -f water.pdb ]; then
    log "Warning: water.pdb not found. Creating a default template."
    # Create water molecule template only if it doesn't exist
    cat > water.pdb << 'EOF_WATER'
ATOM      1  OW  SOL     1       0.000   0.000   0.000  1.00  0.00            
ATOM      2  HW1 SOL     1       0.957   0.000   0.000  1.00  0.00            
ATOM      3  HW2 SOL     1      -0.240   0.927   0.000  1.00  0.00            
END
EOF_WATER
  else
    log "Using existing water.pdb file as template"
  fi
  
  # Copy the water_box.inp file from configs
  cp "${CONFIGS_DIR}/water_box.inp" ./
  
  # Run Packmol
  julia -e 'using Packmol; run_packmol("water_box.inp")' || { log "Packmol failed"; exit 1; }
  log "Water box generated."
"

# Step 1.5: Check available force fields and water models
run_step "check_ff" "Check available force fields and water models" "
  cd "${DATA_DIR}"
  
  # Check available force fields
  log "Checking available force fields..."
  gmx pdb2gmx -h | grep -A 5 "Force fields" | tee ff_list.txt
  
  # Check available water models
  log "Checking available water models..."
  gmx pdb2gmx -h | grep -A 10 "Water models" | tee water_models.txt
  
  # List the force field directories in GMXLIB
  log "Force field directories in GMXLIB=${GMXLIB}:"
  ls -la "${GMXLIB}"/*.ff 2>/dev/null || echo "No .ff directories found in GMXLIB"
  
  # Check if oplsaa.ff exists
  if [ -d "${GMXLIB}/oplsaa.ff" ]; then
    log "Found oplsaa.ff in GMXLIB"
    ls -la "${GMXLIB}/oplsaa.ff"
  else
    log "oplsaa.ff not found in GMXLIB, checking system locations..."
    find /usr/local/gromacs/share/gromacs/top -name "oplsaa.ff" 2>/dev/null || echo "oplsaa.ff not found in system locations"
  fi
"

# Step 2: Generate topology with GROMACS using TIP4P/Ice water model and OPLS-AA force field
run_step "topology" "Generate topology with TIP4P/Ice water model" "
  cd "${DATA_DIR}"
  
  # First, we need to create a proper topology using GROMACS's built-in tools
  # We'll use a direct approach to create a topology for our water box
  log "Creating topology with TIP4P/Ice water model..."
  
  # Convert the PDB file to GRO format with a larger box size to avoid overlaps
  log "Converting water_box.pdb to GRO format with a larger box size..."
  # Use a much larger box size (15 nm) to avoid overlaps during energy minimization
  if [ -f water_box.pdb ]; then
    gmx editconf -f water_box.pdb -o water_box.gro -box 15 15 15 -center 7.5 7.5 7.5 || { log "Failed to convert PDB to GRO"; exit 1; }
    log "Conversion successful."
  else
    log "Error: water_box.pdb not found!"
    exit 1
  fi
  
  # Count the number of water molecules in the water_box.pdb file
  log "Counting water molecules in water_box.pdb..."
  water_count=$(grep -c "OW  SOL" water_box.pdb || echo 5500)
  if [ -z "$water_count" ] || [ "$water_count" -eq 0 ]; then
    log "Could not count water molecules, using default value of 5500"
    water_count=5500
  fi
  log "Found ${water_count} water molecules"
  
  # Create a topology file with TIP4P/Ice parameters
  log "Creating topology file with TIP4P/Ice parameters..."
  cat > topol.top << EOF_TOPOL
; TIP4P/Ice water topology
; Generated by run_workflow.sh

; Include force field parameters
#define _FF_OPLS
#define _FF_OPLSAA

[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               3               yes             0.5     0.5

; Include the local force field files
#include "./ffnonbonded.itp"
#include "./ffbonded.itp"

; Include TIP4P/Ice water model parameters
[ atomtypes ]
; name     mass      charge   ptype    sigma        epsilon
IW     0             0.000       D   0.0           0.0
OWT4   15.99940      0.000       A   0.31668       0.88211
HW     1.00800       0.000       A   0.00000E+00   0.00000E+00

[ moleculetype ]
; name nrexcl
SOL  1

[ atoms ]
; nr type   resnr  residu atom  cgnr  charge   mass
1     OWT4  1      SOL    OW    1     0        15.9994
2     HW    1      SOL    HW1   1     0.5897   1.008
3     HW    1      SOL    HW2   1     0.5897   1.008
4     IW    1      SOL    MW    1    -1.1794   0.0

[ constraints ]
; i   j   funct   length
1       2       1       0.09572
1       3       1       0.09572
2       3       1       0.15139

[ exclusions ]
1       2       3       4
2       1       3       4
3       1       2       4
4       1       2       3

[ dummies3 ]
; Dummy from    funct   a         b
4       1       2       3       1       0.13458   0.13458

[ system ]
TIP4P/Ice Water System

[ molecules ]
; Compound        nmols
SOL               ${water_count}
EOF_TOPOL
  
  # Verify the topology
  log "Verifying topology:"
  grep -A 5 "[ system ]" topol.top
  grep -A 2 "[ molecules ]" topol.top
  
  # Verify that the required files exist
  if [ ! -f water_box.gro ]; then
    log "Error: water_box.gro file was not created!"
    exit 1
  fi
  
  if [ ! -f topol.top ]; then
    log "Error: topol.top file was not created!"
    exit 1
  fi
  
  log "Topology and coordinate files created successfully with TIP4P/Ice water model."
"

# Step 3: Perform energy minimization
run_step "energy_minimization" "Run energy minimization" "
  cd "${DATA_DIR}"
  
  # Create an improved em.mdp file with PME electrostatics
  log "Creating improved energy minimization configuration file..."
  cat > em_improved.mdp << 'EOF_EM'
; Energy minimization parameters for water box (TIP4P/Ice model)

; Run parameters
integrator               = steep     ; Steepest descent energy minimization
emtol                    = 1000.0    ; Stop when max force < 1000 kJ/mol/nm
emstep                   = 0.00001   ; Very small initial step-size for extremely gentle minimization
nsteps                   = 100000    ; Maximum number of steps
nstxout                  = 500       ; Write coordinates every 500 steps

; Neighbor searching
cutoff-scheme           = Verlet    ; Neighbor search method
nstlist                 = 20        ; Update neighbor list frequency
ns_type                 = grid      ; Method to determine neighbor list
pbc                     = xyz       ; Periodic Boundary Condition in all directions
rlist                   = 1.0       ; Cut-off distance for the short-range neighbor list

; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.0       ; Short-range electrostatic cut-off (nm)
pme_order               = 4         ; Cubic interpolation
fourierspacing          = 0.16      ; Grid spacing for FFT

; VdW
vdwtype                 = Cut-off
rvdw                    = 1.0       ; Short-range van der Waals cut-off (nm)
DispCorr                = EnerPres  ; Apply long range dispersion corrections

; Temperature and Pressure
tcoupl                  = no        ; No temperature coupling during minimization
pcoupl                  = no        ; No pressure coupling during minimization

; Constraints - turning off constraints for minimization
constraints             = none      ; No constraints during minimization to allow adjustments in bond lengths
constraint_algorithm    = Lincs     ; Will use LINCS after minimization

; COM motion removal
comm-mode               = Linear    ; Remove center of mass translation
nstcomm                 = 100       ; Frequency for center of mass motion removal
EOF_EM
  
  # Check if the required files exist
  if [ ! -f water_box.gro ]; then
    log "Error: water_box.gro file not found! Running topology step again..."
    remove_checkpoint "topology"
    exit 1
  fi
  
  if [ ! -f topol.top ]; then
    log "Error: topol.top file not found! Running topology step again..."
    remove_checkpoint "topology"
    exit 1
  fi
  
  # Run grompp with increased warning tolerance
  log "Running grompp for energy minimization..."
  gmx grompp -f em_improved.mdp -c water_box.gro -p topol.top -o em.tpr -maxwarn 10 || { log "grompp for em failed"; exit 1; }
  
  # Run mdrun with appropriate parallelization
  log "Running energy minimization..."
  gmx mdrun -v -deffnm em -ntmpi 1 -ntomp 6 || { log "mdrun for em failed"; exit 1; }
  
  log "Energy minimization completed."
"

# Step 4: Perform NVT equilibration
run_step "nvt_equilibration" "Run NVT equilibration" "
  cd "${DATA_DIR}"
  
  # Create an improved nvt.mdp file with PME electrostatics
  log "Creating improved NVT equilibration configuration file..."
  cat > nvt_improved.mdp << 'EOF_NVT'
; NVT equilibration parameters for TIP4P/Ice water at 273K
integrator               = md        ; leap-frog integrator
nsteps                   = 50000     ; 25 ps with 0.5 fs timestep
dt                       = 0.0005    ; 0.5 fs - very small timestep for stability
nstxout                  = 5000      ; save coordinates every 2.5 ps
nstvout                  = 5000      ; save velocities every 2.5 ps
nstenergy                = 5000      ; save energies every 2.5 ps
nstlog                   = 5000      ; update log file every 2.5 ps

; Bond parameters
continuation             = no        ; first dynamics run
constraint_algorithm     = lincs     ; holonomic constraints 
constraints              = h-bonds   ; constrain bonds with H atoms
lincs_iter               = 4         ; increased accuracy of LINCS
lincs_order              = 8         ; increased accuracy

; Neighbor searching
cutoff-scheme           = Verlet    ; Verlet cutoff scheme (efficient on modern hardware)
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; update neighbor list every 20 steps
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)

; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT

; Temperature coupling
tcoupl                  = V-rescale ; velocity rescaling thermostat (more accurate than Berendsen)
tc-grps                 = System    ; couple entire system as one group
tau_t                   = 1.0       ; longer time constant for gentler coupling
ref_t                   = 273       ; reference temperature 273K (0°C)

; Pressure coupling
pcoupl                  = no        ; no pressure coupling in NVT

; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC

; Dispersion correction
DispCorr                = EnerPres  ; account for cut-off vdW scheme (improves energy and pressure accuracy)

; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 100       ; lower temperature for initial velocities
gen_seed                = -1        ; generate a random seed
EOF_NVT
  
  # Run grompp with increased warning tolerance
  log "Running grompp for NVT equilibration..."
  gmx grompp -f nvt_improved.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 10 || { log "grompp for nvt failed"; exit 1; }
  
  # Run mdrun with appropriate parallelization
  log "Running NVT equilibration..."
  gmx mdrun -v -deffnm nvt -ntmpi 1 -ntomp 6 || { log "mdrun for nvt failed"; exit 1; }
  
  log "NVT equilibration completed."
"

# Step 5: Perform NPT equilibration
run_step "npt_equilibration" "Run NPT equilibration" "
  cd "${DATA_DIR}"
  
  # Create an improved npt.mdp file with PME electrostatics
  log "Creating improved NPT equilibration configuration file..."
  cat > npt_improved.mdp << 'EOF_NPT'
; NPT equilibration parameters for TIP4P/Ice water at 273K
integrator               = md        ; leap-frog integrator
nsteps                   = 50000     ; 25 ps with 0.5 fs timestep
dt                       = 0.0005    ; 0.5 fs - very small timestep for stability
nstxout                  = 5000      ; save coordinates every 2.5 ps
nstvout                  = 5000      ; save velocities every 2.5 ps
nstenergy                = 5000      ; save energies every 2.5 ps
nstlog                   = 5000      ; update log file every 2.5 ps

; Bond parameters
continuation             = yes       ; continuing from NVT
constraint_algorithm     = lincs     ; holonomic constraints 
constraints              = h-bonds   ; constrain bonds with H atoms
lincs_iter               = 4         ; increased accuracy of LINCS
lincs_order              = 8         ; increased accuracy

; Neighbor searching
cutoff-scheme           = Verlet    ; Verlet cutoff scheme (efficient on modern hardware)
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; update neighbor list every 20 steps
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)

; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT

; Temperature coupling
tcoupl                  = V-rescale ; velocity rescaling thermostat (more accurate than Berendsen)
tc-grps                 = System    ; couple entire system as one group
tau_t                   = 1.0       ; longer time constant for gentler coupling
ref_t                   = 273       ; reference temperature 273K (0°C)

; Pressure coupling
pcoupl                  = Parrinello-Rahman ; Pressure coupling on in NPT
pcoupltype              = isotropic ; uniform scaling of box vectors
tau_p                   = 2.0       ; time constant, in ps
ref_p                   = 1.0       ; reference pressure, in bar
compressibility         = 4.5e-5    ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com       ; Scale COM of reference coordinates

; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC

; Dispersion correction
DispCorr                = EnerPres  ; account for cut-off vdW scheme (improves energy and pressure accuracy)

; Velocity generation
gen_vel                 = no        ; Velocity generation is off
EOF_NPT
  
  # Run grompp with increased warning tolerance
  log "Running grompp for NPT equilibration..."
  gmx grompp -f npt_improved.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 10 || { log "grompp for npt failed"; exit 1; }
  
  # Run mdrun with appropriate parallelization
  log "Running NPT equilibration..."
  gmx mdrun -v -deffnm npt -ntmpi 1 -ntomp 6 || { log "mdrun for npt failed"; exit 1; }
  
  log "NPT equilibration completed."
"

# Step 6: Perform production MD
run_step "production_md" "Run production MD" "
  cd "${DATA_DIR}"
  
  # Create an improved md.mdp file with PME electrostatics
  log "Creating improved production MD configuration file..."
  cat > md_improved.mdp << 'EOF_MD'
; Production MD parameters for TIP4P/Ice water at 273K
integrator               = md        ; leap-frog integrator
nsteps                   = 1000000   ; 1 ns with 1 fs timestep
dt                       = 0.001     ; 1 fs - small timestep for stability
nstxout                  = 5000      ; save coordinates every 5 ps
nstvout                  = 5000      ; save velocities every 5 ps
nstenergy                = 5000      ; save energies every 5 ps
nstlog                   = 5000      ; update log file every 5 ps
nstxout-compressed       = 5000      ; save compressed coordinates every 5 ps
compressed-x-grps        = System    ; save the whole system

; Bond parameters
continuation             = yes       ; continuing from NPT
constraint_algorithm     = lincs     ; holonomic constraints 
constraints              = h-bonds   ; constrain bonds with H atoms
lincs_iter               = 4         ; increased accuracy of LINCS
lincs_order              = 8         ; increased accuracy

; Neighbor searching
cutoff-scheme           = Verlet    ; Verlet cutoff scheme (efficient on modern hardware)
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; update neighbor list every 20 steps
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)

; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT

; Temperature coupling
tcoupl                  = V-rescale ; velocity rescaling thermostat (more accurate than Berendsen)
tc-grps                 = System    ; couple entire system as one group
tau_t                   = 1.0       ; longer time constant for gentler coupling
ref_t                   = 273       ; reference temperature 273K (0°C)

; Pressure coupling
pcoupl                  = Parrinello-Rahman ; Pressure coupling on in NPT
pcoupltype              = isotropic ; uniform scaling of box vectors
tau_p                   = 2.0       ; time constant, in ps
ref_p                   = 1.0       ; reference pressure, in bar
compressibility         = 4.5e-5    ; isothermal compressibility of water, bar^-1

; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC

; Dispersion correction
DispCorr                = EnerPres  ; account for cut-off vdW scheme (improves energy and pressure accuracy)

; Velocity generation
gen_vel                 = no        ; Velocity generation is off
EOF_MD
  
  # Run grompp with increased warning tolerance
  log "Running grompp for production MD..."
  gmx grompp -f md_improved.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 10 || { log "grompp for md failed"; exit 1; }
  
  # Run mdrun with appropriate parallelization
  log "Running production MD..."
  gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 6 || { log "mdrun for md failed"; exit 1; }
  
  log "Production MD completed."
"

# Step 7: Run analysis
log "Step 7: Running analysis..."

# Calculate Radial Distribution Functions (RDFs)
run_step "rdf_analysis_OO" "Calculate Oxygen-Oxygen RDF" "rdf_analysis_OO"
run_step "rdf_analysis_OH" "Calculate Oxygen-Hydrogen RDF" "rdf_analysis_OH"
run_step "rdf_analysis_HH" "Calculate Hydrogen-Hydrogen RDF" "rdf_analysis_HH"
run_step "plot_rdf" "Plot RDF data" "plot_rdf"

# Calculate Mean Square Displacement (MSD) for diffusion coefficient
run_step "msd_analysis" "Calculate Mean Square Displacement for diffusion coefficient" "
  cd "${DATA_DIR}"
  gmx msd -f md.xtc -s md.tpr -o "${ANALYSIS_DIR}/data/msd.xvg" -beginfit 1000 -endfit 2000 -sel "name OW" -dt 50 -trestart 100 || { log "MSD analysis failed"; exit 1; }
"

# Additional analysis
run_step "additional_analysis" "Perform additional analysis" "
  cd "${DATA_DIR}"
  log "Performing additional analysis..."
  
  # Remove existing XPM file if it exists
  if [ -f "${ANALYSIS_DIR}/data/density_map.xpm" ]; then
      rm "${ANALYSIS_DIR}/data/density_map.xpm"
  fi
  
  # Generate density map
  log "Generating density map..."
  echo "1" | gmx densmap -f md.xtc -s md.tpr -o "${ANALYSIS_DIR}/data/density_map.xpm" -bin 0.05 || {
      log "Standard densmap approach failed, trying alternative parameters..."
      echo "SOL" | gmx densmap -f md.xtc -s md.tpr -o "${ANALYSIS_DIR}/data/density_map.xpm" -bin 0.1 -unit nm-3 || {
          log "All densmap approaches failed."
          exit 1
      }
  }
  
  # Verify the XPM file was created
  if [ -f "${ANALYSIS_DIR}/data/density_map.xpm" ]; then
    log "Density map XPM file created successfully."
  else
    log "Error: Failed to create density map XPM file."
    exit 1
  fi
  
  log "Additional data files generated successfully."
"

# Step 12: Hydrogen bond analysis
run_step "hbond_analysis" "Perform hydrogen bond analysis" "
  cd "${DATA_DIR}"
  
  # Debug: Check if required files exist
  log "Checking for required files..."
  ls -l md.xtc md.tpr
  
  # Create a proper index file with correct atom selections for TIP4P/Ice
  log "Creating index file with correct TIP4P/Ice atom selections..."
  
  # First, let's check what atom names are actually in the topology
  echo q | gmx dump -s md.tpr 2>&1 | grep -A 20 "atoms" | head -20
  
  # Create index file with proper TIP4P/Ice groups
  cat > make_hbond_ndx.txt << 'EOF_NDX'
a H*
name 3 WAT_H
a OW
name 4 WAT_O
q
EOF_NDX
  
  gmx make_ndx -f md.tpr -o hbond.ndx < make_hbond_ndx.txt || { 
    log "Failed to create index file"; 
    exit 1; 
  }
  
  # Debug: Show index file contents
  log "Index file contents:"
  cat hbond.ndx
  
  # Create input file for hbond analysis
  cat > hbond.input << 'EOF_HBOND'
3
4
EOF_HBOND
  
  # Run hbond analysis with the input file
  log "Running hydrogen bond analysis for TIP4P/Ice..."
  
  # Run the analysis
  gmx hbond -f md.xtc -s md.tpr -n hbond.ndx \
       -num "${ANALYSIS_DIR}/data/hbnum.xvg" \
       -dist "${ANALYSIS_DIR}/data/hbdist.xvg" \
       -ang "${ANALYSIS_DIR}/data/hbang.xvg" < hbond.input || {
    
    log "Standard hbond analysis failed, trying with just number of hydrogen bonds..."
    
    gmx hbond -f md.xtc -s md.tpr -n hbond.ndx \
         -num "${ANALYSIS_DIR}/data/hbnum.xvg" < hbond.input || {
      
      log "All hydrogen bond analysis approaches failed."
      exit 1
    }
  }
  
  log "Hydrogen bond analysis completed."
"

# Step 13: Velocity autocorrelation function (VACF) analysis
run_step "vacf_analysis" "Perform velocity autocorrelation function analysis" "
  cd "${DATA_DIR}"
  
  # Check if md.trr exists (velocities are stored in .trr files)
  if [ -f "md.trr" ]; then
    # Create index file if it doesn't exist
    if [ ! -f "index.ndx" ]; then
      echo "q" | gmx make_ndx -f md.tpr -o index.ndx || { log "Failed to create index file for VACF analysis"; exit 1; }
    fi
    
    # Run VACF analysis using the correct command (velacc)
    echo "SOL" | gmx velacc -f md.trr -s md.tpr -n index.ndx -o "${ANALYSIS_DIR}/data/vacf.xvg" -os "${ANALYSIS_DIR}/data/vacf_spectrum.xvg" -acflen 1000 -nonormalize || { log "VACF analysis failed"; exit 1; }
    log "VACF analysis completed."
  else
    log "Warning: md.trr file not found. VACF analysis requires velocities which are stored in .trr files."
    
    # Generate velocities using a more robust approach
    log "Generating velocities from positions..."
    gmx grompp -f "${CONFIGS_DIR}/md.mdp" -c md.gro -p topol.top -o md_with_vel.tpr -maxwarn 2 || { log "TPR generation failed"; exit 1; }
    
    # Then run a short simulation to generate velocities
    gmx mdrun -s md_with_vel.tpr -rerun md.xtc -o md_with_vel.trr -v || { log "TRR generation failed"; exit 1; }
    
    # Now try VACF analysis with the generated TRR file
    if [ -f "md_with_vel.trr" ]; then
      # Create index file for VACF analysis
      echo "q" | gmx make_ndx -f md_with_vel.tpr -o index.ndx || { log "Failed to create index file for VACF analysis"; exit 1; }
      
      # Run VACF analysis with the index file
      echo "SOL" | gmx velacc -f md_with_vel.trr -s md_with_vel.tpr -n index.ndx -o "${ANALYSIS_DIR}/data/vacf.xvg" -os "${ANALYSIS_DIR}/data/vacf_spectrum.xvg" -acflen 1000 -nonormalize || { log "VACF analysis failed"; exit 1; }
      log "VACF analysis completed after generating velocities."
    else
      log "Failed to generate velocity data. VACF analysis cannot be performed."
      exit 1
    fi
  fi
"

# Step 14: Generate plots
run_step "generate_plots" "Generate plots" "
  cd "${ANALYSIS_DIR}"
  log "Generating plots using Python scripts..."
  
  # Use the Python plotting scripts instead of gnuplot
  python3 "${ANALYSIS_DIR}/plotting_scripts/run_all_plots.py" --analysis-dir "${ANALYSIS_DIR}" --plots-dir "${ANALYSIS_DIR}/plots" --verbose
  
  log "Plots generated successfully."
"

# Clean up temporary files
cleanup

log "Workflow completed successfully at $(date)"

# List all completed steps
list_checkpoints

log "======================================"
log "Workflow completed at $(date)!"
log "Results are in ${DATA_DIR} and ${ANALYSIS_DIR}"
log "Plots and summary report are in ${ANALYSIS_DIR}/plots"
log "Summary report is available at ${ANALYSIS_DIR}/plots/tip4pice_water_analysis_summary.png"
log "Text summary is available at ${ANALYSIS_DIR}/plots/tip4pice_water_analysis_summary.txt"
log "Log file saved to ${LOG_FILE}"

# Print a helpful message about how to rerun specific steps
log ""
log "To rerun a specific step, use: ./run_workflow.sh --rerun STEP_NAME"
log "To run only a specific step, use: ./run_workflow.sh --only STEP_NAME"
log ""
log "For example, to regenerate only the plots:"
log "./run_workflow.sh --only generate_plots"

# Add cleanup to the end of the script
cleanup

# RDF Analysis - Oxygen-Oxygen
rdf_analysis_OO() {
    echo "Running RDF analysis for Oxygen-Oxygen..."
    
    # Create index file for atom selections
    echo -e "a OW\nname 1 OW\na HW*\nname 2 HW\nq" | gmx make_ndx -f "${DATA_DIR}/md.tpr" -o "${DATA_DIR}/index.ndx"
    
    # Run RDF analysis with exclusion of intramolecular contributions
    echo -e "3\n3" | gmx rdf -f "${DATA_DIR}/md.xtc" -s "${DATA_DIR}/md.tpr" -n "${DATA_DIR}/index.ndx" -o "${ANALYSIS_DIR}/data/rdf_OO.xvg" -excl
    
    touch "${CHECKPOINT_DIR}/rdf_analysis_OO.done"
}

# RDF Analysis - Oxygen-Hydrogen
rdf_analysis_OH() {
    echo "Running RDF analysis for Oxygen-Hydrogen..."
    
    # Ensure index file exists
    if [ ! -f "${DATA_DIR}/index.ndx" ]; then
        echo -e "a OW\nname 1 OW\na HW*\nname 2 HW\nq" | gmx make_ndx -f "${DATA_DIR}/md.tpr" -o "${DATA_DIR}/index.ndx"
    fi
    
    # Run RDF analysis with exclusion of intramolecular contributions
    echo -e "3\n4" | gmx rdf -f "${DATA_DIR}/md.xtc" -s "${DATA_DIR}/md.tpr" -n "${DATA_DIR}/index.ndx" -o "${ANALYSIS_DIR}/data/rdf_OH.xvg" -excl
    
    touch "${CHECKPOINT_DIR}/rdf_analysis_OH.done"
}

# RDF Analysis - Hydrogen-Hydrogen
rdf_analysis_HH() {
    echo "Running RDF analysis for Hydrogen-Hydrogen..."
    
    # Ensure index file exists
    if [ ! -f "${DATA_DIR}/index.ndx" ]; then
        echo -e "a OW\nname 1 OW\na HW*\nname 2 HW\nq" | gmx make_ndx -f "${DATA_DIR}/md.tpr" -o "${DATA_DIR}/index.ndx"
    fi
    
    # Run RDF analysis with exclusion of intramolecular contributions
    echo -e "4\n4" | gmx rdf -f "${DATA_DIR}/md.xtc" -s "${DATA_DIR}/md.tpr" -n "${DATA_DIR}/index.ndx" -o "${ANALYSIS_DIR}/data/rdf_HH.xvg" -excl
    
    touch "${CHECKPOINT_DIR}/rdf_analysis_HH.done"
}

# Plot RDF data
plot_rdf() {
    echo "Plotting RDF data..."
    
    # Use env -i to clear environment variables that might interfere with Python
    env -i PATH=$PATH HOME=$HOME python3 "${ANALYSIS_DIR}/plotting_scripts/plot_rdf.py" "${ANALYSIS_DIR}" "${ANALYSIS_DIR}/plots"
    
    touch "${CHECKPOINT_DIR}/plot_rdf.done"
}
