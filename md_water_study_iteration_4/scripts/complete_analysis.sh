#!/bin/bash
# Script to complete all analysis for md_water_study_iteration_4
# This script will ensure all the analyses from iteration 3 are also performed in iteration 4
# Now with checkpoints to avoid repeating successful steps

set -e  # Exit on error

# Define directories
ITERATION_DIR="/home/konrad_guest/Documents/research/cursor/md_water_study_iteration_4"
DATA_DIR="${ITERATION_DIR}/data"
ANALYSIS_DIR="${ITERATION_DIR}/analysis"
PLOTS_DIR="${ANALYSIS_DIR}/plots"
CHECKPOINT_DIR="${ITERATION_DIR}/checkpoints"

# Make sure directories exist
mkdir -p "${ANALYSIS_DIR}"
mkdir -p "${PLOTS_DIR}"
mkdir -p "${CHECKPOINT_DIR}"

# Checkpoint functions
# Function to create a checkpoint
create_checkpoint() {
    local step_name="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S')" > "${CHECKPOINT_DIR}/${step_name}.done"
    echo "✓ Checkpoint created: ${step_name}"
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

# Function to run a step with checkpoint checking
run_step() {
    local step_name="$1"
    local step_description="$2"
    local command="$3"
    
    if check_checkpoint "${step_name}"; then
        echo "✓ Step '${step_description}' already completed. Skipping..."
    else
        echo "▶ Running step: ${step_description}..."
        if eval "${command}"; then
            create_checkpoint "${step_name}"
            echo "✓ Step '${step_description}' completed successfully."
        else
            echo "✗ Step '${step_description}' failed. No checkpoint created."
            return 1
        fi
    fi
    return 0
}

# List all completed steps
list_checkpoints() {
    echo "=== Completed analysis steps ==="
    for checkpoint in "${CHECKPOINT_DIR}"/*.done; do
        if [ -f "$checkpoint" ]; then
            step_name=$(basename "$checkpoint" .done)
            completion_time=$(cat "$checkpoint")
            echo "  ✓ ${step_name} (completed at ${completion_time})"
        fi
    done
    echo "=============================="
}

echo "Running complete analysis for TIP4P water model with checkpoints..."
cd "${DATA_DIR}"

# 1. RDF Analyses
run_step "rdf_OO" "Calculate Oxygen-Oxygen RDF" "
    if [ ! -f \"${ANALYSIS_DIR}/rdf_OO.xvg\" ]; then
        gmx rdf -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/rdf_OO.xvg\" -ref \"name OW\" -sel \"name OW\" -seltype atom -selrpos atom
    else
        echo \"O-O RDF file already exists, skipping calculation.\"
    fi
"

run_step "rdf_OH" "Calculate Oxygen-Hydrogen RDF" "
    gmx rdf -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/rdf_OH.xvg\" -ref \"name OW\" -sel \"name HW1 or name HW2\" -seltype atom -selrpos atom -excl
"

run_step "rdf_HH" "Calculate Hydrogen-Hydrogen RDF" "
    gmx rdf -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/rdf_HH.xvg\" -ref \"name HW1 or name HW2\" -sel \"name HW1 or name HW2\" -seltype atom -selrpos atom -excl
"

# 2. MSD and Diffusion Coefficient
run_step "msd_analysis" "Calculate Mean Square Displacement" "
    if [ ! -f \"${ANALYSIS_DIR}/msd.xvg\" ]; then
        gmx msd -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/msd.xvg\" -beginfit 1000 -endfit 2000 -sel \"name OW\"
    else
        echo \"MSD file already exists, skipping calculation.\"
    fi
"

# 3. Hydrogen Bond Analysis
# Fix for gmx hbond - check for available options in GROMACS 2024.3
# The -life parameter seems to be unsupported, let's split this into separate commands
run_step "hbond_number" "Calculate number of hydrogen bonds" "
    # For TIP4P water, use 'SOL' as the selection group, with selections on separate lines
    echo -e \"SOL\nSOL\" | gmx hbond -f md.xtc -s md.tpr -num \"${ANALYSIS_DIR}/hbnum.xvg\"
"

run_step "hbond_dist" "Calculate hydrogen bond distances" "
    # For TIP4P water, use 'SOL' as the selection group, with selections on separate lines
    echo -e \"SOL\nSOL\" | gmx hbond -f md.xtc -s md.tpr -dist \"${ANALYSIS_DIR}/hbdist.xvg\"
"

run_step "hbond_ang" "Calculate hydrogen bond angles" "
    # For TIP4P water, use 'SOL' as the selection group, with selections on separate lines
    echo -e \"SOL\nSOL\" | gmx hbond -f md.xtc -s md.tpr -ang \"${ANALYSIS_DIR}/hbang.xvg\"
"

# Add hydrogen bond lifetime analysis using the legacy command
run_step "hbond_lifetime" "Calculate hydrogen bond lifetimes" \
  "cd ${DATA_DIR} && \
  echo -e \"SOL\nSOL\" | gmx hbond-legacy -f md.xtc -s md.tpr -b 0 -e 100 -life \"${ANALYSIS_DIR}/hblife.xvg\" -ac \"${ANALYSIS_DIR}/hbac.xvg\" || \
  echo \"Hydrogen bond lifetime analysis failed. This may be due to box shrinkage or other issues.\""

# 4. Density Profiles
run_step "density_profile" "Calculate density profile" "
    if [ ! -f \"${ANALYSIS_DIR}/density.xvg\" ]; then
        gmx density -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/density.xvg\" -dens mass -d Z
    else
        echo \"Density profile file already exists, skipping calculation.\"
    fi
"

run_step "density_radial" "Calculate radial density profile" "
    echo \"System\" | gmx densmap -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/density_radial.xpm\" -bin 0.05
"

# 5. VACF Analysis with fallback for missing md.trr
run_step "vacf_analysis" "Perform velocity autocorrelation function analysis" "
    if [ -f \"md.trr\" ]; then
        # Create an index file using here-document to feed commands to the interactive tool
        # 'q' means 'save and quit'
        gmx make_ndx -f md.tpr -o vacf.ndx << EOF
q
EOF
        
        # Check if the index file was created successfully
        if [ ! -f \"vacf.ndx\" ]; then
            echo \"Failed to create index file. Aborting VACF analysis.\"
            exit 1
        fi
        
        # Display the index file contents for debugging
        echo \"Index file created with the following groups:\"
        grep -A 10 \"\\[ \" vacf.ndx | head -20
        
        # Locate SOL or Water group number from the index file
        SOL_GROUP=\$(grep -A 1 \"\\[ SOL \\]\" vacf.ndx | tail -1 | grep -o \"[0-9]*\" | head -1)
        if [ -z \"\$SOL_GROUP\" ]; then
            # If SOL not found, try Water group
            SOL_GROUP=\$(grep -A 1 \"\\[ Water \\]\" vacf.ndx | tail -1 | grep -o \"[0-9]*\" | head -1)
            if [ -z \"\$SOL_GROUP\" ]; then
                # If Water not found either, default to group 1
                SOL_GROUP=1
                echo \"Neither SOL nor Water group found in index, defaulting to group \$SOL_GROUP\"
            else
                echo \"Found Water as group \$SOL_GROUP\"
            fi
        else
            echo \"Found SOL as group \$SOL_GROUP\"
        fi
        
        # Run velacc with the index file, selecting the appropriate group number
        # IMPORTANT: We're using a here-document to feed the group number to the interactive command
        echo \"Running gmx velacc with group \$SOL_GROUP from index file\"
        gmx velacc -f md.trr -s md.tpr -n vacf.ndx -o \"${ANALYSIS_DIR}/vacf.xvg\" -os \"${ANALYSIS_DIR}/vacf_spectrum.xvg\" -acflen 1000 -nonormalize << EOF
\$SOL_GROUP
EOF
    elif [ -f \"npt.gro\" ] && [ -f \"npt.cpt\" ]; then
        # Create a modified mdp file for short production with velocity output
        cp \"${CONFIGS_DIR}/md.mdp\" ./md_velout.mdp
        # Update the mdp file to save velocities
        sed -i 's/nstvout\s*=\s*0/nstvout = 100/g' md_velout.mdp
        if ! grep -q \"nstvout\" md_velout.mdp; then
            echo \"nstvout = 100\" >> md_velout.mdp
        fi
        
        # Run a short production (200 ps) with velocity output
        echo \"Running short production MD with velocity output...\"
        gmx grompp -f md_velout.mdp -c npt.gro -t npt.cpt -p topol.top -o md_velout.tpr -maxwarn 1
        gmx mdrun -deffnm md_velout -v -nsteps 100000  # 200 ps with 2 fs timestep
        
        # Create an index file with proper here-document
        gmx make_ndx -f md_velout.tpr -o vacf.ndx << EOF
q
EOF
        
        # Check if the index file was created successfully
        if [ ! -f \"vacf.ndx\" ]; then
            echo \"Failed to create index file. Aborting VACF analysis.\"
            exit 1
        fi
        
        # Display the index file contents for debugging
        echo \"Index file created with the following groups:\"
        grep -A 10 \"\\[ \" vacf.ndx | head -20
        
        # Locate SOL or Water group number from the index file
        SOL_GROUP=\$(grep -A 1 \"\\[ SOL \\]\" vacf.ndx | tail -1 | grep -o \"[0-9]*\" | head -1)
        if [ -z \"\$SOL_GROUP\" ]; then
            # If SOL not found, try Water group
            SOL_GROUP=\$(grep -A 1 \"\\[ Water \\]\" vacf.ndx | tail -1 | grep -o \"[0-9]*\" | head -1)
            if [ -z \"\$SOL_GROUP\" ]; then
                # If Water not found either, default to group 1
                SOL_GROUP=1
                echo \"Neither SOL nor Water group found in index, defaulting to group \$SOL_GROUP\"
            else
                echo \"Found Water as group \$SOL_GROUP\"
            fi
        else
            echo \"Found SOL as group \$SOL_GROUP\"
        fi
        
        # Run VACF analysis with the index file, selecting the appropriate group
        echo \"Running gmx velacc with group \$SOL_GROUP from index file\"
        gmx velacc -f md_velout.trr -s md_velout.tpr -n vacf.ndx -o \"${ANALYSIS_DIR}/vacf.xvg\" -os \"${ANALYSIS_DIR}/vacf_spectrum.xvg\" -acflen 500 -nonormalize << EOF
\$SOL_GROUP
EOF
    else
        echo \"Cannot perform VACF analysis - required files not found.\"
        exit 1
    fi
"

# 6. Extract thermodynamic properties
run_step "temperature_analysis" "Extract temperature data" "
    echo \"Temperature\" | gmx energy -f md.edr -o \"${ANALYSIS_DIR}/temperature.xvg\"
"

run_step "pressure_analysis" "Extract pressure data" "
    echo \"Pressure\" | gmx energy -f md.edr -o \"${ANALYSIS_DIR}/pressure.xvg\"
"

run_step "energy_analysis" "Extract energy components" "
    echo \"Potential Kinetic-En. Total-Energy\" | gmx energy -f md.edr -o \"${ANALYSIS_DIR}/energy.xvg\"
"

run_step "potential_energy" "Extract potential energy" "
    echo \"Potential\" | gmx energy -f md.edr -o \"${ANALYSIS_DIR}/potential.xvg\"
"

run_step "energy_terms" "Extract detailed energy terms" "
    echo \"Angle Proper-Dih. LJ-14 Coulomb-14 LJ-(SR) Coulomb-(SR) Coul.-recip.\" | gmx energy -f md.edr -o \"${ANALYSIS_DIR}/energy_terms.xvg\"
"

# 7. RMSD calculation
run_step "rmsd_analysis" "Calculate RMSD" "
    # Provide both inputs: one for least-squares fit and one for RMSD calculation
    echo -e \"System\nSystem\" | gmx rms -s md.tpr -f md.xtc -o \"${ANALYSIS_DIR}/rmsd.xvg\" -tu ns
"

# 8. Generate plots (only if all analyses are done)
run_step "generate_plots" "Generate plots for all analyses" "
    # Ensure plotting scripts directory exists
    mkdir -p \"${ANALYSIS_DIR}/plotting_scripts\"
    
    # Make sure all plotting scripts are executable
    chmod +x \"${ANALYSIS_DIR}/plotting_scripts\"/*.py
    
    # Fix potential Python environment issues
    export PYTHONHOME=/usr
    export PYTHONPATH=/usr/lib/python3/dist-packages:/usr/lib/python3
    
    # Run each plotting script individually
    echo \"Running RDF plotting script...\"
    /usr/bin/python3 \"${ANALYSIS_DIR}/plotting_scripts/plot_rdf.py\" \"${ANALYSIS_DIR}\" \"${PLOTS_DIR}\"
    
    echo \"Running thermodynamic properties plotting script...\"
    /usr/bin/python3 \"${ANALYSIS_DIR}/plotting_scripts/plot_temperature.py\" \"${ANALYSIS_DIR}\" \"${PLOTS_DIR}\"
    
    echo \"Running RMSD plotting script...\"
    /usr/bin/python3 \"${ANALYSIS_DIR}/plotting_scripts/plot_rmsd.py\" \"${ANALYSIS_DIR}\" \"${PLOTS_DIR}\"
    
    echo \"Running MSD plotting script...\"
    /usr/bin/python3 \"${ANALYSIS_DIR}/plotting_scripts/plot_msd.py\" \"${ANALYSIS_DIR}\" \"${PLOTS_DIR}\"
    
    echo \"Running hydrogen bond plotting script...\"
    /usr/bin/python3 \"${ANALYSIS_DIR}/plotting_scripts/plot_hbond.py\" \"${ANALYSIS_DIR}\" \"${PLOTS_DIR}\"
    
    echo \"Running density plotting script...\"
    /usr/bin/python3 \"${ANALYSIS_DIR}/plotting_scripts/plot_density.py\" \"${ANALYSIS_DIR}\" \"${PLOTS_DIR}\"
    
    echo \"Running VACF plotting script...\"
    /usr/bin/python3 \"${ANALYSIS_DIR}/plotting_scripts/plot_vacf.py\" \"${ANALYSIS_DIR}\" \"${PLOTS_DIR}\"
    
    echo \"Generating summary report...\"
    /usr/bin/python3 \"${ANALYSIS_DIR}/plotting_scripts/generate_summary_report.py\" \"${ANALYSIS_DIR}\" \"${PLOTS_DIR}\"
    
    echo \"All plots and summary report completed.\"
"

# List all completed steps
list_checkpoints

echo "Analysis complete! All results are in ${ANALYSIS_DIR} with plots in ${PLOTS_DIR}" 