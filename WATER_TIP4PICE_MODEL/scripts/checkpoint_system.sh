#!/bin/bash
# Checkpoint system for MD water study workflow
# This script provides functions to manage checkpoints in the MD workflow

# Define the checkpoint directory
CHECKPOINT_DIR="${ITERATION_DIR}/checkpoints"
mkdir -p "${CHECKPOINT_DIR}"

# Function to create a checkpoint
# Usage: create_checkpoint "step_name"
create_checkpoint() {
    local step_name="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S')" > "${CHECKPOINT_DIR}/${step_name}.done"
    echo "Checkpoint created: ${step_name}"
}

# Function to check if a checkpoint exists
# Usage: if check_checkpoint "step_name"; then echo "Step already done"; fi
check_checkpoint() {
    local step_name="$1"
    if [ -f "${CHECKPOINT_DIR}/${step_name}.done" ]; then
        return 0  # Checkpoint exists
    else
        return 1  # Checkpoint doesn't exist
    fi
}

# Function to remove a checkpoint (e.g., if you want to force a step to run again)
# Usage: remove_checkpoint "step_name"
remove_checkpoint() {
    local step_name="$1"
    if [ -f "${CHECKPOINT_DIR}/${step_name}.done" ]; then
        rm "${CHECKPOINT_DIR}/${step_name}.done"
        echo "Checkpoint removed: ${step_name}"
    else
        echo "No checkpoint found for: ${step_name}"
    fi
}

# Function to list all checkpoints
# Usage: list_checkpoints
list_checkpoints() {
    echo "Completed steps:"
    for checkpoint in "${CHECKPOINT_DIR}"/*.done; do
        if [ -f "$checkpoint" ]; then
            step_name=$(basename "$checkpoint" .done)
            completion_time=$(cat "$checkpoint")
            echo "  - ${step_name} (completed at ${completion_time})"
        fi
    done
}

# Function to run a step with checkpoint checking
# Usage: run_step "step_name" "step_description" command_to_run
run_step() {
    local step_name="$1"
    local step_description="$2"
    local command="$3"
    
    if check_checkpoint "${step_name}"; then
        echo "Step '${step_description}' already completed. Skipping..."
    else
        echo "Running step: ${step_description}..."
        eval "${command}"
        if [ $? -eq 0 ]; then
            create_checkpoint "${step_name}"
            echo "Step '${step_description}' completed successfully."
        else
            echo "Step '${step_description}' failed. No checkpoint created."
            return 1
        fi
    fi
    return 0
}

# Example usage in workflow script:
# 
# source "${SCRIPTS_DIR}/checkpoint_system.sh"
# 
# # Generate water box
# run_step "water_box" "Generate water box" "
#   cd \"${DATA_DIR}\"
#   julia -e 'using Packmol; run_packmol(\"water_box.inp\")'
# "
# 
# # Generate topology
# run_step "topology" "Generate topology" "
#   cd \"${DATA_DIR}\"
#   echo \"15\" | gmx pdb2gmx -f water_box.pdb -o water_box.gro -water tip4p -p topol.top -ff oplsaa
# "
#
# # And so on for other steps...
#
# # Run analysis steps
# run_step "rdf_analysis" "Calculate RDFs" "
#   cd \"${DATA_DIR}\"
#   gmx rdf -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/rdf_OO.xvg\" -ref \"name OW\" -sel \"name OW\"
# "
#
# run_step "msd_analysis" "Calculate MSD" "
#   cd \"${DATA_DIR}\"
#   gmx msd -f md.xtc -s md.tpr -o \"${ANALYSIS_DIR}/msd.xvg\" -beginfit 1000 -endfit 2000 -sel \"name OW\"
# " 