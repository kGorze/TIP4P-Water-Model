#!/bin/bash
# Script to generate plots and summary report for md_water_study_iteration_4
# This script can be run independently of the main workflow to regenerate plots

# Define directories
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ITERATION_DIR="${SCRIPT_DIR}"
ANALYSIS_DIR="${ITERATION_DIR}/analysis"
DATA_DIR="${ANALYSIS_DIR}/data"
PLOTS_DIR="${ANALYSIS_DIR}/plots"
LOGS_DIR="${ITERATION_DIR}/logs"

# Create log file with timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOGS_DIR}/plotting_${TIMESTAMP}.log"

# Function to log messages to both console and log file
log() {
    echo "$@" | tee -a "${LOG_FILE}"
}

# Redirect all output to log file while also displaying on console
exec > >(tee -a "${LOG_FILE}") 2>&1

log "Starting plot generation at $(date)"
log "======================================"

# Create directories if they don't exist
log "Creating necessary directories..."
mkdir -p "${DATA_DIR}"
mkdir -p "${PLOTS_DIR}"
mkdir -p "${LOGS_DIR}"

# Check if the plotting scripts directory exists
if [ ! -d "${ANALYSIS_DIR}/plotting_scripts" ]; then
    log "Error: Plotting scripts directory not found at ${ANALYSIS_DIR}/plotting_scripts"
    log "Please make sure the plotting scripts are in the correct location."
    exit 1
fi

# Check if data files exist
if [ ! "$(ls -A ${DATA_DIR}/*.xvg 2>/dev/null)" ]; then
    log "Warning: No data files found in ${DATA_DIR}"
    log "Please make sure the data files are in the correct location."
fi

# Run all plotting scripts using the central coordinator script
log "Running all plotting scripts..."
python3 "${ANALYSIS_DIR}/plotting_scripts/run_all_plots.py" --analysis-dir "${ANALYSIS_DIR}" --data-dir "${DATA_DIR}" --plots-dir "${PLOTS_DIR}" --verbose

# Check if the plots were generated successfully
if [ $? -eq 0 ]; then
    log "Plots and summary report generated successfully!"
    log "All plots are available in ${PLOTS_DIR}"
    log "Summary report is available at ${PLOTS_DIR}/tip4pice_water_analysis_summary.png"
    log "Text summary is available at ${PLOTS_DIR}/tip4pice_water_analysis_summary.txt"
else
    log "Error: Failed to generate plots and summary report."
    log "Please check the log file for details."
    exit 1
fi

log "======================================"
log "Plot generation finished at $(date)!"
log "All plots are in ${PLOTS_DIR}"
log "Log file saved to ${LOG_FILE}" 