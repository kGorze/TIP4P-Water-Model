#!/bin/bash
# Simple script to run just the plotting part of the workflow

# Get the directory of the script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PLOTS_DIR="${SCRIPT_DIR}/analysis/plots"

# Run the fixed plotting script
"${SCRIPT_DIR}/run_plots_fixed.sh"

echo "Plotting completed! Results are in ${PLOTS_DIR} directory."
echo "Summary report is available at ${PLOTS_DIR}/tip4pice_water_analysis_summary.png"
echo "Text summary is available at ${PLOTS_DIR}/tip4pice_water_analysis_summary.txt" 