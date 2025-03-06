#!/bin/bash
# Script to run the plotting scripts with the correct Python environment

# Get the directory of the script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ANALYSIS_DIR="${SCRIPT_DIR}/analysis"
PLOTS_DIR="${ANALYSIS_DIR}/plots"
DATA_DIR="${ANALYSIS_DIR}/data"

# Set up environment variables
export PYTHONPATH=""
unset PYTHONHOME

# Create plots directory if it doesn't exist
mkdir -p "${PLOTS_DIR}"

# Run each plotting script individually
echo "Running RDF plotting script..."
python3 "${ANALYSIS_DIR}/plotting_scripts/plot_rdf.py" "${ANALYSIS_DIR}" "${PLOTS_DIR}"

echo "Running density plotting script..."
python3 "${ANALYSIS_DIR}/plotting_scripts/plot_density.py" "${ANALYSIS_DIR}" "${PLOTS_DIR}"

echo "Running MSD plotting script..."
python3 "${ANALYSIS_DIR}/plotting_scripts/plot_msd.py" "${ANALYSIS_DIR}" "${PLOTS_DIR}"

echo "Running hydrogen bond plotting script..."
python3 "${ANALYSIS_DIR}/plotting_scripts/plot_hbond.py" "${ANALYSIS_DIR}" "${PLOTS_DIR}"

echo "Running temperature plotting script..."
python3 "${ANALYSIS_DIR}/plotting_scripts/plot_temperature.py" "${ANALYSIS_DIR}" "${PLOTS_DIR}"

echo "Running energy analysis plotting script..."
python3 "${ANALYSIS_DIR}/plotting_scripts/plot_energy_analysis.py" "${ANALYSIS_DIR}" "${PLOTS_DIR}"

# Add RMSD plotting
echo "Running RMSD plotting script..."
python3 "${ANALYSIS_DIR}/plotting_scripts/plot_rmsd.py" "${ANALYSIS_DIR}" "${PLOTS_DIR}"

# Check if VACF data exists
if [ -f "${DATA_DIR}/vacf.xvg" ]; then
    echo "Running VACF plotting script..."
    python3 "${ANALYSIS_DIR}/plotting_scripts/plot_vacf.py" "${ANALYSIS_DIR}" "${PLOTS_DIR}"
else
    echo "VACF data not found, skipping VACF plots"
fi

# Generate summary report
echo "Generating summary report..."
python3 "${ANALYSIS_DIR}/plotting_scripts/generate_summary_report.py" "${ANALYSIS_DIR}" "${PLOTS_DIR}"

echo "Plotting completed! Results are in ${PLOTS_DIR} directory."
echo "Summary report is available at ${PLOTS_DIR}/tip4pice_water_analysis_summary.png"
echo "Text summary is available at ${PLOTS_DIR}/tip4pice_water_analysis_summary.txt" 