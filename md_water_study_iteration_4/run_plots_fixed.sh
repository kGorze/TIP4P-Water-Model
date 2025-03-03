#!/bin/bash
# Script to run the plotting scripts with the correct Python environment

# Set up environment variables
export PYTHONPATH=""
unset PYTHONHOME

# Create plots directory if it doesn't exist
mkdir -p analysis/plots

# Run each plotting script individually
echo "Running RDF plotting script..."
python3 analysis/plotting_scripts/plot_rdf.py analysis analysis/plots

echo "Running density plotting script..."
python3 analysis/plotting_scripts/plot_density.py analysis analysis/plots

echo "Running MSD plotting script..."
python3 analysis/plotting_scripts/plot_msd.py analysis analysis/plots

echo "Running hydrogen bond plotting script..."
python3 analysis/plotting_scripts/plot_hbond.py analysis analysis/plots

echo "Running temperature plotting script..."
python3 analysis/plotting_scripts/plot_temperature.py analysis analysis/plots

echo "Running energy analysis plotting script..."
python3 analysis/plotting_scripts/plot_energy_analysis.py analysis analysis/plots

# Add RMSD plotting
echo "Running RMSD plotting script..."
python3 analysis/plotting_scripts/plot_rmsd.py analysis analysis/plots

# Check if VACF data exists
if [ -f "analysis/data/vacf.xvg" ]; then
    echo "Running VACF plotting script..."
    python3 analysis/plotting_scripts/plot_vacf.py analysis analysis/plots
else
    echo "VACF data not found, skipping VACF plots"
fi

# Generate summary report
echo "Generating summary report..."
python3 analysis/plotting_scripts/generate_summary_report.py analysis analysis/plots

echo "Plotting completed! Results are in analysis/plots directory."
echo "Summary report is available at analysis/plots/tip4p_water_analysis_summary.png"
echo "Text summary is available at analysis/plots/tip4p_water_analysis_summary.txt" 