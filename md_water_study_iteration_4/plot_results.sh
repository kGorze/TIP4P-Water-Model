#!/bin/bash
# Simple script to run just the plotting part of the workflow

# Run the fixed plotting script
./run_plots_fixed.sh

echo "Plotting completed! Results are in analysis/plots directory."
echo "Summary report is available at analysis/plots/tip4p_water_analysis_summary.png"
echo "Text summary is available at analysis/plots/tip4p_water_analysis_summary.txt" 