#!/bin/bash
# Script to fix the failed MSD analysis step

set -e  # Exit on error

# Define directories
ITERATION_DIR="/home/konrad_guest/Documents/research/cursor/md_water_study_iteration_4"
DATA_DIR="${ITERATION_DIR}/data"
ANALYSIS_DIR="${ITERATION_DIR}/analysis"

# Make sure analysis directory exists
mkdir -p "${ANALYSIS_DIR}"

echo "Running fixed MSD analysis..."
cd "${DATA_DIR}"

# Calculate Mean Square Displacement (MSD) for diffusion coefficient with proper selection
echo "Calculating Mean Square Displacement for diffusion coefficient..."
# We select "name OW" to track oxygen atoms of water molecules for diffusion calculation
gmx msd -f md.xtc -s md.tpr -o "${ANALYSIS_DIR}/msd.xvg" -beginfit 1000 -endfit 2000 -sel "name OW" || { echo "MSD analysis failed"; exit 1; }

echo "MSD analysis completed successfully!" 