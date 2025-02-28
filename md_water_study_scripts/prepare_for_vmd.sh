#!/bin/bash

# Script to prepare GROMACS files for VMD visualization
# This script converts GROMACS files to formats that VMD can read more easily

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$( dirname "$SCRIPT_DIR" )"

# Check if GROMACS is installed
if ! command -v gmx &> /dev/null; then
    echo "Error: GROMACS is not installed or not in PATH."
    exit 1
fi

# Set paths
MODEL="tip4p"
TEMP="273"
DATA_DIR="$BASE_DIR/data/$MODEL/${TEMP}K"
VMD_DIR="$DATA_DIR/vmd"

# Create VMD directory if it doesn't exist
mkdir -p "$VMD_DIR"

echo "==================================================================="
echo "          Preparing GROMACS files for VMD visualization"
echo "==================================================================="

# Step 1: Convert the final structure to PDB format
echo "Converting final structure to PDB format..."
gmx editconf -f "$DATA_DIR/md.gro" -o "$VMD_DIR/md_final.pdb"

# Step 2: Convert the trajectory to a VMD-friendly format
echo "Converting trajectory to VMD-friendly format..."
echo "0" | gmx trjconv -s "$DATA_DIR/md.tpr" -f "$DATA_DIR/md.xtc" -o "$VMD_DIR/md_traj.pdb" -pbc mol -dt 100

echo "Files prepared for VMD:"
echo "- Structure: $VMD_DIR/md_final.pdb"
echo "- Trajectory: $VMD_DIR/md_traj.pdb"
echo ""
echo "To view in VMD, run:"
echo "vmd $VMD_DIR/md_final.pdb $VMD_DIR/md_traj.pdb"
echo ""
echo "Would you like to launch VMD now? (y/n)"
read -r response
if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
    vmd "$VMD_DIR/md_final.pdb" "$VMD_DIR/md_traj.pdb"
fi 