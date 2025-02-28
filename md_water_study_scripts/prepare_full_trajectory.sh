#!/bin/bash

# Script to prepare a full trajectory for VMD visualization
# This script creates a trajectory with ALL frames (dt=1) for the smoothest possible visualization

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
VMD_DIR="$DATA_DIR/vmd_full"

# Create VMD directory if it doesn't exist
mkdir -p "$VMD_DIR"

echo "==================================================================="
echo "      Preparing full trajectory for VMD visualization"
echo "==================================================================="
echo "This will create a trajectory with ALL frames for maximum smoothness."
echo "Note: This may create a large file and take longer to process."
echo ""

# Step 1: Convert the final structure to PDB format
echo "Converting final structure to PDB format..."
gmx editconf -f "$DATA_DIR/md.gro" -o "$VMD_DIR/md_final.pdb"

# Step 2: First, create a centered trajectory with whole molecules
echo "Creating centered trajectory with whole molecules..."
# Use System group for centering and keep molecules whole
echo "0 0" | gmx trjconv -s "$DATA_DIR/md.tpr" -f "$DATA_DIR/md.xtc" -o "$VMD_DIR/md_centered.xtc" -pbc whole -center

# Step 3: Convert the centered trajectory to a VMD-friendly format with ALL frames
echo "Converting to VMD-friendly format with ALL frames (dt=1)..."
# Use dt=1 to keep ALL frames from the original trajectory
echo "0" | gmx trjconv -s "$DATA_DIR/md.tpr" -f "$VMD_DIR/md_centered.xtc" -o "$VMD_DIR/md_traj.pdb" -dt 1

# Step 4: Clean up intermediate files
echo "Cleaning up intermediate files..."
rm "$VMD_DIR/md_centered.xtc"

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