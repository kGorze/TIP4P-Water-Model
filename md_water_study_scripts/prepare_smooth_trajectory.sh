#!/bin/bash

# Script to prepare a smooth trajectory for VMD visualization
# This script converts GROMACS trajectory files to PDB format for better VMD visualization

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$( dirname "$SCRIPT_DIR" )"

# Check if GROMACS is installed (for trjconv)
if ! command -v gmx &> /dev/null; then
    echo "Error: GROMACS is not installed or not in PATH."
    echo "Please install GROMACS to convert trajectory files."
    exit 1
fi

# Set default values
MODEL="tip4p"
TEMP="273"
STRIDE=1  # Take every frame for maximum detail

# Parse command line options
while [[ $# -gt 0 ]]; do
    case $1 in
        --model)
            MODEL="$2"
            shift 2
            ;;
        --temp)
            TEMP="$2"
            shift 2
            ;;
        --stride)
            STRIDE="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--model MODEL] [--temp TEMP] [--stride N]"
            exit 1
            ;;
    esac
done

# Set paths
TRAJECTORY="${BASE_DIR}/data/${MODEL}/${TEMP}K/outputs/md.xtc"
STRUCTURE="${BASE_DIR}/data/${MODEL}/${TEMP}K/outputs/md.gro"
TPR_FILE="${BASE_DIR}/data/${MODEL}/${TEMP}K/outputs/md.tpr"
VMD_DIR="${BASE_DIR}/data/${MODEL}/${TEMP}K/vmd_smooth"

# Check if files exist
if [ ! -f "$TRAJECTORY" ]; then
    echo "Error: Trajectory file not found: $TRAJECTORY"
    exit 1
fi

if [ ! -f "$STRUCTURE" ]; then
    STRUCTURE="${BASE_DIR}/data/${MODEL}/${TEMP}K/outputs/em.gro"
    if [ ! -f "$STRUCTURE" ]; then
        echo "Error: Structure file not found: ${BASE_DIR}/data/${MODEL}/${TEMP}K/outputs/md.gro or em.gro"
        exit 1
    fi
fi

if [ ! -f "$TPR_FILE" ]; then
    echo "Error: TPR file not found: $TPR_FILE"
    exit 1
fi

# Create output directory
mkdir -p "$VMD_DIR"

echo "Preparing smoothed trajectory for VMD visualization..."
echo "This may take some time depending on the size of your trajectory."
echo "Structure file: $STRUCTURE"
echo "Trajectory file: $TRAJECTORY"
echo "TPR file: $TPR_FILE"
echo "Stride: Every $STRIDE frames (maximum detail)"

# Convert the final structure to PDB format
echo "Converting structure to PDB format..."
# Use echo to provide input for group selection (0 for System, 0 for output)
echo -e "0\n0" | gmx trjconv -s "$TPR_FILE" -f "$STRUCTURE" -o "${VMD_DIR}/md_final.pdb" -pbc mol -center

# Convert the trajectory to PDB format with reduced frames
echo "Converting trajectory to smoothed PDB format (stride=$STRIDE)..."
# Use echo to provide input for group selection (0 for System, 0 for output)
echo -e "0\n0" | gmx trjconv -s "$TPR_FILE" -f "$TRAJECTORY" -o "${VMD_DIR}/md_traj.pdb" -pbc mol -center -dt $STRIDE

echo "Conversion complete!"
echo "Structure file: ${VMD_DIR}/md_final.pdb"
echo "Trajectory file: ${VMD_DIR}/md_traj.pdb"
echo ""
echo "To visualize, run: ${SCRIPT_DIR}/launch_smooth_vmd.sh" 