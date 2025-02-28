#!/bin/bash

# Script to view water simulation in VMD

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$( dirname "$SCRIPT_DIR" )"

# Check if VMD is installed
if ! command -v vmd &> /dev/null; then
    echo "Error: VMD is not installed or not in PATH."
    echo "Please install VMD from https://www.ks.uiuc.edu/Research/vmd/"
    exit 1
fi

# Check if the trajectory file exists
if [ ! -f "$BASE_DIR/data/tip4p/273K/md.xtc" ]; then
    echo "Error: Trajectory file not found: $BASE_DIR/data/tip4p/273K/md.xtc"
    exit 1
fi

# Check if the topology file exists
if [ ! -f "$BASE_DIR/data/tip4p/273K/md.tpr" ]; then
    echo "Error: Topology file not found: $BASE_DIR/data/tip4p/273K/md.tpr"
    exit 1
fi

# Create movie directory if it doesn't exist
mkdir -p "$BASE_DIR/data/tip4p/273K/movie"

# Display instructions
echo "==================================================================="
echo "                 Water Simulation Viewer"
echo "==================================================================="
echo "This script will open VMD with your water simulation loaded."
echo ""
echo "To create a movie of your simulation:"
echo "1. Go to Extensions -> Visualization -> Movie Maker"
echo "2. Set the movie settings:"
echo "   - Trajectory: 0 to [last frame]"
echo "   - Step Size: 1 or higher (higher = faster movie, smaller file)"
echo "   - Format: MPEG or AVI"
echo "   - Renderer: Snapshot"
echo "   - Working Directory: $BASE_DIR/data/tip4p/273K/movie"
echo "   - Filename: water_simulation"
echo "3. Click 'Make Movie'"
echo ""
echo "For more advanced visualization, you can copy and paste commands"
echo "from the file: $SCRIPT_DIR/vmd_commands.txt"
echo "==================================================================="
echo ""
echo "Press Enter to launch VMD..."
read

# Launch VMD with the simulation
cd "$BASE_DIR"
vmd data/tip4p/273K/md.tpr data/tip4p/273K/md.xtc 