#!/bin/bash

# Script to launch VMD with the full trajectory (all frames)

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$( dirname "$SCRIPT_DIR" )"

# Set paths
MODEL="tip4p"
TEMP="273"
VMD_DIR="$BASE_DIR/data/$MODEL/${TEMP}K/vmd_full"
MOVIE_DIR="$BASE_DIR/data/$MODEL/${TEMP}K/movie_full"

# Check if VMD is installed
if ! command -v vmd &> /dev/null; then
    echo "Error: VMD is not installed or not in PATH."
    echo "Please install VMD from https://www.ks.uiuc.edu/Research/vmd/"
    exit 1
fi

# Check if the converted files exist
if [ ! -f "$VMD_DIR/md_final.pdb" ]; then
    echo "Error: Structure file not found: $VMD_DIR/md_final.pdb"
    echo "Please run the prepare_full_trajectory.sh script first."
    exit 1
fi

if [ ! -f "$VMD_DIR/md_traj.pdb" ]; then
    echo "Error: Trajectory file not found: $VMD_DIR/md_traj.pdb"
    echo "Please run the prepare_full_trajectory.sh script first."
    exit 1
fi

# Create movie directory if it doesn't exist
mkdir -p "$MOVIE_DIR"

# Display instructions
echo "==================================================================="
echo "             Full Water Simulation Viewer"
echo "==================================================================="
echo "This script will open VMD with your water simulation loaded with ALL frames."
echo "This will give you the smoothest possible visualization."
echo ""
echo "To create a movie of your simulation:"
echo "1. Go to Extensions -> Visualization -> Movie Maker"
echo "2. Set the movie settings:"
echo "   - Working Directory: $MOVIE_DIR"
echo "   - Filename: water_simulation_full"
echo "   - Format: MPEG or AVI"
echo "   - Trajectory: 0 to [last frame]"
echo "   - Step Size: 1 or higher (higher = faster movie, smaller file)"
echo "   - Renderer: Snapshot"
echo "3. Click 'Make Movie'"
echo ""
echo "For better visualization, once VMD is open:"
echo "1. Go to Graphics -> Representations"
echo "2. Delete the default representation"
echo "3. Create a new representation with:"
echo "   - Drawing Method: CPK"
echo "   - Selection: all"
echo "   - Coloring Method: Name"
echo "4. Create another representation with:"
echo "   - Drawing Method: HBonds"
echo "   - Selection: name O"
echo "   - Coloring Method: ColorID 4 (blue)"
echo "==================================================================="
echo ""
echo "Press Enter to launch VMD..."
read

# Launch VMD with the converted files
cd "$BASE_DIR"
vmd "$VMD_DIR/md_final.pdb" "$VMD_DIR/md_traj.pdb" 