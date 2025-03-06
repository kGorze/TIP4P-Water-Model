#!/bin/bash

# Script to create a high-quality movie of TIP4P-ICE water simulation using VMD
# This approach renders individual frames and then combines them with ffmpeg
# Note: This script handles the case where the trajectory has limited frames
# by using frame interpolation to create a smoother movie

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DATA_DIR="${SCRIPT_DIR}/data"

# Check if VMD is installed
if ! command -v vmd &> /dev/null; then
    echo "Error: VMD is not installed or not in PATH."
    echo "Please install VMD from https://www.ks.uiuc.edu/Research/vmd/"
    exit 1
fi

# Check if ffmpeg is installed
if ! command -v ffmpeg &> /dev/null; then
    echo "Error: ffmpeg is not installed or not in PATH."
    echo "Please install ffmpeg to create the movie."
    exit 1
fi

# Clean up any previous movie files
rm -rf "${DATA_DIR}/movie"
mkdir -p "${DATA_DIR}/movie/frames"

# Find the best trajectory file - prioritize the largest one for more frames
# Check file sizes and use the largest one
NPT_SIZE=0
NVT_SIZE=0
MD_SIZE=0

if [ -f "${DATA_DIR}/npt.trr" ]; then
    NPT_SIZE=$(stat -c%s "${DATA_DIR}/npt.trr")
fi

if [ -f "${DATA_DIR}/nvt.trr" ]; then
    NVT_SIZE=$(stat -c%s "${DATA_DIR}/nvt.trr")
fi

if [ -f "${DATA_DIR}/md.trr" ]; then
    MD_SIZE=$(stat -c%s "${DATA_DIR}/md.trr")
fi

# Compare sizes and select the largest file
if [ $NPT_SIZE -ge $NVT_SIZE ] && [ $NPT_SIZE -ge $MD_SIZE ] && [ $NPT_SIZE -gt 0 ]; then
    TRAJECTORY="${DATA_DIR}/npt.trr"
    echo "Using NPT trajectory file (largest at $(( NPT_SIZE / 1024 / 1024 )) MB): $TRAJECTORY"
elif [ $NVT_SIZE -ge $NPT_SIZE ] && [ $NVT_SIZE -ge $MD_SIZE ] && [ $NVT_SIZE -gt 0 ]; then
    TRAJECTORY="${DATA_DIR}/nvt.trr"
    echo "Using NVT trajectory file (largest at $(( NVT_SIZE / 1024 / 1024 )) MB): $TRAJECTORY"
elif [ $MD_SIZE -gt 0 ]; then
    TRAJECTORY="${DATA_DIR}/md.trr"
    echo "Using MD trajectory file (largest at $(( MD_SIZE / 1024 / 1024 )) MB): $TRAJECTORY"
else
    echo "Error: No suitable trajectory file found."
    exit 1
fi

# Find the best structure file
if [ -f "${DATA_DIR}/md.gro" ]; then
    STRUCTURE="${DATA_DIR}/md.gro"
elif [ -f "${DATA_DIR}/npt.gro" ]; then
    STRUCTURE="${DATA_DIR}/npt.gro"
elif [ -f "${DATA_DIR}/nvt.gro" ]; then
    STRUCTURE="${DATA_DIR}/nvt.gro"
elif [ -f "${DATA_DIR}/em.gro" ]; then
    STRUCTURE="${DATA_DIR}/em.gro"
else
    echo "Error: No suitable structure file found."
    exit 1
fi

echo "Using structure file: $STRUCTURE"

echo "Starting VMD to render individual frames..."
echo "This may take some time depending on the size of your trajectory."

# Run VMD with the TCL script to render frames
cd "${DATA_DIR}"
vmd "$STRUCTURE" "$TRAJECTORY" -e "render_frames.tcl"

# Check if the frames were created
if [ -d "${DATA_DIR}/movie/frames" ] && [ "$(ls -A "${DATA_DIR}/movie/frames")" ]; then
    echo "Frames rendered successfully."
    
    # Run the ffmpeg script to create the movie
    echo "Creating movie from frames using ffmpeg..."
    bash "${DATA_DIR}/movie/create_movie.sh"
    
    # Check if the movie was created
    if [ -f "${DATA_DIR}/movie/tip4pice_water.mp4" ]; then
        echo "Movie created successfully: ${DATA_DIR}/movie/tip4pice_water.mp4"
        
        # Create an analysis directory if it doesn't exist
        mkdir -p "${SCRIPT_DIR}/analysis"
        
        # Copy the movie to the analysis directory
        cp "${DATA_DIR}/movie/tip4pice_water.mp4" "${SCRIPT_DIR}/analysis/"
        
        echo "Movie copied to analysis directory: ${SCRIPT_DIR}/analysis/tip4pice_water.mp4"
        
        # Try to play the movie if a player is available
        if command -v xdg-open &> /dev/null; then
            echo "Opening movie..."
            xdg-open "${SCRIPT_DIR}/analysis/tip4pice_water.mp4"
        else
            echo "You can view the movie at: ${SCRIPT_DIR}/analysis/tip4pice_water.mp4"
        fi
        
        echo ""
        echo "NOTE: The trajectory files only contain a limited number of frames."
        echo "To get more frames in your movie, you would need to:"
        echo "1. Check your GROMACS simulation parameters to ensure frames are saved more frequently"
        echo "2. Rerun the simulation with proper frame-saving parameters"
        echo "3. This movie uses frame interpolation to create a smoother appearance despite limited frames"
    else
        echo "Error: Movie creation failed."
        echo "You can try to create the movie manually by running:"
        echo "bash ${DATA_DIR}/movie/create_movie.sh"
    fi
else
    echo "Error: Frame rendering failed."
    echo "Check the VMD output for errors."
fi 