#!/bin/bash

# Script to create a movie of water simulation using VMD
# This script runs VMD in text mode to render frames, then uses ffmpeg to create a movie

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$( dirname "$SCRIPT_DIR" )"

# Check if VMD is installed
if ! command -v vmd &> /dev/null; then
    echo "Error: VMD is not installed or not in PATH."
    echo "Please install VMD from https://www.ks.uiuc.edu/Research/vmd/"
    exit 1
fi

# Check if ffmpeg is installed
if ! command -v ffmpeg &> /dev/null; then
    echo "Warning: ffmpeg is not installed or not in PATH."
    echo "You will need to install ffmpeg to create the final movie."
    echo "The script will still render individual frames."
fi

# Set default values
MODEL="tip4p"
TEMP="273"
RENDER_MODE="text"  # text or gui

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
        --gui)
            RENDER_MODE="gui"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--model MODEL] [--temp TEMP] [--gui]"
            exit 1
            ;;
    esac
done

# Check if the trajectory file exists
TRAJECTORY="${BASE_DIR}/data/${MODEL}/${TEMP}K/outputs/md.xtc"
if [ ! -f "$TRAJECTORY" ]; then
    echo "Error: Trajectory file not found: $TRAJECTORY"
    exit 1
fi

# Check if structure files exists - try md.gro, then em.gro
STRUCTURE="${BASE_DIR}/data/${MODEL}/${TEMP}K/outputs/md.gro"
if [ ! -f "$STRUCTURE" ]; then
    STRUCTURE="${BASE_DIR}/data/${MODEL}/${TEMP}K/outputs/em.gro"
    if [ ! -f "$STRUCTURE" ]; then
        echo "Error: Structure file not found: ${BASE_DIR}/data/${MODEL}/${TEMP}K/outputs/md.gro or em.gro"
        exit 1
    fi
fi

# Create output directory
mkdir -p "${BASE_DIR}/data/${MODEL}/${TEMP}K/movie/frames"

echo "Starting VMD to render frames..."
echo "This may take some time depending on the size of your trajectory."
echo "Structure file: $STRUCTURE"
echo "Trajectory file: $TRAJECTORY"

# Run VMD with the render script
cd "$BASE_DIR"
if [ "$RENDER_MODE" = "text" ]; then
    # Run in text mode (no GUI)
    vmd -dispdev text "$STRUCTURE" "$TRAJECTORY" -e "$SCRIPT_DIR/create_vmd_movie.tcl"
else
    # Run in GUI mode
    vmd "$STRUCTURE" "$TRAJECTORY" -e "$SCRIPT_DIR/create_vmd_movie.tcl"
fi

# Check if the movie creation script was created
if [ -f "${BASE_DIR}/data/${MODEL}/${TEMP}K/movie/create_movie.sh" ]; then
    echo "Frame rendering complete."
    echo "Creating movie using ffmpeg..."
    
    # Run the ffmpeg script
    bash "${BASE_DIR}/data/${MODEL}/${TEMP}K/movie/create_movie.sh"
    
    # Check if the movie was created
    if [ -f "${BASE_DIR}/data/${MODEL}/${TEMP}K/movie/water_simulation.mp4" ]; then
        echo "Movie created successfully: ${BASE_DIR}/data/${MODEL}/${TEMP}K/movie/water_simulation.mp4"
        
        # Try to play the movie if a player is available
        if command -v xdg-open &> /dev/null; then
            echo "Opening movie..."
            xdg-open "${BASE_DIR}/data/${MODEL}/${TEMP}K/movie/water_simulation.mp4"
        else
            echo "You can view the movie at: ${BASE_DIR}/data/${MODEL}/${TEMP}K/movie/water_simulation.mp4"
        fi
    else
        echo "Error: Movie creation failed."
        echo "You can try to create the movie manually by running:"
        echo "bash ${BASE_DIR}/data/${MODEL}/${TEMP}K/movie/create_movie.sh"
    fi
else
    echo "Error: Frame rendering failed."
    echo "Check the VMD output for errors."
fi 