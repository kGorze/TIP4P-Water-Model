#!/bin/bash

# Script to create a movie of TIP4P-ICE water simulation using VMD
# This script runs VMD to render frames, then uses ffmpeg to create a movie

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
    echo "Warning: ffmpeg is not installed or not in PATH."
    echo "You will need to install ffmpeg to create the final movie."
    echo "The script will still render individual frames."
fi

# Check if the trajectory files exist
# We'll use the largest trajectory file for the best movie quality
if [ -f "${DATA_DIR}/npt.trr" ]; then
    TRAJECTORY="${DATA_DIR}/npt.trr"
    echo "Using NPT trajectory file (largest): $TRAJECTORY"
elif [ -f "${DATA_DIR}/nvt.trr" ]; then
    TRAJECTORY="${DATA_DIR}/nvt.trr"
    echo "Using NVT trajectory file: $TRAJECTORY"
elif [ -f "${DATA_DIR}/md.trr" ]; then
    TRAJECTORY="${DATA_DIR}/md.trr"
    echo "Using MD trajectory file: $TRAJECTORY"
else
    echo "Error: No suitable trajectory file found."
    exit 1
fi

# Check if structure files exists - try md.gro, then npt.gro, then nvt.gro, then em.gro
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

# Create output directory
mkdir -p "${DATA_DIR}/movie/frames"

# Create the VMD script for rendering
cat > "${DATA_DIR}/render_tip4pice_movie.tcl" << 'EOF'
# VMD script for rendering frames of the TIP4P-ICE water simulation

# Get the current molecule ID
set mol_id [molinfo top]
puts "Rendering molecule ID: $mol_id"

# Get the number of frames
set num_frames [molinfo $mol_id get numframes]
puts "Number of frames: $num_frames"

# Create output directory for frames
set out_dir [file join [pwd] "data" "movie" "frames"]
file mkdir $out_dir
puts "Output directory: $out_dir"

# Create directory for the final movie
set movie_dir [file join [pwd] "data" "movie"]
file mkdir $movie_dir

# Set up nice visualization
mol delrep 0 $mol_id

# Display water molecules as CPK
mol representation CPK 1.0 0.3 18.0 16.0
mol selection {all}
mol color Name
mol addrep $mol_id

# Add hydrogen bonds
mol representation HBonds 3.0 30.0 5.0
mol selection {name OW}
mol color ColorID 4  ;# blue for hydrogen bonds
mol addrep $mol_id

# Set up display - use even numbers for width and height for ffmpeg compatibility
display resize 1280 720
display projection Orthographic
display depthcue off
display nearclip set 0.01
color Display Background white

# Center the view
mol showperiodic $mol_id 1
molinfo $mol_id set center_matrix "{{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}"

# Rotate to a nice viewing angle
rotate x by -60
rotate y by 30
rotate z by 15
display resetview

# Render frames
set frame_dir $out_dir
set render_params "-size 1280 720 -format TGA -renderer TachyonInternal"

puts "Rendering $num_frames frames..."

# Use a small step size to get the most frames possible
# but ensure we don't create too many frames if the trajectory is very large
set step_size 1
if {$num_frames > 2000} {
    set step_size [expr $num_frames / 1000]
    puts "Trajectory is very large, using step size of $step_size to limit to ~1000 frames"
} elseif {$num_frames > 500} {
    set step_size [expr $num_frames / 500]
    puts "Trajectory is large, using step size of $step_size to limit to ~500 frames"
}

# Create a shell script to combine the frames into a movie
set fp [open [file join $movie_dir "create_movie.sh"] w]
puts $fp "#!/bin/bash\n"
puts $fp "# Shell script to create a movie from rendered frames\n"
puts $fp "cd [pwd]\n"
puts $fp "# Create movie with ffmpeg"
puts $fp "ffmpeg -framerate 30 -pattern_type glob -i '$out_dir/water_*.tga' \\"
puts $fp "       -vf 'scale=1280:720:force_original_aspect_ratio=decrease,pad=1280:720:(ow-iw)/2:(oh-ih)/2' \\"
puts $fp "       -c:v libx264 -pix_fmt yuv420p -crf 18 \\"
puts $fp "       '$movie_dir/tip4pice_water_simulation.mp4'\n"
puts $fp "echo \"Movie created: $movie_dir/tip4pice_water_simulation.mp4\"\n"
close $fp

# Make the shell script executable
file attributes [file join $movie_dir "create_movie.sh"] -permissions 0755

# Loop through frames and render
for {set i 0} {$i < $num_frames} {incr i $step_size} {
    # Update to the current frame
    animate goto $i
    # Force update of the display
    display update
    # Create the file name with leading zeros for proper sorting
    set filename [format "water_%04d.tga" $i]
    set filepath [file join $frame_dir $filename]
    # Render the frame
    render TachyonInternal $filepath
    puts "Rendered frame $i/$num_frames: $filepath"
}

puts "Rendering complete. Frames saved to: $out_dir"
puts "To create the final movie, run: $movie_dir/create_movie.sh"

# Exit VMD when done
quit
EOF

echo "Starting VMD to render frames..."
echo "This may take some time depending on the size of your trajectory."
echo "Structure file: $STRUCTURE"
echo "Trajectory file: $TRAJECTORY"

# Run VMD with the render script
cd "$SCRIPT_DIR"
vmd "$STRUCTURE" "$TRAJECTORY" -e "${DATA_DIR}/render_tip4pice_movie.tcl"

# Check if the movie creation script was created
if [ -f "${DATA_DIR}/movie/create_movie.sh" ]; then
    echo "Frame rendering complete."
    echo "Creating movie using ffmpeg..."
    
    # Run the ffmpeg script
    bash "${DATA_DIR}/movie/create_movie.sh"
    
    # Check if the movie was created
    if [ -f "${DATA_DIR}/movie/tip4pice_water_simulation.mp4" ]; then
        echo "Movie created successfully: ${DATA_DIR}/movie/tip4pice_water_simulation.mp4"
        
        # Create an analysis directory if it doesn't exist
        mkdir -p "${SCRIPT_DIR}/analysis"
        
        # Copy the movie to the analysis directory
        cp "${DATA_DIR}/movie/tip4pice_water_simulation.mp4" "${SCRIPT_DIR}/analysis/"
        
        echo "Movie copied to analysis directory: ${SCRIPT_DIR}/analysis/tip4pice_water_simulation.mp4"
        
        # Try to play the movie if a player is available
        if command -v xdg-open &> /dev/null; then
            echo "Opening movie..."
            xdg-open "${SCRIPT_DIR}/analysis/tip4pice_water_simulation.mp4"
        else
            echo "You can view the movie at: ${SCRIPT_DIR}/analysis/tip4pice_water_simulation.mp4"
        fi
    else
        echo "Error: Movie creation failed."
        echo "You can try to create the movie manually by running:"
        echo "bash ${DATA_DIR}/movie/create_movie.sh"
    fi
else
    echo "Error: Frame rendering failed."
    echo "Check the VMD output for errors."
fi 