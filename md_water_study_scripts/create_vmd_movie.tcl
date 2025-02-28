# VMD script for rendering frames of the TIP4P water simulation

# Check if mol exists
if {[molinfo num] == 0} {
    puts "Error: No molecule loaded."
    puts "Please load a trajectory before running this script."
    exit
}

# Get the current molecule ID
set mol_id [molinfo top]
puts "Rendering molecule ID: $mol_id"

# Get the number of frames
set num_frames [molinfo $mol_id get numframes]
puts "Number of frames: $num_frames"

# Create output directory for frames
set out_dir [file join [pwd] "data" "tip4p" "273K" "movie" "frames"]
file mkdir $out_dir
puts "Output directory: $out_dir"

# Create directory for the final movie
set movie_dir [file join [pwd] "data" "tip4p" "273K" "movie"]
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
mol selection {name O}
mol color ColorID 4  ;# blue for hydrogen bonds
mol addrep $mol_id

# Set up display
display resize 800 600
display projection Orthographic
display depthcue off
display nearclip set 0.01

# Center the view
mol showperiodic $mol_id 1
molinfo $mol_id set center_matrix "{{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}"

# Rotate to a nice viewing angle
rotate x by 20
rotate y by 30
rotate z by 10
display resetview

# Render frames
set frame_dir $out_dir
set render_params "-size 800 600 -format RGB -renderer snapshot"

puts "Rendering $num_frames frames..."

# Only render a subset of frames if the trajectory is very large
set step_size 1
if {$num_frames > 1000} {
    set step_size [expr $num_frames / 200]
    puts "Trajectory is large, using step size of $step_size"
}

# Create a shell script to combine the frames into a movie
set fp [open [file join $movie_dir "create_movie.sh"] w]
puts $fp "#!/bin/bash\n"
puts $fp "# Shell script to create a movie from rendered frames\n"
puts $fp "cd [pwd]\n"
puts $fp "# Create movie with ffmpeg"
puts $fp "ffmpeg -framerate 30 -pattern_type glob -i '$out_dir/water_*.rgb' \\"
puts $fp "       -c:v libx264 -pix_fmt yuv420p -crf 23 \\"
puts $fp "       '$movie_dir/water_simulation.mp4'\n"
puts $fp "echo \"Movie created: $movie_dir/water_simulation.mp4\"\n"
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
    set filename [format "water_%04d.rgb" $i]
    set filepath [file join $frame_dir $filename]
    # Render the frame
    render snapshot $filepath
    puts "Rendered frame $i/$num_frames: $filepath"
}

puts "Rendering complete. Frames saved to: $out_dir"
puts "To create the final movie, run: $movie_dir/create_movie.sh"

# Don't exit VMD automatically - let the user decide
puts "Script finished. You can explore the visualization or exit VMD." 