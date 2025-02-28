# VMD script to visualize water simulation and render a movie
# Usage: vmd -dispdev text -e render_movie.tcl

# Set the input files
set tpr_file "data/tip4p/273K/md.tpr"
set xtc_file "data/tip4p/273K/md.xtc"
set output_dir "data/tip4p/273K/movie"

# Create output directory if it doesn't exist
file mkdir $output_dir

# Load the simulation
mol new $tpr_file waitfor all
mol addfile $xtc_file waitfor all

# Set the visualization style
mol delrep 0 top
mol representation CPK 0.8 0.3 12.0 12.0
mol color Name
mol selection "all"
mol material Opaque
mol addrep top

# Add a representation for hydrogen bonds
mol representation HBonds 3.0 30.0
mol color ColorID 4
mol selection "name OW"
mol material Opaque
mol addrep top

# Set the background color to white
color Display Background white

# Center the view on the system
molinfo top set center_matrix "{{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}"
display resetview

# Rotate to a nice angle
rotate x by -60
rotate y by 30
rotate z by 15

# Set the number of frames to render
set num_frames [molinfo top get numframes]
set start_frame 0
set end_frame [expr $num_frames - 1]
set step 10  # Render every 10th frame to reduce file size

# Render frames as PNG images
for {set i $start_frame} {$i <= $end_frame} {incr i $step} {
    # Go to the frame
    animate goto $i
    
    # Update the display
    display update
    
    # Render the image
    render TachyonInternal $output_dir/frame_[format "%04d" $i].png
    
    # Print progress
    puts "Rendered frame $i of $end_frame"
}

# Create a script to convert the frames to a movie using ffmpeg
set ffmpeg_script "$output_dir/create_movie.sh"
set fp [open $ffmpeg_script w]
puts $fp "#!/bin/bash"
puts $fp "cd $output_dir"
puts $fp "ffmpeg -framerate 30 -pattern_type glob -i 'frame_*.png' -c:v libx264 -pix_fmt yuv420p -crf 23 water_simulation.mp4"
puts $fp "echo 'Movie created: water_simulation.mp4'"
close $fp

# Make the script executable
exec chmod +x $ffmpeg_script

# Exit VMD
puts "Rendering complete. Run $ffmpeg_script to create the movie."
quit 