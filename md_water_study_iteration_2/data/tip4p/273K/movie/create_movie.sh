#!/bin/bash
cd data/tip4p/273K/movie
ffmpeg -framerate 30 -pattern_type glob -i 'frame_*.png' -c:v libx264 -pix_fmt yuv420p -crf 23 water_simulation.mp4
echo 'Movie created: water_simulation.mp4'
