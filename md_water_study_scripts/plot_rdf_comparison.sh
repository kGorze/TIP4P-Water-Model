#!/bin/bash

# This script runs the plot_results.py script with the --compare option
# to generate an RDF plot that's comparable to Figure 2 in the article

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Check if model and temperature are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <model> <temperature>"
    echo "Example: $0 tip4p-ice 298"
    exit 1
fi

MODEL=$1
TEMP=$2

# Run the plotting script with the compare option
echo "Generating RDF comparison plot for $MODEL at ${TEMP}K..."
python3 "$SCRIPT_DIR/plot_results.py" --model "$MODEL" --temp "$TEMP" --compare

# Check if the plot was generated successfully
DATA_DIR="$SCRIPT_DIR/../data/$MODEL/${TEMP}K"
PLOT_FILE="$DATA_DIR/rdf_plot_comparison.png"

if [ -f "$PLOT_FILE" ]; then
    echo "Plot generated successfully: $PLOT_FILE"
    
    # Try to display the plot if a display is available
    if command -v display &> /dev/null; then
        echo "Displaying plot..."
        display "$PLOT_FILE"
    elif command -v xdg-open &> /dev/null; then
        echo "Opening plot..."
        xdg-open "$PLOT_FILE"
    else
        echo "Plot saved to: $PLOT_FILE"
        echo "Install ImageMagick or use a file browser to view the plot."
    fi
else
    echo "Error: Plot was not generated."
    exit 1
fi 