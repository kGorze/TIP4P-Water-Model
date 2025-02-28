#!/bin/bash

# Unset Python-related environment variables
unset PYTHONHOME
unset PYTHONPATH

# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate

# Install required packages with specific versions
pip install 'numpy<2' matplotlib seaborn pandas

# Run the plotting script with python3
python3 plot_all_results.py --model tip4p --temp 273

# Deactivate virtual environment
deactivate
