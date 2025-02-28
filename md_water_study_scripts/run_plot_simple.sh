#!/bin/bash

# Unset Python-related environment variables
unset PYTHONHOME
unset PYTHONPATH

# Run the plotting script with the system Python
/usr/bin/python3 -c "
import sys
print('Python executable:', sys.executable)
print('Python version:', sys.version)
"

/usr/bin/python3 plot_all_results.py --model tip4p --temp 273 