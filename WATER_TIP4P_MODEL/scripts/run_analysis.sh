#!/bin/bash
# Script to run analysis on the TIP4P water simulation results

set -e  # Exit on error

# Create analysis directory
mkdir -p ../analysis

# Move to the data directory
cd ../data

# Check that the simulation data exists
if [ ! -f md.xtc ] || [ ! -f md.tpr ]; then
    echo "Error: Simulation files (md.xtc and/or md.tpr) not found. Please run the simulation first."
    exit 1
fi

echo "Starting analysis of TIP4P water simulation at 273K..."

# 1. Calculate Radial Distribution Functions (RDFs)
echo "1. Calculating RDFs..."
# O-O RDF
echo "   - Oxygen-Oxygen RDF..."
gmx rdf -f md.xtc -s md.tpr -o ../analysis/rdf_OO.xvg -ref "name OW" -sel "name OW" << EOF
EOF

# O-H RDF
echo "   - Oxygen-Hydrogen RDF..."
gmx rdf -f md.xtc -s md.tpr -o ../analysis/rdf_OH.xvg -ref "name OW" -sel "name HW" << EOF
EOF

# H-H RDF
echo "   - Hydrogen-Hydrogen RDF..."
gmx rdf -f md.xtc -s md.tpr -o ../analysis/rdf_HH.xvg -ref "name HW" -sel "name HW" << EOF
EOF

# 2. Calculate Mean Square Displacement (MSD) for diffusion coefficient
echo "2. Calculating Mean Square Displacement for diffusion coefficient..."
gmx msd -f md.xtc -s md.tpr -o ../analysis/msd.xvg -beginfit 1000 -endfit 2000 << EOF
0
EOF

# 3. Calculate hydrogen bonding statistics
echo "3. Calculating hydrogen bonding statistics..."
gmx hbond -f md.xtc -s md.tpr -num ../analysis/hbonds.xvg -life ../analysis/hblife.xvg << EOF
0
0
EOF

# 4. Calculate density profile (checking for uniformity)
echo "4. Calculating density profile..."
gmx density -f md.xtc -s md.tpr -o ../analysis/density.xvg -dens mass -d Z << EOF
0
EOF

# 5. Calculate velocity autocorrelation function (if velocities are available)
if [ -f md.trr ]; then
    echo "5. Calculating velocity autocorrelation function..."
    gmx velacc -f md.trr -s md.tpr -o ../analysis/vacf.xvg -acflen 1000 << EOF
0
EOF
else
    echo "5. Skipping VACF calculation (md.trr not found)"
fi

# 6. Calculate thermodynamic properties
echo "6. Extracting thermodynamic properties..."
# Temperature
gmx energy -f md.edr -o ../analysis/temperature.xvg << EOF
Temperature
EOF

# Pressure
gmx energy -f md.edr -o ../analysis/pressure.xvg << EOF
Pressure
EOF

# Density
gmx energy -f md.edr -o ../analysis/box_density.xvg << EOF
Density
EOF

# Total Energy
gmx energy -f md.edr -o ../analysis/energy.xvg << EOF
Potential
Kinetic-En.
Total-Energy
EOF

echo "Analysis completed! Results saved in the analysis directory."
echo "You may want to visualize the .xvg files using a plotting program like Grace, matplotlib, or Excel." 