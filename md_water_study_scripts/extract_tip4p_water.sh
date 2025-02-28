#!/bin/bash

# Create a directory for temporary files
mkdir -p temp_water
cd temp_water

# Create a small box of TIP4P water molecules
echo -e "1\n" | gmx solvate -cs tip4p -box 1 1 1 -o tip4p_box.gro

# Extract a single water molecule
echo -e "0 & r 1\n" | gmx trjconv -f tip4p_box.gro -s tip4p_box.gro -o single_tip4p.gro -pbc mol

# Convert to PDB format
gmx editconf -f single_tip4p.gro -o single_tip4p.pdb

# Move the files to the parent directory
mv single_tip4p.gro single_tip4p.pdb ..

# Clean up
cd ..
rm -rf temp_water

echo "Extracted a single TIP4P water molecule from GROMACS" 