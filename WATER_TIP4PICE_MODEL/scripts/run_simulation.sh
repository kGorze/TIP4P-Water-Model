#!/bin/bash
# Script to run the complete TIP4P water simulation workflow

set -e  # Exit on error

# Output directories
mkdir -p ../data

# Move to the data directory
cd ../data

# 1. Generate water box using PACKMOL
echo "Step 1: Generating water box with PACKMOL..."
# Copy water.pdb if available
if [ ! -f water.pdb ]; then
    if [ -f ../../md_water_study_iteration_3/data/water.pdb ]; then
        cp ../../md_water_study_iteration_3/data/water.pdb .
    else
        echo "Error: water.pdb not found. Please place a water.pdb file in the data directory."
        exit 1
    fi
fi

packmol < ../configs/water_box.inp
echo "Water box generated."

# 2. Generate topology with GROMACS using TIP4P water model
echo "Step 2: Generating topology with TIP4P water model..."
gmx pdb2gmx -f water_box.pdb -o water_box.gro -water tip4p -p topol.top || { echo "pdb2gmx failed"; exit 1; }
echo "Topology generated."

# 3. Perform energy minimization
echo "Step 3: Running energy minimization..."
gmx grompp -f ../configs/em.mdp -c water_box.gro -p topol.top -o em.tpr || { echo "grompp for em failed"; exit 1; }
gmx mdrun -deffnm em || { echo "mdrun for em failed"; exit 1; }
echo "Energy minimization completed."

# 4. Perform NVT equilibration
echo "Step 4: Running NVT equilibration..."
gmx grompp -f ../configs/nvt.mdp -c em.gro -p topol.top -o nvt.tpr || { echo "grompp for nvt failed"; exit 1; }
gmx mdrun -deffnm nvt || { echo "mdrun for nvt failed"; exit 1; }
echo "NVT equilibration completed."

# 5. Perform NPT equilibration
echo "Step 5: Running NPT equilibration..."
gmx grompp -f ../configs/npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr || { echo "grompp for npt failed"; exit 1; }
gmx mdrun -deffnm npt || { echo "mdrun for npt failed"; exit 1; }
echo "NPT equilibration completed."

# 6. Perform production MD
echo "Step 6: Running production MD (2 ns)..."
gmx grompp -f ../configs/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr || { echo "grompp for md failed"; exit 1; }
gmx mdrun -deffnm md || { echo "mdrun for md failed"; exit 1; }
echo "Production MD completed."

echo "All simulation steps completed successfully!"
echo "You can now run the analysis scripts from the scripts directory." 