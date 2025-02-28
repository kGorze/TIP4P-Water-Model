# Water MD Study Iteration 4: TIP4P Water at 273K

This directory contains the configuration files and scripts for running a molecular dynamics simulation of TIP4P water at 273K. This is the 4th iteration of the water study, with significant improvements over previous iterations:

## Improvements in Iteration 4

1. **TIP4P Water Model**: Using the TIP4P water model with proper virtual site (M-site) handling for better accuracy near the freezing point.
2. **Anisotropic Pressure Coupling**: Implemented anisotropic pressure coupling to allow the simulation box to change shape independently in each dimension, which is important near the freezing point where water may form structured regions.
3. **Optimized Minimization Protocol**: Using constraints=none during minimization to allow better relaxation of bad contacts.
4. **Extended Production Run**: The production MD is now 2 ns (1,000,000 steps) for better sampling of water properties.
5. **Improved Barostat Selection**: Using Berendsen/c-rescale barostat for equilibration and Parrinello-Rahman for production.

## Directory Structure

- `configs/`: Contains all GROMACS configuration files and the PACKMOL input script
- `scripts/`: Analysis scripts for post-processing
- `data/`: Will contain simulation outputs
- `analysis/`: Will contain analysis results and plots
- `docs/`: Documentation

## Simulation Protocol

### System Setup
- 5500 TIP4P water molecules
- Initial density: ~0.93 g/cm³ (lower initial density for minimization)
- Target box size: ~5.35 nm per side at equilibrium

### Workflow
1. **System Generation**: Use PACKMOL to create an initial configuration in a 7.8 nm box
2. **Topology Creation**: Generate topology with GROMACS pdb2gmx using the TIP4P water model
3. **Energy Minimization**: Remove bad contacts with steepest descent minimization
4. **NVT Equilibration**: Stabilize temperature at 273K (100 ps)
5. **NPT Equilibration**: Adjust volume to reach target density at 1 bar (100 ps)
6. **Production Run**: 2 ns MD simulation in the NPT ensemble

### Analysis
After the simulation, the following analyses will be performed:
- Radial Distribution Functions (RDFs) for structural characterization
- Self-diffusion coefficient from Mean Square Displacement (MSD)
- Hydrogen bonding analysis (number and lifetime)
- Radial density profile to ensure uniform density
- Velocity Autocorrelation Function (VACF) analysis (if velocity data is sufficient)

## Running the Simulation

The simulation workflow can be executed using the following steps:

```bash
# 1. Generate the initial configuration with PACKMOL
packmol < configs/water_box.inp

# 2. Generate the topology with GROMACS (using TIP4P water model)
gmx pdb2gmx -f water_box.pdb -o water_box.gro -water tip4p -p topol.top

# 3. Perform energy minimization
gmx grompp -f configs/em.mdp -c water_box.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em

# 4. Perform NVT equilibration
gmx grompp -f configs/nvt.mdp -c em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

# 5. Perform NPT equilibration
gmx grompp -f configs/npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

# 6. Perform production MD
gmx grompp -f configs/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md
```

## Analysis Commands

Examples of analysis commands (to be run after the simulation):

```bash
# RDF Analysis (O-O)
gmx rdf -f md.xtc -s md.tpr -o rdf_OO.xvg -ref "name OW" -sel "name OW"

# MSD for diffusion coefficient
gmx msd -f md.xtc -s md.tpr -o msd.xvg -beginfit 1000 -endfit 2000

# Hydrogen bond analysis
gmx hbond -f md.xtc -s md.tpr -num hbonds.xvg -life hblife.xvg

# Velocity autocorrelation (requires velocity data)
gmx velacc -f md.trr -s md.tpr -o vacf.xvg -acflen 1000
```

For detailed scripts and additional analyses, see the `scripts/` directory. 