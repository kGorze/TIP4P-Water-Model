# Molecular Dynamics Water Study

This project is an investigation of water properties using molecular dynamics simulations with the TIP4P water model at different temperatures. It provides a comprehensive analysis system for comparing water properties across different conditions.

## Project Overview

The MD Water Study project aims to:

1. **Analyze water properties** at different temperatures (273K, 298K)
2. **Provide reproducible workflows** for water simulations
3. **Generate quality analyses** and visualizations
4. **Establish temperature-dependent** water properties

## Project Structure

The project consists of:

1. `tip4p_273K/` - TIP4P water model at 273K
2. `tip4p_298K/` - TIP4P water model at 298K
3. `analysis/` - Analysis scripts and results
4. `visualization_scripts/` - Visualization tools
5. `visualization_outputs/` - Generated visualizations
6. `comparison/` - Comparative analysis results

Each simulation directory follows a consistent structure:

```
simulation_directory/
├── analysis/                # Analysis results
├── scripts/                 # Analysis scripts
├── md.mdp                   # Production MD parameters
├── nvt.mdp                  # NVT equilibration parameters
├── em.mdp                   # Energy minimization parameters
├── topol.top               # System topology
├── box_5500.gro            # Initial system coordinates
├── md.gro                   # Final system coordinates
├── md.xtc                   # Production trajectory
├── md.trr                   # Production trajectory (full precision)
├── md.edr                   # Energy data
├── md.log                   # MD log file
└── various checkpoint files
```

## Water Model Details

### TIP4P Model
- **Description**: 4-site model with negative charge on dummy atom (M-site)
- **Key Parameters**:
  - Oxygen sigma: 0.31644 nm
  - Oxygen epsilon: 0.7749 kJ/mol
  - Hydrogen charges: +0.52e
  - M-site charge: -1.04e

## Simulation Methodology

All simulations follow a consistent protocol:

1. **System Creation**:
   - Generate a water box with 5500 water molecules
   - Create topology with TIP4P water model parameters

2. **Energy Minimization**:
   - Steepest descent algorithm to remove bad contacts
   - Careful minimization with small step size

3. **Equilibration**:
   - NVT equilibration to stabilize temperature
   - NPT equilibration to stabilize pressure (1 bar) and density

4. **Production**:
   - Production run with 1-2 fs timestep
   - Periodic boundary conditions
   - PME electrostatics with 1.0 nm cutoff
   - V-rescale thermostat and Parrinello-Rahman barostat

## Analysis Capabilities

The project includes comprehensive analysis tools:

1. **Structural Analysis**:
   - RMSD analysis and drift visualization
   - Radial distribution functions
   - Coordination numbers

2. **Dynamic Properties**:
   - Mean square displacement
   - Diffusion coefficients
   - Velocity autocorrelation functions

3. **Thermodynamic Properties**:
   - Temperature and pressure profiles
   - Energy components
   - Density calculations

4. **Comparative Analysis**:
   - Temperature-dependent comparisons
   - Pressure analysis
   - Energy comparisons
   - RMSD comparisons

## Visualization Tools

The project includes advanced visualization capabilities:

1. **Trajectory Visualization**:
   - RMSD drift visualization
   - Combined RMSD plots
   - Temperature and pressure profiles

2. **Analysis Visualization**:
   - Energy component plots
   - Structural property plots
   - Dynamic property plots

## Dependencies

- **GROMACS**: For simulation and analysis
- **Python 3** with:
  - NumPy
  - Matplotlib
  - SciPy
  - Seaborn

## Usage

### Running Simulations

```bash
cd tip4p_273K  # or tip4p_298K
gmx mdrun -v -deffnm md
```

### Running Analysis

```bash
python analysis/rmsd_drift_visual.py
python analysis/energy_comparison.py
python analysis/temperature_comparison.py
python analysis/pressure_comparison.py
```

### Generating Visualizations

```bash
python visualization_scripts/plot_combined_rmsd.py
python visualization_scripts/rmsd_comparison.py
```

## Future Work

Planned extensions to this project include:

1. Additional temperature points
2. Pressure dependence studies
3. Interface and solvation studies
4. Machine learning analysis of water structure and dynamics
