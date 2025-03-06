# Molecular Dynamics Water Study

This project is a comprehensive investigation of water properties using molecular dynamics simulations with different water models. It provides a centralized scripts system to ensure consistent methodology across multiple iterations, allowing for systematic comparison of different water models and simulation conditions.

## Project Overview

The MD Water Study project aims to:

1. **Compare different water models** (TIP3P, TIP4P, TIP4P/Ice, SPC/E, etc.)
2. **Analyze water properties** at different temperatures and pressures
3. **Provide reproducible workflows** for water simulations
4. **Generate publication-quality analyses** and visualizations
5. **Establish best practices** for water model selection based on the properties of interest

## Project Structure

The project consists of:

1. `md_water_study_scripts/` - Central directory containing all shared scripts
2. `md_water_study_iteration_X/` - Individual iteration directories for different water models/conditions
3. `md_water_study_iteration_result/` - Directory for combined results from all iterations

Current iterations include:

- **WATER_TIP4P_MODEL/** - Standard TIP4P water model at 273K
- **WATER_TIP4PICE_MODEL/** - TIP4P/Ice water model at 273K (optimized for ice phases)

Each iteration directory follows a consistent structure:

```
iteration_directory/
├── analysis/                # Analysis results
│   ├── data/                # Raw analysis data (.xvg, .xpm)
│   ├── plots/               # Generated plots
│   └── plotting_scripts/    # Python scripts for plotting
├── configs/                 # Template configuration files
├── data/                    # Simulation data
│   ├── Coordinate files (.gro, .pdb)
│   ├── Trajectory files (.trr, .xtc)
│   ├── Energy files (.edr)
│   ├── Log files (.log)
│   └── Checkpoint files (.cpt)
├── docs/                    # Documentation
├── logs/                    # Log files
├── checkpoints/             # Checkpoint files for resuming simulations
├── scripts/                 # Iteration-specific scripts
└── README.md                # Iteration-specific documentation
```

## Water Models Comparison

The project currently includes detailed studies of the following water models:

### TIP4P Model
- **Description**: 4-site model with negative charge on dummy atom (M-site)
- **Key Parameters**:
  - Oxygen sigma: 0.31644 nm
  - Oxygen epsilon: 0.7749 kJ/mol
  - Hydrogen charges: +0.52e
  - M-site charge: -1.04e
- **Strengths**: Good for liquid water at ambient conditions
- **Results Highlights**:
  - Density: 1003.0 ± 14 kg/m³ (+0.83% from reference)
  - Diffusion: 4.869 × 10⁻⁹ m²/s (+111.70% from reference)
  - Energy per molecule: -43.08 ± 0.1 kJ/mol (-2.10% from reference)

### TIP4P/Ice Model
- **Description**: Modified TIP4P optimized for ice phases and low temperatures
- **Key Parameters**:
  - Oxygen sigma: 0.31668 nm
  - Oxygen epsilon: 0.88211 kJ/mol
  - Hydrogen charges: +0.5897e
  - M-site charge: -1.1794e
- **Strengths**: Better for ice phases and water at low temperatures
- **Results Highlights**:
  - Density: 966.6 ± 5 kg/m³ (-2.66% from reference)
  - Diffusion: 9.081 × 10⁻¹⁰ m²/s (-17.45% from reference)
  - Energy per molecule: -55.54 ± 2 kJ/mol (+18.17% from reference)

## Simulation Methodology

All iterations follow a consistent simulation protocol:

1. **System Creation**:
   - Generate a water box with 5500 water molecules using Packmol
   - Create topology with appropriate water model parameters

2. **Energy Minimization**:
   - Steepest descent algorithm to remove bad contacts
   - Careful minimization with small step size

3. **Equilibration**:
   - NVT equilibration to stabilize temperature (273K)
   - NPT equilibration to stabilize pressure (1 bar) and density

4. **Production**:
   - 1-2 ns production run with 1-2 fs timestep
   - Periodic boundary conditions
   - PME electrostatics with 1.0 nm cutoff
   - V-rescale thermostat and Parrinello-Rahman barostat

5. **Analysis**:
   - Structural properties (RDF, coordination numbers)
   - Dynamic properties (MSD, diffusion coefficient)
   - Thermodynamic properties (energy, temperature, pressure)
   - Hydrogen bonding analysis
   - Vibrational analysis (VACF)

## Using the Centralized Scripts

### Running Scripts

The central scripts directory provides a utility to run any script on any iteration:

```bash
cd /home/konrad_guest/Documents/research/cursor/md_water_study_scripts
./run_on_iteration.py script_name.py iteration_name [additional arguments]
```

For example:
```bash
# Run the workflow on TIP4P model
./run_on_iteration.py run_workflow.py WATER_TIP4P_MODEL --model tip4p --temp 273

# Run analysis on TIP4P/Ice model
./run_on_iteration.py analyze_results.py WATER_TIP4PICE_MODEL --model tip4pice --temp 273
```

### Listing Available Scripts and Iterations

To list all available scripts:
```bash
./run_on_iteration.py list any_iteration
```

To list all available iterations:
```bash
./run_on_iteration.py any_script list
```

### Running a Script on All Iterations

To run a script on all iterations:
```bash
./run_on_iteration.py script_name.py all [additional arguments]
```

## Individual Iteration Workflows

Each iteration also includes its own workflow script that can be run independently:

```bash
cd /home/konrad_guest/Documents/research/cursor/WATER_TIP4P_MODEL
./run_workflow.sh

# Or for specific steps
./run_workflow.sh --only generate_plots
./run_workflow.sh --rerun msd_analysis
```

## Visualization Capabilities

The project includes advanced visualization capabilities:

- **Publication-quality plots** for all analyses
- **High-quality movies** of water dynamics (TIP4P/Ice model)
- **Comprehensive summary reports** with statistical analysis
- **Radial density maps** and other specialized visualizations

## Key Findings

Comparative analysis of the water models reveals:

1. **Density**: 
   - TIP4P: Slightly overestimates density at 273K (+0.83%)
   - TIP4P/Ice: Slightly underestimates density at 273K (-2.66%)

2. **Diffusion**:
   - TIP4P: Significantly overestimates diffusion (+111.70%)
   - TIP4P/Ice: Moderately underestimates diffusion (-17.45%)

3. **Structure**:
   - Both models capture the tetrahedral coordination of water
   - TIP4P/Ice better reproduces the first minimum in the O-O RDF

4. **Energy**:
   - TIP4P: Energy per molecule close to reference (-2.10%)
   - TIP4P/Ice: Higher energy per molecule (+18.17%)

5. **Hydrogen Bonding**:
   - TIP4P/Ice provides better hydrogen bond geometry and lifetimes
   - TIP4P underestimates the number of hydrogen bonds per molecule

## Dependencies

- **GROMACS**: For simulation and analysis
- **Julia with Packmol**: For generating water boxes
- **VMD**: For visualization (especially for TIP4P/Ice)
- **ffmpeg**: For movie creation
- **Python 3** with:
  - NumPy
  - Matplotlib
  - SciPy
  - Seaborn

## Future Work

Planned extensions to this project include:

1. Additional water models (SPC/E, OPC, TIP5P)
2. Temperature dependence studies (273K to 373K)
3. Pressure dependence studies
4. Interface and solvation studies
5. Machine learning analysis of water structure and dynamics
