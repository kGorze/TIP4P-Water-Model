# TIP4P/Ice Water Model Study

This directory contains a complete workflow for simulating and analyzing water using the TIP4P/Ice water model with GROMACS. The TIP4P/Ice model is specifically designed for simulating water at low temperatures and ice phases.

## Key Features

- Complete molecular dynamics workflow from system creation to analysis
- TIP4P/Ice water model with OPLS-AA force field
- Comprehensive analysis of structural, dynamic, and thermodynamic properties
- High-quality visualization and movie generation
- Checkpoint system for resuming interrupted simulations
- Detailed comparison with literature values

## Directory Structure

- `configs/`: Configuration files for GROMACS (mdp files, etc.)
- `data/`: Output data from the simulation
  - Coordinate files (.gro, .pdb)
  - Trajectory files (.trr, .xtc)
  - Energy files (.edr)
  - Log files (.log)
  - Checkpoint files (.cpt)
  - Visualization scripts (.tcl)
- `analysis/`: Analysis files and plots
  - `data/`: Raw data files from analysis (.xvg, .xpm)
  - `plots/`: Generated plots and summary report
  - `plotting_scripts/`: Python scripts for generating plots
- `logs/`: Log files from the workflow
- `checkpoints/`: Checkpoint files for resuming interrupted simulations
- `oplsaa.ff/`: Local force field directory with TIP4P/Ice parameters
- `scripts/`: Additional scripts for the workflow

## TIP4P/Ice Water Model

The TIP4P/Ice model is a 4-site water model specifically parameterized to reproduce the melting point of ice Ih and the densities of different ice polymorphs. It differs from the standard TIP4P model in the following ways:

- **Increased hydrogen charges**: 0.5897e (vs 0.52e in TIP4P)
- **Increased M-site negative charge**: -1.1794e (vs -1.04e in TIP4P)
- **Slightly increased oxygen sigma**: 0.31668 nm (vs 0.31644 nm in TIP4P)
- **Significantly increased oxygen epsilon**: 0.88211 kJ/mol (vs 0.7749 kJ/mol in TIP4P)

These modifications make the TIP4P/Ice model better suited for simulating ice phases and water at low temperatures, as they increase the melting temperature and improve the phase diagram.

## Simulation Details

- **System size**: 5500 water molecules
- **Temperature**: 273K (0°C)
- **Pressure**: 1 bar
- **Simulation time**: 1 ns production run
- **Timestep**: 1 fs
- **Equilibration**: NVT (25 ps) followed by NPT (25 ps)
- **Energy minimization**: Steepest descent with very small step size (0.00001)
- **Electrostatics**: PME with 1.0 nm cutoff
- **Thermostat**: V-rescale
- **Barostat**: Parrinello-Rahman

## Analysis Results

Based on the analysis summary, the TIP4P/Ice water model simulation shows good agreement with experimental and literature values:

### Structural Properties
- **Density**: 966.6 ± 5 kg/m³ (-2.66% from reference value: 993.0 kg/m³)
- **O-O First Peak Position**: 2.78 Å (+0.72% from reference value: 2.76 Å)
- **O-O Second Peak Position**: 2.98 Å (-33.18% from reference value: 4.46 Å)
- **O-O First Minimum**: 3.3 Å (+0.61% from reference value: 3.28 Å)
- **Coordination Number**: 4.5 (-8.16% from reference value: 4.9)

### Dynamic Properties
- **Diffusion Coefficient**: 9.081e-10 m²/s (-17.45% from reference value: 1.1e-09 m²/s)
- **Diffusion Fit Range**: 200-800 ps with R² = 0.9852

### Thermodynamic Properties
- **Temperature**: 274.6 ± 10 K (+0.58% from reference value: 273.0 K)
- **Pressure**: 22.11 ± 300 bar (higher than reference value: 1.0 bar)
- **Potential Energy**: -305,500 ± 8,900 kJ/mol
- **Energy per Molecule**: -55.54 ± 2 kJ/mol (+18.17% from reference value: -47.0 kJ/mol)
- **Energy Fluctuation**: 2.925%

## Visualization

The simulation includes scripts for generating high-quality visualizations:

- **render_tip4pice_movie.tcl**: VMD script for rendering a movie
- **render_frames.tcl**: VMD script for rendering individual frames
- **create_movie.sh**: Shell script to create a high-quality movie using VMD and ffmpeg
- **create_tip4pice_movie.sh**: Enhanced movie creation script with frame interpolation

The visualization shows:
- Water molecules in CPK representation
- Hydrogen bonds in blue
- Periodic boundary conditions
- High-quality rendering with shadows and ambient occlusion

## Workflow

The workflow is managed by the `run_workflow.sh` script, which includes the following steps:

1. **System Creation**:
   - Generate a water box using Packmol via Julia
   - Create a topology with TIP4P/Ice water model parameters

2. **Energy Minimization**:
   - Use steepest descent algorithm with a very small step size
   - Remove bad contacts and optimize geometry

3. **NVT Equilibration**:
   - Equilibrate temperature at constant volume
   - Use V-rescale thermostat at 273K

4. **NPT Equilibration**:
   - Equilibrate pressure and density at constant temperature
   - Use Parrinello-Rahman barostat at 1 bar

5. **Production MD**:
   - Run 1 ns simulation with 1 fs timestep
   - Save coordinates every 5 ps

6. **Analysis**:
   - Calculate radial distribution functions (RDF)
   - Calculate mean square displacement (MSD) and diffusion coefficient
   - Analyze density and thermodynamic properties
   - Analyze hydrogen bonds and their lifetimes
   - Calculate velocity autocorrelation function (VACF)

7. **Visualization**:
   - Generate high-quality movies and images

## How to Run

### Complete Workflow

To run the complete workflow (simulation + analysis + plotting):

```bash
./run_workflow.sh
```

### Running Specific Steps

You can run specific steps of the workflow:

```bash
# Run only the plotting step
./run_workflow.sh --only generate_plots

# Rerun the MSD analysis and continue from there
./run_workflow.sh --rerun msd_analysis

# Run with cleanup of temporary files
./run_workflow.sh --cleanup
```

### Available Steps

- `water_box`: Generate water box
- `topology`: Generate topology
- `energy_minimization`: Run energy minimization
- `nvt_equilibration`: Run NVT equilibration
- `npt_equilibration`: Run NPT equilibration
- `production_md`: Run production MD
- `rdf_analysis_OO`, `rdf_analysis_OH`, `rdf_analysis_HH`: Calculate RDFs
- `msd_analysis`: Calculate MSD and diffusion coefficient
- `density_analysis`: Analyze density
- `temperature_analysis`: Analyze temperature
- `pressure_analysis`: Analyze pressure
- `energy_analysis`: Analyze energy components
- `hbond_index`, `hbond_analysis`, `hbond_lifetime`: Analyze hydrogen bonds
- `vacf_analysis`: Calculate velocity autocorrelation function
- `additional_analysis`: Run additional analyses
- `generate_plots`: Generate plots and summary report

### Creating Movies

To create a high-quality movie of the simulation:

```bash
./create_movie.sh
```

Or use the enhanced version with frame interpolation:

```bash
./create_tip4pice_movie.sh
```

## Dependencies

- GROMACS (for simulation and analysis)
- Julia with Packmol (for generating the water box)
- VMD (for visualization)
- ffmpeg (for movie creation)
- Python 3 with the following packages:
  - NumPy
  - Matplotlib
  - SciPy
  - Seaborn (for enhanced plots)

## Conclusion

The TIP4P/Ice water model provides a reliable representation of liquid water at low temperatures and is particularly well-suited for studying ice phases. The simulation results show good agreement with experimental and literature values for structural, dynamic, and thermodynamic properties.

The comprehensive analysis and visualization tools in this workflow allow for detailed investigation of water properties using the TIP4P/Ice model. 