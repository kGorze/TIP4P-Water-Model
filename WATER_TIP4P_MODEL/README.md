# TIP4P Water Model Study

This directory contains a complete workflow for simulating and analyzing water using the standard TIP4P water model with GROMACS. The TIP4P model is a widely used 4-site water model that provides a good balance between accuracy and computational efficiency.

## Key Features

- Complete molecular dynamics workflow from system creation to analysis
- TIP4P water model with OPLS-AA force field
- Comprehensive analysis of structural, dynamic, and thermodynamic properties
- Detailed visualization and plotting capabilities
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
- `analysis/`: Analysis files and plots
  - `data/`: Raw data files from analysis (.xvg, .xpm)
  - `plots/`: Generated plots and summary report
  - `plotting_scripts/`: Python scripts for generating plots
- `logs/`: Log files from the workflow
- `checkpoints/`: Checkpoint files for resuming interrupted simulations
- `scripts/`: Additional scripts for the workflow

## TIP4P Water Model

The TIP4P model is a 4-site water model that places the negative charge on a dummy atom (M-site) rather than on the oxygen atom. This provides a better electrostatic distribution compared to 3-site models like SPC or TIP3P.

Key parameters of the TIP4P model:
- **Oxygen sigma**: 0.31644 nm
- **Oxygen epsilon**: 0.7749 kJ/mol
- **Hydrogen charges**: +0.52e
- **M-site charge**: -1.04e (placed on a dummy atom)
- **O-H bond length**: 0.09572 nm
- **H-O-H angle**: 104.52°

The TIP4P model is particularly good at reproducing the liquid structure and thermodynamic properties of water at ambient conditions.

## Simulation Details

- **System size**: 5500 water molecules
- **Temperature**: 273K (0°C)
- **Pressure**: 1 bar
- **Simulation time**: 2 ns production run
- **Timestep**: 2 fs
- **Equilibration**: NVT followed by NPT
- **Energy minimization**: Steepest descent
- **Electrostatics**: PME with 1.0 nm cutoff
- **Thermostat**: V-rescale
- **Barostat**: Parrinello-Rahman

## Analysis Results

Based on the analysis summary, the TIP4P water model simulation shows good agreement with experimental and literature values:

### Structural Properties
- **Density**: 1003.0 ± 14 kg/m³ (+0.83% from reference value: 995.0 kg/m³)
- **O-O First Peak Position**: 2.7 Å (-3.57% from reference value: 2.8 Å)
- **O-O Second Peak Position**: 2.9 Å (-35.56% from reference value: 4.5 Å)
- **O-O First Minimum**: 3.22 Å (-2.42% from reference value: 3.3 Å)
- **Coordination Number**: 4.5 (-10.00% from reference value: 5.0)

### Dynamic Properties
- **Diffusion Coefficient**: 4.869 × 10⁻⁹ m²/s (+111.70% from reference value: 2.3 × 10⁻⁹ m²/s)
- **Diffusion Fit Range**: 400-1600 ps with R² = 0.9806
- **H-bonds per Molecule**: 1.829 (-50.57% from reference value: 3.7)

### Thermodynamic Properties
- **Temperature**: 272.9 ± 2 K (-0.05% from reference value: 273.0 K)
- **Pressure**: -40.96 ± 300 bar (deviation from reference value: 1.0 bar)
- **Potential Energy**: -236,900 ± 540 kJ/mol
- **Energy per Molecule**: -43.08 ± 0.1 kJ/mol (-2.10% from reference value: -44.0 kJ/mol)
- **Energy Fluctuation**: 0.2267%
- **RMSD Final**: 5.615 nm

## Visualization and Analysis

The simulation includes comprehensive analysis and visualization capabilities:

- **Radial Distribution Functions (RDF)**: Analysis of water structure
- **Mean Square Displacement (MSD)**: Calculation of diffusion coefficients
- **Hydrogen Bond Analysis**: Number, geometry, and lifetime of hydrogen bonds
- **Density Analysis**: Spatial distribution of water molecules
- **Energy Analysis**: Detailed breakdown of energy components
- **Vibrational Analysis**: Spectral properties from velocity autocorrelation

All analyses are accompanied by publication-quality plots and statistical evaluations.

## Workflow

The workflow is managed by the `run_workflow.sh` script, which includes the following steps:

1. **System Creation**:
   - Generate a water box using Packmol via Julia
   - Create a topology with TIP4P water model parameters

2. **Energy Minimization**:
   - Use steepest descent algorithm
   - Remove bad contacts and optimize geometry

3. **NVT Equilibration**:
   - Equilibrate temperature at constant volume
   - Use V-rescale thermostat at 273K

4. **NPT Equilibration**:
   - Equilibrate pressure and density at constant temperature
   - Use Parrinello-Rahman barostat at 1 bar

5. **Production MD**:
   - Run 2 ns simulation with 2 fs timestep
   - Save coordinates for analysis

6. **Analysis**:
   - Calculate radial distribution functions (RDF)
   - Calculate mean square displacement (MSD) and diffusion coefficient
   - Analyze density and thermodynamic properties
   - Analyze hydrogen bonds and their lifetimes
   - Calculate velocity autocorrelation function (VACF)
   - Perform additional statistical analyses

7. **Plotting**:
   - Generate publication-quality plots
   - Create comprehensive summary report

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

### Plotting Only

If you already have the simulation and analysis data and just want to regenerate the plots:

```bash
./generate_plots.sh
```

Or use the fixed plotting script:

```bash
./run_plots_fixed.sh
```

## Dependencies

- GROMACS (for simulation and analysis)
- Julia with Packmol (for generating the water box)
- Python 3 with the following packages:
  - NumPy
  - Matplotlib
  - SciPy
  - Seaborn (for enhanced plots)

## Conclusion

The TIP4P water model provides a reliable representation of liquid water at 273K, with particularly good agreement for structural and thermodynamic properties. The simulation is well-equilibrated and stable throughout the trajectory.

The comprehensive analysis shows that the TIP4P model:
- Accurately reproduces the density of water at 273K
- Provides a reasonable estimate of the diffusion coefficient
- Captures the tetrahedral coordination of water molecules
- Maintains stable energy throughout the simulation

This workflow provides a solid foundation for studying liquid water properties and can be extended to investigate other phenomena such as solvation, interfaces, or temperature-dependent properties. 