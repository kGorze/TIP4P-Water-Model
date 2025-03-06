# Molecular Dynamics Water Study - Iteration 2

This directory contains a comprehensive framework for studying water models using molecular dynamics simulations with GROMACS and PackMol. This iteration focuses on the TIP4P water model with a larger system size and longer simulation time compared to previous iterations.

## Key Features

- Comprehensive simulation framework for water models
- Modular workflow for setup, simulation, and analysis
- Support for multiple water models (TIP4P, TIP4P-ICE)
- Larger system size (4000 water molecules)
- Extended production run (2 ns)
- Detailed analysis of water properties

## Directory Structure

The project follows a standardized organization:

```
md_water_study_iteration_2/
├── analysis/                # Analysis results
│   ├── tip4p_simulation_analysis.md             # Detailed analysis report
│   └── tip4p_simulation_analysis_enhanced.md    # Enhanced analysis report
├── configs/                 # Template configuration files
│   ├── em.mdp               # Energy minimization config
│   ├── md.mdp               # Molecular dynamics config
│   ├── nvt.mdp              # NVT equilibration config
│   └── water_box.inp        # Water box generation config
├── data/                    # Simulation data
│   └── tip4p/
│       └── 273K/
│           ├── inputs/      # Input files (.mdp, .top, .inp, .pdb)
│           ├── outputs/     # Output files (.gro, .xtc, .trr, .edr, .tpr)
│           └── logs/        # Log files (.log)
├── docs/                    # Documentation
├── scripts/                 # Python scripts for the workflow
└── README.md                # Project README
```

## Water Models

This iteration supports the following water models:

### TIP4P Model
- **Description**: 4-site model with negative charge on dummy atom (M-site)
- **Key Parameters**:
  - Oxygen sigma: 0.31644 nm
  - Oxygen epsilon: 0.7749 kJ/mol
  - Hydrogen charges: +0.52e
  - M-site charge: -1.04e
  - O-H bond length: 0.09572 nm
  - H-O-H angle: 104.52°

### TIP4P-ICE Model
- **Description**: Modified TIP4P optimized for ice simulations
- **Key Parameters**:
  - Oxygen sigma: 0.31668 nm
  - Oxygen epsilon: 0.88211 kJ/mol
  - Hydrogen charges: +0.5897e
  - M-site charge: -1.1794e

## Simulation Details

This iteration uses the following simulation parameters:

- **System size**: 4000 water molecules
- **Box size**: 70.0 Å
- **Density**: 0.93 g/cm³
- **Temperature**: 273K (0°C)
- **Pressure**: 1 bar
- **Energy minimization**: 100,000 steps
- **NVT equilibration**: 100 ps
- **Production run**: 2 ns (1,000,000 steps with 2 fs timestep)
- **Output frequency**: Every 10 ps
- **Electrostatics**: PME with 1.0 nm cutoff
- **Thermostat**: V-rescale
- **Barostat**: Parrinello-Rahman

## Analysis Approach

The analysis in this iteration is particularly comprehensive, with detailed reports available in the analysis directory:

- **Energy analysis**: Potential, kinetic, and total energy evolution
- **Temperature analysis**: Temperature stability and fluctuations
- **Diffusion analysis**: Mean square displacement and diffusion coefficients
- **Thermodynamic analysis**: Pressure, volume, and density
- **Structural analysis**: Radial distribution functions (O-O, O-H, H-H)
- **Comparison with experimental data**: Validation against known water properties

The analysis is presented in both standard and enhanced formats, with the enhanced version including additional statistical analysis and visualizations.

## Workflow

The typical workflow for this project is:

1. Generate a water box using PackMol
2. Create topology with the selected water model
3. Perform energy minimization
4. Run NVT equilibration
5. Run production MD
6. Analyze the results

## Using the Scripts

To run a complete simulation with the TIP4P water model at 273K:

```bash
cd scripts
./run_workflow.py --model tip4p --temp 273
```

For individual steps:

```bash
# Generate water box only
./generate_water_box.py --model tip4p --temp 273

# Run energy minimization only
./run_em.py --model tip4p --temp 273

# Run analysis only
./analyze_results.py --model tip4p --temp 273
```

## Requirements

- GROMACS (version 2020 or newer)
- PackMol (with Julia)
- Python 3.6+
- Required Python packages:
  - argparse
  - pathlib
  - subprocess
  - os
  - shutil
  - time

## Key Improvements from Previous Iterations

This iteration includes several improvements over previous versions:

1. **Larger system size**: 4000 water molecules vs. 1000 in previous iterations
2. **Longer simulation time**: 2 ns vs. 1 ns in previous iterations
3. **More detailed analysis**: Comprehensive physical analysis in multiple formats
4. **Improved organization**: Clearer directory structure and file naming
5. **Support for multiple water models**: Framework extended to support TIP4P-ICE

## Troubleshooting

If you encounter Python environment issues, try running with a clean environment:

```bash
env -i HOME=$HOME PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/path/to/julia/bin:/path/to/gromacs/bin python3 scripts/run_workflow.py --model tip4p --temp 273
```

## References

1. Abascal, J. L., Sanz, E., García Fernández, R., & Vega, C. (2005). A potential model for the study of ices and amorphous water: TIP4P/Ice. The Journal of chemical physics, 122(23), 234511.
