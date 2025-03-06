# TIP4P Water Model Study - Iteration 3

This directory contains a focused experimental framework for studying the TIP4P water model at 273K using molecular dynamics simulations with GROMACS. This iteration serves as a template and reference implementation for the centralized scripts system.

## Key Features

- Specialized simulation of TIP4P water model at freezing temperature (273K)
- Experimental framework for parameter optimization
- Organized directory structure for reproducible research
- Comprehensive analysis of water properties at low temperature
- Template for future water model studies

## Directory Structure

The project follows a standardized organization:

```
md_water_study_iteration_3/
├── analysis/                # Analysis results
│   └── tip4p/
│       └── 273K/
│           ├── plots/       # Generated plots
│           └── *.xvg        # Analysis data files
├── configs/                 # Template configuration files
│   ├── em.mdp               # Energy minimization config
│   ├── md.mdp               # Molecular dynamics config
│   ├── npt.mdp              # NPT equilibration config
│   ├── nvt.mdp              # NVT equilibration config
│   └── water_box.inp        # Water box generation config
├── data/                    # Simulation data
│   └── tip4p/
│       └── 273K/
│           ├── inputs/      # Input files (.mdp, .top, .inp, .pdb)
│           ├── outputs/     # Output files (.gro, .xtc, .trr, .edr, .tpr)
│           └── logs/        # Log files (.log)
├── docs/                    # Documentation
├── scripts/                 # Active scripts
├── unused_scripts/          # Deprecated scripts
├── unused_analysis/         # Deprecated analysis
└── README.md                # Project README
```

## TIP4P Water Model

The TIP4P model is a 4-site water model that places the negative charge on a dummy atom (M-site) rather than on the oxygen atom. This provides a better electrostatic distribution compared to 3-site models like SPC or TIP3P.

Key parameters of the TIP4P model:
- **Oxygen sigma**: 0.31644 nm
- **Oxygen epsilon**: 0.7749 kJ/mol
- **Hydrogen charges**: +0.52e
- **M-site charge**: -1.04e (placed on a dummy atom)
- **O-H bond length**: 0.09572 nm
- **H-O-H angle**: 104.52°

## Simulation Details

This iteration uses the following simulation parameters:

- **System size**: 1000 water molecules
- **Box size**: 60.0 Å
- **Temperature**: 273K (0°C)
- **Pressure**: 1 bar
- **Energy minimization**: 100,000 steps
- **NVT equilibration**: 100 ps
- **NPT equilibration**: 100 ps
- **Production run**: 500 ps (250,000 steps with 2 fs timestep)
- **Electrostatics**: PME with 1.0 nm cutoff
- **Thermostat**: V-rescale
- **Barostat**: Parrinello-Rahman

## Analysis Focus

For the TIP4P model at 273K, the analysis focuses on:

- **Structural properties**:
  - Radial distribution functions (O-O, O-H, H-H)
  - Coordination numbers
  - Density profiles

- **Dynamic properties**:
  - Diffusion coefficients
  - Mean square displacement

- **Thermodynamic properties**:
  - Energy components
  - Temperature stability
  - Pressure stability

- **Ice structure stability**:
  - Hydrogen bond network
  - Tetrahedral order parameter

## Workflow

The typical workflow for this project is:

1. Run the simulation using `run_workflow.py`
2. Analyze the results using `analyze_results.py`
3. Generate plots using `plot_all_results.py`
4. Clean up temporary files using `clean_temp_files.py`
5. Organize the project using `organize_project.py`

## Using the Centralized Scripts

This iteration can be used with the centralized scripts system:

```bash
cd /home/konrad_guest/Documents/research/cursor/md_water_study_scripts
./run_on_iteration.py run_workflow.py md_water_study_iteration_3 --model tip4p --temp 273
```

## Requirements

- GROMACS (version 2020 or newer)
- Python 3.6+
- Required Python packages:
  - argparse
  - pathlib
  - subprocess
  - os
  - shutil
  - time
  - numpy
  - matplotlib

## Experimental Notes

This iteration is focused on experimentation and serves as a template for future iterations. Key experimental aspects include:

1. **Parameter optimization**: Testing different simulation parameters to achieve optimal results for the TIP4P model at freezing temperature
2. **Directory structure**: Establishing a standardized directory structure for all water model studies
3. **Analysis pipeline**: Developing a comprehensive analysis pipeline for water properties
4. **Organization scripts**: Creating scripts to maintain a clean and organized project structure

## References

1. Abascal, J. L., Sanz, E., García Fernández, R., & Vega, C. (2005). A potential model for the study of ices and amorphous water: TIP4P/Ice. The Journal of chemical physics, 122(23), 234511.
