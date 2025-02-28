# Molecular Dynamics Water Study

A comprehensive framework for studying water models using molecular dynamics simulations with GROMACS and PackMol.

## Project Overview

This project provides a modular workflow for setting up, running, and analyzing molecular dynamics simulations of water using different water models (TIP4P, TIP4P-ICE) at various temperatures. The workflow is designed to be flexible, allowing users to run the entire pipeline or individual steps as needed.

## Directory Structure

- **scripts/**: Contains all Python scripts for the simulation workflow
- **configs/**: Contains configuration files for GROMACS and PackMol
- **data/**: Stores simulation data organized by water model and temperature
- **analysis/**: Stores analysis results organized by water model

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

## Installation

1. Ensure GROMACS is installed and in your PATH
2. Ensure PackMol is installed and accessible
3. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/md_water_study.git
   cd md_water_study
   ```
4. Make all scripts executable:
   ```bash
   chmod +x scripts/*.py
   ```

## Quick Start

To run a complete simulation with the TIP4P water model at 273K:

```bash
cd scripts
./run_workflow.py --model tip4p --temp 273
```

For more detailed usage instructions, see the [Scripts README](scripts/README.md).

## Water Models

This project supports the following water models:

- **TIP4P**: A 4-site water model with a massless charge site (M-site)
- **TIP4P-ICE**: A variant of TIP4P optimized for ice simulations

## Simulation Parameters

The default simulation parameters are:
- Number of water molecules: 1000
- Box size: 60.0 Å
- Tolerance: 3.0 Å
- Energy minimization steps: 100,000
- NVT equilibration time: 100 ps
- Production MD time: 1 ns

These parameters can be adjusted in the configuration files or via command-line arguments.

## Analysis

The analysis script generates:
- Energy plots (potential, kinetic, total)
- RMSD plots
- Radial distribution functions (O-O, O-H, H-H)

Results are stored in the `analysis/<model>/` directory.

## Troubleshooting

If you encounter Python environment issues, try running with a clean environment:

```bash
env -i HOME=$HOME PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/path/to/julia/bin:/path/to/gromacs/bin python3 scripts/run_workflow.py --model tip4p --temp 273
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- GROMACS development team
- PackMol development team 