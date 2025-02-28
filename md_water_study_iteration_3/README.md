# TIP4P Water Model Study at 273K

A focused experimental framework for studying the TIP4P water model at 273K using molecular dynamics simulations with GROMACS.

## Project Overview

This project is a specialized version of the water model simulations focused exclusively on the TIP4P model at 273K. This iteration is designed to experiment with and refine the simulation parameters for optimal results.

## Directory Structure

- **scripts/**: Contains essential Python scripts for the TIP4P simulation workflow
- **configs/**: Contains configuration files for GROMACS
- **data/**: Stores simulation data for TIP4P at 273K
- **analysis/**: Stores analysis results for TIP4P at 273K

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

## Installation

1. Ensure GROMACS is installed and in your PATH
2. Make all scripts executable:
   ```bash
   chmod +x scripts/*.py
   ```

## Quick Start

To run a TIP4P water model simulation at 273K:

```bash
cd scripts
./run_simulation.py
```

## Water Model Details

The TIP4P model is a 4-site water model with:
- An oxygen atom (O) which is the center of the Lennard-Jones interaction
- Two hydrogen atoms (H) carrying partial positive charges
- A massless M-site carrying a negative charge, positioned near the oxygen atom

## Simulation Parameters

The default simulation parameters are:
- Temperature: 273K
- Number of water molecules: 1000
- Box size: 60.0 Å
- Tolerance: 3.0 Å
- Energy minimization steps: 100,000
- NVT equilibration time: 100 ps
- Production MD time: 1 ns

## Analysis Focus

For the TIP4P model at 273K, the analysis focuses on:
- Stability of the ice structure
- Radial distribution functions (O-O, O-H, H-H)
- Diffusion coefficients
- Density and structural properties

Results are stored in the `analysis/` directory.

## Experimental Notes

This iteration is focused on experimentation. Parameters may be adjusted during the process to achieve optimal results for the TIP4P model at freezing temperature (273K).

## References

1. Jorgensen, W. L., Chandrasekhar, J., Madura, J. D., Impey, R. W., & Klein, M. L. (1983). Comparison of simple potential functions for simulating liquid water. The Journal of chemical physics, 79(2), 926-935.
2. Abascal, J. L., Sanz, E., García Fernández, R., & Vega, C. (2005). A potential model for the study of ices and amorphous water: TIP4P/Ice. The Journal of chemical physics, 122(23), 234511. 