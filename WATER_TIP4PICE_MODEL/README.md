# TIP4P/Ice Water Model Study - Iteration 4

This directory contains a complete workflow for simulating and analyzing TIP4P/Ice water using GROMACS.

## Directory Structure

- `configs/`: Configuration files for GROMACS (mdp files, etc.)
- `data/`: Output data from the simulation
- `analysis/`: Analysis files and plots
  - `data/`: Raw data files from analysis (.xvg, .xpm)
  - `plotting_scripts/`: Python scripts for generating plots
  - `plots/`: Generated plots and summary report
- `logs/`: Log files from the workflow
- `scripts/`: Additional scripts for the workflow
- `oplsaa.ff/`: Local force field directory with TIP4P/Ice parameters

## Implementation Note

This workflow uses the TIP4P/Ice water model with the OPLS-AA force field. The TIP4P/Ice parameters are directly incorporated into the topology file during the workflow execution, ensuring that all simulations use the correct water model parameters.

The TIP4P/Ice model differs from the standard TIP4P model in the following ways:
- Increased hydrogen charges (0.5897 vs 0.52)
- Increased M-site negative charge (-1.1794 vs -1.04)
- Slightly increased oxygen sigma (0.31668 nm vs 0.31644 nm)
- Significantly increased oxygen epsilon (0.88211 kJ/mol vs 0.7749 kJ/mol)

These modifications make the TIP4P/Ice model better suited for simulating ice phases and water at low temperatures, as they increase the melting temperature and improve the phase diagram.

## Workflow Overview

The workflow consists of the following steps:

1. Generate a water box using Packmol via Julia
2. Generate topology with TIP4P/Ice water model and OPLS-AA force field
3. Perform energy minimization
4. Perform NVT equilibration
5. Perform NPT equilibration
6. Perform production MD
7. Run analysis (RDF, MSD, density, thermodynamic properties, hydrogen bonds, etc.)
8. Generate plots and summary report

## How to Run

### Complete Workflow

To run the complete workflow (simulation + analysis + plotting):

```bash
./run_workflow.sh
```

This will run the entire workflow from start to finish, including the simulation, analysis, and plotting.

### Plotting Only

If you already have the simulation and analysis data and just want to regenerate the plots:

```bash
./generate_plots.sh
```

This will run only the plotting scripts and generate the summary report.

### Running Individual Plotting Scripts

You can also run individual plotting scripts directly:

```bash
cd analysis/plotting_scripts
python3 plot_rdf.py ../path/to/analysis ../path/to/plots
```

Or use the coordinator script with specific directories:

```bash
cd analysis/plotting_scripts
python3 run_all_plots.py --analysis-dir ../path/to/analysis --plots-dir ../path/to/plots --verbose
```

## Output

The main outputs of the workflow are:

- Simulation data in the `data/` directory
- Analysis files in the `analysis/` directory
- Plots in the `analysis/plots/` directory
- Summary report in `analysis/plots/tip4pice_water_analysis_summary.png` and `analysis/plots/tip4pice_water_analysis_summary.txt`
- Log files in the `logs/` directory

## Dependencies

- GROMACS (for simulation and analysis)
- Julia with Packmol (for generating the water box)
- Python 3 with the following packages:
  - NumPy
  - Matplotlib
  - SciPy
  - Seaborn (optional, for enhanced plots)

## Notes

- The simulation uses the TIP4P/Ice water model with the OPLS-AA force field
- The production run is set to 2 ns by default
- The analysis includes RDF, MSD, density, thermodynamic properties, hydrogen bonds, and more
- The plotting scripts generate publication-quality plots with statistical analysis
- The summary report provides a comprehensive overview of all analyses 