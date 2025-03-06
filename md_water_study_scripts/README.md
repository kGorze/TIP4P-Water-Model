# MD Water Study Centralized Scripts

This directory contains the centralized scripts for the Molecular Dynamics Water Study project. These scripts provide a consistent framework for simulating and analyzing different water models across all project iterations.

## Overview

The centralized scripts system offers several key advantages:

1. **Consistency**: All iterations use the same script versions, ensuring reproducible results
2. **Maintainability**: Bug fixes and improvements only need to be made in one place
3. **Efficiency**: No need to recreate scripts for each new iteration
4. **Flexibility**: Can specify which iteration to work with for any script

## Script Categories

The scripts are organized into several categories:

### Core Workflow Scripts

- **run_workflow.py**: Complete workflow from system creation to analysis and plotting
- **run_simulation.py**: System creation and simulation only
- **analyze_results.py**: Analysis of simulation results
- **plot_results.py**: Generation of plots from analysis data

### Utility Scripts

- **run_on_iteration.py**: Utility to run any script on any iteration
- **setup_iteration.py**: Create a new iteration directory with the correct structure
- **compare_iterations.py**: Compare results across multiple iterations

### Analysis Scripts

- **rdf_analysis.py**: Radial distribution function analysis
- **msd_analysis.py**: Mean square displacement and diffusion coefficient
- **density_analysis.py**: Density profile and distribution
- **energy_analysis.py**: Energy components and fluctuations
- **hbond_analysis.py**: Hydrogen bond analysis
- **vacf_analysis.py**: Velocity autocorrelation function

## Using the Scripts

### The run_on_iteration.py Utility

The recommended way to use these scripts is through the `run_on_iteration.py` utility, which handles the correct passing of parameters to the underlying scripts:

```bash
./run_on_iteration.py <script_name> <iteration_name> [additional arguments]
```

Examples:

```bash
# Run the complete workflow on the TIP4P model iteration
./run_on_iteration.py run_workflow.py WATER_TIP4P_MODEL --model tip4p --temp 273

# Run only the analysis on the TIP4P/Ice model iteration
./run_on_iteration.py analyze_results.py WATER_TIP4PICE_MODEL --model tip4pice --temp 273

# Generate plots for the TIP4P model iteration
./run_on_iteration.py plot_results.py WATER_TIP4P_MODEL --model tip4p --temp 273
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
./run_on_iteration.py analyze_results.py all --model tip4p --temp 273
```

### Direct Script Usage

You can also use the scripts directly by specifying the iteration directory:

```bash
python run_workflow.py --model tip4p --temp 273 --iteration-dir /path/to/WATER_TIP4P_MODEL
```

## Expected Directory Structure

The scripts expect each iteration directory to have the following structure:

```
iteration_directory/
├── analysis/                # Analysis results
│   ├── data/                # Raw analysis data (.xvg, .xpm)
│   ├── plots/               # Generated plots
│   └── plotting_scripts/    # Python scripts for plotting
├── configs/                 # Template configuration files
│   ├── em.mdp               # Energy minimization config
│   ├── md.mdp               # Molecular dynamics config
│   ├── npt.mdp              # NPT equilibration config
│   ├── nvt.mdp              # NVT equilibration config
│   └── water_box.inp        # Water box generation config
├── data/                    # Simulation data
│   ├── Coordinate files (.gro, .pdb)
│   ├── Trajectory files (.trr, .xtc)
│   ├── Energy files (.edr)
│   ├── Log files (.log)
│   └── Checkpoint files (.cpt)
├── logs/                    # Log files
├── checkpoints/             # Checkpoint files for resuming simulations
└── scripts/                 # Iteration-specific scripts
```

## Script Parameters

Most scripts accept the following common parameters:

- `--model`: Water model to use (e.g., tip3p, tip4p, tip4pice, spce)
- `--temp`: Temperature in Kelvin (e.g., 273, 298, 310)
- `--pressure`: Pressure in bar (default: 1.0)
- `--iteration-dir`: Path to the iteration directory
- `--verbose`: Enable verbose output
- `--help`: Show help message

## Creating a New Iteration

To create a new iteration directory with the correct structure:

```bash
./setup_iteration.py create --name NEW_ITERATION_NAME --base-model tip4p
```

This will create a new iteration directory with all necessary subdirectories and template files.

## Extending the Scripts

The centralized scripts are designed to be extensible. To add a new analysis type:

1. Create a new script in this directory (e.g., `new_analysis.py`)
2. Follow the pattern of existing analysis scripts
3. Use the common parameter handling functions
4. Register the script in `run_on_iteration.py` if needed

## Troubleshooting

If you encounter issues with the scripts:

1. Check that the iteration directory has the expected structure
2. Ensure GROMACS and other dependencies are properly installed
3. Check the log files in the iteration's `logs/` directory
4. Try running with the `--verbose` flag for more detailed output

## Dependencies

These scripts depend on:

- Python 3.6+
- NumPy, Matplotlib, SciPy
- GROMACS (accessible in PATH)
- Julia with Packmol (for system creation)

## Contributing

When contributing to these scripts, please:

1. Follow the existing code style
2. Add comprehensive docstrings
3. Test your changes on multiple iterations
4. Update this README if you add new functionality 