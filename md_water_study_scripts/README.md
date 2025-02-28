# Centralized MD Water Study Scripts

This directory contains the centralized scripts for the MD Water Study project. The scripts are designed to work with any iteration directory and provide consistent functionality across all iterations.

## Usage

The scripts in this directory can be used with any iteration of the MD Water Study by passing the `--iteration-dir` parameter, which specifies the path to the iteration directory to work with.

### Running a Complete Workflow

To run a complete workflow on a specific iteration:

```bash
python run_workflow.py --model tip4p --temp 273 --iteration-dir /path/to/md_water_study_iteration_X
```

### Running Only the Simulation

To run only the simulation part:

```bash
python run_simulation.py --model tip4p --temp 273 --iteration-dir /path/to/md_water_study_iteration_X
```

### Running Analysis Only

To run only the analysis part:

```bash
python analyze_results.py --model tip4p --temp 273 --iteration-dir /path/to/md_water_study_iteration_X
```

### Generating Plots

To generate plots from analysis results:

```bash
python plot_results.py --model tip4p --temp 273 --iteration-dir /path/to/md_water_study_iteration_X
```

## Directory Structure

The scripts expect each iteration directory to have the following structure:

```
md_water_study_iteration_X/
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
└── unused_analysis/         # Deprecated analysis
└── unused_scripts/          # Deprecated scripts
```

## Important Note

These scripts always use the centralized script directory for executing commands, but work on the specified iteration directory for data and configuration. This ensures consistent script behavior across all iterations. 