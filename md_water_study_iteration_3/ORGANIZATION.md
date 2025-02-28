# Project Organization Guide

This document explains the organization of the TIP4P Water Model Study project and provides instructions for maintaining a clean and organized project structure.

## Directory Structure

The project is organized as follows:

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
│   ├── analyze_results.py   # Analysis script
│   ├── check_environment.py # Environment check script
│   ├── cleanup.py           # Cleanup script
│   ├── clean_temp_files.py  # Temporary file cleanup script
│   ├── organize_project.py  # Project organization script
│   ├── plot_all_results.py  # Plotting script
│   ├── run_simulation.py    # Simulation script
│   └── run_workflow.py      # Workflow script
├── unused_scripts/          # Deprecated scripts
├── unused_analysis/         # Deprecated analysis
└── README.md                # Project README
```

## Organization Scripts

Two scripts have been created to help maintain the project organization:

### 1. organize_project.py

This script organizes the project directory by:
- Creating organized directories if they don't exist
- Moving files to appropriate directories based on file type

Usage:
```bash
cd scripts
python organize_project.py --dry-run  # Show what would be done without making changes
python organize_project.py            # Actually perform the organization
```

### 2. clean_temp_files.py

This script cleans up temporary files from simulation directories:
- Removes intermediate step files
- Removes GROMACS backup files
- Removes checkpoint files (except the latest ones)
- Removes generated mdp files
- Removes PackMol input files
- Removes template water molecule files

Usage:
```bash
cd scripts
python clean_temp_files.py --model tip4p --temp 273 --dry-run --verbose  # Show what would be deleted
python clean_temp_files.py --model tip4p --temp 273                      # Actually delete files
```

## File Categories

Files are categorized as follows:

### Input Files
- `.mdp`: GROMACS parameter files
- `.top`: Topology files
- `.inp`: Input files for external programs
- `.pdb`: Protein Data Bank files (initial structures)

### Output Files
- `.gro`: GROMACS structure files
- `.xtc`: Compressed trajectory files
- `.trr`: Full precision trajectory files
- `.edr`: Energy files
- `.tpr`: Run input files
- `.cpt`: Checkpoint files

### Log Files
- `.log`: Log files from simulations

### Analysis Files
- `.xvg`: XVG data files for plotting

## Important Files to Keep

The following files are considered important and should not be deleted:

- `md.edr`: Final energy file
- `md.xtc`: Final trajectory
- `md.trr`: Final trajectory (uncompressed)
- `md.gro`: Final structure
- `md.log`: Final log
- `md.cpt`: Final checkpoint
- `md_prev.cpt`: Previous checkpoint
- `nvt.log`: NVT log
- `em.log`: Energy minimization log
- `topol.top`: Topology
- `processed.gro`: Processed structure
- All `.mdp` files: Configuration files
- All `.xvg` files: Analysis output
- `water_box.pdb`: Initial water box

## Workflow

The typical workflow for this project is:

1. Run the simulation using `run_workflow.py`
2. Analyze the results using `analyze_results.py`
3. Generate plots using `plot_all_results.py`
4. Clean up temporary files using `clean_temp_files.py`
5. Organize the project using `organize_project.py`

## Best Practices

1. Always run scripts with `--dry-run` first to see what changes will be made
2. Keep backups of important files before making major changes
3. Document any changes to the project structure
4. Use the organization scripts regularly to maintain a clean project 