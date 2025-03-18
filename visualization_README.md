# Molecular Dynamics Visualization Scripts

This directory contains scripts to visualize key concepts in molecular dynamics simulations, specifically focusing on water models and simulation parameters.

## Overview

These visualization scripts illustrate four important concepts in molecular dynamics simulations:

1. **Water Molecule Constraints** - Visualizes the fixed bond lengths and angles in water models
2. **Non-Bonded Interactions Cutoffs** - Illustrates the cutoff radius for non-bonded interactions and PME
3. **Integration Time Step** - Demonstrates how constraints enable larger time steps
4. **Sampling Methods for Water Models** - Shows different approaches to explore the free energy landscape

## Scripts

### 1. Water Molecule Constraint (`water_molecule_constraint.py`)

This script creates 2D and 3D visualizations of a water molecule with constrained bond lengths and angles. It shows:

- Fixed OH bond lengths (0.9572 Å)
- Fixed HOH angle (104.52°)
- Both 2D schematic and 3D model representations
- Explanation of how constraints benefit MD simulations

### 2. Non-Bonded Interactions Cutoffs (`minimal_nonbonded_cutoffs.py`)

This script visualizes the concept of interaction cutoffs in MD simulations:

- Lennard-Jones and Coulomb potential energy curves with cutoff
- 2D representation of the cutoff sphere (typically 1.0-1.2 nm)
- PME grid representation for long-range electrostatics
- Explanation of how cutoffs and PME work together

The script creates three separate visualizations:
- `nonbonded_potentials.png` - Shows the potential energy curves
- `cutoff_2d.png` - Shows the 2D representation of the cutoff sphere
- `pme_grid.png` - Shows the PME grid for long-range electrostatics

### 3. Integration Time Step Illustration (`timestep_visualization.py`)

This script demonstrates how constraints enable larger time steps:

- Comparison of different time steps (0.5 fs, 1.0 fs, 2.0 fs)
- Timeline visualization showing computational efficiency gains
- Animations of constrained vs. unconstrained water molecules
- Explanation of the Nyquist limit and high-frequency bond vibrations

### 4. Sampling Methods for Water Models (`sampling_methods_visualization.py`)

This script visualizes different sampling approaches for water models in MD simulations:

- 2D free energy landscape representation of water configurations
- Comparison between standard MD sampling and enhanced sampling techniques
- Animation showing how different methods explore the energy landscape over time
- Explanation of various enhanced sampling methods (Replica Exchange, Umbrella Sampling, etc.)
- Demonstration of how multiple short simulations can provide better sampling than a single long one

The script creates:
- `sampling_methods.png` - Static visualization of different sampling approaches
- `sampling_comparison.gif` - Animation comparing standard vs. enhanced sampling over time

## Running the Scripts

You can run all visualizations at once using the provided runner script:

```bash
unset PYTHONHOME && unset PYTHONPATH && python3 run_visualizations.py
```

Or run individual scripts:

```bash
unset PYTHONHOME && unset PYTHONPATH && python3 water_molecule_constraint.py
unset PYTHONHOME && unset PYTHONPATH && python3 minimal_nonbonded_cutoffs.py
unset PYTHONHOME && unset PYTHONPATH && python3 timestep_visualization.py
unset PYTHONHOME && unset PYTHONPATH && python3 sampling_methods_visualization.py
```

## Output

All visualizations are saved to the `visualization_outputs` directory in both PNG and PDF formats. The time step and sampling methods visualizations also create GIF animations.

## Requirements

These scripts require the following Python packages:
- numpy
- matplotlib
- pillow (for GIF animations)

## Notes

- The animations may take a moment to generate, especially the time step and sampling visualizations
- For best results, view the PDF outputs for high-quality vector graphics
- The GIF animations illustrate dynamic concepts that are difficult to show in static images 