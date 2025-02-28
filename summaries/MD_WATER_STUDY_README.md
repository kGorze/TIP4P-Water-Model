# MD Water Study Centralized Scripts System

This project organizes multiple iterations of molecular dynamics water simulations, providing a centralized scripts directory to avoid code duplication and ensure consistent behavior across all iterations.

## Project Structure

The project consists of:

1. `md_water_study_scripts/` - The central directory containing all scripts
2. `md_water_study_iteration_X/` - Individual iteration directories containing data, configs, and analysis results
3. `md_water_study_iteration_result/` - Directory for combined results from all iterations

Each iteration directory follows a consistent structure based on the model in iteration 3:

```
md_water_study_iteration_X/
├── analysis/                # Analysis results
├── configs/                 # Template configuration files
├── data/                    # Simulation data
├── docs/                    # Documentation
├── unused_analysis/         # Deprecated analysis
├── unused_scripts/          # Deprecated scripts
└── README.md                # Iteration-specific documentation
```

## Using the Centralized Scripts

### Running Scripts

The central scripts directory provides a convenient utility to run any script on any iteration:

```bash
cd /home/konrad_guest/Documents/research/cursor/md_water_study_scripts
./run_on_iteration.py script_name.py iteration_name [additional arguments]
```

For example:
```bash
# Run the workflow on iteration 4
./run_on_iteration.py run_workflow.py md_water_study_iteration_4 --model tip4p --temp 273

# Run analysis on iteration 3
./run_on_iteration.py analyze_results.py md_water_study_iteration_3 --model tip4p --temp 273
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
./run_on_iteration.py script_name.py all [additional arguments]
```

For example:
```bash
# Run analysis on all iterations
./run_on_iteration.py analyze_results.py all --model tip4p --temp 273
```

## Direct Script Usage

You can also use the scripts directly by specifying the iteration directory:

```bash
cd /home/konrad_guest/Documents/research/cursor/md_water_study_scripts
python run_workflow.py --model tip4p --temp 273 --iteration-dir /home/konrad_guest/Documents/research/cursor/md_water_study_iteration_3
```

## Workflow Steps

The typical workflow includes:

1. **Setup and Simulation**: Run simulations on water molecules
2. **Analysis**: Calculate RDF, density, MSD, energy components, etc.
3. **Visualization**: Generate plots of simulation results

## Benefits of Centralized Scripts

1. **Consistency**: All iterations use the same script versions
2. **Maintainability**: Bug fixes only need to be made in one place
3. **Efficiency**: No need to recreate scripts for each new iteration
4. **Flexibility**: Can specify which iteration to work with for any script 