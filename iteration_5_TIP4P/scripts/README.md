# TIP4P Water Analysis Scripts

This directory contains Python scripts for analyzing GROMACS molecular dynamics simulations of TIP4P water. The scripts extract structural, dynamic, and thermodynamic properties from GROMACS output files.

## Directory Structure

The expected directory structure is:
```
manual_tip4p/
├── md.tpr, md.xtc, md.edr, etc. (GROMACS files in parent directory)
├── analysis/  (created automatically - output directory)
└── scripts/   (this directory - contains analysis scripts)
```

The scripts will automatically look for GROMACS files in the parent directory if they are not found in the current directory.

## Requirements

- Python 3.6 or higher
- Required Python packages:
  - MDAnalysis (recommended version: 2.0.0)
  - NumPy (recommended version: 1.24.3)
  - Matplotlib (recommended version: 3.7.2)
  - Pandas (recommended version: 2.0.3)
  - SciPy (recommended version: 1.10.1)
  - Seaborn (recommended version: 0.12.2)

You can install the required packages using pip:

```bash
pip install mdanalysis==2.0.0 numpy==1.24.3 matplotlib==3.7.2 pandas==2.0.3 scipy==1.10.1 seaborn==0.12.2
```

## GROMACS Files

The scripts analyze the following GROMACS output files:
- `.tpr` files (binary run input)
- `.xtc` or `.trr` files (trajectory)
- `.edr` files (energy data)
- `.gro` files (structure)

## Compatibility and Troubleshooting

### Package Version Compatibility

These scripts are sensitive to package versions. Common issues include:

1. **NumPy-SciPy compatibility**: Use NumPy < 1.25.0 for best compatibility with SciPy
2. **MDAnalysis module changes**: The `diffusion` module was replaced by `msd` in newer versions
3. **EinsteinMSD API changes**: Newer MDAnalysis versions (>=2.0.0) store data in `results.timeseries` instead of directly in `timeseries`
4. **Matplotlib Axes3D warnings**: These are harmless warnings related to 3D plotting

### Troubleshooting

If you encounter errors:

1. **ImportError for diffusion module**:
   - Solution: Install MDAnalysis 2.0.0 which uses the msd module instead
   - Alternative: Install an older version like MDAnalysis 0.20.0

2. **AttributeError with EinsteinMSD**:
   - Error: `'EinsteinMSD' object has no attribute 'timeseries'`
   - Solution: This script has been updated to handle both older and newer MDAnalysis APIs
   - If errors persist, ensure you're using MDAnalysis 2.0.0 or newer:
   ```bash
   pip install MDAnalysis==2.0.0
   ```

3. **IndexError or TypeError with EinsteinMSD results**:
   - Errors include:
     - `too many indices for array: array is 1-dimensional, but 2 were indexed`
     - `unhashable type: 'slice'`
   - These occur due to differences in how MSD data is structured across MDAnalysis versions
   - The script now has extensive debug output and automatic handling for different data structures
   - If you're having issues, run with `./diffusion_analysis.py` to see detailed debug information
   - In modern MDAnalysis versions (>2.0.0), the structure is often:
     ```
     msd.results.timeseries = [time_array, msd_values_array]
     ```
   - The script now automatically detects and handles this structure

4. **TypeError with HydrogenBondAnalysis**:
   - Error: `AnalysisBase.run() got an unexpected keyword argument 'n_jobs'`
   - This occurs because the `n_jobs` parameter was added in newer versions of MDAnalysis
   - The script now automatically detects if your version supports this parameter
   - If not, it falls back to the version without parallel processing

5. **DeprecationWarning and array size mismatch with HydrogenBondAnalysis**:
   - Warning: `The 'hbonds' attribute was deprecated in MDAnalysis 2.0.0 and will be removed in MDAnalysis 3.0.0`
   - Error: `ValueError: all the input array dimensions except for the concatenation axis must match exactly`
   - This occurs because the API for accessing hydrogen bond data changed in MDAnalysis 2.0.0
   - The script now handles both older and newer APIs and ensures array sizes match

6. **SciPy compatibility errors with NumPy**:
   - Solution: Downgrade NumPy to 1.24.3 or earlier
   ```bash
   pip install numpy==1.24.3
   ```

7. **Matplotlib warnings**:
   - These warnings are harmless and can be ignored
   - If they persist, try reinstalling matplotlib:
   ```bash
   pip install matplotlib==3.7.2
   ```

8. **Water selection errors**:
   - The scripts now try multiple common water residue names (SOL, TIP4, WAT)
   - If selection fails, it falls back to selecting atoms by name (OW, O)

The `run_all_analysis.py` script includes automatic version checking and will warn you about potential compatibility issues.

## Available Analysis Scripts

1. `rdf_analysis.py` - Calculates radial distribution functions (RDF) for O-O, O-H, and H-H pairs
2. `density_analysis.py` - Analyzes density distributions, profiles, and maps
3. `diffusion_analysis.py` - Calculates mean squared displacement (MSD) and diffusion coefficients
4. `hbond_analysis.py` - Analyzes hydrogen bonds, including counts, distributions, and lifetimes
5. `rmsd_analysis.py` - Calculates root mean square deviation (RMSD) and identifies equilibration
6. `thermodynamics_analysis.py` - Analyzes thermodynamic properties (energy, temperature, pressure)

## Utility Scripts

- `utils.py` - Common utility functions used by other scripts
- `run_all_analysis.py` - Runs all analysis scripts and creates a summary

## Usage

### Running Individual Scripts

Each script can be run independently. For example:

```bash
python rdf_analysis.py md.tpr md.xtc
python thermodynamics_analysis.py md.edr
```

The default filenames are `md.tpr`, `md.xtc`, and `md.edr`, so you can run the scripts without arguments if your files have these names.

#### Debugging Diffusion Analysis

The diffusion_analysis.py script has been updated to better handle different MDAnalysis API versions. If you encounter issues:

1. Make the script executable:
   ```bash
   chmod +x diffusion_analysis.py
   ```

2. Run it directly for more detailed output:
   ```bash
   ./diffusion_analysis.py
   ```

The script will print debugging information about the available attributes and data structures to help diagnose any issues.

### Running All Analyses

To run all analysis scripts in sequence and generate a summary:

```bash
python run_all_analysis.py [tpr_file] [trajectory_file] [edr_file]
```

Example:
```bash
python run_all_analysis.py md.tpr md.xtc md.edr
```

## Output

Analysis results are saved to the `../analysis/` directory as:
- Data files (.dat)
- Plot images (.png)
- A summary PDF and image

## Examples of Generated Plots

- RDF plots: `rdf_OO_plot.png`, `rdf_OH_plot.png`, `rdf_HH_plot.png`, `combined_rdf_plot.png`
- Density plots: `density_histogram.png`, `density_profile_plot.png`, `density_map_*.png`
- Diffusion plots: `msd_plot.png`, `diffusion_analysis_plot.png`, `vacf_plot.png`
- H-bond plots: `hbond_count_plot.png`, `hblife_plot.png`, `combined_hbond_plot.png`
- RMSD plots: `rmsd_plot.png`, `rmsd_equilibrated_plot.png`, `rmsf_plot.png`
- Thermodynamic plots: `temperature_plot.png`, `pressure_plot.png`, `energy_plot.png`, etc.

## Summary Report

After running all analyses, a summary report is generated:
- `tip4p_water_analysis_summary.pdf` - Multi-page PDF with key plots
- `tip4p_water_analysis_summary.png` - Single image with key plots

## Customization

You can modify parameters in each script to customize the analysis:
- Change cutoff distances for RDF and hydrogen bond analysis
- Adjust bin sizes for histograms
- Select different atom selections
- Modify plot appearance 