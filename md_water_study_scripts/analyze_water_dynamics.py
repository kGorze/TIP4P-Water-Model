#!/usr/bin/env python3

import os
import subprocess
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib as mpl
import pandas as pd

# Set publication-quality plot style
sns.set_theme(style="whitegrid")
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['font.size'] = 10
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['axes.grid'] = True
mpl.rcParams['grid.alpha'] = 0.3
mpl.rcParams['figure.figsize'] = [8, 6]

def run_command(cmd, cwd=None, input_text=None):
    """Run a shell command and handle errors"""
    try:
        if input_text:
            result = subprocess.run(cmd, check=True, shell=True, cwd=cwd, 
                                   input=input_text.encode(), 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
        else:
            result = subprocess.run(cmd, check=True, shell=True, cwd=cwd,
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
        return True, result.stdout.decode()
    except subprocess.CalledProcessError as e:
        return False, f"Error running command: {cmd}\nError: {e}\nOutput: {e.stdout.decode()}\nError: {e.stderr.decode()}"

def read_xvg(file_path):
    """Read XVG file and return time and data arrays with metadata"""
    time = []
    data = []
    title = ""
    xlabel = "Time (ps)"
    ylabel = ""
    legends = []
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                if "title" in line:
                    title = line.split('"')[1]
                elif "xaxis label" in line:
                    xlabel = line.split('"')[1]
                elif "yaxis label" in line:
                    ylabel = line.split('"')[1]
                elif "legend" in line and '"' in line:
                    legends.append(line.split('"')[1])
                continue
            
            values = [float(val) for val in line.strip().split()]
            if len(values) > 0:
                if len(time) == 0 or values[0] > time[-1]:  # Ensure time is increasing
                    time.append(values[0])
                    data.append(values[1:])
    
    return np.array(time), np.array(data), title, xlabel, ylabel, legends

def analyze_hbonds(output_dir, analysis_dir):
    """Analyze hydrogen bonds in water"""
    print("Analyzing hydrogen bonds...")
    
    # Check if md.xtc exists, if not use the path from output_dir
    if not (output_dir / "md.xtc").exists():
        if (output_dir / "md.trr").exists():
            print("Found md.trr, converting to md.xtc...")
            run_command(
                "echo '0' | gmx trjconv -f md.trr -s md.tpr -o md.xtc",
                cwd=output_dir
            )
        else:
            return False, "Error: Neither md.xtc nor md.trr found. Run production_md.py first."
    
    # Run h-bond analysis (updated for GROMACS 2024.3)
    success, output = run_command(
        "echo '1 1' | gmx hbond -f md.xtc -s md.tpr -num hbnum.xvg -dist hbdist.xvg -ang hbang.xvg",
        cwd=output_dir
    )
    
    if not success:
        print(output)
        return False, "Failed to analyze hydrogen bonds"
    
    # Copy results to analysis directory
    for file in ["hbnum.xvg", "hbdist.xvg", "hbang.xvg"]:
        if (output_dir / file).exists():
            import shutil
            shutil.copy2(output_dir / file, analysis_dir / file)
    
    return True, "Hydrogen bond analysis completed successfully"

def analyze_diffusion(output_dir, analysis_dir):
    """Analyze diffusion coefficient using MSD data"""
    print("Analyzing diffusion coefficient...")
    
    # Check if md.xtc exists
    if not (output_dir / "md.xtc").exists():
        if (output_dir / "md.trr").exists():
            print("Found md.trr, converting to md.xtc...")
            success, output = run_command(
                "echo '0' | gmx trjconv -f md.trr -s md.tpr -o md.xtc",
                cwd=output_dir
            )
            if not success:
                print(output)
                return False, "Failed to convert md.trr to md.xtc"
        else:
            return False, "Error: Neither md.xtc nor md.trr found. Run production_md.py first."
    
    # Run MSD analysis (updated for GROMACS 2024.3)
    success, output = run_command(
        "echo '2' | gmx msd -f md.xtc -s md.tpr -o msd.xvg -mol",
        cwd=output_dir
    )
    
    if not success:
        print(output)
        return False, "Failed to analyze MSD"
    
    # Copy the MSD file to analysis directory
    if (output_dir / "msd.xvg").exists():
        import shutil
        shutil.copy2(output_dir / "msd.xvg", analysis_dir / "msd.xvg")
    
    return True, "Diffusion coefficient analysis completed successfully"

def linear_func(x, a, b):
    """Linear function for MSD fitting: MSD = 6*D*t + b"""
    return a * x + b

def calculate_diffusion_coefficient(msd_file):
    """Calculate diffusion coefficient from MSD data"""
    time, data, _, _, _, _ = read_xvg(msd_file)
    
    # For water, use the first column of data (usually total MSD)
    msd_data = data[:, 0]
    
    # Find the linear region (typically from 10% to 90% of trajectory)
    start_idx = int(len(time) * 0.1)
    end_idx = int(len(time) * 0.9)
    
    # Fit a linear function to the MSD data in the linear region
    # MSD = 6*D*t + b, where D is the diffusion coefficient
    popt, pcov = curve_fit(linear_func, time[start_idx:end_idx], msd_data[start_idx:end_idx])
    
    # Calculate diffusion coefficient (in nm²/ps)
    slope = popt[0]
    diffusion_coeff_nm2_ps = slope / 6.0
    
    # Convert to standard units (cm²/s)
    # 1 nm²/ps = 10⁻⁷ cm²/s
    diffusion_coeff_cm2_s = diffusion_coeff_nm2_ps * 1e-7
    
    # Calculate error from covariance matrix
    perr = np.sqrt(np.diag(pcov))
    error_nm2_ps = perr[0] / 6.0
    error_cm2_s = error_nm2_ps * 1e-7
    
    return diffusion_coeff_nm2_ps, diffusion_coeff_cm2_s, error_nm2_ps, error_cm2_s, time, msd_data, popt

def plot_msd(time, msd_data, popt, diffusion_coeff_cm2_s, error_cm2_s, output_dir):
    """Create detailed MSD plot with diffusion coefficient"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot MSD data
    ax.plot(time, msd_data, 'o', markersize=2, alpha=0.5, label='MSD Data')
    
    # Plot fitted line
    fit_x = np.linspace(time[0], time[-1], 100)
    ax.plot(fit_x, linear_func(fit_x, *popt), 'r-', linewidth=2, 
            label=f'Linear Fit (Slope = {popt[0]:.4g})')
    
    # Add diffusion coefficient annotation
    ax.text(0.05, 0.95, 
            f'Diffusion Coefficient:\n{diffusion_coeff_cm2_s:.4g} ± {error_cm2_s:.2g} cm²/s\n'
            f'({diffusion_coeff_cm2_s*1e5:.4g} ± {error_cm2_s*1e5:.2g} ×10⁻⁵ cm²/s)', 
            transform=ax.transAxes, fontsize=12, va='top', 
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Add reference value for experimental water diffusion at this temperature
    # For water at 273K (0°C), D ≈ 1.1 × 10⁻⁵ cm²/s
    exp_diff = 1.1e-5  # cm²/s
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.text(0.65, 0.15, f'Experimental value at 273K:\n{exp_diff:.2g} cm²/s', 
            transform=ax.transAxes, fontsize=10, 
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('MSD (nm²)')
    ax.set_title('Mean Square Displacement and Diffusion Coefficient')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / "diffusion_coefficient.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_hydrogen_bonds(num_file, dist_file, ang_file, output_dir):
    """Create comprehensive hydrogen bond analysis plots"""
    # Create a 2x2 figure for different h-bond properties
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot number of hydrogen bonds
    if num_file.exists():
        time, data, _, _, _, _ = read_xvg(num_file)
        axes[0, 0].plot(time, data[:, 0], 'b-', alpha=0.7)
        
        # Calculate average number of h-bonds
        avg_hbonds = np.mean(data[:, 0])
        std_hbonds = np.std(data[:, 0])
        
        axes[0, 0].axhline(y=avg_hbonds, color='r', linestyle='--')
        axes[0, 0].text(0.05, 0.95, f'Average: {avg_hbonds:.2f} ± {std_hbonds:.2f} H-bonds', 
                      transform=axes[0, 0].transAxes, fontsize=10, va='top', 
                      bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        axes[0, 0].set_xlabel('Time (ps)')
        axes[0, 0].set_ylabel('Number of H-bonds')
        axes[0, 0].set_title('Hydrogen Bond Count')
    else:
        axes[0, 0].text(0.5, 0.5, "Hydrogen bond number data not found", 
                      ha='center', va='center', transform=axes[0, 0].transAxes)
    
    # Plot hydrogen bond distance distribution
    if dist_file.exists():
        dist, data, _, _, _, _ = read_xvg(dist_file)
        axes[0, 1].plot(dist, data[:, 0], 'g-', linewidth=2)
        
        # Mark the peak (most common h-bond distance)
        peak_idx = np.argmax(data[:, 0])
        peak_dist = dist[peak_idx]
        
        axes[0, 1].axvline(x=peak_dist, color='r', linestyle='--')
        axes[0, 1].text(0.05, 0.95, f'Most common distance: {peak_dist:.3f} nm', 
                      transform=axes[0, 1].transAxes, fontsize=10, va='top', 
                      bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        axes[0, 1].set_xlabel('Distance (nm)')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('H-bond Distance Distribution')
    else:
        axes[0, 1].text(0.5, 0.5, "Hydrogen bond distance data not found", 
                      ha='center', va='center', transform=axes[0, 1].transAxes)
    
    # Plot hydrogen bond angle distribution
    if ang_file.exists():
        angle, data, _, _, _, _ = read_xvg(ang_file)
        axes[1, 0].plot(angle, data[:, 0], 'purple', linewidth=2)
        
        # Mark the peak (most common h-bond angle)
        peak_idx = np.argmax(data[:, 0])
        peak_angle = angle[peak_idx]
        
        axes[1, 0].axvline(x=peak_angle, color='r', linestyle='--')
        axes[1, 0].text(0.05, 0.95, f'Most common angle: {peak_angle:.1f}°', 
                      transform=axes[1, 0].transAxes, fontsize=10, va='top', 
                      bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        axes[1, 0].set_xlabel('Angle (degrees)')
        axes[1, 0].set_ylabel('Frequency')
        axes[1, 0].set_title('H-bond Angle Distribution')
    else:
        axes[1, 0].text(0.5, 0.5, "Hydrogen bond angle data not found", 
                      ha='center', va='center', transform=axes[1, 0].transAxes)
    
    # We removed the lifetime plot since that option was not available
    axes[1, 1].text(0.5, 0.5, "Additional hydrogen bond data:\n\nH-bond cutoffs:\nDistance < 0.35 nm\nAngle < 30°", 
                  ha='center', va='center', transform=axes[1, 1].transAxes,
                  fontsize=12, bbox=dict(boxstyle='round', facecolor='whitesmoke', alpha=0.8))
    axes[1, 1].set_title('H-bond Parameters')
    axes[1, 1].axis('off')
    
    plt.tight_layout()
    plt.savefig(output_dir / "hydrogen_bonds.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_water_dynamics_report(diffusion_coeff, hbond_data, output_dir):
    """Create a comprehensive report on water dynamics"""
    # Extract diffusion coefficient data
    D_nm2_ps, D_cm2_s, err_nm2_ps, err_cm2_s = diffusion_coeff[:4]
    
    # Extract hydrogen bond data if available
    avg_hbonds = hbond_data.get('avg_hbonds', 'N/A')
    std_hbonds = hbond_data.get('std_hbonds', 'N/A') 
    peak_dist = hbond_data.get('peak_dist', 'N/A')
    peak_angle = hbond_data.get('peak_angle', 'N/A')
    
    # Reference experimental values
    exp_diff_273K = 1.1e-5  # cm²/s at 273K
    exp_hbonds = "3.5-4.0"  # per water molecule
    exp_lifetime = "1-2"    # ps
    
    # Calculate percentage difference from experimental value
    if D_cm2_s > 0:
        diff_percent = (D_cm2_s - exp_diff_273K) / exp_diff_273K * 100
    else:
        diff_percent = "N/A"
    
    # Create report text
    report_text = f"""# Water Dynamics Analysis Report

## Diffusion Coefficient Results

- **Diffusion coefficient (simulation)**: {D_cm2_s:.6g} ± {err_cm2_s:.2g} cm²/s
- **Diffusion coefficient (×10⁻⁵)**: {D_cm2_s*1e5:.4g} ± {err_cm2_s*1e5:.2g} ×10⁻⁵ cm²/s
- **Experimental value at 273K**: {exp_diff_273K:.2g} cm²/s (1.1 ×10⁻⁵ cm²/s)
- **Difference from experiment**: {diff_percent}%

## Hydrogen Bond Analysis

- **Average number of H-bonds per molecule**: {avg_hbonds}
- **Standard deviation**: {std_hbonds}
- **Experimental reference value**: {exp_hbonds} H-bonds per molecule
- **Most common H-bond distance**: {peak_dist} nm
- **Most common H-bond angle**: {peak_angle}°
- **Typical H-bond lifetime**: {exp_lifetime} ps (experimental reference)

## Water Structure and Dynamics Interpretation

The diffusion coefficient measures how quickly water molecules move through the system. A value close to the 
experimental reference (1.1×10⁻⁵ cm²/s at 273K) indicates the simulation correctly captures water dynamics.

Hydrogen bonding statistics reveal the structural network of water. With approximately 3.5-4 hydrogen 
bonds per molecule in liquid water, the hydrogen bond structure is fundamental to water's unique properties.

The hydrogen bond lifetime is typically around 1-2 ps, highlighting the dynamic nature of the hydrogen 
bond network in liquid water, with bonds constantly breaking and reforming.

## Conclusions

This analysis provides insight into how well the simulation reproduces the correct behavior of water
at the simulated temperature. Significant deviations in diffusion coefficient may indicate issues with:
- Simulation parameters (timestep, cutoffs)
- Water model parameterization
- Equilibration issues
- Temperature or pressure control problems

Similarly, hydrogen bond statistics that differ substantially from expected values can indicate structural
problems with the water model or simulation conditions.
"""
    
    # Write report to file
    with open(output_dir / "water_dynamics_report.md", "w") as f:
        f.write(report_text)
    
    print(f"Created water dynamics report: {output_dir / 'water_dynamics_report.md'}")
    return True

def main():
    parser = argparse.ArgumentParser(description="Analyze water dynamics")
    parser.add_argument("--model", type=str, default="tip4p", help="Water model (tip4p)")
    parser.add_argument("--temp", type=int, default=273, help="Temperature in K")
    args = parser.parse_args()
    
    base_dir = Path(__file__).resolve().parent.parent
    data_dir = base_dir / "data" / args.model / f"{args.temp}K" / "outputs"
    analysis_dir = base_dir / "analysis" / args.model / f"{args.temp}K"
    plots_dir = analysis_dir / "plots"
    
    if not data_dir.exists():
        print(f"Data directory not found: {data_dir}")
        return 1
    
    # Create analysis and plots directories
    analysis_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Analyze hydrogen bonds
    success, message = analyze_hbonds(data_dir, analysis_dir)
    if success:
        print(message)
    else:
        print(f"Warning: {message}")
    
    # Analyze diffusion
    success, message = analyze_diffusion(data_dir, analysis_dir)
    if success:
        print(message)
    else:
        print(f"Warning: {message}")
    
    # Calculate diffusion coefficient from MSD data
    if (analysis_dir / "msd.xvg").exists():
        diffusion_results = calculate_diffusion_coefficient(analysis_dir / "msd.xvg")
        D_nm2_ps, D_cm2_s, err_nm2_ps, err_cm2_s, time, msd_data, popt = diffusion_results
        
        print(f"Diffusion coefficient: {D_cm2_s:.6g} ± {err_cm2_s:.2g} cm²/s")
        print(f"Diffusion coefficient: {D_cm2_s*1e5:.4g} ± {err_cm2_s*1e5:.2g} ×10⁻⁵ cm²/s")
        
        # Plot MSD and diffusion coefficient
        plot_msd(time, msd_data, popt, D_cm2_s, err_cm2_s, plots_dir)
        print("Generated diffusion coefficient plot")
    else:
        print("MSD data not found, skipping diffusion coefficient calculation")
        diffusion_results = (0, 0, 0, 0, None, None, None)
    
    # Create hydrogen bond plots
    hbond_data = {}
    if any([(analysis_dir / f).exists() for f in ["hbnum.xvg", "hbdist.xvg", "hbang.xvg"]]):
        # Extract hydrogen bond statistics for the report
        if (analysis_dir / "hbnum.xvg").exists():
            time, data, _, _, _, _ = read_xvg(analysis_dir / "hbnum.xvg")
            hbond_data['avg_hbonds'] = np.mean(data[:, 0])
            hbond_data['std_hbonds'] = np.std(data[:, 0])
        
        if (analysis_dir / "hbdist.xvg").exists():
            dist, data, _, _, _, _ = read_xvg(analysis_dir / "hbdist.xvg")
            peak_idx = np.argmax(data[:, 0])
            hbond_data['peak_dist'] = dist[peak_idx]
        
        if (analysis_dir / "hbang.xvg").exists():
            angle, data, _, _, _, _ = read_xvg(analysis_dir / "hbang.xvg")
            peak_idx = np.argmax(data[:, 0])
            hbond_data['peak_angle'] = angle[peak_idx]
        
        plot_hydrogen_bonds(
            analysis_dir / "hbnum.xvg", 
            analysis_dir / "hbdist.xvg", 
            analysis_dir / "hbang.xvg", 
            plots_dir
        )
        print("Generated hydrogen bond plots")
    else:
        print("Hydrogen bond data not found, skipping plots")
    
    # Create comprehensive report
    create_water_dynamics_report(diffusion_results, hbond_data, plots_dir)
    
    print("Water dynamics analysis completed!")
    return 0

if __name__ == "__main__":
    exit(main()) 