#!/usr/bin/env python3

import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def plot_results(model, temp, iteration_dir=None):
    # Set base directory
    if iteration_dir:
        base_dir = Path(iteration_dir).resolve()
    else:
        base_dir = Path(__file__).resolve().parent.parent
    
    print(f"Plotting results for {model} at {temp}K")
    print(f"Working with iteration directory: {base_dir}")
    
    # Define directories
    analysis_dir = base_dir / "analysis" / model / f"{temp}K"
    plots_dir = analysis_dir / "plots"
    
    # Create plots directory if it doesn't exist
    os.makedirs(plots_dir, exist_ok=True)
    
    # Plot RDF
    try:
        rdf_file = analysis_dir / "rdf_OO.xvg"
        if os.path.exists(rdf_file):
            data = np.loadtxt(rdf_file, comments=["#", "@"])
            r = data[:, 0]  # Distance in nm
            g_r = data[:, 1]  # RDF g(r)
            
            plt.figure(figsize=(10, 6))
            plt.plot(r, g_r, 'b-', linewidth=2)
            plt.xlabel('r (nm)')
            plt.ylabel('g(r)')
            plt.title(f'Radial Distribution Function of {model.upper()} Water at {temp}K')
            plt.grid(True, alpha=0.3)
            plt.savefig(plots_dir / "rdf.png", dpi=300)
            plt.close()
            print(f"Generated RDF plot: {plots_dir}/rdf.png")
    except Exception as e:
        print(f"Error plotting RDF: {e}")
    
    # Plot MSD
    try:
        msd_file = analysis_dir / "msd.xvg"
        if os.path.exists(msd_file):
            data = np.loadtxt(msd_file, comments=["#", "@"])
            time = data[:, 0]  # Time in ps
            msd = data[:, 1]  # MSD in nm²
            
            plt.figure(figsize=(10, 6))
            plt.plot(time, msd, 'r-', linewidth=2)
            plt.xlabel('Time (ps)')
            plt.ylabel('MSD (nm²)')
            plt.title(f'Mean Square Displacement of {model.upper()} Water at {temp}K')
            plt.grid(True, alpha=0.3)
            plt.savefig(plots_dir / "msd.png", dpi=300)
            plt.close()
            print(f"Generated MSD plot: {plots_dir}/msd.png")
            
            # Calculate diffusion coefficient
            # D = slope / (2 * dimensions)
            if len(time) > 100:
                # Use the latter half of the data for linear fit
                half_idx = len(time) // 2
                slope, _ = np.polyfit(time[half_idx:], msd[half_idx:], 1)
                D = slope / 6  # 3D diffusion
                
                with open(plots_dir / "diffusion_coefficient.txt", "w") as f:
                    f.write(f"Diffusion coefficient (D): {D:.6f} cm²/s\n")
                print(f"Calculated diffusion coefficient: {D:.6f} cm²/s")
    except Exception as e:
        print(f"Error plotting MSD: {e}")
    
    # Plot Energy
    try:
        energy_file = analysis_dir / "energy.xvg"
        if os.path.exists(energy_file):
            data = np.loadtxt(energy_file, comments=["#", "@"])
            time = data[:, 0]  # Time in ps
            potential = data[:, 1]  # Potential energy
            kinetic = data[:, 2]  # Kinetic energy
            total = data[:, 3]  # Total energy
            
            plt.figure(figsize=(12, 8))
            plt.plot(time, potential, 'b-', linewidth=1.5, label='Potential')
            plt.plot(time, kinetic, 'r-', linewidth=1.5, label='Kinetic')
            plt.plot(time, total, 'g-', linewidth=2, label='Total')
            plt.xlabel('Time (ps)')
            plt.ylabel('Energy (kJ/mol)')
            plt.title(f'Energy Components of {model.upper()} Water at {temp}K')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(plots_dir / "energy.png", dpi=300)
            plt.close()
            print(f"Generated energy plot: {plots_dir}/energy.png")
    except Exception as e:
        print(f"Error plotting energy: {e}")
    
    # Plot Temperature
    try:
        temp_file = analysis_dir / "temperature.xvg"
        if os.path.exists(temp_file):
            data = np.loadtxt(temp_file, comments=["#", "@"])
            time = data[:, 0]  # Time in ps
            temperature = data[:, 1]  # Temperature in K
            
            plt.figure(figsize=(10, 6))
            plt.plot(time, temperature, 'r-', linewidth=1.5)
            plt.axhline(y=temp, color='k', linestyle='--', alpha=0.7, label=f'Target: {temp}K')
            plt.xlabel('Time (ps)')
            plt.ylabel('Temperature (K)')
            plt.title(f'Temperature of {model.upper()} Water Simulation')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(plots_dir / "temperature.png", dpi=300)
            plt.close()
            print(f"Generated temperature plot: {plots_dir}/temperature.png")
            
            # Calculate average temperature and fluctuation
            avg_temp = np.mean(temperature)
            std_temp = np.std(temperature)
            with open(plots_dir / "temperature_stats.txt", "w") as f:
                f.write(f"Average temperature: {avg_temp:.2f} K\n")
                f.write(f"Standard deviation: {std_temp:.2f} K\n")
                f.write(f"Target temperature: {temp} K\n")
            print(f"Average temperature: {avg_temp:.2f} K ± {std_temp:.2f} K")
    except Exception as e:
        print(f"Error plotting temperature: {e}")
    
    # Plot Pressure
    try:
        pressure_file = analysis_dir / "pressure.xvg"
        if os.path.exists(pressure_file):
            data = np.loadtxt(pressure_file, comments=["#", "@"])
            time = data[:, 0]  # Time in ps
            pressure = data[:, 1]  # Pressure in bar
            
            plt.figure(figsize=(10, 6))
            plt.plot(time, pressure, 'g-', linewidth=1.5)
            plt.axhline(y=1.0, color='k', linestyle='--', alpha=0.7, label='Target: 1 bar')
            plt.xlabel('Time (ps)')
            plt.ylabel('Pressure (bar)')
            plt.title(f'Pressure of {model.upper()} Water at {temp}K')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(plots_dir / "pressure.png", dpi=300)
            plt.close()
            print(f"Generated pressure plot: {plots_dir}/pressure.png")
            
            # Calculate average pressure and fluctuation
            avg_pressure = np.mean(pressure)
            std_pressure = np.std(pressure)
            with open(plots_dir / "pressure_stats.txt", "w") as f:
                f.write(f"Average pressure: {avg_pressure:.2f} bar\n")
                f.write(f"Standard deviation: {std_pressure:.2f} bar\n")
                f.write(f"Target pressure: 1.0 bar\n")
            print(f"Average pressure: {avg_pressure:.2f} bar ± {std_pressure:.2f} bar")
    except Exception as e:
        print(f"Error plotting pressure: {e}")
    
    # Plot density
    try:
        density_file = analysis_dir / "density.xvg"
        if os.path.exists(density_file):
            data = np.loadtxt(density_file, comments=["#", "@"])
            z = data[:, 0]  # z-coordinate
            density = data[:, 1]  # Density
            
            plt.figure(figsize=(10, 6))
            plt.plot(z, density, 'b-', linewidth=1.5)
            plt.xlabel('z (nm)')
            plt.ylabel('Density (kg/m³)')
            plt.title(f'Density Profile of {model.upper()} Water at {temp}K')
            plt.grid(True, alpha=0.3)
            plt.savefig(plots_dir / "density.png", dpi=300)
            plt.close()
            print(f"Generated density plot: {plots_dir}/density.png")
            
            # Calculate average density
            avg_density = np.mean(density)
            with open(plots_dir / "density_stats.txt", "w") as f:
                f.write(f"Average density: {avg_density:.2f} kg/m³\n")
            print(f"Average density: {avg_density:.2f} kg/m³")
    except Exception as e:
        print(f"Error plotting density: {e}")
    
    # Plot Hydrogen Bonds
    try:
        hbond_file = analysis_dir / "hbnum.xvg"
        if os.path.exists(hbond_file):
            data = np.loadtxt(hbond_file, comments=["#", "@"])
            time = data[:, 0]  # Time in ps
            hbonds = data[:, 1]  # Number of hydrogen bonds
            
            plt.figure(figsize=(10, 6))
            plt.plot(time, hbonds, 'b-', linewidth=1.5)
            plt.xlabel('Time (ps)')
            plt.ylabel('Number of Hydrogen Bonds')
            plt.title(f'Hydrogen Bonds in {model.upper()} Water at {temp}K')
            plt.grid(True, alpha=0.3)
            plt.savefig(plots_dir / "hydrogen_bonds.png", dpi=300)
            plt.close()
            print(f"Generated hydrogen bonds plot: {plots_dir}/hydrogen_bonds.png")
            
            # Calculate average number of hydrogen bonds
            avg_hbonds = np.mean(hbonds)
            std_hbonds = np.std(hbonds)
            with open(plots_dir / "hbond_stats.txt", "w") as f:
                f.write(f"Average number of hydrogen bonds: {avg_hbonds:.2f}\n")
                f.write(f"Standard deviation: {std_hbonds:.2f}\n")
            print(f"Average number of hydrogen bonds: {avg_hbonds:.2f} ± {std_hbonds:.2f}")
            
            # Plot hydrogen bond distribution if available
            hbdist_file = analysis_dir / "hbdist.xvg"
            if os.path.exists(hbdist_file):
                data = np.loadtxt(hbdist_file, comments=["#", "@"])
                r = data[:, 0]  # Distance in nm
                dist = data[:, 1]  # Distribution
                
                plt.figure(figsize=(10, 6))
                plt.plot(r, dist, 'r-', linewidth=1.5)
                plt.xlabel('Distance (nm)')
                plt.ylabel('Count')
                plt.title(f'Hydrogen Bond Distance Distribution in {model.upper()} Water at {temp}K')
                plt.grid(True, alpha=0.3)
                plt.savefig(plots_dir / "hbond_distance.png", dpi=300)
                plt.close()
                print(f"Generated hydrogen bond distance plot: {plots_dir}/hbond_distance.png")
            
            # Plot hydrogen bond angle distribution if available
            hbang_file = analysis_dir / "hbang.xvg"
            if os.path.exists(hbang_file):
                data = np.loadtxt(hbang_file, comments=["#", "@"])
                angle = data[:, 0]  # Angle in degrees
                dist = data[:, 1]  # Distribution
                
                plt.figure(figsize=(10, 6))
                plt.plot(angle, dist, 'g-', linewidth=1.5)
                plt.xlabel('Angle (degrees)')
                plt.ylabel('Count')
                plt.title(f'Hydrogen Bond Angle Distribution in {model.upper()} Water at {temp}K')
                plt.grid(True, alpha=0.3)
                plt.savefig(plots_dir / "hbond_angle.png", dpi=300)
                plt.close()
                print(f"Generated hydrogen bond angle plot: {plots_dir}/hbond_angle.png")
    except Exception as e:
        print(f"Error plotting hydrogen bonds: {e}")
    
    # Plot Velocity Autocorrelation Function (VACF)
    try:
        vacf_file = analysis_dir / "vacf.xvg"
        if os.path.exists(vacf_file):
            data = np.loadtxt(vacf_file, comments=["#", "@"])
            time = data[:, 0]  # Time in ps
            vacf = data[:, 1]  # VACF
            
            plt.figure(figsize=(10, 6))
            plt.plot(time, vacf, 'b-', linewidth=1.5)
            plt.xlabel('Time (ps)')
            plt.ylabel('Velocity Autocorrelation')
            plt.title(f'Velocity Autocorrelation Function of {model.upper()} Water at {temp}K')
            plt.grid(True, alpha=0.3)
            plt.savefig(plots_dir / "vacf.png", dpi=300)
            plt.close()
            print(f"Generated VACF plot: {plots_dir}/vacf.png")
            
            # Plot VACF spectrum if available
            spectrum_file = analysis_dir / "vacf_spectrum.xvg"
            if os.path.exists(spectrum_file):
                data = np.loadtxt(spectrum_file, comments=["#", "@"])
                freq = data[:, 0]  # Frequency in 1/ps or cm^-1
                intensity = data[:, 1]  # Spectral intensity
                
                plt.figure(figsize=(10, 6))
                plt.plot(freq, intensity, 'r-', linewidth=1.5)
                plt.xlabel('Frequency (cm$^{-1}$)')
                plt.ylabel('Intensity')
                plt.title(f'Vibrational Spectrum of {model.upper()} Water at {temp}K')
                plt.grid(True, alpha=0.3)
                plt.savefig(plots_dir / "vacf_spectrum.png", dpi=300)
                plt.close()
                print(f"Generated VACF spectrum plot: {plots_dir}/vacf_spectrum.png")
    except Exception as e:
        print(f"Error plotting VACF: {e}")
    
    print(f"All plots have been saved to {plots_dir}")
    return True

def main():
    parser = argparse.ArgumentParser(description="Plot water simulation results")
    parser.add_argument("--model", type=str, default="tip4p", help="Water model (default: tip4p)")
    parser.add_argument("--temp", type=int, default=273, help="Temperature in K (default: 273)")
    parser.add_argument("--iteration-dir", type=str, help="Path to the iteration directory (default: parent directory of this script)")
    args = parser.parse_args()
    
    plot_results(args.model, args.temp, args.iteration_dir)

if __name__ == "__main__":
    main() 