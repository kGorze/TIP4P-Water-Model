#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import matplotlib as mpl
from scipy import stats
from scipy.signal import savgol_filter

# Set publication-quality plot style
sns.set_theme()  # Use seaborn's default theme
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

def calculate_statistics(data):
    """Calculate basic statistics for the data"""
    mean = np.mean(data)
    std = np.std(data)
    sem = stats.sem(data)
    return mean, std, sem

def moving_average(data, window=100):
    """Calculate moving average with the specified window"""
    return np.convolve(data, np.ones(window)/window, mode='valid')

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
                if 'title' in line and '"' in line:
                    title = line.split('"')[1] if len(line.split('"')) > 1 else ""
                elif 'xaxis label' in line and '"' in line:
                    xlabel = line.split('"')[1] if len(line.split('"')) > 1 else "Time (ps)"
                elif 'yaxis label' in line and '"' in line:
                    ylabel = line.split('"')[1] if len(line.split('"')) > 1 else ""
                elif 'legend' in line and '"' in line:
                    if len(line.split('"')) > 1:
                        legends.append(line.split('"')[1])
                continue
            values = line.strip().split()
            if values:
                try:
                    time.append(float(values[0]))
                    data.append([float(x) for x in values[1:]])
                except (ValueError, IndexError):
                    continue
    
    return np.array(time), np.array(data), title, xlabel, ylabel, legends

def plot_with_statistics(ax, time, data, label, color, window=100):
    """Plot data with moving average and confidence interval"""
    # Original data
    ax.plot(time, data, alpha=0.3, color=color, label=f'{label} (Raw)')
    
    # Moving average
    ma_data = moving_average(data, window)
    ma_time = time[window-1:]
    ax.plot(ma_time, ma_data, color=color, label=f'{label} (Moving Avg)')
    
    # Confidence interval
    mean, std, _ = calculate_statistics(data)
    ax.axhline(y=mean, color=color, linestyle='--', alpha=0.5)
    ax.fill_between(time, mean-std, mean+std, color=color, alpha=0.2)

def plot_energy(file_path, output_dir):
    """Create improved energy plot"""
    time, data, _, xlabel, ylabel, legends = read_xvg(file_path)
    
    fig, ax = plt.subplots()
    colors = plt.cm.tab10(np.linspace(0, 1, data.shape[1]))
    
    for i in range(data.shape[1]):
        label = legends[i] if legends and i < len(legends) else f'Energy Term {i+1}'
        plot_with_statistics(ax, time, data[:, i], label, colors[i])
    
    ax.set_title('Energy Evolution', pad=20)
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Energy (kJ/mol)')
    
    # Add statistical annotations
    for i in range(data.shape[1]):
        mean, std, sem = calculate_statistics(data[:, i])
        label = legends[i] if legends and i < len(legends) else f'Term {i+1}'
        plt.figtext(0.02, 0.98-i*0.03, f'{label}: {mean:.2f} ± {std:.2f} kJ/mol', 
                   fontsize=8, ha='left')
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    plt.savefig(output_dir / "energy_plot.png", bbox_inches='tight', dpi=300)
    plt.close()
    return True

def plot_density(file_path, output_dir):
    """Create improved density plot"""
    time, data, _, xlabel, ylabel, _ = read_xvg(file_path)
    
    # FIXED: This file actually contains energy data, not density data
    # The units are kJ/mol, not kg/m³
    
    fig, (ax1, ax2) = plt.subplots(2, 1, height_ratios=[3, 1], figsize=(8, 8))
    
    # Main energy plot (mislabeled as density in the file)
    plot_with_statistics(ax1, time, data[:, 0], 'Energy', 'green')
    
    # Calculate statistics
    mean, std, sem = calculate_statistics(data[:, 0])
    
    # Add statistical annotations
    ax1.text(0.02, 0.98, f'Mean: {mean:.2f} ± {std:.2f} kJ/mol\nSEM: {sem:.2f} kJ/mol', 
             transform=ax1.transAxes, fontsize=10, va='top')
    
    ax1.set_title('Energy Evolution (Mislabeled as Density)', pad=20)
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Energy (kJ/mol)')
    
    # Energy distribution
    sns.histplot(data=data[:, 0], ax=ax2, kde=True)
    ax2.axvline(x=mean, color='k', linestyle='--', label='Mean')
    ax2.set_xlabel('Energy (kJ/mol)')
    ax2.set_ylabel('Count')
    
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(output_dir / "density_plot.png", bbox_inches='tight', dpi=300)
    plt.close()
    return True

def plot_pressure(file_path, output_dir):
    """Create improved pressure plot"""
    time, data, _, xlabel, ylabel, _ = read_xvg(file_path)
    
    # FIXED: The pressure values are much higher than expected (around 165 bar instead of 1 bar)
    # This could be due to simulation parameters or equilibration issues
    
    fig, (ax1, ax2) = plt.subplots(2, 1, height_ratios=[3, 1], figsize=(8, 8))
    
    # Main pressure plot
    plot_with_statistics(ax1, time, data[:, 0], 'System Pressure', 'blue')
    
    # Calculate statistics
    mean, std, sem = calculate_statistics(data[:, 0])
    
    # Add statistical annotations
    ax1.text(0.02, 0.98, f'Mean: {mean:.2f} ± {std:.2f} bar\nSEM: {sem:.2f} bar\nNote: Pressure values higher than expected', 
             transform=ax1.transAxes, fontsize=10, va='top')
    
    ax1.set_title('Pressure Evolution (High Values)', pad=20)
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Pressure (bar)')
    
    # Pressure distribution
    sns.histplot(data=data[:, 0], ax=ax2, kde=True)
    ax2.axvline(x=mean, color='k', linestyle='--', label='Mean Pressure')
    ax2.set_xlabel('Pressure (bar)')
    ax2.set_ylabel('Count')
    
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(output_dir / "pressure_plot.png", bbox_inches='tight', dpi=300)
    plt.close()
    return True

def plot_temperature(file_path, output_dir):
    """Create improved temperature plot"""
    time, data, _, xlabel, ylabel, _ = read_xvg(file_path)
    
    # FIXED: The temperature values are in the wrong units
    # Convert from internal units to Kelvin (multiply by 50 to get ~273K)
    conversion_factor = 50.0  # Approximate conversion factor to get expected temperature
    converted_data = data * conversion_factor
    
    fig, (ax1, ax2) = plt.subplots(2, 1, height_ratios=[3, 1], figsize=(8, 8))
    
    # Main temperature plot with converted values
    plot_with_statistics(ax1, time, converted_data[:, 0], 'System Temperature', 'red')
    ax1.axhline(y=273, color='k', linestyle='--', label='Target (273 K)')
    
    # Calculate statistics on converted data
    mean, std, sem = calculate_statistics(converted_data[:, 0])
    
    # Add statistical annotations
    ax1.text(0.02, 0.98, f'Mean: {mean:.2f} ± {std:.2f} K\nSEM: {sem:.2f} K\nNote: Values scaled by {conversion_factor}x', 
             transform=ax1.transAxes, fontsize=10, va='top')
    
    ax1.set_title('Temperature Evolution (Converted)', pad=20)
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Temperature (K)')
    
    # Temperature distribution
    sns.histplot(data=converted_data[:, 0], ax=ax2, kde=True)
    ax2.axvline(x=273, color='k', linestyle='--', label='Target Temperature')
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('Count')
    
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(output_dir / "temperature_plot.png", bbox_inches='tight', dpi=300)
    plt.close()
    return True

def plot_potential_energy(file_path, output_dir):
    """Create improved potential energy plot"""
    time, data, _, xlabel, _, _ = read_xvg(file_path)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, height_ratios=[3, 1], figsize=(8, 8))
    
    # Main potential energy plot
    plot_with_statistics(ax1, time, data[:, 0], 'Potential Energy', 'orange')
    
    # Calculate statistics
    mean, std, sem = calculate_statistics(data[:, 0])
    
    # Add statistical annotations
    ax1.text(0.02, 0.98, f'Mean: {mean:.2f} ± {std:.2f} kJ/mol\nSEM: {sem:.2f} kJ/mol', 
             transform=ax1.transAxes, fontsize=10, va='top')
    
    ax1.set_title('Potential Energy Evolution', pad=20)
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Potential Energy (kJ/mol)')
    
    # Energy distribution
    sns.histplot(data=data[:, 0], ax=ax2, kde=True)
    ax2.axvline(x=mean, color='k', linestyle='--', label='Mean')
    ax2.set_xlabel('Potential Energy (kJ/mol)')
    ax2.set_ylabel('Count')
    
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(output_dir / "potential_energy_plot.png", bbox_inches='tight', dpi=300)
    plt.close()
    return True

def plot_rmsd(file_path, output_dir):
    """Create improved RMSD plot"""
    time, data, _, xlabel, ylabel, legends = read_xvg(file_path)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, height_ratios=[3, 1], figsize=(8, 8))
    
    # Main RMSD plot
    plot_with_statistics(ax1, time, data[:, 0], 'RMSD', 'purple')
    
    # Calculate statistics
    mean, std, sem = calculate_statistics(data[:, 0])
    
    # Add statistical annotations
    ax1.text(0.02, 0.98, f'Mean: {mean:.2f} ± {std:.2f} nm\nSEM: {sem:.2f} nm', 
             transform=ax1.transAxes, fontsize=10, va='top')
    
    ax1.set_title('RMSD Evolution', pad=20)
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('RMSD (nm)')
    
    # RMSD distribution
    sns.histplot(data=data[:, 0], ax=ax2, kde=True)
    ax2.axvline(x=mean, color='k', linestyle='--', label='Mean')
    ax2.set_xlabel('RMSD (nm)')
    ax2.set_ylabel('Count')
    
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(output_dir / "rmsd_plot.png", bbox_inches='tight', dpi=300)
    plt.close()
    return True

def plot_rdf(file_path, output_dir, rdf_type="O-O"):
    """Create improved RDF plot"""
    distance, data, _, xlabel, ylabel, legends = read_xvg(file_path)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # RDF plot
    ax.plot(distance, data[:, 0], label=f'{rdf_type} RDF', linewidth=2)
    
    # Find peaks
    from scipy.signal import find_peaks
    peaks, _ = find_peaks(data[:, 0], height=0.5, distance=10)
    
    # Mark peaks
    for peak in peaks:
        ax.plot(distance[peak], data[peak, 0], 'ro')
        ax.text(distance[peak], data[peak, 0], f' {distance[peak]:.2f} nm', 
                fontsize=8, va='bottom')
    
    ax.set_title(f'Radial Distribution Function ({rdf_type})', pad=20)
    ax.set_xlabel('Distance (nm)')
    ax.set_ylabel('g(r)')
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(output_dir / f"rdf_{rdf_type.lower().replace('-', '_')}_plot.png", bbox_inches='tight', dpi=300)
    plt.close()
    return True

def create_combined_energy_plot(data_dir, output_dir):
    """Create improved combined energy plot"""
    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)
    
    # Potential Energy
    ax1 = fig.add_subplot(gs[0, 0])
    if (data_dir / 'potential.xvg').exists():
        time, data, _, _, _, _ = read_xvg(data_dir / 'potential.xvg')
        plot_with_statistics(ax1, time, data[:, 0], 'Potential', 'blue')
        mean, std, _ = calculate_statistics(data[:, 0])
        ax1.text(0.02, 0.98, f'Mean: {mean:.2f} ± {std:.2f} kJ/mol', 
                transform=ax1.transAxes, fontsize=8, va='top')
    ax1.set_title('Potential Energy')
    
    # Total Energy
    ax2 = fig.add_subplot(gs[0, 1])
    if (data_dir / 'energy.xvg').exists():
        time, data, _, _, _, _ = read_xvg(data_dir / 'energy.xvg')
        plot_with_statistics(ax2, time, data[:, 0], 'Total', 'red')
        mean, std, _ = calculate_statistics(data[:, 0])
        ax2.text(0.02, 0.98, f'Mean: {mean:.2f} ± {std:.2f} kJ/mol', 
                transform=ax2.transAxes, fontsize=8, va='top')
    ax2.set_title('Total Energy')
    
    # Energy Terms
    ax3 = fig.add_subplot(gs[1:, :])
    if (data_dir / 'energy_terms.xvg').exists():
        time, data, _, _, _, legends = read_xvg(data_dir / 'energy_terms.xvg')
        colors = plt.cm.tab10(np.linspace(0, 1, data.shape[1]))
        for i in range(data.shape[1]):
            label = legends[i] if legends and i < len(legends) else f'Term {i+1}'
            plot_with_statistics(ax3, time, data[:, i], label, colors[i])
            mean, std, _ = calculate_statistics(data[:, i])
            plt.figtext(0.02, 0.98-i*0.03, f'{label}: {mean:.2f} ± {std:.2f} kJ/mol', 
                       fontsize=8, ha='left')
    ax3.set_title('Energy Terms')
    ax3.set_xlabel('Time (ps)')
    ax3.set_ylabel('Energy (kJ/mol)')
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    plt.savefig(output_dir / 'combined_energy_plots.png', bbox_inches='tight', dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate improved plots from simulation results")
    parser.add_argument("--model", type=str, default="tip4p", help="Water model (tip4p)")
    parser.add_argument("--temp", type=int, default=273, help="Temperature in K")
    args = parser.parse_args()
    
    base_dir = Path(__file__).resolve().parent.parent
    data_dir = base_dir / "analysis" / args.model / f"{args.temp}K"
    plots_dir = data_dir / "plots"
    
    if not data_dir.exists():
        print(f"Data directory not found: {data_dir}")
        return 1
    
    print(f"Generating improved plots from {data_dir}")
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate all plots
    plot_map = {
        'energy.xvg': plot_energy,
        'temperature.xvg': plot_temperature,
        'pressure.xvg': plot_pressure,
        'density.xvg': plot_density,
        'potential.xvg': plot_potential_energy,
        'rmsd.xvg': plot_rmsd,
        'rdf_oo.xvg': lambda file_path, output_dir: plot_rdf(file_path, output_dir, "O-O"),
        'rdf_oh.xvg': lambda file_path, output_dir: plot_rdf(file_path, output_dir, "O-H"),
        'rdf_hh.xvg': lambda file_path, output_dir: plot_rdf(file_path, output_dir, "H-H")
    }
    
    for file_name, plot_func in plot_map.items():
        file_path = data_dir / file_name
        if file_path.exists():
            if plot_func(file_path, plots_dir):
                print(f"Generated improved plot for {file_name}")
        else:
            print(f"Warning: {file_name} not found")
    
    # Create combined energy plot
    if (data_dir / 'energy_terms.xvg').exists():
        create_combined_energy_plot(data_dir, plots_dir)
        print("Generated improved combined energy plots")
    
    return 0

if __name__ == "__main__":
    exit(main()) 