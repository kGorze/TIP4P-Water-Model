#!/usr/bin/env python3
"""
Diffusion and Mean Squared Displacement (MSD) Analysis for TIP4P Water

This script calculates the MSD and diffusion coefficient of water molecules
from GROMACS trajectory files.

Usage:
    python diffusion_analysis.py [tpr_file] [trajectory_file]

Default:
    Uses md.tpr and md.xtc in the current directory
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from utils import load_universe, save_plot
import MDAnalysis as mda

# Import the appropriate MSD analysis module
# Modern MDAnalysis versions use the msd module directly
try:
    from MDAnalysis.analysis.msd import EinsteinMSD
except ImportError:
    try:
        # Older versions might have it in diffusion
        from MDAnalysis.analysis.diffusion import EinsteinMSD
    except ImportError:
        print("Error: Could not import EinsteinMSD from either msd or diffusion modules")
        print("Please install a compatible version of MDAnalysis")
        sys.exit(1)

def calculate_msd_and_diffusion(universe, selection='name OW', msd_dimension='xyz', n_frames=None):
    """
    Calculate mean squared displacement (MSD) and diffusion coefficient
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory
    selection : str
        Atom selection string
    msd_dimension : str
        Dimension for MSD calculation ('xyz', 'xy', 'z', etc.)
    n_frames : int
        Number of frames to analyze (if None, use all frames)
        
    Returns:
    --------
    tuple
        (time_array, msd_array, D, r_squared)
    """
    print("Calculating mean squared displacement (MSD)...")
    
    # Select atoms
    atoms = universe.select_atoms(selection)
    if len(atoms) == 0:
        print(f"Warning: Selection '{selection}' returned no atoms")
        # Try alternative selections for water
        for alt_sel in ['name O', 'resname SOL and name O*', 'resname TIP4 and name O*', 'resname WAT and name O*']:
            atoms = universe.select_atoms(alt_sel)
            if len(atoms) > 0:
                print(f"Using alternative selection: '{alt_sel}'")
                break
        
        if len(atoms) == 0:
            print("Error: Could not find water oxygen atoms with any selection")
            return None, None, None, None
    
    print(f"Selected {len(atoms)} atoms for MSD calculation")
    
    # Limit the number of frames if specified
    if n_frames is not None:
        n_frames = min(n_frames, len(universe.trajectory))
        frames = np.linspace(0, len(universe.trajectory)-1, n_frames, dtype=int)
    else:
        frames = None
    
    # Try to use the newer MSD module (MDAnalysis >= 2.0.0)
    try:
        from MDAnalysis.analysis import msd
        print("Using newer MDAnalysis.analysis.msd module")
        
        # Create MSD analyzer
        msd_analyzer = msd.EinsteinMSD(atoms, msd_type=msd_dimension, fft=True)
        
        # Run analysis
        if frames is not None:
            msd_analyzer.run(frames=frames)
        else:
            msd_analyzer.run()
        
        # Get results - handle different API versions
        if hasattr(msd_analyzer, 'results') and hasattr(msd_analyzer.results, 'timeseries'):
            print("Using msd_analyzer.results.timeseries")
            msd_array = msd_analyzer.results.timeseries
            
            # Try to get time array
            if hasattr(msd_analyzer, 'times'):
                print("Using msd_analyzer.times for time")
                time_array = msd_analyzer.times
            elif hasattr(msd_analyzer.results, 'times'):
                print("Using msd_analyzer.results.times for time")
                time_array = msd_analyzer.results.times
            else:
                # Generate time array based on trajectory
                print("Generating time array from trajectory")
                dt = universe.trajectory.dt
                time_array = np.arange(len(msd_array)) * dt
        
        elif hasattr(msd_analyzer, 'timeseries'):
            print("Using msd_analyzer.timeseries")
            msd_array = msd_analyzer.timeseries
            
            # Try to get time array
            if hasattr(msd_analyzer, 'times'):
                print("Using msd_analyzer.times for time")
                time_array = msd_analyzer.times
            else:
                # Generate time array based on trajectory
                print("Generating time array from trajectory")
                dt = universe.trajectory.dt
                time_array = np.arange(len(msd_array)) * dt
        
        else:
            # Last resort: check if we can access the MSD data directly
            print("Could not find standard attributes, trying alternative approaches")
            
            # Try to access the MSD data directly from the universe
            dt = universe.trajectory.dt
            n_frames_actual = len(universe.trajectory) if frames is None else len(frames)
            
            # Generate time array
            time_array = np.arange(n_frames_actual) * dt
            
            # Calculate MSD manually if needed
            print("Calculating MSD manually")
            msd_array = np.zeros(n_frames_actual)
            
            # Get reference positions
            universe.trajectory[0]
            ref_pos = atoms.positions.copy()
            
            # Calculate MSD for each frame
            for i, ts in enumerate(universe.trajectory[:n_frames_actual]):
                if i % 10 == 0:
                    print(f"  Frame {i}/{n_frames_actual}")
                
                # Calculate displacement
                disp = atoms.positions - ref_pos
                
                # Calculate MSD based on selected dimensions
                if msd_dimension == 'xyz':
                    msd_array[i] = np.mean(np.sum(disp**2, axis=1))
                elif msd_dimension == 'xy':
                    msd_array[i] = np.mean(np.sum(disp[:, :2]**2, axis=1))
                elif msd_dimension == 'z':
                    msd_array[i] = np.mean(disp[:, 2]**2)
    
    # Fall back to the older diffusion module (MDAnalysis < 2.0.0)
    except ImportError:
        try:
            from MDAnalysis.analysis import diffusion
            print("Using older MDAnalysis.analysis.diffusion module")
            
            # Create MSD analyzer
            msd_analyzer = diffusion.EinsteinMSD(atoms, msd_type=msd_dimension, fft=True)
            
            # Run analysis
            if frames is not None:
                msd_analyzer.run(frames=frames)
            else:
                msd_analyzer.run()
            
            # Get results
            if hasattr(msd_analyzer, 'timeseries'):
                msd_array = msd_analyzer.timeseries
            else:
                print("Error: Could not find MSD data in the analyzer")
                return None, None, None, None
            
            # Try to get time array
            if hasattr(msd_analyzer, 'times'):
                time_array = msd_analyzer.times
            else:
                # Generate time array based on trajectory
                dt = universe.trajectory.dt
                time_array = np.arange(len(msd_array)) * dt
        
        except ImportError:
            print("Error: Could not import either msd or diffusion module")
            return None, None, None, None
    
    # Ensure time_array and msd_array have the same length
    min_len = min(len(time_array), len(msd_array))
    time_array = time_array[:min_len]
    msd_array = msd_array[:min_len]
    
    print(f"Final time_array shape: {time_array.shape}")
    print(f"Final msd_array shape: {msd_array.shape}")
    
    # Calculate diffusion coefficient using Einstein relation: MSD = 2*d*D*t
    # where d is the dimensionality (1, 2, or 3)
    dimensionality = len(msd_dimension)
    
    # Filter out zero or negative values for log plot
    valid_indices = (time_array > 0) & (msd_array > 0)
    log_time = np.log10(time_array[valid_indices])
    log_msd = np.log10(msd_array[valid_indices])
    
    # Determine appropriate time range for fitting
    # For water, the diffusive regime typically starts after ~5-10 ps
    # and continues until ~80-90% of the trajectory
    
    # First, try to identify the diffusive regime using log-log slope analysis
    # In the diffusive regime, log(MSD) vs log(time) should have a slope of ~1
    
    # Calculate slopes using a sliding window
    window_size = max(5, len(log_time) // 20)  # Use at least 5 points or 5% of data
    slopes = []
    midpoints = []
    
    for i in range(len(log_time) - window_size):
        x = log_time[i:i+window_size]
        y = log_msd[i:i+window_size]
        slope, _, r_value, _, _ = linregress(x, y)
        slopes.append(slope)
        midpoints.append(np.mean(x))
    
    slopes = np.array(slopes)
    midpoints = np.array(midpoints)
    
    # Find regions where slope is close to 1 (diffusive regime)
    # and the R² value is good (indicating a good linear fit)
    diffusive_indices = np.where((np.abs(slopes - 1.0) < 0.15))[0]
    
    # If we found diffusive regions, use the longest continuous one
    if len(diffusive_indices) > 0:
        # Find continuous regions
        regions = []
        current_region = [diffusive_indices[0]]
        
        for i in range(1, len(diffusive_indices)):
            if diffusive_indices[i] == diffusive_indices[i-1] + 1:
                current_region.append(diffusive_indices[i])
            else:
                regions.append(current_region)
                current_region = [diffusive_indices[i]]
        
        regions.append(current_region)
        
        # Find the longest region
        longest_region = max(regions, key=len)
        
        # Convert log-space indices back to original time indices
        start_idx = longest_region[0]
        end_idx = longest_region[-1] + window_size  # Add window_size to get the end of the window
        
        # Map back to original time array indices
        start_time = 10**midpoints[start_idx]
        end_time = 10**midpoints[min(end_idx, len(midpoints)-1)]
        
        # Find the closest indices in the original time array
        fit_start_idx = np.argmin(np.abs(time_array - start_time))
        fit_end_idx = np.argmin(np.abs(time_array - end_time))
        
        # Ensure we have enough points for fitting
        if fit_end_idx - fit_start_idx < 10:
            print("Warning: Automatically detected diffusive region is too short")
            print("Falling back to default range (10-80% of trajectory)")
            # Use a reasonable default range: 10-80% of the trajectory
            fit_start_idx = int(len(time_array) * 0.1)
            fit_end_idx = int(len(time_array) * 0.8)
    else:
        print("Warning: Could not automatically detect diffusive regime")
        print("Falling back to default range (10-80% of trajectory)")
        # Use a reasonable default range: 10-80% of the trajectory
        fit_start_idx = int(len(time_array) * 0.1)
        fit_end_idx = int(len(time_array) * 0.8)
    
    # Extract the time and MSD data for fitting
    fit_time = time_array[fit_start_idx:fit_end_idx+1]
    fit_msd = msd_array[fit_start_idx:fit_end_idx+1]
    
    # Perform linear regression on the selected region
    slope, intercept, r_value, p_value, std_err = linregress(fit_time, fit_msd)
    
    # Calculate diffusion coefficient (in Å²/ps)
    D = slope / (2 * dimensionality)
    
    # Convert to standard units (10⁻⁵ cm²/s)
    D_standard = D * 0.1  # Convert from Å²/ps to 10⁻⁵ cm²/s
    
    # Calculate R-squared value
    r_squared = r_value**2
    
    print(f"Diffusion coefficient: {D_standard:.6f} × 10⁻⁵ cm²/s (R² = {r_squared:.4f})")
    print(f"Fitted region: {fit_time[0]:.1f}-{fit_time[-1]:.1f} ps")
    
    # Save data
    np.savetxt('../analysis/msd_data.dat', 
               np.column_stack((time_array, msd_array)),
               header='Time (ps)\tMSD (Å²)', 
               comments='# ')
    
    # Save slopes data for reference
    if len(slopes) > 0:
        np.savetxt('../analysis/msd_slopes.dat',
                   np.column_stack((10**midpoints, slopes)),
                   header='Time (ps)\tSlope', 
                   comments='# ')
    
    # Plot MSD with improved styling
    plt.figure(figsize=(12, 8))
    plt.plot(time_array, msd_array, 'b-', linewidth=2, label='MSD')
    
    # Plot the linear fit over the fitted region
    fit_line = intercept + slope * fit_time
    plt.plot(fit_time, fit_line, 'r--', linewidth=2,
             label=f'Linear fit ({fit_time[0]:.1f}-{fit_time[-1]:.1f} ps)\nD = {D_standard:.6f} × 10⁻⁵ cm²/s\nR² = {r_squared:.4f}')
    
    # Highlight the fitted region
    plt.axvspan(fit_time[0], fit_time[-1], alpha=0.2, color='gray', label=f'Fitted region')
    
    plt.xlabel('Time (ps)', fontsize=12, fontweight='bold')
    plt.ylabel('MSD (Å²)', fontsize=12, fontweight='bold')
    plt.title('Mean Squared Displacement vs. Time', fontsize=14, fontweight='bold')
    plt.legend(fontsize=10, framealpha=0.9)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    plt.savefig('../analysis/msd_plot.png', dpi=300, bbox_inches='tight')
    
    # Create a more detailed analysis plot
    plt.figure(figsize=(12, 12))
    
    # Plot 1: MSD with linear fit
    plt.subplot(2, 1, 1)
    plt.plot(time_array, msd_array, 'b-', linewidth=2, label='MSD')
    plt.plot(fit_time, fit_line, 'r--', linewidth=2,
             label=f'Linear fit ({fit_time[0]:.1f}-{fit_time[-1]:.1f} ps)\nD = {D_standard:.6f} × 10⁻⁵ cm²/s\nR² = {r_squared:.4f}')
    
    # Highlight the fitted region
    plt.axvspan(fit_time[0], fit_time[-1], alpha=0.2, color='gray')
    
    plt.xlabel('Time (ps)', fontsize=12, fontweight='bold')
    plt.ylabel('MSD (Å²)', fontsize=12, fontweight='bold')
    plt.title('Mean Squared Displacement vs. Time', fontsize=14, fontweight='bold')
    plt.legend(fontsize=10, framealpha=0.9)
    plt.grid(True, alpha=0.3)
    
    # Plot 2: Log-log plot to identify diffusive regime
    plt.subplot(2, 1, 2)
    
    # Plot the log-log data
    plt.scatter(log_time, log_msd, s=30, alpha=0.7, c='blue', label='MSD data')
    
    # Add reference lines for different regimes
    # Slope = 2 (ballistic regime)
    x_range = np.linspace(min(log_time), max(log_time), 100)
    plt.plot(x_range, 2*x_range + log_msd[0] - 2*log_time[0], 'g--', linewidth=2, label='Slope = 2 (Ballistic)')
    
    # Slope = 0.5 (sub-diffusive regime)
    plt.plot(x_range, 0.5*x_range + log_msd[0] - 0.5*log_time[0], 'y--', linewidth=2, label='Slope = 0.5 (Sub-diffusive)')
    
    # Slope = 1 (diffusive regime)
    plt.plot(x_range, 1*x_range + log_msd[-1] - 1*log_time[-1] - 1, 'r--', linewidth=2, label='Slope = 1 (Diffusive)')
    
    # Highlight the fitted region in log scale
    log_fit_start = np.log10(max(fit_time[0], min(time_array[valid_indices])))
    log_fit_end = np.log10(min(fit_time[-1], max(time_array[valid_indices])))
    plt.axvspan(log_fit_start, log_fit_end, alpha=0.2, color='gray', label=f'Fitted region')
    
    # Plot slopes if available
    if len(slopes) > 0:
        ax2 = plt.gca().twinx()
        ax2.plot(midpoints, slopes, 'k-', alpha=0.5, linewidth=1.5, label='Local slope')
        ax2.axhline(1.0, color='r', linestyle=':', alpha=0.5, linewidth=1.5)
        ax2.set_ylabel('Local Slope', fontsize=12, fontweight='bold')
        ax2.set_ylim(0, 2.5)
        
        # Add legend for the second y-axis
        lines, labels = plt.gca().get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc='lower right', fontsize=10, framealpha=0.9)
    else:
        plt.legend(fontsize=10, framealpha=0.9)
    
    plt.xlabel('Log Time (ps)', fontsize=12, fontweight='bold')
    plt.ylabel('Log MSD (Å²)', fontsize=12, fontweight='bold')
    plt.title('Log-Log Plot of MSD vs. Time (to identify diffusive regime)', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('../analysis/diffusion_analysis_plot.png', dpi=300, bbox_inches='tight')
    
    # Create a third plot showing the different regimes more clearly
    plt.figure(figsize=(12, 8))
    
    # Plot MSD on log-log scale with regime annotations
    plt.loglog(time_array[valid_indices], msd_array[valid_indices], 'b-', linewidth=2, label='MSD')
    
    # Identify approximate regime boundaries
    # Typically, ballistic regime is at very short times
    ballistic_end = min(time_array[valid_indices]) * 5  # First 5% of time range
    # Diffusive regime starts around the fit_time[0]
    diffusive_start = fit_time[0]
    
    # Highlight regimes
    plt.axvspan(min(time_array[valid_indices]), ballistic_end, alpha=0.2, color='green', label='Ballistic regime')
    plt.axvspan(ballistic_end, diffusive_start, alpha=0.2, color='yellow', label='Sub-diffusive regime')
    plt.axvspan(diffusive_start, max(time_array[valid_indices]), alpha=0.2, color='red', label='Diffusive regime')
    
    # Add annotations
    plt.annotate('Ballistic\nSlope ≈ 2', 
                xy=(min(time_array[valid_indices])*2, min(msd_array[valid_indices])*4),
                xytext=(min(time_array[valid_indices])*2, min(msd_array[valid_indices])*10),
                arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
                fontsize=12, fontweight='bold', ha='center')
    
    plt.annotate('Sub-diffusive\nSlope < 1', 
                xy=((ballistic_end + diffusive_start)/2, np.interp((ballistic_end + diffusive_start)/2, 
                                                                  time_array[valid_indices], msd_array[valid_indices])),
                xytext=((ballistic_end + diffusive_start)/2, np.interp((ballistic_end + diffusive_start)/2, 
                                                                      time_array[valid_indices], msd_array[valid_indices])*3),
                arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
                fontsize=12, fontweight='bold', ha='center')
    
    plt.annotate('Diffusive\nSlope = 1', 
                xy=(diffusive_start*2, np.interp(diffusive_start*2, time_array[valid_indices], msd_array[valid_indices])),
                xytext=(diffusive_start*2, np.interp(diffusive_start*2, time_array[valid_indices], msd_array[valid_indices])/3),
                arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
                fontsize=12, fontweight='bold', ha='center')
    
    plt.xlabel('Time (ps)', fontsize=12, fontweight='bold')
    plt.ylabel('MSD (Å²)', fontsize=12, fontweight='bold')
    plt.title('MSD Regimes in Log-Log Scale', fontsize=14, fontweight='bold')
    plt.legend(fontsize=10, framealpha=0.9)
    plt.grid(True, alpha=0.3, which='both')
    plt.tight_layout()
    
    plt.savefig('../analysis/msd_regimes_plot.png', dpi=300, bbox_inches='tight')
    
    # Print the diffusion coefficient in standard units
    print(f"Diffusion coefficient: {D_standard:.6e} cm²/s (R² = {r_squared:.4f})")
    
    return time_array, msd_array, D_standard, r_squared

def calculate_velocity_autocorrelation(universe, selection='name OW', max_frames=1000):
    """
    Calculate the velocity autocorrelation function (VACF) for selected atoms
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory
    selection : str
        Selection string for atoms to analyze (default: oxygen atoms)
    max_frames : int
        Maximum number of frames to analyze
    
    Returns:
    --------
    tuple
        (time, vacf_array)
    """
    # Select atoms
    atoms = universe.select_atoms(selection)
    print(f"Calculating VACF for {len(atoms)} atoms")
    
    # Limit frames if needed
    n_frames = min(len(universe.trajectory), max_frames)
    
    # Store velocities
    velocities = np.zeros((n_frames, len(atoms), 3))
    
    # Collect velocities from trajectory
    for i, ts in enumerate(universe.trajectory[:n_frames]):
        if i % 100 == 0:
            print(f"Reading frame {i}/{n_frames}")
        
        # If velocities are available in the trajectory, use them
        if hasattr(ts, 'has_velocities') and ts.has_velocities:
            velocities[i] = atoms.velocities
        # Otherwise, calculate approximate velocities from positions
        elif i > 0:
            dt = ts.time - universe.trajectory[i-1].time
            if dt > 0:
                velocities[i] = (atoms.positions - atoms.positions_old) / dt
    
    # Calculate time step
    dt = universe.trajectory[1].time - universe.trajectory[0].time if n_frames > 1 else 1.0
    
    # Calculate VACF
    # Maximum lag to consider (half the trajectory)
    max_lag = n_frames // 2
    
    # Initialize arrays
    vacf = np.zeros(max_lag)
    time = np.arange(max_lag) * dt
    
    # For each lag
    for lag in range(max_lag):
        if lag % 10 == 0:
            print(f"Computing lag {lag}/{max_lag}")
        
        # For each frame where we can calculate the product
        for i in range(n_frames - lag):
            # Dot product of velocities at time t and t+lag
            v_t = velocities[i]
            v_t_lag = velocities[i + lag]
            
            # Sum over all atoms and dimensions
            dot_products = np.sum(v_t * v_t_lag, axis=1)
            vacf[lag] += np.mean(dot_products)
    
    # Normalize by the zero-lag value
    vacf = vacf / vacf[0]
    
    # Save VACF data
    np.savetxt('../analysis/vacf_data.dat', 
               np.column_stack((time, vacf)),
               header='Time (ps)\tVACF (normalized)', 
               comments='# ')
    
    # Plot VACF
    plt.figure(figsize=(10, 6))
    plt.plot(time, vacf, 'b-')
    plt.xlabel('Time lag (ps)')
    plt.ylabel('Velocity Autocorrelation')
    plt.title('Velocity Autocorrelation Function for TIP4P Water')
    plt.grid(True, alpha=0.3)
    save_plot(plt, 'vacf_plot.png')
    
    # Calculate diffusion coefficient from VACF using Green-Kubo relation
    # D = (1/3) ∫ <v(0)·v(t)> dt
    # Integrate using trapezoidal rule
    # Note: velocities should be in Å/ps for consistency with units
    
    # Check if we have proper velocities before calculating D
    if np.abs(vacf[0] - 1.0) < 1e-6:  # If normalized correctly
        vacf_integral = np.trapz(vacf, x=time)
        diffusion_coeff_vacf = vacf_integral / 3.0  # in Å²/ps
        diffusion_coeff_vacf_cm2_s = diffusion_coeff_vacf * 1e-4  # Convert to cm²/s
        
        with open('../analysis/diffusion_coefficient_vacf.dat', 'w') as f:
            f.write(f"# Diffusion coefficient from VACF for {selection}\n")
            f.write(f"# VACF Integral: {vacf_integral:.6f} Å²/ps\n")
            f.write(f"# Diffusion coefficient: {diffusion_coeff_vacf:.6f} Å²/ps\n")
            f.write(f"# Diffusion coefficient: {diffusion_coeff_vacf_cm2_s:.6e} cm²/s\n")
        
        print(f"Diffusion coefficient from VACF: {diffusion_coeff_vacf_cm2_s:.6e} cm²/s")
    
    return time, vacf

def main():
    # Get command line arguments if provided
    if len(sys.argv) > 2:
        tpr_file = sys.argv[1]
        trajectory_file = sys.argv[2]
    else:
        tpr_file = 'md.tpr'
        trajectory_file = 'md.xtc'
    
    # Load the trajectory
    universe = load_universe(tpr_file, trajectory_file)
    
    # Calculate MSD and diffusion coefficient for water oxygens
    print("Calculating MSD and diffusion coefficient...")
    try:
        time, msd, diffusion_coeff, r_squared = calculate_msd_and_diffusion(
            universe, 
            selection='name OW',
            msd_dimension='xyz',
            n_frames=None  # Use all frames
        )
        
        print(f"Diffusion coefficient: {diffusion_coeff:.6e} cm²/s (R² = {r_squared:.4f})")
        
        # Calculate velocity autocorrelation function if velocities are available
        try:
            print("Calculating velocity autocorrelation function...")
            calculate_velocity_autocorrelation(universe, selection='name OW', max_frames=1000)
        except Exception as e:
            print(f"Could not calculate VACF: {e}")
            print("This is non-critical, continuing with other analyses.")
        
        print("Diffusion analysis complete!")
        return True
    except Exception as e:
        print(f"ERROR in diffusion analysis: {e}")
        print("Consider using an alternative method or checking MDAnalysis version compatibility.")
        import traceback
        traceback.print_exc()
        return False

if __name__ == '__main__':
    main() 