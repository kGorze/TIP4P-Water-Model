#!/usr/bin/env python3
"""
RMSD Analysis for TIP4P Water

This script calculates the Root Mean Square Deviation (RMSD) of water molecules 
from GROMACS trajectory files and identifies equilibration periods.

Usage:
    python rmsd_analysis.py [tpr_file] [trajectory_file]

Default:
    Uses md.tpr and md.xtc in the current directory
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from utils import load_universe, save_plot
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import align

def calculate_rmsd(universe, reference_frame=0, selection='all', n_frames=None):
    """
    Calculate RMSD of a trajectory with respect to a reference frame
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory
    reference_frame : int
        Frame to use as reference
    selection : str
        Atom selection string
    n_frames : int
        Number of frames to analyze (if None, use all frames)
        
    Returns:
    --------
    tuple
        (time_array, rmsd_array)
    """
    print("Calculating RMSD...")
    
    # Select atoms
    try:
        atoms = universe.select_atoms(selection)
        if len(atoms) == 0:
            raise ValueError(f"No atoms selected with '{selection}'")
    except Exception as e:
        print(f"Error with selection '{selection}': {e}")
        print("Trying alternative selections...")
        
        # Try common water oxygen selections for pure water systems
        for alt_sel in ['name OW', 'name O', 'resname SOL and name O*', 'resname TIP* and name O*', 'resname WAT and name O*']:
            try:
                atoms = universe.select_atoms(alt_sel)
                if len(atoms) > 0:
                    print(f"Using alternative selection: '{alt_sel}' selected {len(atoms)} atoms")
                    break
            except:
                continue
        
        if len(atoms) == 0:
            print("Could not find suitable atoms, using all atoms")
            atoms = universe.atoms
    
    print(f"Selected {len(atoms)} atoms for RMSD calculation")
    
    # Set reference coordinates
    universe.trajectory[reference_frame]
    reference_coords = atoms.positions.copy()
    
    # Calculate RMSD for each frame
    time_array = []
    rmsd_array = []
    
    # Limit the number of frames if specified
    if n_frames is not None:
        n_frames = min(n_frames, len(universe.trajectory))
        frames = np.linspace(0, len(universe.trajectory)-1, n_frames, dtype=int).tolist()
    else:
        frames = range(len(universe.trajectory))
    
    for frame_idx in frames:
        universe.trajectory[frame_idx]
        time_array.append(universe.trajectory.time)
        rmsd_array.append(rms.rmsd(atoms.positions, reference_coords, superposition=True))
    
    time_array = np.array(time_array)
    rmsd_array = np.array(rmsd_array)
    
    # Save data
    np.savetxt('../analysis/rmsd_data.dat', 
               np.column_stack((time_array, rmsd_array)),
               header='Time (ps)\tRMSD (Å)', 
               comments='# ')
    
    # Plot RMSD
    plt.figure(figsize=(12, 8))
    plt.plot(time_array, rmsd_array, 'b-', linewidth=1.5)
    
    # Calculate statistics
    mean_rmsd = np.mean(rmsd_array)
    std_rmsd = np.std(rmsd_array)
    
    # Add horizontal line for the mean
    plt.axhline(y=mean_rmsd, color='r', linestyle='--',
               label=f'Mean: {mean_rmsd:.2f} Å')
    
    # Add horizontal lines for ±1 standard deviation
    plt.axhline(y=mean_rmsd + std_rmsd, color='g', linestyle=':',
               label=f'±1σ: {std_rmsd:.2f} Å')
    plt.axhline(y=mean_rmsd - std_rmsd, color='g', linestyle=':', alpha=0.7)
    
    # Add shaded region for ±1 standard deviation
    plt.fill_between(time_array, 
                    mean_rmsd - std_rmsd, 
                    mean_rmsd + std_rmsd, 
                    color='b', alpha=0.1)
    
    plt.xlabel('Time (ps)')
    plt.ylabel('RMSD (Å)')
    plt.title('Root Mean Square Deviation (RMSD) vs. Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add text box with statistics
    textstr = '\n'.join((
        f'Mean RMSD: {mean_rmsd:.2f} Å',
        f'Std Dev: {std_rmsd:.2f} Å',
        f'Min RMSD: {np.min(rmsd_array):.2f} Å',
        f'Max RMSD: {np.max(rmsd_array):.2f} Å'))
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0.02, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    # Add annotation about RMSD for water systems
    plt.annotate('For pure water, high RMSD values (30-40 Å) are normal\nas molecules diffuse away from initial positions', 
                 xy=(0.5, 0.05), xycoords='axes fraction',
                 bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="orange", alpha=0.8),
                 ha='center')
    
    plt.tight_layout()
    plt.savefig('../analysis/rmsd_plot.png', dpi=300, bbox_inches='tight')
    
    return time_array, rmsd_array

def detect_equilibration(time_array, rmsd_array, window_size=100):
    """
    Detect equilibration point in RMSD data
    
    Parameters:
    -----------
    time_array : numpy.ndarray
        Array of time points
    rmsd_array : numpy.ndarray
        Array of RMSD values
    window_size : int
        Size of the window for moving average
        
    Returns:
    --------
    tuple
        (equilibration_time, equilibration_index)
    """
    print("Detecting equilibration point...")
    
    # Ensure arrays are numpy arrays
    time_array = np.array(time_array)
    rmsd_array = np.array(rmsd_array)
    
    # Calculate moving average
    if len(rmsd_array) <= window_size:
        print(f"Warning: RMSD array length ({len(rmsd_array)}) is less than window size ({window_size})")
        window_size = max(10, len(rmsd_array) // 10)
        print(f"Using window size of {window_size} instead")
    
    # Calculate moving average and standard deviation
    moving_avg = np.convolve(rmsd_array, np.ones(window_size)/window_size, mode='valid')
    
    # Pad the beginning of the moving average to match the original array length
    padding = len(rmsd_array) - len(moving_avg)
    moving_avg = np.pad(moving_avg, (padding, 0), 'edge')
    
    # Calculate the derivative of the moving average
    derivative = np.gradient(moving_avg)
    
    # Find where the derivative becomes small and stays small
    # This indicates that the RMSD has plateaued
    derivative_threshold = 0.01 * np.max(np.abs(derivative))  # 1% of max derivative
    
    # Find regions where the derivative is below the threshold
    stable_regions = np.abs(derivative) < derivative_threshold
    
    # Find the first point where the system is stable for at least window_size frames
    equilibration_index = None
    for i in range(len(stable_regions) - window_size):
        if np.all(stable_regions[i:i+window_size]):
            equilibration_index = i
            break
    
    # If no stable region found, use the point where the derivative is smallest
    if equilibration_index is None:
        print("Warning: Could not find a stable region, using point of minimum derivative")
        equilibration_index = np.argmin(np.abs(derivative))
    
    equilibration_time = time_array[equilibration_index]
    
    print(f"Detected equilibration at {equilibration_time:.2f} ps (frame {equilibration_index})")
    
    # Calculate statistics for the equilibrated region
    equilibrated_rmsd = rmsd_array[equilibration_index:]
    mean_rmsd = np.mean(equilibrated_rmsd)
    std_rmsd = np.std(equilibrated_rmsd)
    
    print(f"Equilibrated RMSD: {mean_rmsd:.2f} ± {std_rmsd:.2f} Å")
    
    # Create a plot showing the equilibration detection
    plt.figure(figsize=(12, 8))
    
    # Plot the RMSD data
    plt.plot(time_array, rmsd_array, 'b-', alpha=0.5, label='RMSD')
    
    # Plot the moving average
    plt.plot(time_array, moving_avg, 'r-', linewidth=2, label=f'Moving average (window={window_size})')
    
    # Mark the equilibration point
    plt.axvline(x=equilibration_time, color='g', linestyle='--', 
               label=f'Equilibration: {equilibration_time:.2f} ps')
    
    # Add horizontal line for the mean equilibrated RMSD
    plt.axhline(y=mean_rmsd, color='k', linestyle=':', 
               label=f'Mean equilibrated RMSD: {mean_rmsd:.2f} Å')
    
    # Shade the equilibrated region
    plt.axvspan(equilibration_time, time_array[-1], alpha=0.2, color='g', 
               label='Equilibrated region')
    
    plt.xlabel('Time (ps)')
    plt.ylabel('RMSD (Å)')
    plt.title('RMSD Equilibration Detection')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add text box with statistics
    textstr = '\n'.join((
        f'Equilibration time: {equilibration_time:.2f} ps',
        f'Equilibrated RMSD: {mean_rmsd:.2f} ± {std_rmsd:.2f} Å',
        f'Min RMSD: {np.min(equilibrated_rmsd):.2f} Å',
        f'Max RMSD: {np.max(equilibrated_rmsd):.2f} Å'))
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0.02, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig('../analysis/rmsd_equilibration_plot.png', dpi=300, bbox_inches='tight')
    
    # Create a more detailed analysis plot
    plt.figure(figsize=(12, 12))
    
    # Plot 1: RMSD with equilibration
    plt.subplot(2, 1, 1)
    plt.plot(time_array, rmsd_array, 'b-', alpha=0.5, label='RMSD')
    plt.plot(time_array, moving_avg, 'r-', linewidth=2, label=f'Moving average')
    plt.axvline(x=equilibration_time, color='g', linestyle='--', 
               label=f'Equilibration: {equilibration_time:.2f} ps')
    plt.axhline(y=mean_rmsd, color='k', linestyle=':', 
               label=f'Mean: {mean_rmsd:.2f} Å')
    plt.axvspan(equilibration_time, time_array[-1], alpha=0.2, color='g')
    
    plt.xlabel('Time (ps)')
    plt.ylabel('RMSD (Å)')
    plt.title('RMSD vs. Time with Equilibration Detection')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot 2: Derivative of moving average
    plt.subplot(2, 1, 2)
    plt.plot(time_array, derivative, 'b-', label='Derivative of moving average')
    plt.axvline(x=equilibration_time, color='g', linestyle='--', 
               label=f'Equilibration: {equilibration_time:.2f} ps')
    plt.axhline(y=derivative_threshold, color='r', linestyle=':', 
               label=f'Threshold: {derivative_threshold:.4f}')
    plt.axhline(y=-derivative_threshold, color='r', linestyle=':', alpha=0.5)
    plt.axhspan(-derivative_threshold, derivative_threshold, alpha=0.2, color='r',
               label='Stable region')
    
    plt.xlabel('Time (ps)')
    plt.ylabel('RMSD derivative (Å/ps)')
    plt.title('RMSD Rate of Change')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('../analysis/rmsd_derivative_plot.png', dpi=300, bbox_inches='tight')
    
    # Create a histogram of equilibrated RMSD values
    plt.figure(figsize=(10, 6))
    
    # Plot histogram
    n, bins, patches = plt.hist(equilibrated_rmsd, bins=30, alpha=0.7, 
                               color='blue', edgecolor='black', density=True)
    
    # Add vertical lines for statistics
    plt.axvline(mean_rmsd, color='r', linestyle='--',
               label=f'Mean: {mean_rmsd:.2f} Å')
    plt.axvline(mean_rmsd + std_rmsd, color='g', linestyle=':',
               label=f'±1σ: {std_rmsd:.2f} Å')
    plt.axvline(mean_rmsd - std_rmsd, color='g', linestyle=':', alpha=0.7)
    
    # Add a fitted normal distribution
    from scipy.stats import norm
    x = np.linspace(mean_rmsd - 4*std_rmsd, mean_rmsd + 4*std_rmsd, 1000)
    plt.plot(x, norm.pdf(x, mean_rmsd, std_rmsd), 'r-', 
             label='Normal distribution')
    
    plt.xlabel('RMSD (Å)')
    plt.ylabel('Probability density')
    plt.title('Distribution of Equilibrated RMSD Values')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add text about RMSD interpretation for water
    textstr = '\n'.join((
        'RMSD interpretation for water box:',
        '• High RMSD is expected as molecules diffuse',
        '• Plateau indicates structural equilibration',
        '• For pure water, absolute value depends on box size',
        '• Stability of RMSD is more important than its value'))
    
    props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5)
    plt.text(0.02, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig('../analysis/rmsd_distribution_plot.png', dpi=300, bbox_inches='tight')
    
    return equilibration_time, equilibration_index

def calculate_rmsf(universe, selection='name OW', start_frame=None, n_frames=None):
    """
    Calculate Root Mean Square Fluctuations (RMSF) for selected atoms
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory
    selection : str
        Selection string for atoms to include in RMSF calculation
    start_frame : int or None
        First frame to use (None = first frame)
    n_frames : int or None
        Number of frames to analyze (None = all frames)
    """
    from MDAnalysis.analysis import align
    
    # Select atoms
    selected_atoms = universe.select_atoms(selection)
    print(f"Calculating RMSF for {len(selected_atoms)} atoms")
    
    # Determine frame range
    if start_frame is None:
        start_frame = 0
    
    if n_frames is None:
        end_frame = len(universe.trajectory)
    else:
        end_frame = min(start_frame + n_frames, len(universe.trajectory))
    
    # Create a reference structure (first frame)
    universe.trajectory[start_frame]
    reference = universe.select_atoms(selection)
    ref_coordinates = reference.positions.copy()
    
    # Calculate average structure through iterative alignment
    avg_pos = np.zeros((len(selected_atoms), 3))
    n_frames_used = 0
    
    # First pass: align to first frame and accumulate positions
    for ts in universe.trajectory[start_frame:end_frame]:
        # Align current frame to reference
        align.alignto(selected_atoms, reference, weights="mass")
        avg_pos += selected_atoms.positions
        n_frames_used += 1
    
    # Calculate average positions
    avg_pos /= n_frames_used
    
    # Second pass: align to average and calculate RMSF
    rmsf = np.zeros(len(selected_atoms))
    
    # Reset to start frame
    universe.trajectory[start_frame]
    
    for ts in universe.trajectory[start_frame:end_frame]:
        # Align to reference again
        align.alignto(selected_atoms, reference, weights="mass")
        # Calculate distance from average position
        delta = selected_atoms.positions - avg_pos
        rmsf += np.sum(delta**2, axis=1)
    
    rmsf = np.sqrt(rmsf / n_frames_used)
    
    # Save data
    np.savetxt('../analysis/rmsf_data.dat', 
               np.column_stack((np.arange(len(rmsf)), rmsf)),
               header='Atom index\tRMSF (Å)', 
               comments='# ')
    
    # Plot RMSF
    plt.figure(figsize=(10, 6))
    plt.plot(rmsf, 'b-')
    plt.axhline(np.mean(rmsf), color='r', linestyle='--', 
                label=f'Mean: {np.mean(rmsf):.4f} Å')
    plt.xlabel('Atom Index')
    plt.ylabel('RMSF (Å)')
    plt.title(f'Root Mean Square Fluctuations ({selection})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    save_plot(plt, 'rmsf_plot.png')
    
    return rmsf

def calculate_rmsd_multi(universes, phase_frames=None, phase_names=None):
    """
    Calculate RMSD for multiple trajectories
    
    Parameters:
    -----------
    universes : list of (MDAnalysis.Universe, int, int)
        List of (universe, start_frame, end_frame) tuples for each trajectory
    phase_frames : list of int
        Number of frames in each phase
    phase_names : list of str
        Names of each phase
        
    Returns:
    --------
    tuple
        (all_times, all_rmsd) - arrays containing time and RMSD values
    """
    # Initialize arrays to store all RMSD values and times
    all_rmsd = []
    all_times = []
    
    # Track cumulative time
    cumulative_time = 0.0
    
    # Define phase durations in ns (adjust these values as needed)
    phase_durations = [0.1, 3.0, 3.0, 5.0]  # EM, NVT, NPT, MD in ns
    
    # Process each universe
    for i, (universe, start_frame, end_frame) in enumerate(universes):
        phase_name = phase_names[i] if phase_names and i < len(phase_names) else f"Phase {i+1}"
        print(f"Processing {phase_name} trajectory...")
        
        # Calculate RMSD for this phase
        time_array, rmsd_array = calculate_rmsd(
            universe=universe,
            reference_frame=0,  # Use first frame of each phase as reference
            selection='name OW'  # Use water oxygens
        )
        
        # Get phase duration
        if i < len(phase_durations):
            phase_duration = phase_durations[i]
        else:
            phase_duration = 1.0  # Default 1 ns if not specified
            
        # Adjust times to be continuous across phases
        n_frames = len(time_array)
        if n_frames > 1:
            phase_times = np.linspace(cumulative_time, 
                                    cumulative_time + phase_duration, 
                                    n_frames)
        else:
            phase_times = [cumulative_time]
            
        # Update cumulative time for next phase
        cumulative_time += phase_duration
        
        # Add to overall arrays
        all_rmsd.extend(rmsd_array)
        all_times.extend(phase_times)
        
        print(f"  Added {len(rmsd_array)} frames")
    
    # Convert to numpy arrays
    all_times = np.array(all_times)
    all_rmsd = np.array(all_rmsd)
    
    # Calculate statistics
    mean_rmsd = np.mean(all_rmsd)
    std_rmsd = np.std(all_rmsd)
    
    # Create enhanced plot
    plt.figure(figsize=(12, 7))
    
    # Plot RMSD
    plt.plot(all_times, all_rmsd, color='blue', linewidth=1.5, label='RMSD')
    
    # Add mean line
    plt.axhline(mean_rmsd, color='red', linestyle='dashed', 
                linewidth=2, label=f'Mean: {mean_rmsd:.2f} Å')
    
    # Add standard deviation bands
    plt.axhline(mean_rmsd + std_rmsd, color='red', linestyle='dotted', linewidth=1)
    plt.axhline(mean_rmsd - std_rmsd, color='red', linestyle='dotted', linewidth=1)
    plt.fill_between(all_times, mean_rmsd - std_rmsd, mean_rmsd + std_rmsd, 
                    color='red', alpha=0.1, label=f'±1σ: {std_rmsd:.2f} Å')
    
    # Add phase boundaries if provided
    if phase_names:
        current_time = 0
        for i, phase_name in enumerate(phase_names):
            if i < len(phase_durations):
                phase_duration = phase_durations[i]
                # Add vertical line at phase boundary
                if i > 0:  # Don't add line before first phase
                    plt.axvline(current_time, color='black', linestyle='-', alpha=0.5)
                
                # Add phase label
                plt.text(current_time + phase_duration/2, plt.ylim()[1] * 0.95,
                        phase_name, horizontalalignment='center',
                        bbox=dict(facecolor='white', alpha=0.7))
                
                current_time += phase_duration
    
    # Add statistics text box
    stats_text = (
        f"Mean RMSD: {mean_rmsd:.2f} Å\n"
        f"Std Dev: {std_rmsd:.2f} Å\n"
        f"Min RMSD: {np.min(all_rmsd):.2f} Å\n"
        f"Max RMSD: {np.max(all_rmsd):.2f} Å"
    )
    
    plt.text(0.97, 0.97, stats_text, transform=plt.gca().transAxes,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(facecolor='white', edgecolor='gray', alpha=0.9, 
                     boxstyle='round,pad=0.4'),
            fontsize=10)
    
    # Improve axis labels and title
    plt.xlabel('Simulation Time (ns)', fontsize=12)
    plt.ylabel('RMSD (Å)', fontsize=12)
    plt.title('TIP4P Water RMSD vs. Time', fontsize=14)
    
    # Add grid and legend
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='upper left')
    
    # Save plot
    save_plot(plt, 'rmsd_combined_plot.png')
    
    return all_times, all_rmsd

def main():
    # Get command line arguments if provided
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze RMSD from GROMACS trajectories')
    parser.add_argument('--tpr', nargs='+', default=['em.tpr', 'nvt.tpr', 'npt.tpr', 'md.tpr'],
                        help='TPR files for each simulation phase')
    parser.add_argument('--trj', nargs='+', default=['em.trr', 'nvt.trr', 'npt.trr', 'md.trr'],
                        help='Trajectory files for each simulation phase')
    parser.add_argument('--phases', nargs='+', default=['EM', 'NVT', 'NPT', 'MD'],
                        help='Names of simulation phases')
    parser.add_argument('--start-frames', nargs='+', type=int, 
                        help='Starting frame for each phase')
    parser.add_argument('--end-frames', nargs='+', type=int,
                        help='Ending frame for each phase')
    parser.add_argument('--single-phase', action='store_true',
                        help='Analyze only a single phase (md.tpr and md.xtc)')
    parser.add_argument('--frames', type=int, default=None,
                        help='Number of frames to analyze (None = all frames)')
    
    args = parser.parse_args()
    
    # If single-phase flag is set or using traditional command line format
    if args.single_phase or len(sys.argv) == 3:
        if len(sys.argv) == 3:
            tpr_file = sys.argv[1]
            trajectory_file = sys.argv[2]
        else:
            tpr_file = 'md.tpr'
            trajectory_file = 'md.xtc'
        
        print(f"Analyzing single phase with {tpr_file} and {trajectory_file}")
        
        # Load the trajectory
        universe = load_universe(tpr_file, trajectory_file)
        
        # Calculate RMSD
        print("Calculating RMSD...")
        time_array, rmsd_array = calculate_rmsd(
            universe, 
            reference_frame=0,
            selection='name OW'  # Only use oxygen atoms
        )
        
        # Detect equilibration
        print("Detecting equilibration...")
        equilibration_time, equilibration_index = detect_equilibration(
            time_array, 
            rmsd_array,
            window_size=100
        )
        
        print(f"System equilibrated at {equilibration_time:.2f} ps (frame {equilibration_index})")
        
        # Calculate RMSF using the equilibrated part of the trajectory
        print("Calculating RMSF...")
        rmsf = calculate_rmsf(
            universe,
            selection='name OW',
            start_frame=equilibration_index
        )
    else:
        # Multi-phase analysis
        # Ensure the number of TPR files, trajectory files, and phase names match
        if len(args.tpr) != len(args.trj) or (args.phases and len(args.tpr) != len(args.phases)):
            print("Error: Number of TPR files, trajectory files, and phase names must match")
            sys.exit(1)
        
        # Convert start and end frames to lists of integers if provided
        start_frames = args.start_frames if args.start_frames else [0] * len(args.tpr)
        end_frames = args.end_frames if args.end_frames else [None] * len(args.tpr)
        
        print(f"Analyzing {len(args.tpr)} phases: {', '.join(args.phases)}")
        
        # Load trajectories
        from utils import load_combined_universe
        universes, frame_counts = load_combined_universe(
            tpr_files=args.tpr,
            trajectory_files=args.trj,
            start_frames=start_frames,
            end_frames=end_frames
        )
        
        # Calculate RMSD for all phases
        print("Calculating RMSD across all phases...")
        all_times, all_rmsd = calculate_rmsd_multi(
            universes=universes,
            phase_frames=frame_counts,
            phase_names=args.phases
        )
        
        # For RMSF analysis, use only the last phase (MD)
        print("Calculating RMSF using MD phase...")
        universe, start_frame, end_frame = universes[-1]
        rmsf = calculate_rmsf(
            universe,
            selection='name OW',
            start_frame=start_frame,
            n_frames=args.frames
        )
    
    print("RMSD analysis complete!")

if __name__ == '__main__':
    main() 