#!/usr/bin/env python3
"""
Check the number of frames in trajectory files
"""

import os
import sys
import MDAnalysis as mda

def count_frames(tpr_file, trajectory_file):
    """Count the number of frames in a trajectory file"""
    try:
        u = mda.Universe(tpr_file, trajectory_file)
        n_frames = len(u.trajectory)
        dt = u.trajectory.dt  # Time step in ps
        total_time = n_frames * dt  # Total time in ps
        
        print(f"File: {trajectory_file}")
        print(f"  Number of frames: {n_frames}")
        print(f"  Time step: {dt} ps")
        print(f"  Total time: {total_time} ps = {total_time/1000:.2f} ns")
        print(f"  First frame time: {u.trajectory[0].time} ps")
        print(f"  Last frame time: {u.trajectory[-1].time} ps")
        
        return n_frames, dt, total_time
    except Exception as e:
        print(f"Error processing {trajectory_file}: {e}")
        return 0, 0, 0

def main():
    parent_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
    
    # List of trajectory files to check
    tpr_files = ['em.tpr', 'nvt.tpr', 'npt.tpr', 'md.tpr']
    trajectory_files = ['em.trr', 'nvt.trr', 'npt.trr', 'md.xtc']
    
    total_frames = 0
    total_time = 0
    
    for tpr, traj in zip(tpr_files, trajectory_files):
        tpr_path = os.path.join(parent_dir, tpr)
        traj_path = os.path.join(parent_dir, traj)
        
        if os.path.exists(tpr_path) and os.path.exists(traj_path):
            n_frames, dt, time = count_frames(tpr_path, traj_path)
            total_frames += n_frames
            total_time += time
        else:
            if not os.path.exists(tpr_path):
                print(f"Warning: {tpr_path} does not exist")
            if not os.path.exists(traj_path):
                print(f"Warning: {traj_path} does not exist")
    
    print("\nSummary:")
    print(f"Total frames across all trajectories: {total_frames}")
    print(f"Total simulation time: {total_time} ps = {total_time/1000:.2f} ns")

if __name__ == "__main__":
    main() 