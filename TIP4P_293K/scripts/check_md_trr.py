#!/usr/bin/env python3
"""
Check the number of frames in md.trr file
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
    
    # Check md.trr file
    tpr_path = os.path.join(parent_dir, 'md.tpr')
    traj_path = os.path.join(parent_dir, 'md.trr')
    
    if os.path.exists(tpr_path) and os.path.exists(traj_path):
        count_frames(tpr_path, traj_path)
    else:
        if not os.path.exists(tpr_path):
            print(f"Warning: {tpr_path} does not exist")
        if not os.path.exists(traj_path):
            print(f"Warning: {traj_path} does not exist")

if __name__ == "__main__":
    main() 