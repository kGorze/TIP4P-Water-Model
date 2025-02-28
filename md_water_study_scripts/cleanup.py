#!/usr/bin/env python3

import os
import glob
import argparse
from pathlib import Path

def cleanup_simulation_dir(directory):
    """Clean up non-essential files from a simulation directory"""
    print(f"Cleaning up directory: {directory}")
    
    # Files to delete
    patterns_to_delete = [
        "step*.pdb",       # Intermediate step files
        "#*#",             # GROMACS backup files
        "*.trr",           # Trajectory files (can be large)
        "temp_water",      # Temporary directory
        "*.edr",           # Energy files (keep only final ones)
        "*.xtc",           # Compressed trajectory files (keep only final ones)
        "*.cpt",           # Checkpoint files
        "*.top.*",         # Topology backup files
        "posre.*",         # Position restraint files
        "mdout.mdp",       # Generated mdp files
        "pack.inp",        # PackMol input file
        "water.pdb"        # Template water molecule
    ]
    
    # Files to keep
    files_to_keep = [
        "md.edr",          # Final energy file
        "md.xtc",          # Final trajectory
        "md.gro",          # Final structure
        "md.log",          # Final log
        "nvt.log",         # NVT log
        "em.log",          # Energy minimization log
        "topol.top",       # Topology
        "processed.gro",   # Processed structure
        "*.mdp",           # Configuration files
        "*.xvg",           # Analysis output
        "water_box.pdb"    # Initial water box
    ]
    
    # Count files before cleanup
    total_files_before = sum(len(files) for _, _, files in os.walk(directory))
    
    # Delete files matching patterns
    deleted_count = 0
    for pattern in patterns_to_delete:
        for file_path in glob.glob(os.path.join(directory, pattern)):
            if os.path.isfile(file_path):
                # Check if this file should be kept
                filename = os.path.basename(file_path)
                keep_file = False
                for keep_pattern in files_to_keep:
                    if glob.fnmatch.fnmatch(filename, keep_pattern):
                        keep_file = True
                        break
                
                if not keep_file:
                    try:
                        os.remove(file_path)
                        deleted_count += 1
                        print(f"Deleted: {file_path}")
                    except Exception as e:
                        print(f"Error deleting {file_path}: {e}")
            elif os.path.isdir(file_path):
                try:
                    import shutil
                    shutil.rmtree(file_path)
                    deleted_count += 1
                    print(f"Deleted directory: {file_path}")
                except Exception as e:
                    print(f"Error deleting directory {file_path}: {e}")
    
    # Count files after cleanup
    total_files_after = sum(len(files) for _, _, files in os.walk(directory))
    
    print(f"Cleanup complete. Deleted {deleted_count} files/directories.")
    print(f"Files before: {total_files_before}, Files after: {total_files_after}")
    print(f"Freed approximately {(total_files_before - total_files_after) * 0.5:.2f} MB")

def main():
    parser = argparse.ArgumentParser(description="Clean up non-essential files from simulation directories")
    parser.add_argument("--model", type=str, choices=["tip4p", "tip4p-ice", "all"], default="all", 
                        help="Water model to clean up")
    parser.add_argument("--temp", type=int, choices=[150, 200, 273, 298, 0], default=0,
                        help="Temperature to clean up (0 for all)")
    args = parser.parse_args()
    
    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / "data"
    
    if args.model == "all":
        models = ["tip4p", "tip4p-ice"]
    else:
        models = [args.model]
    
    for model in models:
        model_dir = data_dir / model
        if not model_dir.exists():
            print(f"Model directory {model_dir} does not exist. Skipping.")
            continue
        
        if args.temp == 0:
            # Clean all temperature directories
            for temp_dir in model_dir.glob("*K"):
                if temp_dir.is_dir():
                    cleanup_simulation_dir(temp_dir)
        else:
            # Clean specific temperature directory
            temp_dir = model_dir / f"{args.temp}K"
            if temp_dir.exists():
                cleanup_simulation_dir(temp_dir)
            else:
                print(f"Temperature directory {temp_dir} does not exist. Skipping.")

if __name__ == "__main__":
    main() 