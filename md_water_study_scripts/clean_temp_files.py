#!/usr/bin/env python3

import os
import glob
import argparse
import shutil
from pathlib import Path

def clean_temp_files(directory, dry_run=False, verbose=False):
    """Clean up temporary files from a simulation directory"""
    print(f"Cleaning up directory: {directory}")
    
    # Files to delete
    patterns_to_delete = [
        "step*.pdb",       # Intermediate step files
        "#*#",             # GROMACS backup files
        "*.cpt",           # Checkpoint files (except the latest ones)
        "mdout.mdp",       # Generated mdp files
        "pack.inp",        # PackMol input file
        "water.pdb"        # Template water molecule
    ]
    
    # Files to keep
    files_to_keep = [
        "md.edr",          # Final energy file
        "md.xtc",          # Final trajectory
        "md.trr",          # Final trajectory (uncompressed)
        "md.gro",          # Final structure
        "md.log",          # Final log
        "md.cpt",          # Final checkpoint
        "md_prev.cpt",     # Previous checkpoint
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
    total_size_before = sum(os.path.getsize(os.path.join(root, file)) 
                           for root, _, files in os.walk(directory) 
                           for file in files if os.path.isfile(os.path.join(root, file)))
    
    # Delete files matching patterns
    deleted_count = 0
    deleted_size = 0
    
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
                        file_size = os.path.getsize(file_path)
                        if verbose:
                            print(f"Would delete: {file_path} ({file_size/1024/1024:.2f} MB)")
                        
                        if not dry_run:
                            os.remove(file_path)
                            deleted_count += 1
                            deleted_size += file_size
                            if verbose:
                                print(f"Deleted: {file_path} ({file_size/1024/1024:.2f} MB)")
                    except Exception as e:
                        print(f"Error processing {file_path}: {e}")
            elif os.path.isdir(file_path):
                try:
                    dir_size = sum(os.path.getsize(os.path.join(root, file)) 
                                  for root, _, files in os.walk(file_path) 
                                  for file in files if os.path.isfile(os.path.join(root, file)))
                    
                    if verbose:
                        print(f"Would delete directory: {file_path} ({dir_size/1024/1024:.2f} MB)")
                    
                    if not dry_run:
                        shutil.rmtree(file_path)
                        deleted_count += 1
                        deleted_size += dir_size
                        if verbose:
                            print(f"Deleted directory: {file_path} ({dir_size/1024/1024:.2f} MB)")
                except Exception as e:
                    print(f"Error deleting directory {file_path}: {e}")
    
    # Count files after cleanup
    if not dry_run:
        total_files_after = sum(len(files) for _, _, files in os.walk(directory))
        total_size_after = sum(os.path.getsize(os.path.join(root, file)) 
                              for root, _, files in os.walk(directory) 
                              for file in files if os.path.isfile(os.path.join(root, file)))
        
        print(f"Cleanup complete. Deleted {deleted_count} files/directories.")
        print(f"Files before: {total_files_before}, Files after: {total_files_after}")
        print(f"Size before: {total_size_before/1024/1024:.2f} MB, Size after: {total_size_after/1024/1024:.2f} MB")
        print(f"Freed approximately {deleted_size/1024/1024:.2f} MB")
    else:
        print(f"Dry run complete. Would delete {deleted_count} files/directories.")
        print(f"Would free approximately {deleted_size/1024/1024:.2f} MB")
        print("Run without --dry-run to actually delete files.")

def main():
    parser = argparse.ArgumentParser(description="Clean up temporary files from simulation directories")
    parser.add_argument("--model", type=str, choices=["tip4p", "tip4p-ice", "all"], default="all", 
                        help="Water model to clean up")
    parser.add_argument("--temp", type=int, choices=[150, 200, 273, 298, 0], default=0,
                        help="Temperature to clean up (0 for all)")
    parser.add_argument("--dry-run", action="store_true", 
                        help="Show what would be deleted without actually deleting")
    parser.add_argument("--verbose", action="store_true", 
                        help="Show detailed information about files being deleted")
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
                    clean_temp_files(temp_dir, dry_run=args.dry_run, verbose=args.verbose)
        else:
            # Clean specific temperature directory
            temp_dir = model_dir / f"{args.temp}K"
            if temp_dir.exists():
                clean_temp_files(temp_dir, dry_run=args.dry_run, verbose=args.verbose)
            else:
                print(f"Temperature directory {temp_dir} does not exist. Skipping.")

if __name__ == "__main__":
    main() 