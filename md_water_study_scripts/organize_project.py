#!/usr/bin/env python3

import os
import glob
import shutil
import argparse
from pathlib import Path
import re

def organize_project(dry_run=False):
    """
    Organize the project directory by:
    1. Creating organized directories if they don't exist
    2. Moving files to appropriate directories
    """
    base_dir = Path(__file__).parent.parent
    print(f"Base directory: {base_dir}")
    
    # 1. Organize data files
    data_dir = base_dir / "data"
    if data_dir.exists():
        for model_dir in data_dir.glob("*"):
            if not model_dir.is_dir():
                continue
                
            for temp_dir in model_dir.glob("*"):
                if not temp_dir.is_dir():
                    continue
                
                # Organize by file type
                file_categories = {
                    "inputs": [".mdp", ".top", ".inp", ".pdb"],
                    "outputs": [".gro", ".xtc", ".trr", ".edr", ".tpr", ".cpt"],
                    "logs": [".log"]
                }
                
                for category, extensions in file_categories.items():
                    for ext in extensions:
                        files = list(temp_dir.glob(f"*{ext}"))
                        if not files:
                            continue
                            
                        # Create category directory if it doesn't exist
                        category_dir = temp_dir / category
                        if not category_dir.exists() and not dry_run:
                            print(f"Creating directory: {category_dir}")
                            category_dir.mkdir(exist_ok=True)
                        
                        for file_path in files:
                            # Skip if file is already in the right directory
                            if category in file_path.parts:
                                continue
                                
                            dest_path = category_dir / file_path.name
                            print(f"Moving file: {file_path} -> {dest_path}")
                            if not dry_run:
                                # Use move instead of copy to save disk space
                                shutil.move(file_path, dest_path)
    
    # 2. Organize analysis files
    analysis_dir = base_dir / "analysis"
    if analysis_dir.exists():
        for model_dir in analysis_dir.glob("*"):
            if not model_dir.is_dir():
                continue
                
            for temp_dir in model_dir.glob("*"):
                if not temp_dir.is_dir():
                    continue
                
                # Make sure plots directory exists
                plots_dir = temp_dir / "plots"
                if not plots_dir.exists() and not dry_run:
                    print(f"Creating directory: {plots_dir}")
                    plots_dir.mkdir(exist_ok=True)
                
                # Move plot files to plots directory
                for plot_file in temp_dir.glob("*.png"):
                    if "plots" not in plot_file.parts:
                        dest_path = plots_dir / plot_file.name
                        print(f"Moving plot: {plot_file} -> {dest_path}")
                        if not dry_run:
                            shutil.move(plot_file, dest_path)
    
    print("Organization complete!")
    if dry_run:
        print("This was a dry run. No files were actually moved. Run without --dry-run to actually perform the operations.")

def main():
    parser = argparse.ArgumentParser(description="Organize the project directory")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be done without actually doing it")
    args = parser.parse_args()
    
    organize_project(dry_run=args.dry_run)

if __name__ == "__main__":
    main() 