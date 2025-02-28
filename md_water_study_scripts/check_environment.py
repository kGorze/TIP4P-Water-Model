#!/usr/bin/env python3

import os
import sys
import subprocess
import shutil
import argparse
from pathlib import Path

def run_command(command, check=True):
    """Run a shell command and return the output."""
    try:
        result = subprocess.run(
            command,
            shell=True,
            check=check,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.CalledProcessError as e:
        return False, "", str(e)

def check_python_version():
    """Check if Python version is 3.6+."""
    version = sys.version_info
    if version.major < 3 or (version.major == 3 and version.minor < 6):
        print(f"❌ Python version {version.major}.{version.minor} detected. Version 3.6+ is required.")
        return False
    else:
        print(f"✅ Python version {version.major}.{version.minor}.{version.micro} detected.")
        return True

def check_python_modules():
    """Check if required Python modules are installed."""
    required_modules = ["argparse", "pathlib", "subprocess", "os", "shutil", "time"]
    missing_modules = []
    
    for module in required_modules:
        try:
            __import__(module)
            print(f"✅ Python module '{module}' is installed.")
        except ImportError:
            print(f"❌ Python module '{module}' is missing.")
            missing_modules.append(module)
    
    return len(missing_modules) == 0

def check_gromacs():
    """Check if GROMACS is installed and accessible."""
    success, stdout, stderr = run_command("gmx --version", check=False)
    
    if success:
        version = stdout.strip().split("\n")[0]
        print(f"✅ GROMACS is installed: {version}")
        return True
    else:
        print("❌ GROMACS is not installed or not in PATH.")
        return False

def check_packmol():
    """Check if PackMol is installed and accessible."""
    packmol_path = shutil.which("packmol")
    
    if packmol_path:
        print(f"✅ PackMol is installed at: {packmol_path}")
        return True
    else:
        print("❌ PackMol is not installed or not in PATH.")
        return False

def check_julia():
    """Check if Julia is installed and accessible."""
    success, stdout, stderr = run_command("julia --version", check=False)
    
    if success:
        version = stdout.strip()
        print(f"✅ Julia is installed: {version}")
        return True
    else:
        print("❌ Julia is not installed or not in PATH.")
        return False

def check_directory_structure():
    """Check if the required directory structure exists."""
    base_dir = Path(__file__).parent.parent
    required_dirs = [
        "scripts",
        "configs",
        "data",
        "analysis"
    ]
    
    missing_dirs = []
    for dir_name in required_dirs:
        dir_path = base_dir / dir_name
        if dir_path.exists() and dir_path.is_dir():
            print(f"✅ Directory '{dir_name}/' exists.")
        else:
            print(f"❌ Directory '{dir_name}/' is missing.")
            missing_dirs.append(dir_name)
    
    return len(missing_dirs) == 0

def check_config_files():
    """Check if required configuration files exist."""
    base_dir = Path(__file__).parent.parent
    config_dir = base_dir / "configs"
    
    if not config_dir.exists():
        print("❌ Config directory does not exist.")
        return False
    
    required_files = [
        "em.mdp",
        "nvt.mdp",
        "md.mdp"
    ]
    
    missing_files = []
    for file_name in required_files:
        file_path = config_dir / file_name
        if file_path.exists() and file_path.is_file():
            print(f"✅ Config file '{file_name}' exists.")
        else:
            print(f"❌ Config file '{file_name}' is missing.")
            missing_files.append(file_name)
    
    return len(missing_files) == 0

def check_script_permissions():
    """Check if all Python scripts are executable."""
    base_dir = Path(__file__).parent
    python_scripts = list(base_dir.glob("*.py"))
    
    non_executable_scripts = []
    for script in python_scripts:
        if os.access(script, os.X_OK):
            print(f"✅ Script '{script.name}' is executable.")
        else:
            print(f"❌ Script '{script.name}' is not executable.")
            non_executable_scripts.append(script.name)
    
    return len(non_executable_scripts) == 0

def main():
    parser = argparse.ArgumentParser(description="Check environment and dependencies for MD water study.")
    parser.add_argument("--fix", action="store_true", help="Attempt to fix issues automatically")
    args = parser.parse_args()
    
    print("\n=== Checking Python Environment ===")
    python_ok = check_python_version()
    modules_ok = check_python_modules()
    
    print("\n=== Checking External Dependencies ===")
    gromacs_ok = check_gromacs()
    packmol_ok = check_packmol()
    julia_ok = check_julia()
    
    print("\n=== Checking Project Structure ===")
    dirs_ok = check_directory_structure()
    configs_ok = check_config_files()
    
    print("\n=== Checking Script Permissions ===")
    permissions_ok = check_script_permissions()
    
    # Fix permissions if requested
    if args.fix and not permissions_ok:
        print("\n=== Fixing Script Permissions ===")
        base_dir = Path(__file__).parent
        success, stdout, stderr = run_command(f"chmod +x {base_dir}/*.py")
        if success:
            print("✅ Fixed script permissions.")
        else:
            print(f"❌ Failed to fix script permissions: {stderr}")
    
    # Create missing directories if requested
    if args.fix and not dirs_ok:
        print("\n=== Creating Missing Directories ===")
        base_dir = Path(__file__).parent.parent
        for dir_name in ["scripts", "configs", "data", "analysis"]:
            dir_path = base_dir / dir_name
            if not dir_path.exists():
                try:
                    dir_path.mkdir(parents=True)
                    print(f"✅ Created directory '{dir_name}/'.")
                except Exception as e:
                    print(f"❌ Failed to create directory '{dir_name}/': {e}")
    
    # Summary
    print("\n=== Environment Check Summary ===")
    all_ok = python_ok and modules_ok and gromacs_ok and packmol_ok and julia_ok and dirs_ok and configs_ok and permissions_ok
    
    if all_ok:
        print("✅ All checks passed! Your environment is ready for MD water simulations.")
    else:
        print("❌ Some checks failed. Please address the issues above before running simulations.")
        
        # Provide specific recommendations
        if not gromacs_ok:
            print("\nTo install GROMACS, follow the instructions at: http://manual.gromacs.org/current/install-guide/index.html")
        
        if not packmol_ok:
            print("\nTo install PackMol, follow the instructions at: http://leandro.iqm.unicamp.br/m3g/packmol/download.shtml")
        
        if not julia_ok:
            print("\nTo install Julia, visit: https://julialang.org/downloads/")
        
        if not permissions_ok:
            print("\nTo make all scripts executable, run: chmod +x scripts/*.py")
    
    return 0 if all_ok else 1

if __name__ == "__main__":
    sys.exit(main()) 