import os
import numpy as np  # Compatible version: >=1.17.3 and <1.25.0 for most scipy operations
import matplotlib.pyplot as plt
import matplotlib as mpl
import MDAnalysis as mda
# Splitting imports to handle potential missing modules
from MDAnalysis.analysis import rdf, hydrogenbonds, rms
# Note: diffusion module is not available in newer MDAnalysis versions
# Import msd module for diffusion analysis instead
try:
    from MDAnalysis.analysis import diffusion  # For older MDAnalysis versions
except ImportError:
    try:
        from MDAnalysis.analysis import msd  # For newer MDAnalysis versions
    except ImportError:
        print("Warning: Neither diffusion nor msd module available in MDAnalysis")
        
import pandas as pd
from scipy import stats
import seaborn as sns

# Set plotting style
plt.style.use('seaborn-v0_8-whitegrid')
mpl.rcParams['figure.figsize'] = (10, 6)
mpl.rcParams['font.size'] = 12

def load_universe(tpr_file='md.tpr', trajectory_file='md.xtc'):
    """
    Load an MDAnalysis Universe from GROMACS files
    
    Parameters:
    -----------
    tpr_file : str
        Path to the .tpr file
    trajectory_file : str
        Path to the trajectory file (.xtc or .trr)
        
    Returns:
    --------
    MDAnalysis.Universe
    """
    # Adjust paths to look in parent directory if files don't exist in current directory
    parent_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
    
    # Check if files exist in current directory
    if not os.path.exists(tpr_file):
        # Try in parent directory
        parent_tpr = os.path.join(parent_dir, os.path.basename(tpr_file))
        if os.path.exists(parent_tpr):
            print(f"File {tpr_file} not found in current directory, using {parent_tpr} instead")
            tpr_file = parent_tpr
        else:
            print(f"Warning: Could not find {tpr_file} in current or parent directory")
    
    if not os.path.exists(trajectory_file):
        # Try in parent directory
        parent_trajectory = os.path.join(parent_dir, os.path.basename(trajectory_file))
        if os.path.exists(parent_trajectory):
            print(f"File {trajectory_file} not found in current directory, using {parent_trajectory} instead")
            trajectory_file = parent_trajectory
        else:
            print(f"Warning: Could not find {trajectory_file} in current or parent directory")
    
    print(f"Loading Universe from {tpr_file} and {trajectory_file}...")
    return mda.Universe(tpr_file, trajectory_file)

def save_plot(plt, filename, analysis_dir='../analysis', dpi=300):
    """
    Save a matplotlib plot to the analysis directory
    
    Parameters:
    -----------
    plt : matplotlib.pyplot
        The pyplot instance
    filename : str
        Filename for the plot
    analysis_dir : str
        Directory where plots should be saved
    dpi : int
        Resolution for saved image
    """
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    
    filepath = os.path.join(analysis_dir, filename)
    plt.tight_layout()
    plt.savefig(filepath, dpi=dpi)
    print(f"Plot saved to {filepath}")
    
def parse_edr_file(edr_file='md.edr', properties=None):
    """
    Parse GROMACS energy file (.edr) by first converting to .xvg using gmx energy
    
    Parameters:
    -----------
    edr_file : str
        Path to the .edr file
    properties : list
        List of properties to extract (if None, extracts temperature, pressure, 
        potential/kinetic/total energy)
        
    Returns:
    --------
    dict : Dictionary of pandas DataFrames with property names as keys
    """
    import subprocess
    import tempfile
    
    # Adjust path to look in parent directory if file doesn't exist in current directory
    parent_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
    
    # Check if file exists in current directory
    if not os.path.exists(edr_file):
        # Try in parent directory
        parent_edr = os.path.join(parent_dir, os.path.basename(edr_file))
        if os.path.exists(parent_edr):
            print(f"File {edr_file} not found in current directory, using {parent_edr} instead")
            edr_file = parent_edr
        else:
            print(f"Warning: Could not find {edr_file} in current or parent directory")
    
    if properties is None:
        properties = ['Temperature', 'Pressure', 'Potential', 'Kinetic', 'Total-Energy']
    
    property_data = {}
    
    for prop in properties:
        # Create temporary file
        with tempfile.NamedTemporaryFile(suffix='.xvg') as tmp:
            tmp_name = tmp.name
            
        # Use gmx energy to extract property
        echo_str = f"echo {prop}"
        command = f"echo {prop} | gmx energy -f {edr_file} -o {tmp_name}"
        
        try:
            subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # Parse the XVG file
            time_col = []
            prop_col = []
            
            with open(tmp_name, 'r') as f:
                for line in f:
                    if line.startswith('#') or line.startswith('@'):
                        continue
                    cols = line.strip().split()
                    if len(cols) >= 2:
                        time_col.append(float(cols[0]))
                        prop_col.append(float(cols[1]))
            
            # Create DataFrame
            df = pd.DataFrame({
                'Time (ps)': time_col,
                prop: prop_col
            })
            
            property_data[prop] = df
            
            # Clean up
            os.remove(tmp_name)
            
        except subprocess.CalledProcessError as e:
            print(f"Error extracting {prop}: {e}")
            property_data[prop] = None
            
    return property_data 