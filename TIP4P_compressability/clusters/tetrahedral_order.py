#!/usr/bin/env python3
"""
Water Structure Analysis Script
================================

This script combines different methods for analyzing water structure in MD trajectories:
1. Tetrahedral parameter (q) - a measure of local water ordering.
2. RDF (Radial Distribution Function) - radial distributions (g(r)) between selected atoms (e.g., O-O, O-H).
3. Hydrogen bond analysis - number of H-bonds, their lifetime, etc.

The code assumes 'md.tpr' and 'md.xtc' files in the '../../tip4p_273K' directory.
Adjust paths and selector names according to your configuration.
"""

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import os
import random
import json
import logging
import time
# Use imageio.v2 to avoid deprecation warnings
import imageio.v2 as imageio
from PIL import Image
import matplotlib as mpl

from mpl_toolkits.mplot3d import Axes3D  # necessary for 3D plots
from matplotlib import cm
from scipy.spatial import cKDTree

# Modules from MDAnalysis for RDF and H-bonds
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis


###############################################################################
# LOGGING CONFIGURATION
###############################################################################
log_dir = "tetrahedral_logs"
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(os.path.join(log_dir, "water_structure_analysis.log")),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("water_structure_analysis")

# Configure matplotlib for better quality 3D plots
mpl.rcParams['figure.figsize'] = [10.0, 8.0]
mpl.rcParams['figure.dpi'] = 100
mpl.rcParams['savefig.dpi'] = 200


###############################################################################
# FILE PATHS
###############################################################################
# Using relative path from the script execution location
base_dir = os.path.join("..", "..", "tip4p_273K")
tpr_file = os.path.join(base_dir, "md.tpr")
xtc_file = os.path.join(base_dir, "md.xtc")

# Check if files exist
if not os.path.exists(tpr_file):
    logger.error(f"Error: File {tpr_file} does not exist!")
    logger.error(f"Current directory: {os.getcwd()}")
    exit(1)
if not os.path.exists(xtc_file):
    logger.error(f"Error: File {xtc_file} does not exist!")
    exit(1)

logger.info(f"Using files: {tpr_file} and {xtc_file}")


###############################################################################
# LOADING TRAJECTORY
###############################################################################
start_time = time.time()
u = mda.Universe(tpr_file, xtc_file)
logger.info(f"Loaded trajectory in {time.time() - start_time:.2f} sec.")


###############################################################################
# WATER ATOM SELECTORS
###############################################################################
# For TIP4P: oxygen `OW`, hydrogens `HW1`, `HW2`
oxygen_atoms = u.select_atoms("name OW or name O*")
water_group = u.select_atoms("resname SOL or resname TIP4")  # all water atoms

logger.info(f"Number of water oxygen atoms: {len(oxygen_atoms)}")


###############################################################################
# OUTPUT DIRECTORY
###############################################################################
output_dir = "water_analysis_results"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


###############################################################################
# PARAMETER SETTINGS
###############################################################################
# 1) Tetrahedral order
q_threshold = 0.78       # threshold above which water is considered strongly tetrahedral
min_cluster_size = 15     # minimum number of molecules in a cluster to report
max_frames_to_analyze = 250  # maximum number of frames to analyze

# 2) RDF
rdf_range = (0.0, 10.0)  # range [Å] for RDF calculations
rdf_nbins = 200          # number of bins

# 3) Hydrogen bond analysis
hbond_distance = 3.5     # [Å]
hbond_angle = 150        # [degrees]

# 4) Other parameters
max_molecules_to_show = 500  # limit on the number of molecules to visualize in clusters
animation_frame_stride = 5   # use every n-th frame for the animation to save memory
max_animation_frames = 100   # maximum number of frames to include in animation


###############################################################################
# HELPER FUNCTIONS
###############################################################################
def calculate_rdf(universe, sel1, sel2, range=(0.0,10.0), nbins=200, start=None, stop=None, step=None, label="g(r)"):
    """
    Calculates RDF between two groups of atoms (sel1 and sel2) in a given trajectory.
    Returns InterRDF object which can be further used for plotting.
    """
    logger.info(f"Calculating {label}...")
    rdf = InterRDF(sel1, sel2, range=range, nbins=nbins, universe=universe,
                   start=start, stop=stop, step=step)
    rdf.run()
    return rdf

def plot_rdf(rdf, title="RDF", filename="rdf.png"):
    """Simple plot for InterRDF object."""
    plt.figure(figsize=(8,6))
    # Note: In newer MDAnalysis .bins and .rdf are in rdf.results
    # If you get deprecation warnings, you can replace rdf.bins -> rdf.results.bins, rdf.rdf -> rdf.results.rdf
    plt.plot(rdf.bins, rdf.rdf, label=title)
    plt.xlabel("r [Å]")
    plt.ylabel("g(r)")
    plt.title(title)
    plt.legend()
    plt.savefig(os.path.join(output_dir, filename), dpi=300)
    plt.close()

def calculate_hbonds(universe, donors_sel, hydrogens_sel, acceptors_sel, dist=3.5, angle=150, start=None, stop=None, step=None):
    """
    Hydrogen bond analysis: count H-bonds over time.
    Note: In newer versions of MDAnalysis we use parameters:
    donors_sel, hydrogens_sel, acceptors_sel, d_a_cutoff, d_h_a_angle_cutoff
    
    Parameters:
    -----------
    donors_sel : str
        Selection string for donors - string compatible with MDAnalysis selection syntax
    hydrogens_sel : str
        Selection string for hydrogens - string compatible with MDAnalysis selection syntax
    acceptors_sel : str
        Selection string for acceptors - string compatible with MDAnalysis selection syntax
    """
    logger.info("Starting hydrogen bond analysis...")

    h = HydrogenBondAnalysis(
        universe=universe,
        donors_sel=donors_sel,
        hydrogens_sel=hydrogens_sel,
        acceptors_sel=acceptors_sel,
        d_a_cutoff=dist,           # Zamiast distance
        d_h_a_angle_cutoff=angle,  # Zamiast angle
    )

    h.run(start=start, stop=stop, step=step)
    return h

def calculate_tetrahedral_order(positions, k=4):
    """
    Calculates tetrahedral parameter for each oxygen atom.
    Returns:
      - q_values: np.ndarray with q values
      - neighbors_list: list of lists with indices of nearest neighbors
    """
    n_atoms = len(positions)
    q_values = np.zeros(n_atoms)
    all_neighbors = []

    kdtree = cKDTree(positions)

    for i in range(n_atoms):
        # k+1 nearest neighbors (first is the atom itself)
        distances, indices = kdtree.query(positions[i], k=k+1)
        neighbors = indices[1:k+1]
        all_neighbors.append(neighbors)

        vectors = positions[neighbors] - positions[i]
        norms = np.linalg.norm(vectors, axis=1)
        vectors = vectors / norms[:, np.newaxis]

        q = 0.0
        for j in range(k):
            for l in range(j+1, k):
                cos_angle = np.dot(vectors[j], vectors[l])
                q += (cos_angle + 1/3)**2
        q = 1 - 3/8 * q
        q_values[i] = q

    return q_values, all_neighbors


def find_connected_molecules(positions, q_values, q_threshold, neighbors_list, distance_cutoff=3.2):
    """
    Finds connected groups (clusters) of oxygen atoms with q >= q_threshold,
    where "connected" means they are within < distance_cutoff.
    Returns list of clusters.
    """
    high_q_indices = np.where(q_values >= q_threshold)[0]
    if len(high_q_indices) == 0:
        return []

    kdtree = cKDTree(positions)
    high_q_neighbors = {}

    for idx in high_q_indices:
        neighbors_indices = kdtree.query_ball_point(positions[idx], distance_cutoff)
        valid_neighbors = [n for n in neighbors_indices if n in high_q_indices and n != idx]
        high_q_neighbors[idx] = valid_neighbors

    visited = set()
    clusters = []

    for idx in high_q_indices:
        if idx not in visited:
            cluster = []
            stack = [idx]
            while stack:
                current = stack.pop()
                if current not in visited:
                    visited.add(current)
                    cluster.append(current)
                    for neigh in high_q_neighbors.get(current, []):
                        if neigh not in visited:
                            stack.append(neigh)
            clusters.append(cluster)

    clusters.sort(key=len, reverse=True)
    return clusters

def visualize_tetrahedral_clusters(positions, clusters, frame_num, q_values, q_threshold,
                                   min_cluster_size=15, max_molecules_show=500):
    """
    Visualizes (3D) largest clusters with q >= q_threshold.
    Saves plot to file in output_dir.
    """
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Background: all oxygen atoms (gray dots)
    ax.scatter(
        positions[:, 0],
        positions[:, 1],
        positions[:, 2],
        color='gray',
        s=1,
        alpha=0.02,
        marker='.'
    )

    shown = 0
    sc = None

    for i, cluster in enumerate(clusters):
        if len(cluster) < min_cluster_size:
            continue
        if shown >= 3:
            break
        shown += 1

        cluster_pos = positions[cluster]
        if len(cluster) > max_molecules_show:
            indices = random.sample(list(cluster), max_molecules_show)
            cluster_pos = positions[indices]

        cvals = q_values[cluster][:len(cluster_pos)]
        sc = ax.scatter(
            cluster_pos[:, 0],
            cluster_pos[:, 1],
            cluster_pos[:, 2],
            c=cvals,
            cmap='viridis',
            s=20,
            alpha=0.8,
            marker='o',
            label=f'Cluster {i+1} ({len(cluster)} atoms)'
        )

    if sc is not None:
        cbar = plt.colorbar(sc)
        cbar.set_label('q Parameter')

    plt.title(f"Frame {frame_num} - Tetrahedral Clusters (q >= {q_threshold})")
    ax.set_xlabel("X [Å]")
    ax.set_ylabel("Y [Å]")
    ax.set_zlabel("Z [Å]")
    ax.legend()

    plt.savefig(os.path.join(output_dir, f"tetra_clusters_frame_{frame_num}.png"), dpi=300)
    plt.close()

def plot_q_distribution(q_values, frame_num, q_threshold, save=True):
    """
    Histogram of tetrahedral parameter distribution.
    """
    plt.figure(figsize=(8,6))
    mean_q = np.mean(q_values)
    median_q = np.median(q_values)
    min_qv = np.min(q_values)
    max_qv = np.max(q_values)
    std_q = np.std(q_values)
    high_q_fraction = np.sum(q_values >= q_threshold) / len(q_values)

    n, bins, patches = plt.hist(q_values, bins=50, alpha=0.7, color='blue', density=True)
    plt.axvline(x=q_threshold, color='red', linestyle='--', label=f'q Threshold={q_threshold}')
    plt.axvline(x=mean_q, color='green', linestyle='-', label=f'Mean = {mean_q:.4f}')
    plt.axvline(x=median_q, color='orange', linestyle='-', label=f'Median = {median_q:.4f}')

    plt.xlabel('Tetrahedral Parameter (q)')
    plt.ylabel('Probability Density')
    plt.title(f'Frame {frame_num} - q Distribution\n'
              f'Min={min_qv:.4f}, Max={max_qv:.4f}, Std={std_q:.4f}, frac q>={q_threshold}={high_q_fraction:.4f}')
    plt.legend()

    if save:
        plt.savefig(os.path.join(output_dir, f"q_distribution_frame_{frame_num}.png"), dpi=300)
    plt.close()

    return {
        'mean': float(mean_q),
        'median': float(median_q),
        'min': float(min_qv),
        'max': float(max_qv),
        'std': float(std_q),
        'high_q_fraction': float(high_q_fraction)
    }

def create_combined_visualization(positions, clusters, frame_num, q_values, q_threshold,
                              min_cluster_size=15, max_molecules_show=500):
    """
    Creates a combined visualization with tetrahedral clusters on top and Q distribution on bottom.
    Saves the combined image for each frame.
    """
    # Create a figure with two subplots - top for clusters, bottom for Q distribution
    # Use smaller figure size and lower DPI for animation frames to save memory
    fig = plt.figure(figsize=(8, 12))
    
    # Top plot - 3D clusters visualization
    ax_top = fig.add_subplot(2, 1, 1, projection='3d')
    
    # Background: all oxygen atoms (gray dots)
    # Reduce the number of background points to improve performance
    background_sample = np.random.choice(len(positions), size=min(1000, len(positions)), replace=False)
    ax_top.scatter(
        positions[background_sample, 0],
        positions[background_sample, 1],
        positions[background_sample, 2],
        color='gray',
        s=1,
        alpha=0.02,
        marker='.'
    )

    shown = 0
    sc = None

    for i, cluster in enumerate(clusters):
        if len(cluster) < min_cluster_size:
            continue
        if shown >= 3:
            break
        shown += 1

        cluster_pos = positions[cluster]
        if len(cluster) > max_molecules_show:
            indices = random.sample(list(cluster), max_molecules_show)
            cluster_pos = positions[indices]

        cvals = q_values[cluster][:len(cluster_pos)]
        sc = ax_top.scatter(
            cluster_pos[:, 0],
            cluster_pos[:, 1],
            cluster_pos[:, 2],
            c=cvals,
            cmap='viridis',
            s=15,  # Smaller point size
            alpha=0.8,
            marker='o',
            label=f'Cluster {i+1} ({len(cluster)} atoms)'
        )

    if sc is not None:
        cbar = plt.colorbar(sc, ax=ax_top)
        cbar.set_label('q Parameter')

    ax_top.set_title(f"Frame {frame_num} - Tetrahedral Clusters (q >= {q_threshold})")
    ax_top.set_xlabel("X [Å]")
    ax_top.set_ylabel("Y [Å]")
    ax_top.set_zlabel("Z [Å]")
    ax_top.legend(fontsize='small')
    
    # Bottom plot - Q distribution histogram
    ax_bottom = fig.add_subplot(2, 1, 2)
    
    mean_q = np.mean(q_values)
    median_q = np.median(q_values)
    min_qv = np.min(q_values)
    max_qv = np.max(q_values)
    std_q = np.std(q_values)
    high_q_fraction = np.sum(q_values >= q_threshold) / len(q_values)

    n, bins, patches = ax_bottom.hist(q_values, bins=50, alpha=0.7, color='blue', density=True)
    ax_bottom.axvline(x=q_threshold, color='red', linestyle='--', label=f'q Threshold={q_threshold}')
    ax_bottom.axvline(x=mean_q, color='green', linestyle='-', label=f'Mean = {mean_q:.4f}')
    ax_bottom.axvline(x=median_q, color='orange', linestyle='-', label=f'Median = {median_q:.4f}')

    ax_bottom.set_xlabel('Tetrahedral Parameter (q)')
    ax_bottom.set_ylabel('Probability Density')
    ax_bottom.set_title(f'Frame {frame_num} - q Distribution\n'
              f'Min={min_qv:.4f}, Max={max_qv:.4f}, Std={std_q:.4f}, q>={q_threshold}={high_q_fraction:.4f}')
    ax_bottom.legend(fontsize='small')
    
    # Adjust layout and save
    plt.tight_layout()
    combined_img_path = os.path.join(output_dir, f"combined_frame_{frame_num:04d}.png")
    plt.savefig(combined_img_path, dpi=100)  # Lower DPI for animation frames
    plt.close()
    
    return combined_img_path

def create_animation_from_frames(frames_dir, output_filename="water_clusters_animation.gif", fps=5, max_frames=max_animation_frames):
    """
    Creates an animated GIF from individual frame images.
    Uses a memory-efficient approach by processing one image at a time.
    
    Parameters:
    -----------
    frames_dir : str
        Directory containing frame images
    output_filename : str
        Output GIF filename
    fps : int
        Frames per second for the animation
    max_frames : int
        Maximum number of frames to include in the animation
    """
    logger.info(f"Creating animation from frames in {frames_dir}...")
    
    # Get list of all frame files sorted by frame number
    frame_files = [f for f in os.listdir(frames_dir) if f.startswith("combined_frame_") and f.endswith(".png")]
    frame_files.sort()
    
    if not frame_files:
        logger.error("No frame files found for animation!")
        return
    
    # Limit number of frames if needed
    if len(frame_files) > max_frames:
        # Sample frames evenly
        logger.info(f"Limiting animation to {max_frames} frames (from {len(frame_files)} available)...")
        indices = np.linspace(0, len(frame_files)-1, max_frames, dtype=int)
        frame_files = [frame_files[i] for i in indices]
    
    # Set the output path
    output_path = os.path.join(frames_dir, output_filename)
    
    try:
        # First, resize images to reduce memory usage
        resized_images_dir = os.path.join(frames_dir, "resized_frames")
        if not os.path.exists(resized_images_dir):
            os.makedirs(resized_images_dir)
            
        logger.info(f"Resizing {len(frame_files)} frames for animation...")
        resized_files = []
        
        for i, filename in enumerate(frame_files):
            if i % 10 == 0:
                logger.info(f"Resizing frame {i}/{len(frame_files)}...")
                
            file_path = os.path.join(frames_dir, filename)
            output_name = os.path.join(resized_images_dir, filename)
            
            # Use PIL for efficient resizing
            with Image.open(file_path) as img:
                # Reduce size by 50%
                new_size = (img.width // 2, img.height // 2)
                img_resized = img.resize(new_size, Image.LANCZOS)
                img_resized.save(output_name, quality=85)
                resized_files.append(output_name)
        
        # Now create the GIF using the resized images
        logger.info(f"Creating GIF with {len(resized_files)} frames at {fps} fps...")
        
        # Get first image to determine parameters
        first_img = imageio.imread(resized_files[0])
        
        # Create writer with first image
        with imageio.get_writer(output_path, mode='I', duration=1000/fps, loop=0) as writer:
            # Add first image
            writer.append_data(first_img)
            
            # Add remaining images one by one
            for i, file_path in enumerate(resized_files[1:], 1):
                if i % 10 == 0:
                    logger.info(f"Adding frame {i}/{len(resized_files)} to GIF...")
                
                try:
                    img = imageio.imread(file_path)
                    writer.append_data(img)
                except Exception as e:
                    logger.error(f"Error processing {file_path}: {str(e)}")
                    
        logger.info(f"Animation created successfully: {output_path}")
        
    except Exception as e:
        logger.error(f"Error creating animation: {str(e)}")
        
    return output_path


###############################################################################
# MAIN PART OF THE SCRIPT
###############################################################################
def main():
    logger.info("Starting comprehensive water structure analysis...")

    # 1) RDF Analysis
    # Calculate O–O and (optionally) O–H RDF
    logger.info("--- [1] RDF Analysis ---")
    sel_O = oxygen_atoms
    sel_H = u.select_atoms("name HW1 or name HW2")

    # O–O RDF
    rdf_oo = calculate_rdf(u, sel_O, sel_O, range=rdf_range, nbins=rdf_nbins, label="O–O")
    plot_rdf(rdf_oo, title="RDF O–O", filename="rdf_OO.png")

    # O–H RDF
    rdf_oh = calculate_rdf(u, sel_O, sel_H, range=rdf_range, nbins=rdf_nbins, label="O–H")
    plot_rdf(rdf_oh, title="RDF O–H", filename="rdf_OH.png")

    # 2) Hydrogen Bond Analysis
    logger.info("--- [2] H-bond Analysis ---")
    # In HydrogenBondAnalysis we need to use selection strings, not AtomGroup objects
    donors_sel = "name OW or name O*"      # Selecion of oxygen as donor
    hydrogens_sel = "name HW1 or name HW2"  # Selecion of hydrogens
    acceptors_sel = "name OW or name O*"    # Selecion of oxygen as acceptor

    # USING MODIFIED FUNCTION WITH CORRECT PARAMETERS
    hba = calculate_hbonds(u, donors_sel, hydrogens_sel, acceptors_sel, dist=hbond_distance, angle=hbond_angle)
    # Save statistics over time
    hbond_stats_file = os.path.join(output_dir, "hbond_stats.csv")
    with open(hbond_stats_file, 'w') as f:
        f.write("frame,nbonds\n")
        # In newer versions of MDAnalysis count_by_time() returns only an array of bond counts
        # without frame indices, so we need to add them ourselves
        nbonds_array = hba.count_by_time()
        for frame_i, nbonds_i in enumerate(nbonds_array):
            f.write(f"{frame_i},{nbonds_i}\n")

    # H-bonds over time plot
    nbonds_hb = hba.count_by_time()
    frames_hb = np.arange(len(nbonds_hb))  # Frame indices
    plt.figure(figsize=(8,6))
    plt.plot(frames_hb, nbonds_hb, label='Number of H-bonds')
    plt.xlabel('Frame')
    plt.ylabel('Number of H-bonds')
    plt.title('H-bonds over time')
    plt.legend()
    plt.savefig(os.path.join(output_dir, "hbond_vs_time.png"), dpi=300)
    plt.close()

    logger.info(f"Saved H-bond statistics to {hbond_stats_file}")

    # 3) Tetrahedral Order + clustering
    logger.info("--- [3] Tetrahedral Order + clustering ---")
    tetra_stats_file = os.path.join(output_dir, "tetrahedral_stats.csv")
    with open(tetra_stats_file, 'w') as f:
        f.write("frame,mean_q,median_q,min_q,max_q,std_q,high_q_fraction,largest_cluster,n_clusters\n")

    frame_counter = 0
    frame_files = []  # To collect paths of combined visualization files
    
    # Count total frames for progress reporting
    total_frames = len(u.trajectory) if max_frames_to_analyze is None else min(len(u.trajectory), max_frames_to_analyze)
    logger.info(f"Starting to process {total_frames} frames for visualization...")
    logger.info(f"Using every {animation_frame_stride}th frame for animation to save memory")

    for ts in u.trajectory:
        if ts.frame >= max_frames_to_analyze:
            break

        positions = oxygen_atoms.positions
        q_values, neighbors_list = calculate_tetrahedral_order(positions)

        mean_q = np.mean(q_values)
        median_q = np.median(q_values)
        min_qv = np.min(q_values)
        max_qv = np.max(q_values)
        std_qv = np.std(q_values)
        high_q_fraction = np.sum(q_values >= q_threshold) / len(q_values)

        clusters = find_connected_molecules(positions, q_values, q_threshold, neighbors_list)
        n_clusters = len(clusters)
        largest_cluster = len(clusters[0]) if clusters else 0

        with open(tetra_stats_file, 'a') as f:
            f.write(f"{ts.frame},{mean_q:.6f},{median_q:.6f},{min_qv:.6f},{max_qv:.6f},{std_qv:.6f},"
                    f"{high_q_fraction:.6f},{largest_cluster},{n_clusters}\n")

        if ts.frame % 10 == 0:
            progress = (frame_counter / total_frames) * 100
            logger.info(f"Frame {ts.frame} ({progress:.1f}%): mean q={mean_q:.4f}, q>={q_threshold} fraction={high_q_fraction:.4f}, "
                        f"n_clusters={n_clusters}, largest={largest_cluster}")

        # Create combined visualization only for selected frames to save memory
        if clusters and ts.frame % animation_frame_stride == 0:
            frame_file = create_combined_visualization(positions, clusters, ts.frame, q_values, q_threshold,
                                        min_cluster_size=min_cluster_size,
                                        max_molecules_show=max_molecules_to_show)
            frame_files.append(frame_file)
        
        frame_counter += 1

    logger.info("Completed tetrahedral analysis.")
    
    # Create animation from frame files
    if frame_files:
        animation_path = create_animation_from_frames(
            output_dir, 
            "water_tetrahedral_animation.gif", 
            fps=5,
            max_frames=max_animation_frames
        )
        logger.info(f"Created animation: {animation_path}")
    else:
        logger.warning("No frames were generated for animation.")
        
    # Load data for making time-evolution summary plot
    data = np.loadtxt(tetra_stats_file, delimiter=',', skiprows=1)
    # Handle case of single frame (data may be 1D)
    if data.ndim == 1:
        data = data[np.newaxis, :]

    frames = data[:, 0]
    mean_q_vals = data[:, 1]
    median_q_vals = data[:, 2]
    high_q_vals = data[:, 6]
    largest_cluster_vals = data[:, 7]
    n_clusters_vals = data[:, 8]

    plt.figure(figsize=(10,8))
    plt.subplot(4,1,1)
    plt.plot(frames, mean_q_vals, label='Mean q')
    plt.plot(frames, median_q_vals, label='Median q')
    plt.ylabel('q')
    plt.title('Tetrahedral Parameter Time Evolution')
    plt.legend()

    plt.subplot(4,1,2)
    plt.plot(frames, high_q_vals, label=f'q>={q_threshold} Fraction')
    plt.ylabel('Fraction')
    plt.legend()

    plt.subplot(4,1,3)
    plt.plot(frames, largest_cluster_vals, label='Largest Cluster')
    plt.ylabel('Size')
    plt.legend()

    plt.subplot(4,1,4)
    plt.plot(frames, n_clusters_vals, label='Number of Clusters')
    plt.ylabel('Clusters')
    plt.xlabel('Frame')
    plt.legend()

    plt.tight_layout()
    summary_plot = os.path.join(output_dir, "tetrahedral_time_evolution.png")
    plt.savefig(summary_plot, dpi=300)
    plt.close()
    logger.info(f"Saved time-evolution summary to {summary_plot}")

    logger.info("Analysis complete!")


if __name__ == "__main__":
    main()
