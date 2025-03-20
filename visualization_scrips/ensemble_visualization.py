import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle, Wedge
from matplotlib.animation import FuncAnimation
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

# Set random seed for reproducibility
np.random.seed(42)

# Create figure with subplots - now with 3 columns instead of 4
fig = plt.figure(figsize=(15, 5))
gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1])

# Common parameters
n_particles = 50
particle_size = 30
box_color = 'black'
box_linewidth = 2

# Color maps for different properties
particle_cmap = plt.cm.viridis
velocity_cmap = plt.cm.plasma
temp_cmap = LinearSegmentedColormap.from_list('temp_cmap', ['blue', 'red'])

# ----- Microcanonical Ensemble (NVE) - fixed N, V, E -----
ax1 = plt.subplot(gs[0])
ax1.set_xlim(0, 100)
ax1.set_ylim(0, 100)
ax1.set_aspect('equal')
ax1.set_title('Microcanonical (NVE)', fontsize=14, weight='bold')
ax1.set_xticks([])
ax1.set_yticks([])

# Fixed box for Microcanonical - thicker to show isolation
box = Rectangle((10, 10), 80, 80, linewidth=box_linewidth*1.5, edgecolor=box_color, facecolor='none')
ax1.add_patch(box)

# Additional outer box to emphasize isolation
isolation_box = Rectangle((8, 8), 84, 84, linewidth=1, edgecolor='black', facecolor='none', linestyle='-')
ax1.add_patch(isolation_box)
ax1.text(50, 95, "Isolated System", ha='center', fontsize=10)

# Initial positions and velocities for Microcanonical
x_nve = np.random.uniform(15, 85, n_particles)
y_nve = np.random.uniform(15, 85, n_particles)

# Velocities with a range of magnitudes
dx_nve = np.random.uniform(-2, 2, n_particles)
dy_nve = np.random.uniform(-2, 2, n_particles)
vel_magnitudes = np.sqrt(dx_nve**2 + dy_nve**2)

# Normalize color mapping to velocity magnitude
norm = Normalize(vmin=0, vmax=np.max(vel_magnitudes)*1.2)
colors_nve = velocity_cmap(norm(vel_magnitudes))

scatter_nve = ax1.scatter(x_nve, y_nve, s=particle_size, c=colors_nve, alpha=0.7)
quiver_nve = ax1.quiver(x_nve, y_nve, dx_nve, dy_nve, color='black', scale=25)

# Add a colorbar for velocity magnitude
sm = ScalarMappable(cmap=velocity_cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax1, shrink=0.7, pad=0.03)
cbar.set_label('Velocity', size=8)
cbar.ax.tick_params(labelsize=6)

# Simplified parameter indicators
ax1.text(50, 97, "Fixed: N, V, E", ha='center', fontsize=12)
ax1.text(50, 3, "T varies", ha='center', fontsize=12, color='orange')

# ----- Canonical Ensemble (NVT) - fixed N, V, T -----
ax2 = plt.subplot(gs[1])
ax2.set_xlim(0, 100)
ax2.set_ylim(0, 100)
ax2.set_aspect('equal')
ax2.set_title('Canonical (NVT)', fontsize=14, weight='bold')
ax2.set_xticks([])
ax2.set_yticks([])

# Add subtle heat bath gradient background
gradient = np.linspace(0, 1, 100).reshape(1, -1)
gradient = np.repeat(gradient, 100, axis=0)
ax2.imshow(gradient, extent=[0, 100, 0, 100], cmap=temp_cmap, alpha=0.05, aspect='auto')

# Fixed box for Canonical
box = Rectangle((10, 10), 80, 80, linewidth=box_linewidth, edgecolor=box_color, facecolor='none')
ax2.add_patch(box)

# Heat bath representation
heat_bath = Rectangle((5, 5), 90, 90, linewidth=1, edgecolor='red', facecolor='none', linestyle='--')
ax2.add_patch(heat_bath)

# Explicit thermal contact - heat transfer symbols
for i in range(8):
    theta = i * 45
    # Small heat transfer circles
    circle = Circle((50 + 45*np.cos(np.radians(theta)), 50 + 45*np.sin(np.radians(theta))), 
                   3, color='red', alpha=0.2)
    ax2.add_patch(circle)
    
    # Heat arrows - show direction of heat transfer
    if i % 2 == 0:  # Alternate directions
        heat_arrow = FancyArrowPatch((50 + 40*np.cos(np.radians(theta)), 50 + 40*np.sin(np.radians(theta))),
                                    (50 + 50*np.cos(np.radians(theta)), 50 + 50*np.sin(np.radians(theta))), 
                                    mutation_scale=8, color='red', alpha=0.6,
                                    arrowstyle='-|>', linewidth=1)
    else:
        heat_arrow = FancyArrowPatch((50 + 50*np.cos(np.radians(theta)), 50 + 50*np.sin(np.radians(theta))),
                                    (50 + 40*np.cos(np.radians(theta)), 50 + 40*np.sin(np.radians(theta))), 
                                    mutation_scale=8, color='red', alpha=0.6,
                                    arrowstyle='-|>', linewidth=1)
    ax2.add_patch(heat_arrow)

ax2.text(50, 95, "Heat Bath", ha='center', color='red', fontsize=10)

# Initial positions for Canonical
x_can = np.random.uniform(15, 85, n_particles)
y_can = np.random.uniform(15, 85, n_particles)
colors_can = particle_cmap(np.random.uniform(0.3, 0.8, n_particles))
sizes_can = particle_size * np.ones(n_particles)  # For size fluctuations
scatter_can = ax2.scatter(x_can, y_can, s=sizes_can, c=colors_can, alpha=0.7)

# Simplified parameter indicators
ax2.text(50, 97, "Fixed: N, V, T", ha='center', fontsize=12)
ax2.text(50, 3, "E varies", ha='center', fontsize=12, color='red')

# ----- NPT Ensemble - fixed N, P, T -----
ax3 = plt.subplot(gs[2])
ax3.set_xlim(0, 100)
ax3.set_ylim(0, 100)
ax3.set_aspect('equal')
ax3.set_title('NPT Ensemble', fontsize=14, weight='bold')
ax3.set_xticks([])
ax3.set_yticks([])

# Add subtle heat bath gradient background
ax3.imshow(gradient, extent=[0, 100, 0, 100], cmap=temp_cmap, alpha=0.05, aspect='auto')

# Box for NPT (will change size)
box_npt_size = 70  # Initial size
box_npt = Rectangle((50-box_npt_size/2, 50-box_npt_size/2), box_npt_size, box_npt_size, 
                   linewidth=box_linewidth, edgecolor=box_color, facecolor='none')
ax3.add_patch(box_npt)

# Volume expansion/contraction indicators (double-headed arrows)
arrow_length = 10
# Horizontal arrows
vol_arrow1 = FancyArrowPatch((50-box_npt_size/2-5, 50), (50-box_npt_size/2, 50), 
                            mutation_scale=10, color='blue', alpha=0.8,
                            arrowstyle='<->', linewidth=1.5)
vol_arrow2 = FancyArrowPatch((50+box_npt_size/2, 50), (50+box_npt_size/2+5, 50), 
                            mutation_scale=10, color='blue', alpha=0.8,
                            arrowstyle='<->', linewidth=1.5)
# Vertical arrows
vol_arrow3 = FancyArrowPatch((50, 50-box_npt_size/2-5), (50, 50-box_npt_size/2), 
                            mutation_scale=10, color='blue', alpha=0.8,
                            arrowstyle='<->', linewidth=1.5)
vol_arrow4 = FancyArrowPatch((50, 50+box_npt_size/2), (50, 50+box_npt_size/2+5), 
                            mutation_scale=10, color='blue', alpha=0.8,
                            arrowstyle='<->', linewidth=1.5)
ax3.add_patch(vol_arrow1)
ax3.add_patch(vol_arrow2)
ax3.add_patch(vol_arrow3)
ax3.add_patch(vol_arrow4)

# Heat bath indicators - similar to Canonical but fewer
for i in range(4):
    theta = i * 90
    heat_arrow = FancyArrowPatch((50 + 45*np.cos(np.radians(theta)), 50 + 45*np.sin(np.radians(theta))),
                                (50 + 50*np.cos(np.radians(theta)), 50 + 50*np.sin(np.radians(theta))), 
                                mutation_scale=8, color='red', alpha=0.6,
                                arrowstyle='-|>', linewidth=1)
    ax3.add_patch(heat_arrow)

# Initial positions for NPT
x_npt = np.random.uniform(50-box_npt_size/2+5, 50+box_npt_size/2-5, n_particles)
y_npt = np.random.uniform(50-box_npt_size/2+5, 50+box_npt_size/2-5, n_particles)
colors_npt = particle_cmap(np.random.uniform(0.3, 0.8, n_particles))
scatter_npt = ax3.scatter(x_npt, y_npt, s=particle_size, c=colors_npt, alpha=0.7)

# Simplified parameter indicators
ax3.text(50, 97, "Fixed: N, P, T", ha='center', fontsize=12)
ax3.text(50, 3, "V varies", ha='center', fontsize=12, color='blue')

# Function to update the animation
def update(frame):
    global box_npt, scatter_npt, vol_arrow1, vol_arrow2, vol_arrow3, vol_arrow4
    
    # ----- Update Microcanonical (NVE) with particle motion -----
    # In a real simulation, this would conserve energy
    x_nve_new = x_nve + dx_nve * 0.5
    y_nve_new = y_nve + dy_nve * 0.5
    
    # Reflect particles off walls to conserve energy
    mask_left = x_nve_new < 15
    mask_right = x_nve_new > 85
    mask_bottom = y_nve_new < 15
    mask_top = y_nve_new > 85
    
    dx_nve[mask_left | mask_right] *= -1
    dy_nve[mask_bottom | mask_top] *= -1
    
    x_nve_new[mask_left] = 15 + (15 - x_nve_new[mask_left])
    x_nve_new[mask_right] = 85 - (x_nve_new[mask_right] - 85)
    y_nve_new[mask_bottom] = 15 + (15 - y_nve_new[mask_bottom])
    y_nve_new[mask_top] = 85 - (y_nve_new[mask_top] - 85)
    
    # Update velocities for color mapping
    vel_magnitudes = np.sqrt(dx_nve**2 + dy_nve**2)
    colors_nve_new = velocity_cmap(norm(vel_magnitudes))
    
    scatter_nve.set_offsets(np.c_[x_nve_new, y_nve_new])
    scatter_nve.set_facecolor(colors_nve_new)
    quiver_nve.set_offsets(np.c_[x_nve_new, y_nve_new])
    quiver_nve.set_UVC(dx_nve, dy_nve)
    
    # Update global variables for next frame
    for i in range(n_particles):
        x_nve[i] = x_nve_new[i]
        y_nve[i] = y_nve_new[i]
        
    # ----- Update Canonical (NVT) -----
    x_can_new = x_can + np.random.uniform(-1, 1, n_particles)
    y_can_new = y_can + np.random.uniform(-1, 1, n_particles)
    
    # Temperature fluctuations - occasionally change particle sizes
    sizes_can_new = sizes_can.copy()
    if frame % 5 == 0:
        # Randomly select a few particles to change size
        idx = np.random.choice(n_particles, 10, replace=False)
        sizes_can_new[idx] = particle_size * np.random.uniform(0.8, 1.4, 10)
    
    # Energy fluctuation represented by different colors
    if frame % 8 == 0:
        idx = np.random.choice(n_particles, 5, replace=False)
        colors_can_new = colors_can.copy()
        for i in idx:
            colors_can_new[i] = particle_cmap(np.random.uniform(0.3, 0.8))
        scatter_can.set_facecolor(colors_can_new)
    
    # Keep particles within boundaries
    x_can_new = np.clip(x_can_new, 15, 85)
    y_can_new = np.clip(y_can_new, 15, 85)
    
    scatter_can.set_offsets(np.c_[x_can_new, y_can_new])
    scatter_can.set_sizes(sizes_can_new)
    
    # ----- Update NPT (changing volume) -----
    box_npt_size = 70 + 10 * np.sin(frame/10)
    box_npt.remove()
    box_npt = Rectangle((50-box_npt_size/2, 50-box_npt_size/2), box_npt_size, box_npt_size, 
                       linewidth=box_linewidth, edgecolor=box_color, facecolor='none')
    ax3.add_patch(box_npt)
    
    # Update volume arrows
    vol_arrow1.remove()
    vol_arrow2.remove()
    vol_arrow3.remove()
    vol_arrow4.remove()
    
    vol_arrow1 = FancyArrowPatch((50-box_npt_size/2-5, 50), (50-box_npt_size/2, 50), 
                                mutation_scale=10, color='blue', alpha=0.8,
                                arrowstyle='<->', linewidth=1.5)
    vol_arrow2 = FancyArrowPatch((50+box_npt_size/2, 50), (50+box_npt_size/2+5, 50), 
                                mutation_scale=10, color='blue', alpha=0.8,
                                arrowstyle='<->', linewidth=1.5)
    vol_arrow3 = FancyArrowPatch((50, 50-box_npt_size/2-5), (50, 50-box_npt_size/2), 
                                mutation_scale=10, color='blue', alpha=0.8,
                                arrowstyle='<->', linewidth=1.5)
    vol_arrow4 = FancyArrowPatch((50, 50+box_npt_size/2), (50, 50+box_npt_size/2+5), 
                                mutation_scale=10, color='blue', alpha=0.8,
                                arrowstyle='<->', linewidth=1.5)
    ax3.add_patch(vol_arrow1)
    ax3.add_patch(vol_arrow2)
    ax3.add_patch(vol_arrow3)
    ax3.add_patch(vol_arrow4)
    
    # Scale particle positions as volume changes
    scale_factor = box_npt_size / 70
    x_npt_new = (x_npt - 50) * scale_factor + 50 + np.random.uniform(-0.5, 0.5, n_particles)
    y_npt_new = (y_npt - 50) * scale_factor + 50 + np.random.uniform(-0.5, 0.5, n_particles)
    
    # Keep particles within boundaries
    x_npt_new = np.clip(x_npt_new, 50-box_npt_size/2+5, 50+box_npt_size/2-5)
    y_npt_new = np.clip(y_npt_new, 50-box_npt_size/2+5, 50+box_npt_size/2-5)
    
    scatter_npt.set_offsets(np.c_[x_npt_new, y_npt_new])
    
    return scatter_nve, quiver_nve, scatter_can, box_npt, scatter_npt

# Create animation
ani = FuncAnimation(fig, update, frames=100, interval=50, blit=False)

plt.tight_layout()

# Save as GIF with higher quality
ani.save('ensemble_visualization_3panels.gif', writer='pillow', fps=20, dpi=150)

plt.show() 