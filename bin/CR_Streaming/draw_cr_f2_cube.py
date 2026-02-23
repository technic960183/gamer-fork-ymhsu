#!/usr/bin/env python3
"""
Script to visualize CR_F2 values from Data_000002 as 2D heatmaps
Shows 16 slices of 16x128 data arranged in a grid
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def load_cr_f2_data(filename):
    """Load CR_F2 data from HDF5 file"""
    with h5py.File(filename, 'r') as f:
        # Load CR_F2 data: shape (64, 8, 8, 8) - (patches, z, y, x)
        cr_f2 = f['GridData']['CR_E'][:]
        
        # Load patch corner positions
        corners = f['Tree']['Corner'][:]
        
        return cr_f2, corners

def reconstruct_grid(cr_f2_patches, corners):
    """
    Reconstruct the full grid from patches using corner coordinates
    cr_f2_patches: shape (n_patches, 8, 8, 8) - patch data
    corners: shape (n_patches, 3) - (x, y, z) corner coordinates in units
    Returns: 3D array with the full grid
    """
    n_patches = cr_f2_patches.shape[0]
    patch_size = 8
    
    # Determine patch spacing from corners
    unique_x = np.unique(corners[:, 0])
    unique_y = np.unique(corners[:, 1])
    unique_z = np.unique(corners[:, 2])
    
    patch_spacing = unique_x[1] - unique_x[0] if len(unique_x) > 1 else 4096
    
    # Calculate full grid dimensions
    nx = len(unique_x) * patch_size
    ny = len(unique_y) * patch_size
    nz = len(unique_z) * patch_size
    
    full_grid = np.zeros((nz, ny, nx))
    
    # Place each patch according to its corner coordinates
    for patch_idx in range(n_patches):
        corner = corners[patch_idx]
        
        # Convert corner coordinates to patch indices
        ix = corner[0] // patch_spacing
        iy = corner[1] // patch_spacing
        iz = corner[2] // patch_spacing
        
        # Convert patch indices to cell indices
        x_start = ix * patch_size
        y_start = iy * patch_size
        z_start = iz * patch_size
        
        # Get patch data (shape: 8, 8, 8)
        patch_data = cr_f2_patches[patch_idx]
        
        # Place patch in full grid
        full_grid[z_start:z_start+patch_size, 
                  y_start:y_start+patch_size, 
                  x_start:x_start+patch_size] = patch_data
    
    return full_grid

def plot_slices(data, output_file='cr_f2_slices.png'):
    """
    Plot ALL slices of the 3D data as 2D heatmaps
    data: shape (nz, ny, nx) - expected (128, 16, 16)
    """
    nz, ny, nx = data.shape
    
    # Show ALL slices
    n_slices = nz  # All 128 slices
    slice_indices = np.arange(nz)
    
    # Create figure with 16x8 grid of subplots (128 total)
    n_cols = 16
    n_rows = n_slices // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*3, n_rows*3))
    fig.suptitle(f'CR_F2 Values - ALL {n_slices} Slices (Grid: {nz}×{ny}×{nx})', fontsize=20, fontweight='bold')
    
    # Find global min/max for consistent colorbar
    vmin, vmax = data.min(), data.max()
    
    # Plot each slice
    for idx, z_idx in enumerate(slice_indices):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]
        
        # Extract slice at this z position
        slice_data = data[z_idx, :, :]  # Shape: (ny, nx) = (16, 16)
        
        # Plot heatmap
        im = ax.imshow(slice_data, cmap='viridis', aspect='auto', 
                      origin='lower', vmin=vmin, vmax=vmax)
        
        # Labels and title
        ax.set_title(f'z={z_idx}', fontsize=6)
        ax.set_xticks([])
        ax.set_yticks([])
        
    # Add a single colorbar for the entire figure
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('CR_F2', rotation=270, labelpad=20, fontsize=14)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Figure saved to {output_file}")
    
    # Also print statistics
    print(f"\nCR_F2 Statistics:")
    print(f"  Shape: {data.shape}")
    print(f"  Min: {vmin:.6e}")
    print(f"  Max: {vmax:.6e}")
    print(f"  Mean: {data.mean():.6e}")
    print(f"  Std: {data.std():.6e}")

def main():
    # Load data
    print("Loading CR_F2 data from Data_000000...")
    cr_f2_patches, corners = load_cr_f2_data('Data_000003')
    print(f"Loaded {cr_f2_patches.shape[0]} patches of shape {cr_f2_patches.shape[1:]}")
    print(f"Corner coordinates shape: {corners.shape}")
    
    # Reconstruct full grid
    print("Reconstructing full grid using corner coordinates...")
    full_grid = reconstruct_grid(cr_f2_patches, corners)
    print(f"Full grid shape: {full_grid.shape}")
    
    # Plot slices
    print("Creating visualization...")
    plot_slices(full_grid)
    
    print("\nDone!")

if __name__ == '__main__':
    main()
