#!/usr/bin/env python3
"""
Script to output statistics for CR_F1, CR_F2, and CR_F3 from all Data_* files
"""

import h5py
import numpy as np
import glob
import os

def reconstruct_grid(patches, corners):
    """
    Reconstruct full grid from patches using corner coordinates
    patches: shape (n_patches, 8, 8, 8)
    corners: shape (n_patches, 3) - (x, y, z) corner coordinates
    Returns: 3D array with the full grid
    """
    n_patches = patches.shape[0]
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
        
        # Place patch in full grid
        full_grid[z_start:z_start+patch_size, 
                  y_start:y_start+patch_size, 
                  x_start:x_start+patch_size] = patches[patch_idx]
    
    return full_grid

def print_stats(data, field_name):
    """Print statistics in the requested format"""
    print(f"{field_name} Statistics:")
    print(f"  Shape: {data.shape}")
    print(f"  Min: {data.min():.6e}")
    print(f"  Max: {data.max():.6e}")
    print(f"  Mean: {data.mean():.6e}")
    print(f"  Std: {data.std():.6e}")
    print()

def main():
    # Find all Data_* files
    data_files = sorted(glob.glob('Data_*'))
    
    if not data_files:
        print("No Data_* files found in current directory")
        return
    
    print(f"Found {len(data_files)} data files\n")
    print("="*60)
    
    # Process each file
    for data_file in data_files:
        print(f"\nFile: {data_file}")
        print("-"*60)
        
        with h5py.File(data_file, 'r') as f:
            # Load corner coordinates
            corners = f['Tree']['Corner'][:]
            
            # Load and process CR_F1, CR_F2, CR_F3
            for field in ['CR_F1', 'CR_F2', 'CR_F3']:
                if field in f['GridData']:
                    patches = f['GridData'][field][:]
                    full_grid = reconstruct_grid(patches, corners)
                    print_stats(full_grid, field)
                else:
                    print(f"{field}: Not found in file")
                    print()

if __name__ == '__main__':
    main()
