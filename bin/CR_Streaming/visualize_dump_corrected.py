#!/usr/bin/env python3
"""
Visualize CR_Streaming dump data from .npy patch files.
Reconstructs full grid with ghost zones and shows patch boundaries.

CORRECTED VERSION: Stacks patches along Z-axis (simulation extends in z-direction)
Grid: 16×16×128 without ghost zones -> 22×22×176 with ghost zones (8 patches along z)
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from collections import defaultdict

# Configuration
DUMP_DIR = "dump/npy/"
OUTPUT_DIR = "fig/"
PATCH_SIZE = 22  # Including ghost zones (16 + 3 + 3)
N_PATCHES = 8    # 8 patches along z-axis


def discover_files(dump_dir):
    """Discover all .npy files and organize by array_type, field, step, patch"""
    pattern = re.compile(r'^(.+)_([^_]+)_step(\d+)_patch(\d+)\.npy$')
    
    files_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    
    if not os.path.isdir(dump_dir):
        raise NotADirectoryError(f"Dump directory not found: {dump_dir}")
    
    for filename in os.listdir(dump_dir):
        if not filename.endswith('.npy'):
            continue
        
        match = pattern.match(filename)
        if not match:
            continue
        
        array_type = match.group(1)
        field = match.group(2)
        step = int(match.group(3))
        patch = int(match.group(4))
        
        files_dict[array_type][field][step][patch] = filename
    
    return files_dict


def load_patches(dump_dir, array_type, field, step, n_patches=8):
    """Load all patches for a given array_type, field, and step"""
    files_dict = discover_files(dump_dir)
    
    if array_type not in files_dict:
        raise ValueError(f"Array type '{array_type}' not found")
    if field not in files_dict[array_type]:
        raise ValueError(f"Field '{field}' not found for array type '{array_type}'")
    if step not in files_dict[array_type][field]:
        raise ValueError(f"Step {step} not found for {array_type}/{field}")
    
    patches = []
    for patch_idx in range(n_patches):
        if patch_idx not in files_dict[array_type][field][step]:
            raise ValueError(f"Patch {patch_idx} missing for {array_type}/{field}/step{step}")
        
        filename = files_dict[array_type][field][step][patch_idx]
        filepath = os.path.join(dump_dir, filename)
        data = np.load(filepath)
        patches.append(data)
    
    return patches


def reconstruct_grid(patches):
    """
    Reconstruct full grid from 8 patches along z-axis.
    Each patch: 22×22×22
    Output: 22×22×176 (8*22 along z-axis)
    """
    n_patches = len(patches)
    patch_size = patches[0].shape[0]  # Should be 22
    
    # Full grid dimensions
    nx = patch_size  # 22
    ny = patch_size  # 22
    nz = n_patches * patch_size  # 176
    
    full_grid = np.zeros((nx, ny, nz))
    
    # Place patches along z-axis
    for patch_idx in range(n_patches):
        z_start = patch_idx * patch_size
        z_end = z_start + patch_size
        full_grid[:, :, z_start:z_end] = patches[patch_idx]
    
    return full_grid


def compute_statistics(data):
    """Compute statistics for a 3D array"""
    stats = {
        'min': float(np.min(data)),
        'max': float(np.max(data)),
        'avg': float(np.mean(data)),
        'n_zero': int(np.sum(data == 0.0)),
        'total': int(data.size)
    }
    stats['pct_zero'] = 100.0 * stats['n_zero'] / stats['total']
    return stats


def visualize_slices_xy(full_grid, array_type, field, step, output_dir):
    """
    Visualize x-y slices at different z positions.
    full_grid: shape (22, 22, 176)
    Shows patch boundaries along z-axis.
    """
    nx, ny, nz = full_grid.shape
    
    # Calculate grid layout for subplots
    # We have 176 slices from 8 patches, each patch contributes 22 slices
    # Arrange as 22 rows × 8 columns (each column = one patch)
    n_cols = 8  # 8 patches
    n_rows = 22  # 22 slices per patch
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*2.5, n_rows*2.5))
    fig.suptitle(f'{array_type} - {field} - Step {step}\nFull Grid: {nx}×{ny}×{nz} (8 patches × 22³) - XY slices', 
                 fontsize=18, fontweight='bold')
    
    # Find global min/max for consistent colorbar
    vmin, vmax = full_grid.min(), full_grid.max()
    
    # Plot each x-y slice at different z positions, organized by patch
    for z_idx in range(nz):
        patch_idx = z_idx // PATCH_SIZE
        slice_in_patch = z_idx % PATCH_SIZE
        
        col = patch_idx
        row = slice_in_patch
        ax = axes[row, col]
        
        # Extract x-y slice at this z position
        slice_data = full_grid[:, :, z_idx]  # Shape: (22, 22)
        
        # Plot heatmap
        im = ax.imshow(slice_data, cmap='viridis', aspect='equal',
                      vmin=vmin, vmax=vmax, origin='lower',
                      interpolation='nearest')
        
        # Add title for first row (patch numbers)
        if row == 0:
            ax.set_title(f'Patch {patch_idx}', fontsize=10, fontweight='bold')
        
        # Add label for first column (slice within patch)
        if col == 0:
            ax.set_ylabel(f'z={z_idx}', fontsize=7)
        
        # Remove ticks
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Highlight patch boundaries with colored border
        if slice_in_patch == 0:
            # First slice in patch - add thick border
            for spine in ax.spines.values():
                spine.set_edgecolor('red')
                spine.set_linewidth(3)
        elif slice_in_patch == PATCH_SIZE - 1:
            # Last slice in patch - add thick border
            for spine in ax.spines.values():
                spine.set_edgecolor('orange')
                spine.set_linewidth(3)
        else:
            # Interior slices
            for spine in ax.spines.values():
                spine.set_edgecolor('gray')
                spine.set_linewidth(0.5)
    
    # Add colorbar
    fig.colorbar(im, ax=axes, orientation='horizontal', 
                 pad=0.01, fraction=0.03, label=field)
    
    # Add legend for boundary colors
    legend_elements = [
        Rectangle((0,0), 1, 1, fc='none', ec='red', linewidth=3, label='Patch start (z mod 22 = 0)'),
        Rectangle((0,0), 1, 1, fc='none', ec='orange', linewidth=3, label='Patch end (z mod 22 = 21)')
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=2, 
               bbox_to_anchor=(0.5, -0.01), fontsize=10)
    
    plt.tight_layout()
    
    # Save figure
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{array_type}_{field}_step{step}_xy_slices.png')
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    plt.close()
    
    return output_file


def visualize_z_profile(full_grid, array_type, field, step, output_dir):
    """
    Visualize z-profile: average value at each z-slice.
    full_grid: shape (22, 22, 176)
    """
    nx, ny, nz = full_grid.shape
    
    # Compute average at each z-slice
    z_profile = np.mean(full_grid, axis=(0, 1))
    z_coords = np.arange(nz)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot profile
    ax.plot(z_coords, z_profile, 'b-', linewidth=2, label='Average value')
    
    # Mark patch boundaries
    for patch_idx in range(N_PATCHES):
        z_boundary = patch_idx * PATCH_SIZE
        if patch_idx > 0:
            ax.axvline(z_boundary, color='red', linestyle='--', alpha=0.5, linewidth=2)
    
    # Add patch labels
    for patch_idx in range(N_PATCHES):
        z_center = (patch_idx + 0.5) * PATCH_SIZE
        ax.text(z_center, ax.get_ylim()[1] * 0.95, f'P{patch_idx}', 
                ha='center', va='top', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.set_xlabel('Z index', fontsize=12)
    ax.set_ylabel(f'{field} (averaged over x-y)', fontsize=12)
    ax.set_title(f'{array_type} - {field} - Step {step}\nZ-profile (Grid: {nx}×{ny}×{nz})', 
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    
    # Save figure
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{array_type}_{field}_step{step}_z_profile.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    return output_file


def visualize_xz_slice(full_grid, array_type, field, step, output_dir, y_slice=11):
    """
    Visualize x-z slice at a fixed y position (default: center at y=11).
    full_grid: shape (22, 22, 176)
    """
    nx, ny, nz = full_grid.shape
    
    # Extract x-z slice at y position
    slice_data = full_grid[:, y_slice, :]  # Shape: (22, 176)
    
    fig, ax = plt.subplots(figsize=(16, 4))
    
    # Plot heatmap
    im = ax.imshow(slice_data.T, cmap='viridis', aspect='auto',
                   origin='lower', interpolation='nearest',
                   extent=[0, nx, 0, nz])
    
    # Add patch boundaries
    for patch_idx in range(1, N_PATCHES):
        z_boundary = patch_idx * PATCH_SIZE
        ax.axhline(z_boundary, color='red', linestyle='--', linewidth=2, alpha=0.7)
    
    # Add patch labels
    for patch_idx in range(N_PATCHES):
        z_center = (patch_idx + 0.5) * PATCH_SIZE
        ax.text(nx * 0.95, z_center, f'Patch {patch_idx}', 
                ha='right', va='center', fontsize=10, fontweight='bold',
                color='white', bbox=dict(boxstyle='round', facecolor='black', alpha=0.5))
    
    ax.set_xlabel('X index', fontsize=12)
    ax.set_ylabel('Z index', fontsize=12)
    ax.set_title(f'{array_type} - {field} - Step {step}\nX-Z slice at y={y_slice} (Grid: {nx}×{ny}×{nz})', 
                 fontsize=14, fontweight='bold')
    
    # Add colorbar
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.1, fraction=0.05)
    cbar.set_label(field, fontsize=12)
    
    plt.tight_layout()
    
    # Save figure
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{array_type}_{field}_step{step}_xz_slice_y{y_slice}.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    return output_file


def process_all(dump_dir=DUMP_DIR, output_dir=OUTPUT_DIR, visualize_all_slices=False):
    """Process all discovered files and generate visualizations + statistics"""
    
    print("=" * 80)
    print("CR_Streaming Data Visualization (CORRECTED - Z-axis stacking)")
    print("=" * 80)
    print(f"Expected grid: 22×22×176 (8 patches of 22³ along z-axis)")
    print(f"Simulation: 16×16×128 without ghost zones")
    print("=" * 80)
    
    # Discover all files
    files_dict = discover_files(dump_dir)
    
    if not files_dict:
        print(f"No .npy files found in {dump_dir}")
        return
    
    # Prepare statistics output
    stats_file = os.path.join(output_dir, 'statistics_corrected.txt')
    os.makedirs(output_dir, exist_ok=True)
    
    with open(stats_file, 'w') as f:
        f.write("CR_Streaming Data Statistics (Corrected Z-axis Reconstruction)\n")
        f.write("=" * 80 + "\n")
        f.write("Grid: 22×22×176 (8 patches of 22³ along z-axis)\n")
        f.write("Simulation: 16×16×128 without ghost zones\n")
        f.write("=" * 80 + "\n\n")
    
    # Process each combination
    for array_type in sorted(files_dict.keys()):
        for field in sorted(files_dict[array_type].keys()):
            for step in sorted(files_dict[array_type][field].keys()):
                print(f"\nProcessing: {array_type} - {field} - Step {step}")
                
                try:
                    # Load and reconstruct
                    patches = load_patches(dump_dir, array_type, field, step, N_PATCHES)
                    full_grid = reconstruct_grid(patches)
                    
                    print(f"  Grid shape: {full_grid.shape}")
                    
                    # Compute statistics
                    stats = compute_statistics(full_grid)
                    
                    # Print to stdout
                    print(f"  Min: {stats['min']:.10e}")
                    print(f"  Max: {stats['max']:.10e}")
                    print(f"  Avg: {stats['avg']:.10e}")
                    print(f"  Number of zero values: {stats['n_zero']} / {stats['total']} ({stats['pct_zero']:.2f}%)")
                    
                    # Write to file
                    with open(stats_file, 'a') as f:
                        f.write(f"ArrayType: {array_type}, Field: {field}, Step: {step}\n")
                        f.write(f"  Grid shape: {full_grid.shape}\n")
                        f.write(f"  Min: {stats['min']:.10e}\n")
                        f.write(f"  Max: {stats['max']:.10e}\n")
                        f.write(f"  Avg: {stats['avg']:.10e}\n")
                        f.write(f"  Number of zero values: {stats['n_zero']} / {stats['total']} ({stats['pct_zero']:.2f}%)\n")
                        f.write("\n")
                    
                    # Visualize z-profile (always)
                    output_file = visualize_z_profile(full_grid, array_type, field, step, output_dir)
                    print(f"  Saved z-profile: {output_file}")
                    
                    # Visualize x-z slice (always)
                    output_file = visualize_xz_slice(full_grid, array_type, field, step, output_dir)
                    print(f"  Saved x-z slice: {output_file}")
                    
                    # Visualize all x-y slices (optional - can be slow)
                    if visualize_all_slices:
                        output_file = visualize_slices_xy(full_grid, array_type, field, step, output_dir)
                        print(f"  Saved x-y slices: {output_file}")
                    
                except Exception as e:
                    print(f"  ERROR: {e}")
                    import traceback
                    traceback.print_exc()
                    continue
    
    print("\n" + "=" * 80)
    print(f"Processing complete!")
    print(f"Figures saved to: {output_dir}")
    print(f"Statistics saved to: {stats_file}")
    print("=" * 80)


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Visualize CR_Streaming dump data (corrected z-axis stacking)')
    parser.add_argument('--all-slices', action='store_true', 
                        help='Generate all x-y slice visualizations (slow, creates large files)')
    args = parser.parse_args()
    
    process_all(visualize_all_slices=args.all_slices)


if __name__ == "__main__":
    main()
