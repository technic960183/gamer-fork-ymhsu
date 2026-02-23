#!/usr/bin/env python3
"""
Visualize CR_Streaming dump data - ALL slices visualization
Dynamically handles different patch sizes (18³, 20³, 22³, etc.)

Uses data from dump/npy/ directory with multiple patches along z-axis
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# Configuration
DUMP_DIR = "dump/npy/"
OUTPUT_DIR = "fig/"
N_PATCHES = 8    # 8 patches along z-axis (typically)


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
    Reconstruct full grid from patches along z-axis.
    Patches can have any size (e.g., 18³, 20³, 22³).
    
    Args:
        patches: List of 3D arrays, each with shape (edge, edge, edge)
    
    Returns:
        3D array with shape (nx, ny, nz) where nz = n_patches * edge
    """
    n_patches = len(patches)
    
    # Get dimensions from first patch
    patch_shape = patches[0].shape
    nx, ny, patch_size_z = patch_shape[0], patch_shape[1], patch_shape[2]
    
    # Verify all patches have the same shape
    for i, patch in enumerate(patches):
        if patch.shape != patch_shape:
            raise ValueError(f"Patch {i} has shape {patch.shape}, expected {patch_shape}")
    
    # Full grid dimensions
    nz = n_patches * patch_size_z
    
    full_grid = np.zeros((nx, ny, nz))
    
    # Place patches along z-axis
    for patch_idx in range(n_patches):
        z_start = patch_idx * patch_size_z
        z_end = z_start + patch_size_z
        full_grid[:, :, z_start:z_end] = patches[patch_idx]
    
    return full_grid


def plot_all_slices(data, array_type, field, step, output_file):
    """
    Plot ALL z-slices of the 3D data as 2D heatmaps
    data: shape (nx, ny, nz) - can be any size
    Shows all nz slices arranged in grid: (slices_per_patch) rows × (n_patches) columns
    """
    nx, ny, nz = data.shape
    
    # Infer number of patches and slices per patch
    # Try common patch counts: 8, 4, 16, etc.
    for n_patches_guess in [8, 4, 16, 2, 1, 32]:
        if nz % n_patches_guess == 0:
            n_patches = n_patches_guess
            slices_per_patch = nz // n_patches
            break
    else:
        # Fallback: treat as single patch
        n_patches = 1
        slices_per_patch = nz
    
    # Create figure with slices_per_patch rows × n_patches columns
    n_cols = n_patches
    n_rows = slices_per_patch
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*2, n_rows*2))
    
    # Handle single patch case (axes won't be 2D)
    if n_patches == 1 and n_rows == 1:
        axes = np.array([[axes]])
    elif n_patches == 1:
        axes = axes.reshape(-1, 1)
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    
    # Title with underscores replaced for better readability
    title_array_type = array_type.replace('_', ' ')
    fig.suptitle(f'{title_array_type} - {field} - Step {step}\n'
                 f'All {nz} Slices (Grid: {nx}×{ny}×{nz}, {n_patches} patches of {slices_per_patch}³)', 
                 fontsize=16, fontweight='bold')
    
    # Find global min/max for consistent colorbar
    vmin, vmax = data.min(), data.max()
    
    # Plot each x-y slice at different z positions
    for z_idx in range(nz):
        patch_idx = z_idx // slices_per_patch
        slice_in_patch = z_idx % slices_per_patch
        
        col = patch_idx
        row = slice_in_patch
        ax = axes[row, col]
        
        # Extract x-y slice at this z position
        slice_data = data[:, :, z_idx]  # Shape: (nx, ny)
        
        # Plot heatmap
        im = ax.imshow(slice_data, cmap='viridis', aspect='equal', 
                      origin='lower', vmin=vmin, vmax=vmax,
                      interpolation='nearest')
        
        # Add title for first row (patch numbers)
        if row == 0:
            ax.set_title(f'P{patch_idx}', fontsize=8, fontweight='bold')
        
        # Add z-index label for first column
        if col == 0:
            ax.set_ylabel(f'z={z_idx}', fontsize=6)
        
        # Remove ticks
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Highlight patch boundaries
        if slice_in_patch == 0:
            # First slice in patch
            for spine in ax.spines.values():
                spine.set_edgecolor('red')
                spine.set_linewidth(2)
        elif slice_in_patch == slices_per_patch - 1:
            # Last slice in patch
            for spine in ax.spines.values():
                spine.set_edgecolor('orange')
                spine.set_linewidth(2)
        else:
            # Interior slices
            for spine in ax.spines.values():
                spine.set_edgecolor('gray')
                spine.set_linewidth(0.3)
    
    # Add a single colorbar for the entire figure
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(field, rotation=270, labelpad=20, fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")
    print(f"  Layout: {n_rows} rows × {n_cols} cols ({slices_per_patch} slices per patch)")


def print_stats(data, array_type, field, step):
    """Print statistics for the data"""
    print(f"\n{array_type} - {field} - Step {step} Statistics:")
    print(f"  Shape: {data.shape}")
    print(f"  Min: {data.min():.6e}")
    print(f"  Max: {data.max():.6e}")
    print(f"  Mean: {data.mean():.6e}")
    print(f"  Std: {data.std():.6e}")


def process_all(dump_dir=DUMP_DIR, output_dir=OUTPUT_DIR, 
                specific_array_type=None, specific_field=None, specific_step=None):
    """
    Process all discovered files and generate visualizations
    
    Parameters:
    - specific_array_type: If provided, only process this array type
    - specific_field: If provided, only process this field
    - specific_step: If provided, only process this step
    """
    print("=" * 80)
    print("CR_Streaming Dump Data - ALL SLICES Visualization")
    print("=" * 80)
    print(f"Dynamically detects patch sizes from data")
    print("=" * 80)
    
    # Discover all files
    files_dict = discover_files(dump_dir)
    
    if not files_dict:
        print(f"No .npy files found in {dump_dir}")
        return
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Count total combinations
    total = 0
    for array_type in files_dict.keys():
        if specific_array_type and array_type != specific_array_type:
            continue
        for field in files_dict[array_type].keys():
            if specific_field and field != specific_field:
                continue
            for step in files_dict[array_type][field].keys():
                if specific_step is not None and step != specific_step:
                    continue
                total += 1
    
    print(f"Found {total} combinations to process\n")
    
    # Process each combination
    count = 0
    for array_type in sorted(files_dict.keys()):
        if specific_array_type and array_type != specific_array_type:
            continue
            
        for field in sorted(files_dict[array_type].keys()):
            if specific_field and field != specific_field:
                continue
                
            for step in sorted(files_dict[array_type][field].keys()):
                if specific_step is not None and step != specific_step:
                    continue
                    
                count += 1
                print(f"\n[{count}/{total}] Processing: {array_type} - {field} - Step {step}")
                
                try:
                    # Load and reconstruct
                    patches = load_patches(dump_dir, array_type, field, step, N_PATCHES)
                    full_grid = reconstruct_grid(patches)
                    
                    print(f"  Reconstructed grid shape: {full_grid.shape}")
                    
                    # Print statistics
                    print_stats(full_grid, array_type, field, step)
                    
                    # Generate output filename
                    output_file = os.path.join(output_dir, 
                                              f'{array_type}_{field}_step{step}_all_slices.png')
                    
                    # Visualize all slices
                    plot_all_slices(full_grid, array_type, field, step, output_file)
                    
                except Exception as e:
                    print(f"  ERROR: {e}")
                    import traceback
                    traceback.print_exc()
                    continue
    
    print("\n" + "=" * 80)
    print(f"Processing complete! {count} visualizations created.")
    print(f"Figures saved to: {output_dir}")
    print("=" * 80)


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Visualize CR_Streaming dump data - ALL slices',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all data
  python visualize_dump_all_slices.py
  
  # Process only CR fields
  python visualize_dump_all_slices.py --field E
  python visualize_dump_all_slices.py --field F1
  
  # Process specific array type
  python visualize_dump_all_slices.py --array-type FC_Flux_1PG_2_CR
  
  # Process specific step
  python visualize_dump_all_slices.py --step 0
  
  # Combine filters
  python visualize_dump_all_slices.py --array-type FC_Flux_1PG_2_CR --field E --step 0
        """)
    
    parser.add_argument('--array-type', type=str, default=None,
                       help='Filter by specific array type (e.g., FC_Flux_1PG_2_CR)')
    parser.add_argument('--field', type=str, default=None,
                       help='Filter by specific field (e.g., E, F1, F2, F3, SIGMA)')
    parser.add_argument('--step', type=int, default=None,
                       help='Filter by specific step number')
    parser.add_argument('--dump-dir', type=str, default=DUMP_DIR,
                       help=f'Dump directory (default: {DUMP_DIR})')
    parser.add_argument('--output-dir', type=str, default=OUTPUT_DIR,
                       help=f'Output directory (default: {OUTPUT_DIR})')
    
    args = parser.parse_args()
    
    process_all(dump_dir=args.dump_dir, 
                output_dir=args.output_dir,
                specific_array_type=args.array_type,
                specific_field=args.field,
                specific_step=args.step)


if __name__ == "__main__":
    main()
