#!/usr/bin/env python3
"""
Visualize g_Flu_Array_In dump as raw double arrays.
Default expected shape: 14 x 22 x 22 x 22 (field, z, y, x).
"""

import argparse
import math
import os

import numpy as np
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description="Plot slices from g_Flu_Array_In[P][8] dump.")
    parser.add_argument(
        "-i",
        "--input",
        default="g_Flu_Array_In.bin",
        help="Path to raw binary file (float64).",
    )
    parser.add_argument(
        "-n",
        "--size",
        type=int,
        default=22,
        help="Cube edge length (default: 22).",
    )
    parser.add_argument(
        "-f",
        "--nfields",
        type=int,
        default=14,
        help="Number of fields in the dump (default: 14).",
    )
    parser.add_argument(
        "--field",
        type=int,
        default=None,
        help="Only plot a single field index (0-based).",
    )
    parser.add_argument(
        "-o",
        "--output-base",
        default="g_Flu_Array_In_slices",
        help="Output image base name (field index is appended).",
    )
    parser.add_argument(
        "--cols",
        type=int,
        default=6,
        help="Number of columns for slice grid.",
    )
    return parser.parse_args()


def load_fields(path, n, nfields):
    expected = nfields * n * n * n
    data = np.fromfile(path, dtype=np.float64)
    if data.size != expected:
        print(
            f"Warning: element count {data.size} != expected {expected} for {nfields}*{n}^3."
        )
    if data.size < expected:
        raise ValueError("Not enough data to reshape to fields.")
    data = data[:expected]
    fields = data.reshape((nfields, n, n, n))
    return fields


def plot_slices(cube, output_file, cols, field_idx=None):
    nz, ny, nx = cube.shape
    rows = math.ceil(nz / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 2.5, rows * 2.5))
    title = f"g_Flu_Array_In Field {field_idx} Slices ({nz}x{ny}x{nx})"
    fig.suptitle(title, fontsize=14)

    vmin = cube.min()
    vmax = cube.max()

    axes = np.atleast_2d(axes)
    for z in range(nz):
        r = z // cols
        c = z % cols
        ax = axes[r, c]
        im = ax.imshow(
            cube[z, :, :],
            cmap="viridis",
            origin="lower",
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_title(f"z={z}", fontsize=7)
        ax.set_xticks([])
        ax.set_yticks([])

    for z in range(nz, rows * cols):
        r = z // cols
        c = z % cols
        axes[r, c].axis("off")

    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("Value", rotation=270, labelpad=12)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    print(f"Saved figure to {output_file}")

    print("\nStatistics:")
    print(f"  Shape: {cube.shape}")
    print(f"  Min: {vmin:.6e}")
    print(f"  Max: {vmax:.6e}")
    print(f"  Mean: {cube.mean():.6e}")
    print(f"  Std: {cube.std():.6e}")


def main():
    args = parse_args()
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")

    fields = load_fields(args.input, args.size, args.nfields)

    if args.field is not None:
        if args.field < 0 or args.field >= args.nfields:
            raise ValueError("Field index out of range.")
        output = f"{args.output_base}_f{args.field:02d}.png"
        plot_slices(fields[args.field], output, args.cols, field_idx=args.field)
        return

    for field_idx in range(args.nfields):
        output = f"{args.output_base}_f{field_idx:02d}.png"
        plot_slices(fields[field_idx], output, args.cols, field_idx=field_idx)


if __name__ == "__main__":
    main()
