import argparse
import os
import re
import sys

import numpy as np

LINE_RE = re.compile(
    r'^0x[0-9a-fA-F]+:\s+'
    r'([+-]?(?:\d+\.\d+|\d+)(?:[eE][+-]?\d+)?)\s+'
    r'([+-]?(?:\d+\.\d+|\d+)(?:[eE][+-]?\d+)?)'
)

# Stride configuration for different array types
# Maps filename patterns to their actual grid stride (not including ghost zones)
# If a pattern matches, use the specified stride; otherwise infer from data size
STRIDE_CONFIG = {
    r'PriVar_Half': 20,      # PriVar_Half arrays are 20×20×20
    r'FC_Flux_1PG_2': 18,    # FC_Flux_1PG_2 arrays are 18×18×18
}

def parse_gdb_log(filename):
    # Match only the two numeric values after each address line.
    data = []
    with open(filename, 'r') as f:
        for line in f:
            match = LINE_RE.match(line.strip())
            if not match:
                continue
            data.append(float(match.group(1)))
            data.append(float(match.group(2)))

    return np.array(data, dtype=np.float64)


def _infer_cube_edge(size):
    edge = int(round(size ** (1.0 / 3.0)))
    return edge if edge ** 3 == size else None


def _get_stride_from_filename(filename):
    """Get the expected stride for a given filename based on STRIDE_CONFIG"""
    basename = os.path.splitext(filename)[0]
    for pattern, stride in STRIDE_CONFIG.items():
        if re.search(pattern, basename):
            return stride
    return None


def reshape_to_cube(data, filename=None):
    warnings = []
    
    # Try to get stride from filename configuration
    expected_stride = None
    if filename:
        expected_stride = _get_stride_from_filename(filename)
    
    # If stride is configured, trim data to that size
    if expected_stride is not None:
        expected_size = expected_stride ** 3
        if data.size >= expected_size:
            if len(np.nonzero(data[expected_size:])[0]) > 0:
                warnings.append(f"{filename}: data has more than {expected_size} elements for stride {expected_stride}, "
                              f"but extra elements are not all zero; check if this is correct")
            data = data[:expected_size]
            edge = expected_stride
        else:
            warnings.append(f"{filename}: expected {expected_size} elements for stride {expected_stride}, "
                          f"but got {data.size}; attempting to infer")
            edge = _infer_cube_edge(data.size)
    else:
        # Infer edge from data size
        edge = _infer_cube_edge(data.size)
    
    if edge is None:
        warnings.append(f"cannot infer cube edge from size {data.size}; saving flat array")
        return data, warnings

    # GAMER's 1D layout: IDX321(i,j,k,Ni,Nj) = (k*Nj + j)*Ni + i
    # This means i (x) varies fastest, j (y) middle, k (z) slowest
    # To create array[x, y, z] (or array[i, j, k]), use Fortran order
    reshaped = data.reshape((edge, edge, edge), order='F')
    
    return reshaped, warnings


def parse_dump_dir(dump_base):
    warnings = []
    dump_dir = dump_base
    raw_dir = os.path.join(dump_base, "raw")
    if os.path.isdir(raw_dir):
        dump_dir = raw_dir

    if not os.path.isdir(dump_dir):
        raise NotADirectoryError(f"dump directory not found: {dump_dir}")

    files = [f for f in os.listdir(dump_dir) if f.endswith(".txt")]
    if not files:
        warnings.append(f"no .txt files found in {dump_dir}")
        return 0, warnings

    saved_count = 0
    for filename in sorted(files):
        file_path = os.path.join(dump_dir, filename)
        data = parse_gdb_log(file_path)

        cube, reshape_warnings = reshape_to_cube(data, filename=filename)
        warnings.extend(reshape_warnings)

        npy_dir = os.path.join(dump_base, "npy")
        if not os.path.exists(npy_dir):
            os.makedirs(npy_dir)
        npy_path = os.path.join(npy_dir, os.path.splitext(filename)[0] + ".npy")
        np.save(npy_path, cube)
        saved_count += 1

    return saved_count, warnings


def main():
    parser = argparse.ArgumentParser(
        description="Parse GDB dump files and save as NumPy cubes."
    )
    parser.add_argument("--dump-dir", default="dump", help="Dump directory to parse")
    args = parser.parse_args()

    try:
        saved_count, warnings = parse_dump_dir(args.dump_dir)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    for warning in warnings:
        print(f"Warning: {warning}", file=sys.stderr)

    print(f"Parsed {saved_count} files into .npy arrays.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())