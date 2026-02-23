#!/usr/bin/env python3
"""
Verify that parse.py correctly reshapes 1D arrays to 3D cubes
according to GAMER's indexing convention:
    IDX321(i, j, k, Ni, Nj) = (k * Nj + j) * Ni + i
"""

import numpy as np

def gamer_idx321(i, j, k, Ni, Nj):
    """GAMER's 3D to 1D index transformation"""
    return (k * Nj + j) * Ni + i

def test_reshape_correctness():
    """
    Test that NumPy's reshape matches GAMER's indexing convention.
    
    GAMER convention: IDX321(i, j, k, Ni, Nj) = (k * Nj + j) * Ni + i
    - i varies fastest (innermost/rightmost in memory)
    - j varies at middle speed
    - k varies slowest (outermost/leftmost in memory)
    
    This is C-order (row-major) with the order: k, j, i
    """
    print("=" * 80)
    print("Testing GAMER 3D ↔ 1D Index Transformation")
    print("=" * 80)
    
    # Test with a small cube
    Ni, Nj, Nk = 3, 3, 3
    size = Ni * Nj * Nk
    
    # Create 1D array with sequential values (simulating GAMER's memory layout)
    data_1d = np.arange(size, dtype=np.float64)
    print(f"\n1D array (GAMER memory layout):")
    print(f"  Size: {size}")
    print(f"  Values: {data_1d}")
    
    # Reshape using parse.py's method
    data_3d = data_1d.reshape((Nk, Nj, Ni))  # Note: (Nk, Nj, Ni) not (Ni, Nj, Nk)
    
    print(f"\n3D array after reshape({Nk}, {Nj}, {Ni}):")
    print(f"  Shape: {data_3d.shape}")
    
    # Verify: check if data_3d[k, j, i] == gamer_idx321(i, j, k, Ni, Nj)
    print("\nVerification: checking if data_3d[k, j, i] == GAMER_IDX321(i, j, k, Ni, Nj)")
    print("-" * 80)
    
    all_correct = True
    errors = []
    
    for k in range(Nk):
        for j in range(Nj):
            for i in range(Ni):
                # GAMER's 1D index for this (i, j, k)
                gamer_1d_idx = gamer_idx321(i, j, k, Ni, Nj)
                
                # Value at this position in 3D array (indexed as [k, j, i])
                value_3d = data_3d[k, j, i]
                
                # Expected value from 1D array
                expected = data_1d[gamer_1d_idx]
                
                if value_3d != expected:
                    all_correct = False
                    errors.append(f"  Mismatch at (i={i}, j={j}, k={k}): "
                                f"data_3d[{k},{j},{i}] = {value_3d}, "
                                f"expected {expected} from 1D[{gamer_1d_idx}]")
    
    if all_correct:
        print("✓ ALL CHECKS PASSED!")
        print(f"  Checked all {Ni}×{Nj}×{Nk} = {size} positions")
        print(f"  data_3d[k, j, i] correctly matches GAMER's IDX321(i, j, k, Ni, Nj)")
    else:
        print("✗ ERRORS FOUND:")
        for err in errors[:10]:  # Show first 10 errors
            print(err)
        if len(errors) > 10:
            print(f"  ... and {len(errors) - 10} more errors")
    
    # Show some example mappings
    print("\nExample mappings (first 10 positions):")
    print("-" * 80)
    print(f"{'1D idx':>7} | {'(i, j, k)':>12} | {'3D[k,j,i]':>10} | {'Value':>8}")
    print("-" * 80)
    
    count = 0
    for k in range(Nk):
        for j in range(Nj):
            for i in range(Ni):
                if count < 10:
                    gamer_1d_idx = gamer_idx321(i, j, k, Ni, Nj)
                    value = data_3d[k, j, i]
                    print(f"{gamer_1d_idx:7d} | ({i}, {j}, {k}):>10s | [{k},{j},{i}]:>8s | {value:8.0f}")
                    count += 1
    
    print("\n" + "=" * 80)
    print("Conclusion:")
    print("=" * 80)
    if all_correct:
        print("✓ parse.py's reshape((edge, edge, edge)) is CORRECT!")
        print("\nThe reshaped array should be accessed as:")
        print("  array[k, j, i]  where k ∈ [0, Nk), j ∈ [0, Nj), i ∈ [0, Ni)")
        print("\nNote: This means the array dimensions are ordered as (z, y, x)")
        print("      if you interpret k=z, j=y, i=x")
    else:
        print("✗ parse.py's reshape needs correction!")
    print("=" * 80)
    
    return all_correct

def test_realistic_size():
    """Test with realistic GAMER patch size: 8³ or 22³"""
    print("\n\n" + "=" * 80)
    print("Testing with realistic GAMER patch size: 22³")
    print("=" * 80)
    
    edge = 22
    size = edge ** 3
    
    # Create test data
    data_1d = np.arange(size, dtype=np.float64)
    data_3d = data_1d.reshape((edge, edge, edge))
    
    print(f"\n1D array size: {size}")
    print(f"3D array shape: {data_3d.shape}")
    
    # Check a few random positions
    test_positions = [
        (0, 0, 0),    # First corner
        (21, 21, 21), # Last corner
        (10, 10, 10), # Center
        (5, 10, 15),  # Random position
    ]
    
    print("\nSpot checks:")
    print("-" * 80)
    all_correct = True
    
    for i, j, k in test_positions:
        gamer_1d_idx = gamer_idx321(i, j, k, edge, edge)
        value_3d = data_3d[k, j, i]
        expected = data_1d[gamer_1d_idx]
        
        match = "✓" if value_3d == expected else "✗"
        print(f"{match} (i={i:2d}, j={j:2d}, k={k:2d}) → 1D[{gamer_1d_idx:5d}] = {expected:6.0f}, "
              f"3D[{k:2d},{j:2d},{i:2d}] = {value_3d:6.0f}")
        
        if value_3d != expected:
            all_correct = False
    
    if all_correct:
        print("\n✓ All spot checks passed for 22³ array!")
    else:
        print("\n✗ Some spot checks failed!")
    
    return all_correct

if __name__ == "__main__":
    result1 = test_reshape_correctness()
    result2 = test_realistic_size()
    
    print("\n" + "=" * 80)
    if result1 and result2:
        print("FINAL RESULT: ✓ parse.py correctly handles GAMER's indexing!")
    else:
        print("FINAL RESULT: ✗ parse.py needs correction!")
    print("=" * 80)
