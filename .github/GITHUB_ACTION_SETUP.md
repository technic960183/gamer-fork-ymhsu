# GitHub Actions GAMER Setup

This repository includes configurations for compiling and running GAMER on GitHub Actions runners.

## Files Created

1. **configs/github_action.config** - Machine configuration file for GitHub Actions runners
   - CXX: `g++`
   - CXX_MPI: `mpicxx`
   - MPI_PATH: `/usr` (where OpenMPI installs on Ubuntu)
   - CUDA_PATH: `/usr` (CUDA toolkit installation path)
   - HDF5_PATH: `/usr/lib/x86_64-linux-gnu/hdf5/serial` (HDF5 library path)
   - Compiler flags: `-g -O2`

2. **.github/install_gamer_env.sh** - Setup script for GAMER environment
   - Installs OpenMPI (openmpi-bin and libopenmpi-dev)
   - Installs HDF5 (libhdf5-dev)
   - Installs CUDA toolkit (nvidia-cuda-toolkit)
   - Configures GAMER to use `github_action` machine settings
   - Verifies installation

## Using the Setup Script

To use the GAMER environment setup script in your workflow:

```yaml
name: My GAMER Workflow

on:
  push:
    branches: [ main ]

jobs:
  compile-and-test:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout Repository
        uses: actions/checkout@v4
      
      - name: Setup GAMER Environment
        run: bash .github/install_gamer_env.sh
      
      - name: Generate Makefile
        run: |
          cd src
          cp ../example/test_problem/Hydro/Riemann/generate_make.sh ./
          bash generate_make.sh --openmp=false --mpi=false --gpu=false
      
      - name: Compile GAMER
        run: |
          cd src
          make clean
          make -j4
      
      - name: Run Test
        run: |
          mkdir -p bin/test
          cd bin/test
          cp -r ../../example/test_problem/Hydro/Riemann/* .
          cp ../../src/gamer .
          ./gamer
          tail -n 3 Record__Note
```

## Manual Setup

To manually set up the environment:

1. Install dependencies (OpenMPI, HDF5, CUDA):
   ```bash
   sudo apt-get update -qq
   sudo apt-get install -y openmpi-bin libopenmpi-dev libhdf5-dev nvidia-cuda-toolkit
   ```

2. Configure GAMER machine settings:
   ```bash
   cd tool/config
   bash set_settings.sh --local --machine=github_action
   ```

3. Generate Makefile and compile:
   ```bash
   cd src
   cp ../example/test_problem/Hydro/Riemann/generate_make.sh ./
   bash generate_make.sh --openmp=false --mpi=false --gpu=false
   make clean
   make -j4
   ```

4. Run test:
   ```bash
   mkdir -p bin/shocktube
   cd bin/shocktube
   cp -r ../../example/test_problem/Hydro/Riemann/* .
   cp ../../src/gamer .
   ./gamer
   tail -n 3 Record__Note
   ```

## Verified Tests

The setup has been verified with multiple configurations:

### 1. CPU-only without OpenMP
- Successfully compiles GAMER
- Runs the 1D Shock Tube (Riemann) test problem
- Completes successfully showing "~ GAME OVER ~"
- Records processing time in `Record__Note`

### 2. CPU-only with OpenMP
- Successfully compiles GAMER with OpenMP support
- Runs the 1D Shock Tube test problem
- Default `OMP_NTHREAD` (-1, auto) detects **4 threads** on GitHub Actions runners
- OpenMP diagnosis shows proper thread distribution across CPU cores

### 3. Hybrid OpenMP/GPU
- Successfully compiles GAMER with both OpenMP and GPU support
- GPU code compilation works with CUDA toolkit
- **Note**: Cannot run tests as GitHub Actions runners don't have physical GPUs

### 4. MPI + OpenMP + HDF5 (3D Blast Wave)
- Successfully compiles with MPI, OpenMP, and HDF5 support
- Runs with 2 MPI processes
- OpenMP shows **2 threads per MPI rank**
- Successfully creates HDF5 output files (`Data_000000`, etc.)
- HDF5 files verified with valid structure (GridData, Info, Tree groups)

## Notes

- OpenMPI 4.1.6 is installed from Ubuntu repositories
- CUDA toolkit 12.0 is available for GPU compilation
- HDF5 library is installed for data output
- GitHub Actions runners have 2-4 CPU cores available
- GPU code compiles but cannot run without physical GPU hardware
- All MPI, OpenMP, and HDF5 features are fully functional
