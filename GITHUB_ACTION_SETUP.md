# GitHub Actions GAMER Setup

This repository includes configurations for compiling and running GAMER on GitHub Actions runners.

## Files Created

1. **configs/github_action.config** - Machine configuration file for GitHub Actions runners
   - CXX: `g++`
   - CXX_MPI: `mpicxx`
   - MPI_PATH: `/usr` (where OpenMPI installs on Ubuntu)
   - Compiler flags: `-g -O2`

2. **.github/workflows/install_gamer_env.yaml** - Reusable workflow for setting up GAMER environment
   - Installs OpenMPI (openmpi-bin and libopenmpi-dev)
   - Configures GAMER to use `github_action` machine settings
   - Verifies installation

## Using the Reusable Workflow

To use the GAMER environment setup in your own workflow:

```yaml
name: My GAMER Workflow

on:
  push:
    branches: [ main ]

jobs:
  setup-environment:
    uses: ./.github/workflows/install_gamer_env.yaml
  
  compile-and-test:
    needs: setup-environment
    runs-on: ubuntu-latest
    steps:
      - name: Checkout Repository
        uses: actions/checkout@v4
      
      - name: Install OpenMPI
        run: |
          sudo apt-get update -qq
          sudo apt-get install -y openmpi-bin libopenmpi-dev
      
      - name: Configure GAMER
        run: |
          cd tool/config
          bash set_settings.sh --local --machine=github_action
      
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

1. Install OpenMPI:
   ```bash
   sudo apt-get update -qq
   sudo apt-get install -y openmpi-bin libopenmpi-dev
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

## Verified Test

The setup has been verified with the 1D Shock Tube test problem:
- Successfully compiles GAMER in CPU-only mode (no GPU, MPI, or OpenMP)
- Runs the Riemann test problem
- Completes successfully showing "~ GAME OVER ~"
- Records processing time in `Record__Note`

## Notes

- The current setup is for CPU-only compilation without GPU support
- MPI is installed but not enabled in the test configuration (use `--mpi=true` to enable)
- OpenMP can be enabled with `--openmp=true`
- GPU support requires CUDA installation and proper GPU hardware
