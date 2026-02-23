#!/bin/bash

# Script to recompile specific CUDA files with -G debug flag and relink
# Usage: sh recompile_with_debug.sh

echo "Recompiling CUDA files with -G flag for debugging..."

# Change to src directory
cd "$(dirname "$0")"

# NVCC compiler path
NVCC="/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/cuda/12.2/bin/nvcc"

# Common flags (extracted from Makefile)
COMMON_FLAGS="-I../include -IModel_Hydro/GPU_Hydro -DMODEL=HYDRO -DNCOMP_PASSIVE_USER=0 -DFLU_SCHEME=MHM_RP -DLR_SCHEME=PLM -DRSOLVER=HLLE -DMHD -DCOSMIC_RAY -DEOS=EOS_COSMIC_RAY -DCR_DIFFUSION -DCR_STREAMING -DNLEVEL=10 -DMAX_PATCH=1000000 -DPATCH_SIZE=8 -DBITWISE_REPRODUCIBILITY -DTIMING -DFLOAT8 -DSUPPORT_HDF5 -DRANDOM_NUMBER=RNG_GNU_EXT -DGPU -DSERIAL -DGPU_COMPUTE_CAPABILITY=800"
NVCC_FLAGS="-O0 -g -I/opt/cray/pe/mpich/8.1.28/ofi/gnu/12.3/include -gencode arch=compute_80,code=\"compute_80,sm_80\""
FLU_FLAGS="-Xptxas -dlcm=ca -prec-div=false -ftz=true --maxrregcount=192"

# Add debug flag
DEBUG_FLAG="-G -Xcompiler -gdwarf-4"

# Recompile CUFLU_FluidSolver_MHM.cu
echo "Compiling CUFLU_FluidSolver_MHM.cu with -G flag..."
$NVCC $COMMON_FLAGS $NVCC_FLAGS $FLU_FLAGS $DEBUG_FLAG \
    -o Object/__gpu__CUFLU_FluidSolver_MHM.o \
    -dc Model_Hydro/GPU_Hydro/CUFLU_FluidSolver_MHM.cu

if [ $? -ne 0 ]; then
    echo "Error: Failed to compile CUFLU_FluidSolver_MHM.cu"
    exit 1
fi

# Recompile CUFLU_CR_TwoMoment.cu
echo "Compiling CUFLU_CR_TwoMoment.cu with -G flag..."
$NVCC $COMMON_FLAGS $NVCC_FLAGS $FLU_FLAGS $DEBUG_FLAG \
    -o Object/__gpu__CUFLU_CR_TwoMoment.o \
    -dc Model_Hydro/GPU_Hydro/CUFLU_CR_TwoMoment.cu

if [ $? -ne 0 ]; then
    echo "Error: Failed to compile CUFLU_CR_TwoMoment.cu"
    exit 1
fi

# Relink GPU objects
echo "Relinking GPU objects..."
$NVCC -o Object/gpu_link.o Object/__gpu__*.o \
    -gencode arch=compute_80,code=\"compute_80,sm_80\" -dlink

if [ $? -ne 0 ]; then
    echo "Error: Failed to link GPU objects"
    exit 1
fi

# Relink final executable
echo "Linking final executable..."
CC -target-accel=nvidia80 -g3 -O0 \
    -o gamer Object/__cpu__*.o Object/__gpu__*.o Object/gpu_link.o \
    -L/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/cuda/12.2/lib64 \
    -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/cuda/12.2/lib64 \
    -lcudart \
    -L/opt/cray/pe/hdf5/default/gnu/12.3/lib \
    -lhdf5 \
    -Wl,-rpath,/opt/cray/pe/hdf5/default/gnu/12.3/lib

if [ $? -ne 0 ]; then
    echo "Error: Failed to link final executable"
    exit 1
fi

# Copy to bin directory
echo "Copying executable to bin/..."
cp gamer ../bin/

if [ $? -ne 0 ]; then
    echo "Error: Failed to copy executable"
    exit 1
fi

echo "Done! CUDA files recompiled with -G flag and executable updated."
