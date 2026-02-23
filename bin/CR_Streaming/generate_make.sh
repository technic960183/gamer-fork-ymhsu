# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --model=HYDRO --flu_scheme=MHM_RP --flux=HLLE --mhd=true --cr_streaming=true --slope=PLM --bitwise_reproducibility=true\
                       --cosmic_ray=true --eos=COSMIC_RAY --cr_diffusion=true --hdf5=true --openmp=false --double=true --gpu=true "$@"
