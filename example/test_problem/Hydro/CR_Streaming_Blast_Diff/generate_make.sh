# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --model=HYDRO --flu_scheme=MHM_RP --flux=HLLE --mhd=true --cr_streaming=true --slope=PLM \
                       --cosmic_ray=true --eos=COSMIC_RAY --hdf5=true --openmp=false --double=true "$@"
