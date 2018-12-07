#!/bin/sh
LIB_BASE="_USER_LIB_BASE_"
export PATH=${LIB_BASE}/mpich-install/bin:$PATH
export PATH=${LIB_BASE}/usr/torc/bin:$PATH

echo "Use mpicc:"
which mpicc

#Number of MPI ranks
M="$1"

mpirun -np $M ./engine_tmcmc

./extractInferenceOutput.sh
