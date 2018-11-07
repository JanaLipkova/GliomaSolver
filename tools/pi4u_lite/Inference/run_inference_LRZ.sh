#!/bin/bash
#SBATCH -o myjob.%j.%N.out
#SBATCH -D .
#SBATCH -J UQ
#SBATCH --clusters=mpp3
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
# mpp3 64, mpp2 28 
#SBATCH --export=NONE
#SBATCH --time=06:00:00

# Modules
source /etc/profile.d/modules.sh
module purge
module load admin/1.0 lrz/default
module load gsl/2.3 gcc/4.9
module list

# Libraries
LIB_BASE=/home/hpc/txh01/di49zin/GliomaAdvance/lib
export LD_LIBRARY_PATH=$LIB_BASE/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LIB_BASE/myVTK/lib/vtk-5.4/:$LD_LIBRARY_PATH
export PATH=$HOME/usr/torc/bin:$PATH
export LD_LIBRARY_PATH=$HOME/usr/torc/bin:$LD_LIBRARY_PATH
export PATH=$PATH:$HOME/pi4u-libs/mpich-install/bin/
export LD_LIBRARY_PATH=$HOME/pi4u-libs/mpich-install/lib/:$LD_LIBRARY_PATH

echo "Using this mpicc:"
which mpicc

export LANG=C
export LC_ALL=C
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

echo "In the directory: $PWD"
echo "Running program on $SLURM_TASKS nodes with $SLURM_CPUS_PER_TASK tasks, each with $SLURM_CPUS_PER_TASK cores."

mpirun -env TORC_WORKERS 1 ./engine_tmcmc
#mpirun -np 64 -env TORC_WORKERS 1 ./engine_tmcmc
