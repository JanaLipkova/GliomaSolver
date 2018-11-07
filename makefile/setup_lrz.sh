# Load modules
module purge
module load admin/1.0 lrz/default mkl/2017
module unload mpi.intel
module load gcc/7  mpi.ompi/2.1/gcc blast
module list

# Set up paths libraries
LIB_BASE=/home/hpc/txh01/di49zin/GliomaAdvance/lib
export LD_LIBRARY_PATH=$LIB_BASE/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LIB_BASE/myVTK/lib/vtk-5.4/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LIB_BASE/hypre-2.10.0b/src/hypre/lib/:$LD_LIBRARY_PATH
export PATH=$HOME/usr/torc/bin:$PATH

# needed on LRZ
export LANG=C
export LC_ALL=C

# NOTE: To use hypre library provided by LRZ do:
# module load hypre
# export LD_LIBRARY_PATH=/lrz/sys/libraries/hypre/2.11.2_impi51/lib/:$LD_LIBRARY_PATH


