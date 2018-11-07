# Load modules
module purge
module load admin/1.0 lrz/default intel/17.0 mkl/2017
module load gsl/2.3 blast gcc/4.9
echo "Currently loaded modules:"

# Set up paths libraries
LIB_BASE=/home/hpc/txh01/di49zin/GliomaAdvance/lib
export LD_LIBRARY_PATH=$LIB_BASE/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LIB_BASE/myVTK/lib/vtk-5.4/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LIB_BASE/hypre-2.10.0b/src/hypre/lib/:$LD_LIBRARY_PATH
export PATH=$HOME/usr/torc/bin:$PATH
export PATH=$PATH:$HOME/pi4u-libs/mpich-install/bin/
export LD_LIBRARY_PATH=$HOME/pi4u-libs/mpich-install/lib/:$LD_LIBRARY_PATH
echo "we use this mpicc:"
which mpicc

#needed on LRZ for MAC login discompatiblity
export LANG=C
export LC_ALL=C

# NOTE: To use hypre library provided by LRZ do:
# module load hypre
# export LD_LIBRARY_PATH=/lrz/sys/libraries/hypre/2.11.2_impi51/lib/:$LD_LIBRARY_PATH


