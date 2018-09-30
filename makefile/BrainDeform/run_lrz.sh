#!/bin/bash
#SBATCH -o myjob.%j.%N.out
#SBATCH -D .
#SBATCH -J DefTest
#SBATCH --get-user-env
#SBATCH --clusters=mpp3
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
# cpus per task max value: mpp1 = 16, mpp2=28, mpp3=64
#SBATCH --export=NONE
#SBATCH --time=00:10:00


# modules:
source /etc/profile.d/modules.sh
module purge
module load admin/1.0 lrz/default mkl/2017
module unload mpi.intel
module load gcc/7  mpi.ompi/2.1/gcc 
module list

# libraries:
export LIB_BASE=$HOME/GliomaAdvance/lib
export LD_LIBRARY_PATH=$LIB_BASE/myVTK/lib/vtk-5.4/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LIB_BASE/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LIB_BASE/hypre-2.10.0b/src/hypre/lib/:$LD_LIBRARY_PATH

# Threads:
N=2
export OMP_NUM_THREADS=$N

# brain simulation set up
program=brain
model=deform
verbose=1
profiler=1
bDumpIC=0
dumpfreq=2
vtk=1
CFL=0.8

# model parameters
tend=10
rho=0.012
Dw=1.3e-03
kCSF=100
kWM=10
kGM=1
bMobility=0
ICtype=1

echo "In the directory: $PWD"
echo "Running program on $SLURM_NODES nodes with total $SLURM_NTASKS MPI tasks, each with $SLURM_CPUS_PER_TASK resp $N threads."

mpirun -np 64 ./$program -nthreads $N -model $model -verbose $verbose -profiler $profiler -PatFileName $PatFileName -vtk $vtk -bDumpIC $bDumpIC -dumpfreq $dumpfreq -CFL $CFL -tend $tend -rho $rho -Dw $Dw -kCSF $kCSF -kWM $kWM -kGM $kGM -bMobility $bMobility -ICtype $ICtype

