#!/bin/bash
#SBATCH -o myjob.%j.%N.out
#SBATCH -D .
#SBATCH -J DefTest
#SBATCH --get-user-env
#SBATCH --clusters=mpp3
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --get-user-env
# cpus per task max value: mpp1 = 16, mpp2=28, mpp3=64
#SBATCH --export=NONE
#SBATCH --time=00:02:00

source /etc/profile.d/modules.sh
module purge
module load admin/1.0 lrz/default
module load gcc intel mkl
module list

# libraries:
export LIB_BASE=$HOME/GliomaAdvance/lib
export LD_LIBRARY_PATH=$LIB_BASE/myVTK/lib/vtk-5.4/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LIB_BASE/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LIB_BASE/hypre-2.10.0b/src/hypre/lib/:$LD_LIBRARY_PATH

# Threads:
export LANG=C
export LC_ALL=C
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# brain simulation set up
program=brain
model=deform
verbose=0
profiler=1
bDumpIC=0
dumpfreq=30
vtk=1
CFL=0.8

# model parameters
tend=300
rho=0.012
Dw=1.3e-03
kCSF=0.2
kWM=0.02
kGM=0.002
bMobility=0
mCSF=0.1
mWM=0.1
mGM=0.1
ICtype=1
icx=0.6
icy=0.7
icz=0.5

# Path to patient data
PatFileName=Anatomy/Patient00/P00_

echo "====================================================="
echo "In the directory: $PWD"
echo "Running program on $SLURM_NODES nodes with total $SLURM_NTASKS MPI tasks, each with $SLURM_CPUS_PER_TASK OpenMP threads."
echo "====================================================="

mpirun -np $SLURM_NTASKS ./$program -nthreads $SLURM_CPUS_PER_TASK -model $model -verbose $verbose -profiler $profiler -PatFileName $PatFileName -vtk $vtk -bDumpIC $bDumpIC -dumpfreq $dumpfreq -CFL $CFL -tend $tend -rho $rho -Dw $Dw -kCSF $kCSF -kWM $kWM -kGM $kGM -bMobility $bMobility -mCSF $mCSF -mWM $mWM -mGM $mGM -ICtype $ICtype -icx $icx -icy $icy -icz $icz

