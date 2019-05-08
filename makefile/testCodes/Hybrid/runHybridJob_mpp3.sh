#!/bin/bash
#SBATCH -o myjob.%j.%N.out
#SBATCH -D .
#SBATCH -J ex1
#SBATCH --clusters=mpp3
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=4
#SBATCH --get-user-env
# cpus per task max value: mpp1 = 16, mpp2=28, mpp3=64
#SBATCH --export=NONE
#SBATCH --time=00:02:00

source /etc/profile.d/modules.sh
module purge
module load admin/1.0 lrz/default 
module load gcc intel mkl
module list


export LANG=C
export LC_ALL=C
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


echo "====================================================="
echo "In the directory: $PWD"
echo "Running program on $SLURM_NODES nodes with total $SLURM_NTASKS MPI tasks, each with $SLURM_CPUS_PER_TASK OpenMP threads."
echo "====================================================="
#srun --cpus-per-task=2 -n 2 ./ex1
#mpiexec --perhost=2 -n 2 ./ex1
#mpirun -n 32 a.out
mpirun -n $SLURM_NTASKS ./a.out
