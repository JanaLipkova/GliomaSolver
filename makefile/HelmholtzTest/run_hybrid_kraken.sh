# Threads: N OMP, M MPI
N=1
M=8
export OMP_NUM_THREADS=$N

# brain simulation set up
program=brain
model=helmholtzTest
verbose=1
profiler=1
vtk=0

echo "In the directory: $PWD"
echo "Running program on $SLURM_NODES nodes with total $SLURM_NTASKS MPI tasks, each with $SLURM_CPUS_PER_TASK resp $N threads."


mpirun -np $M ./$program -nthreads $N -model $model -verbose $verbose -profiler $profiler -vtk $vtk

