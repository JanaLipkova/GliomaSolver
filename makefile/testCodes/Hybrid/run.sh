# Threads: N OMP, M MPI
N=4
M=2
export OMP_NUM_THREADS=$N

echo "Running program with total $M MPI tasks,each with $N OMP threads."
mpirun -np $M ./a.out
