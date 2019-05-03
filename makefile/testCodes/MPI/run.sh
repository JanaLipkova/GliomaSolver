# Threads: M MPI ranks
M=4

# two options, usually both works fine
#mpirun -n $M ./a.out
mpirun -np $M ./a.out
