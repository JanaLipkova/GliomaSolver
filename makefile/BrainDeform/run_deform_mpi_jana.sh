# Threads:
N=2
export OMP_NUM_THREADS=$N

# brain simulation set up
program=brain
model=deform
adaptive=0
verbose=1
profiler=0
pID=0
bDumpIC=0
dumpfreq=10
vtk=1
CFL=0.8
tend=600
rho=0.012
Dw=1.3e-04
kCSF=100
kWM=10
kGM=1
bMobility=0

echo "In the directory: $PWD"
echo "Running program with total $SLURM_NTASKS MPI tasks,  each with $SLURM_CPUS_PER_TASK."


mpirun -np 1 ./$program -n 2 -nthreads $N -model $model -verbose $verbose -profiler $profiler -pID $pID -vtk $vtk -bDumpIC $bDumpIC -dumpfreq $dumpfreq -CFL $CFL -tend $tend -rho $rho -Dw $Dw -kCSF $kCSF -kWM $kWM -kGM $kGM -bMobility $bMobility

