# Threads: N OMP, M MPI
N=2
M=4
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

# Path to patient data
DATA_BASE=/home/jana/Work/GliomaAdvance/GliomaSolver/Anatomy
PatFileName=$DATA_BASE/Patient00/P00

echo "In the directory: $PWD"
echo "Running program with total $M MPI tasks,each with $N OMP threads."


mpirun -np $M ./$program -nthreads $N -model $model -verbose $verbose -profiler $profiler -PatFileName $PatFileName -vtk $vtk -bDumpIC $bDumpIC -dumpfreq $dumpfreq -CFL $CFL -tend $tend -rho $rho -Dw $Dw -kCSF $kCSF -kWM $kWM -kGM $kGM -bMobility $bMobility -ICtype $ICtype


