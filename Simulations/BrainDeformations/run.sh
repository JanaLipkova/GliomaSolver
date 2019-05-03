# Threads: N OMP, M MPI
N=2
M=4
export OMP_NUM_THREADS=$N

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

echo "In the directory: $PWD"
echo "Running program with total $M MPI tasks,each with $N OMP threads."



mpirun -np $M ./$program -nthreads $N -model $model -verbose $verbose -profiler $profiler -PatFileName $PatFileName -vtk $vtk -bDumpIC $bDumpIC -dumpfreq $dumpfreq -CFL $CFL -tend $tend -rho $rho -Dw $Dw -kCSF $kCSF -kWM $kWM -kGM $kGM -bMobility $bMobility -mCSF $mCSF -mWM $mWM -mGM $mGM -ICtype $ICtype -icx $icx -icy $icy -icz $icz


