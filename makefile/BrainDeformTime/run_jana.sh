# Threads:
N=4
export OMP_NUM_THREADS=$N

# brain simulation set up
program=brain
model=deformTime
verbose=1
profiler=1
adaptive=1
pID=100
vtk=1
dumpfreq=20
bDumpIC=3

#model parameters
Dw=0.0013
rho=0.025
Tend=500
betaCSF=1
betaWM=1
betaGM=1
mobCSF=0.1
mobWM=0.1
mobGM=0.1
CFL=0.9


echo "In the directory: $PWD"
echo "Running program on $N nodes."

./$program -nthreads $N -model $model -adaptive $adaptive -verbose $verbose -profiler $profiler -pID $pID -vtk $vtk -profiler $profiler -dumpfreq $dumpfreq -bDumpIC $bDumpIC -Dw $Dw -rho $rho -Tend $Tend -betaCSF $betaCSF -betaWM $betaWM -betaGM $betaGM -mobCSF $mobCSF -mobWM $mobWM -mobGM $mobGM -CFL $CFL
