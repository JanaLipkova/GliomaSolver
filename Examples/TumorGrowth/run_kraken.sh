# Threads:
N=8
export OMP_NUM_THREADS=$N

# brain simulation set up
program=brain
model=RD
verbose=1
profiler=1
adaptive=1
vtk=1
dumpfreq=50

# Path to patient data
#DATA_BASE=/home/jana/Work/GliomaAdvance/GliomaSolver/Anatomy
DATA_BASE=../../Anatomy
PatFileName=$DATA_BASE/Patient00/P00

#tumor parameters
Dw=0.0013
rho=0.025
Tend=300

echo "In the directory: $PWD"
echo "Running program on $N nodes."

./$program -nthreads $N -model $model -profiler $profiler -verbose $verbose -adaptive $adaptive -PatFileName $PatFileName -vtk $vtk -dumpfreq $dumpfreq -Dw $Dw -rho $rho -Tend $Tend 
