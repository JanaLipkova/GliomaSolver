# Threads:
N=4
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
PatFileName=../../Simulations/PatientInferenceGamma/InputData/Atlas/anatomy_dat/

#tumor parameters
Dw=0.0013
rho=0.025
Tend=300
icx=0.28
icy=0.67
icz=0.35

echo "In the directory: $PWD"
echo "Running program on $N nodes."

./$program -nthreads $N -model $model -profiler $profiler -verbose $verbose -adaptive $adaptive -PatFileName $PatFileName -vtk $vtk -dumpfreq $dumpfreq -Dw $Dw -rho $rho -Tend $Tend -icx $icx -icy $icy -icz $icz 
