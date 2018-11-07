#0) Enviroment
#------------------
N=1
export OMP_NUM_THREADS=

# 1) Glioma growth
#------------------
program=brain
model=RD
verbose=1
adaptive=1
vtk=0
UQ=1
dumpfreq=20
PatFileName="/Users/lipkova 1/WORK/GliomaSolver/tools/scripts"

./$program -nthreads $N -model $model -verbose $verbose -adaptive $adaptive -PatFileName $PatFileName -vtk $vtk -dumpfreq $dumpfreq -UQ $UQ

# Likelihood
#------------------
stepPET=3
stepMRI=1
bROI=1

./likelihood -stepPET $stepPET -stepMRI $stepMRI bROI $bROI
