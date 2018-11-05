#!/bin/bash

echo "Need as input path to data"

N=1
export OMP_NUM_THREADS=$N

# 1) Glioma growth
#------------------
program=brain
model=RD
verbose=1
adaptive=1
vtk=0
UQ=1
dumpfreq=20
PatFileName="$1"

./$program -nthreads $N -model $model -verbose $verbose -adaptive $adaptive -PatFileName $PatFileName -vtk $vtk -dumpfreq $dumpfreq -UQ $UQ

stepPET=3
stepMRI=1
bROI=1

# Likelihood
#------------------
./likelihood -stepPET $stepPET -stepMRI $stepMRI bROI $bROI

