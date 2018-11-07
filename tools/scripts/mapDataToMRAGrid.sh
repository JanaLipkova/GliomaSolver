#!/bin/sh

if [ $# != 1 ]
   then
   echo "Input argument with full path to patient data is required"
   exit
fi

# Threads:
N=4
export OMP_NUM_THREADS=$N

# brain simulation set up
program=brain
model=UQprep
verbose=0
vtk=0

# Path to patient data
PatFileName="$1"

./$program -nthreads $N -model $model -verbose $verbose -PatFileName $PatFileName -vtk $vtk
