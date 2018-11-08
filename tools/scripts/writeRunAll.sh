#!/bin/sh

if [ $# != 1 ]
   then
   echo "Input argument with full path to patient data is required"
   exit
fi

InputFile="$1"

scriptname=runAll.sh
cat > $scriptname << EOF
#0) Enviroment
#------------------
N=1
export OMP_NUM_THREADS=\$N

# 1) Glioma growth
#------------------
program=brain
model=RD
verbose=0
adaptive=1
vtk=0
UQ=1
dumpfreq=20
PatFileName="${InputFile}/"

./\$program -nthreads \$N -model \$model -verbose \$verbose -adaptive \$adaptive -PatFileName \$PatFileName -vtk \$vtk -dumpfreq \$dumpfreq -UQ \$UQ

# Likelihood
#------------------
stepPET=3
stepMRI=1
bROI=1

./likelihood -stepPET \$stepPET -stepMRI \$stepMRI -bROI \$bROI
EOF

chmod 777 ${scriptname}
