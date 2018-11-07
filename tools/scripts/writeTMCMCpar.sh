#!/bin/sh

if [ $# != 1 ]
   then
   echo "An input argument with number of TMCMC samples is required"
   exit
fi

Nsamples="$1"

# Read in log prior range for IC
let i=0
while read line; do
  myArray[i]="$line"
  ((++i))
done < logPriorIC.txt


# Write the TMCMC configuration file
scriptname=tmcmc_glioma.par
cat > $scriptname << EOF
# TMCMC and SN-TMCMC configuration file

# problem dimension
Nth             11

# Max number of generations
MaxStages       1000

# Population size
PopSize         ${Nsamples}

# Chain Lengths
MinChainLength  0
MaxChainLength  1
#MinChainLength 2
#MaxChainLength 1e6

# Boundaries
Bdef             0          4
B0              -8.9480   -3.2702
B1              -5.9145   -1.6607
B2               6.5280    9.8875
B3              ${myArray[0]}    ${myArray[1]} 
B4              ${myArray[2]}    ${myArray[3]}
B5              ${myArray[4]}    ${myArray[5]}
B6              -4.0174   -1.3863
B7              -0.5108    0.0198
B8              -0.5108   -0.2231
B9              -2.9957   -0.5108
B10             -2.9957   -2.3026

TolCOV          1  # 1 !!! don't play with fire
bbeta           0.2

#auxil_size     5
#auxil_data     0.1,0.2,0.3,0.4,5

# RNG seed
#seed           0
seed            280675

# Internal options
opt.MaxIter     1000
opt.Tol         1e-12
opt.Display     1
opt.Step        1e-4

# Online plotting of data, useful for 2D data and debugging
iplot           0

# SN-TMCMC specific options
diffstep        1e-4
posdef          4

EOF


