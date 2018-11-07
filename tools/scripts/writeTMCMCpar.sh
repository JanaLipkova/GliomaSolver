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
done < LogPriorRange.txt


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
B0               ${myArray[0]}
B1               ${myArray[1]}
B2                ${myArray[2]}
B3               ${myArray[3]}
B4               ${myArray[4]}
B5               ${myArray[5]}
B6               ${myArray[6]}
B7               ${myArray[7]}
B8               ${myArray[8]}
B9               ${myArray[9]}
B10              ${myArray[10]}

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


