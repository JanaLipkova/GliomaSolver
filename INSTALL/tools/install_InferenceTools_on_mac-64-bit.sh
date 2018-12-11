#!/bin/bash
echo "==================================================="
echo "  Installing libraries for the Bayesian Inference  "
echo "==================================================="
InstallDir=$(pwd)
SolverDir=$(dirname "${InstallDir}")
cd "${SolverDir}"
mkdir -p lib
cd lib
LIB_BASE=$(pwd)

echo " "
echo "--------------------------------------"
echo ">>> Installing libraries   <<<"
echo "--------------------------------------"

#echo " Installing mpich:"
#echo "--------------------------------------"
#sudo port install mpich
#export PATH="${LIB_BASE}"/mpich-install/bin:$PATH
#echo "---------------"
#echo "mpicc is set to"
#which mpicc
#echo "---------------"

echo " Installing gsl:"
echo "--------------------------------------"
sudo port install gsl

echo " Installing torc:"
echo "--------------------------------------"
cd "${SolverDir}"/tools/pi4u_lite/torc_lite/
autoreconf
automake --add-missing
./configure CC=mpicc --prefix="${LIB_BASE}"/usr/torc --with-maxnodes=1024
make
make install
export PATH="${LIB_BASE}"/usr/torc/bin:$PATH


echo " Compiling Inference tools:"
echo "--------------------------------------"
cd ../Inference
UserName=$(hostname -s)
cp tools.make/Makefile .
cp tools.make/setup_mac.sh setup_${UserName}.sh
cp tools.make/run_inference.sh .

sed -i 's|"_USER_LIB_BASE_"/gsl-install/include/|"/opt/local/include/gsl/"|g' Makefile
sed -i 's|"_USER_LIB_BASE_"/gsl-install/lib/|"/opt/local/lib/"|g' Makefile
sed -i 's|_USER_LIB_BASE_|'"${LIB_BASE}"'|g' setup_${UserName}.sh
sed -i 's|_USER_LIB_BASE_|'"${LIB_BASE}"'|g' run_inference.sh

source setup_${UserName}.sh
make clean
make

echo " Installing pandoc:"
echo "--------------------------------------"
sudo port install pandoc


echo "  "
echo "==============================================="
echo ">>>     The Installation is completed       <<<"
echo "==============================================="
echo " "
echo "Please update your .bashrc or similar file to include path to the installed libraries."
echo " To do so, put the text between the following dashed lines into your .bashrc file:"
echo "---------------------------------------"
cat setup_${UserName}.sh
echo "---------------------------------------"
echo "To run the Inference use ./run_inference.sh or ./run_inference_LRZ.sh"
echo "For more info see the user manual"
echo "==============================================="
