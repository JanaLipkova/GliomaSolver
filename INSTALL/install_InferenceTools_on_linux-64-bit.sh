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

echo ">>> Downloading the libraries   <<<"
echo "--------------------------------------"
wget tdo.sk/~janka/GliomaSolverHome/libs/lib-linux-64-bit/inference_libs.tgz
tar -zxf inference_libs.tgz
mv inference_libs/* .
rm inference_libs.tgz
rm -r inference_libs

echo " "
echo "--------------------------------------"
echo ">>> Installing libraries   <<<"
echo "--------------------------------------"
gsl_src=gsl-src
mpich_src=mpich-3.2.1-src

tar -zxf ${gsl_src}.tgz
tar -zxf ${mpich_src}.tgz

mkdir -p gsl-install
mkdir -p mpich-install

echo " Installing mpich:"
echo "--------------------------------------"
cd ${mpich_src}
./configure --prefix="${LIB_BASE}"/mpich-install
make
make isntall
export PATH="${LIB_BASE}"/mpich-install/bin:$PATH
echo "---------------"
echo "mpicc is set to"
which mpicc
echo "---------------"
cd ../
rm ${mpich_src}.tgz
rm -r rm ${mpich_src}


echo " Installing gsl:"
echo "--------------------------------------"
cd ${gsl_src}
./autogen.sh
sudo apt-get install libtool
./configure --enable-maintainer-mode --disable-dynamic --prefix="${LIB_BASE}"/gsl-install
make
make install
cd ../
rm ${gsl_src}.tgz
rm ${gsl_src}

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
cp tools.make/setup_linux.sh setup_${UserName}.sh
cp tools.make/run_inference.sh .

sed -i 's|_USER_LIB_BASE_|'"${LIB_BASE}"'|g' Makefile
sed -i 's|_USER_LIB_BASE_|'"${LIB_BASE}"'|g' setup_${UserName}.sh
sed -i 's|_USER_LIB_BASE_|'"${LIB_BASE}"'|g' run_inference.sh

source setup_${UserName}.sh
make clean
make

echo " Installing pandoc:"
echo "--------------------------------------"
sudo apt-get install pandoc pandoc-citeproc texlive


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
