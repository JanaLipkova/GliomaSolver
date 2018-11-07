#!/bin/sh

# USER INPUT HERE:
#------------------------------------------------------
DataPath="/Users/lipkova 1/WORK/GliomaSolver/Simulations/PatientInference/InputData_dat/"
SolverPath="/Users/lipkova 1/WORK/GliomaSolver/"
Nsamples=1000
#------------------- User input stops here --------------


echo "------------------------------------------------------"
echo "          SETTING-UP INFERENCE ENVIROMENT             "
echo "------------------------------------------------------"
echo " "

# Folders name
PrepFolder=Preprocessing
TMP=Inference/TMPTMP/
LIKE=Inference/Likelihood/makefile/


# Get Inference tools and scripts
cp -r "$SolverPath/tools/pi4u_lite/Inference" .
cp -r "$SolverPath/tools/scripts" .


echo ">>> Data Preprocessing <<<"
echo "---------------------------------------"
# Map data to MRAG grid
mkdir -p $PrepFolder
cp brain $PrepFolder
cp scripts/mapDataToMRAGrid.sh $PrepFolder
cp scripts/writeRunAll.sh $PrepFolder/
cp scripts/writeTMCMCpar.sh $PrepFolder/
cp scripts/GenericPrior.txt $PrepFolder/


cd $PrepFolder
./mapDataToMRAGrid.sh "$DataPath"
echo "---------------------------------------"

echo " "
echo ">>> Setting up Patient-Specific Enviroment <<<"
echo "---------------------------------------"
./writeRunAll.sh "$DataPath"
./writeTMCMCpar.sh $Nsamples

# Copy where they belong
cp brain ../$TMP
cp *dat ../$TMP
cp runAll.sh ../$TMP
cp tmcmc_glioma.par ../Inference/ 

#echo " "
echo ">>> Creating executables <<<"
echo "---------------------------------------"
cd ../$LIKE
make clean && make
cp likelihood ../../TMPTMP
cd ../../
make clean && make

echo "------------------------------------------------------"
echo "              DONE WITH THE SET-UP                    "
echo "------------------------------------------------------"
echo "Run the inference by:"
echo "./run_inference.sh"
echo "For long simulation it is recomend to use tmux"


