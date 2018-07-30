# libraries:
#export MY_BASE=$HOME/GliomaAdvance/lib
#export LD_LIBRARY_PATH=$MY_BASE/myVTK/lib/vtk-5.4/:$MY_BASE/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$LD_LIBRARY_PATH

# Threads:
export OMP_NUM_THREADS=4

# brain simulation set up
program=brain
model=RD
verbose=1
adaptive=1
pID=0
vtk=1
dumpfreq=50


echo "In the directory: $PWD"
echo "Running program on $SLURM_NODES nodes, with $SLURM_CPUS_ON_NODE cores on node, each with $SLURM_CPUS_PER_TASK cores."

./$program -nthreads $SLURM_CPUS_ON_NODE -model $model -verbose $verbose -adaptive $adaptive -pID $pID -vtk $vtk -dumpfreq $dumpfreq 
