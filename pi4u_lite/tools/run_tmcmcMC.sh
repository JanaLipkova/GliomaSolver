# set up module, paths
#	TMCMC stuff
export PATH=$PATH:/cluster/home/mavt/chatzidp/usr/mpich3/bin/
export LD_LIBRARY_PATH=/cluster/home/mavt/chatzidp/usr/mpich3/lib/:$LD_LIBRARY_PATH
export PATH=$HOME/usr/torc/bin:$PATH

module load gcc

#672 Cores -> 14 nodes -> 14x12 = 168 MPI, 4 threads per MPI rank
#480 Cores -> 10 nodes -> 10x12 = 120 MPI, 4 threads per MPI rank
#432 Cores ->  9 nodes ->  9x12 = 108 MPI, 4 threads per MPI rank
#384 Cores ->  8 nodes ->  8x12 =  96 MPI, 4 threads per MPI rank
#336 Cores ->  7 nodes ->  7x12 =  84 MPI, 4 threads per MPI rank
cores=672

bsub -J P31_6KJ1 -W 11:55 -n $cores -R "span[ptile=48]" -R "select[model==Opteron6174]" < runtmcmc.lsf.sh
bsub -J P31_6KJ2 -w "ended(P31_6KJ1)" -W 11:55 -n $cores -R "span[ptile=48]" -R "select[model==Opteron6174]" < runtmcmc.lsf.sh
bsub -J P31_6KJ3 -w "ended(P31_6KJ2)" -W 11:55 -n $cores -R "span[ptile=48]" -R "select[model==Opteron6174]" < runtmcmc.lsf.sh
bsub -J P31_6KJ4 -w "ended(P31_6KJ3)" -W 11:55 -n $cores -R "span[ptile=48]" -R "select[model==Opteron6174]" < runtmcmc.lsf.sh
bsub -J P31_6KJ5 -w "ended(P31_6KJ4)" -W 11:55 -n $cores -R "span[ptile=48]" -R "select[model==Opteron6174]" < runtmcmc.lsf.sh
