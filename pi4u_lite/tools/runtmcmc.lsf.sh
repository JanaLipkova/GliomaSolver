#672 Cores -> 14 nodes -> 14x12 = 168 MPI, 4 threads per MPI rank
#480 Cores -> 10 nodes -> 10x12 = 120 MPI, 4 threads per MPI rank
#432 Cores ->  9 nodes ->  9x12 = 108 MPI, 4 threads per MPI rank
#384 Cores ->  8 nodes ->  8x12 =  96 MPI, 4 threads per MPI rank
#336 Cores ->  7 nodes ->  7x12 =  84 MPI, 4 threads per MPI rank
MPI=168
ppn=12
mpich_run -n $MPI -ppn $ppn -env TORC_WORKERS 1 -launcher ssh -f $LSB_DJOB_HOSTFILE ./engine_tmcmc



