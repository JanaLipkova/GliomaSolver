#include <stdio.h>
#include <sys/time.h>

void empty()
{
}

double MPI_Wtime()
{
	struct timeval t;      
	gettimeofday(&t, NULL);
	return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
}

