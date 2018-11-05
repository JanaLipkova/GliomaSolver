/*
 *  fibo.c
 *  TORC_Lite
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */
#include <torc.h>
#include <math.h>
#include <stdio.h>

#define FIB_NUM		50
#define LEVELS		5

void slave (unsigned long *pn, unsigned long *res)
{
	int i = 0;
	unsigned long res1 = 0;
	unsigned long res2 = 0;
	int nodes;
	int rand_node;
	unsigned long n = *pn;

	if ((n == 0) || (n==1))
		*res = n;
	else {
		unsigned long n_1, n_2;

		n_1 = n-1;
		n_2 = n-2;
		if (n < (FIB_NUM-LEVELS)) {
		   slave(&n_1, &res1);
		   slave(&n_2, &res2);
		}
		else {
			torc_create(-1, slave, 2, 
				1, MPI_LONG, CALL_BY_COP,
				1, MPI_LONG, CALL_BY_RES,
				&n_1, &res1);

		   	//torc_create(-1, slave, 2,
			//	1, MPI_LONG, CALL_BY_COP,
			//	1, MPI_LONG, CALL_BY_RES,
			//	&n_2, &res2);
			slave(&n_2, &res2);


			torc_waitall();
		}
		*res = res1+res2;
	}
}


int main(int argc, char *argv[])
{
	unsigned long res;
	double time1, time2, elapsed;
	int vp;
	unsigned long n = FIB_NUM;

	if (argc == 2) n = atoi(argv[1]);

	torc_register_task(slave);
	torc_init(argc, argv, MODE_MS);
	printf ("Executing in main [%ld]\n", n);
	torc_enable_stealing();
	time1 = torc_gettime();
	slave(&n, &res);
	time2 = torc_gettime();
	elapsed = time2 - time1;
	printf ("fib(%ld) = %ld\n", n, res);
	printf ("Elapsed :  %lf secs \n", elapsed);

	torc_finalize();
	return 0;
}
