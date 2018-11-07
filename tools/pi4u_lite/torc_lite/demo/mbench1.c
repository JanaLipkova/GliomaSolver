/*
 *  mbench1.c
 *  TORC_Lite
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

/* creation overhead */
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <torc.h>
#ifndef WIN32
#include <sys/time.h>
#endif

#define	DEF_NTHREADS	100
unsigned int count;

volatile unsigned long counter = 0;

void slave (void *arg)
{
//	counter++;
//	printf(".");
}

int my_main (void)
{
	unsigned int i = 0, nthreads = count;
	double time1, time2, time3;
	double creation, execution, total;
	int iteration, niterations = 10;
	double sum_avg_total = 0.0;
	double min_avg = 1e6, max_avg = -1e6;
	
	printf ("Executing in my_main \n");
	fflush(NULL);

	for (iteration = 0; iteration < niterations; iteration++)
	{
		printf("ITERATION (%d)\n", iteration); fflush(NULL);
		counter = 0;
		time1 = torc_gettime();
		i = 0;
		while (i< nthreads) {
			torc_create(-1, slave, 1, 1, MPI_LONG, CALL_BY_COP, &i);
			++i;
		}
		time2 = torc_gettime();
		torc_waitall();
		time3 = torc_gettime();

		creation = (time2 - time1)*1.0E3;
		execution = (time3 - time2)*1.0E3;
		total = creation + execution;
		
		counter = nthreads;

		sum_avg_total += total/counter;
		if (iteration == 0) {
			min_avg = max_avg = (total/counter);
		}
		else {
			if (min_avg > (total/counter))
				min_avg =  (total/counter);
			if (max_avg < (total/counter))
				max_avg =  (total/counter);
		}
		fflush(NULL);
	}
	
	printf("===================================================\n");
	printf("-----------> Total average = %.3lf msecs <-----------\n", sum_avg_total/niterations);
	printf("-----------> Total minimum = %.3lf msecs <-----------\n", min_avg);
	printf("-----------> Total maximum = %.3lf msecs <-----------\n", max_avg);
	printf("===================================================\n");
	fflush(NULL);
	return 55;
}

int main (int argc, char * argv [])
{
	char * p;
	unsigned int nthreads;


	if (argc != 1)
		count = atoi(argv[1]);
	else
		count = DEF_NTHREADS;

	counter = 0;
	nthreads = count;

	torc_register_task(slave);	
	torc_init(argc, argv, MODE_MS);

	my_main();

	torc_finalize();

	return 0;
}
