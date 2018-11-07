/*
 *  masterslave.c
 *  TORC_Lite
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <torc.h>
#include <sys/wait.h>
int times = 0;

void slave(double *pin, double *out)
{
	double in;
	sleep(1);
	in = *pin;
	*out = sqrt(in);
	printf("slave in = %f, *out = %f\n", in, *out);

}

int main(int argc, char *argv[])
{
	int cnt = 4;
	double di;
	double *result;
	double *ii;
	int i;
	double t0, t1;

	torc_register_task(slave);

	printf("address(slave)=%p\n", slave);
	torc_init(argc, argv, MODE_MS);

	result = (double *)malloc(cnt*sizeof(double));
	ii = (double *)malloc(cnt*sizeof(double));

	//torc_taskinit();
	torc_enable_stealing();
	t0 = torc_gettime();	
	for (i=0; i<cnt; i++) {
		di = (double) (i+1);
		result[i] = 100 + i;
		torc_create(-1, slave, 2,
			1, MPI_DOUBLE, CALL_BY_COP,
			1, MPI_DOUBLE, CALL_BY_RES,
			&di, &result[i]);
	}
	torc_waitall();
	t1 = torc_gettime();

	for (i = 0; i < cnt; i++) {
		printf("Received: sqrt(%6.3f)=%6.3f\n",(double) (i+1), result[i]);
	}

	printf("Elapsed time: %.2lf seconds\n", t1-t0);
	torc_finalize();
	return 0;
}
