/*
 *  zerolength.c
 *  TORC_Lite
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/17.
 *  Copyright 2017 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <torc.h>
#include <sys/wait.h>
int times = 0;

void slave(double *pin, double *out, double *x, int *pn)
{
	int n = *pn;
	double in = *pin;

	if (n > 0) {
		int j;
		for (j = 0; j < n; j++) x[j] = (double) torc_worker_id(); 
	}

	*out = sqrt(in);
	printf("B slave in = %f, *out = %f, [x = %p, n = %d]\n", in, *out, x, n); fflush(0);
}

int main(int argc, char *argv[])
{
	int cnt = 4;
	double di;
	double *result;
	double *ii;
	int i;
	double t0, t1;

	int sz = 0;

	if (argc == 2) sz = atoi(argv[1]);

	srand48(33);

	torc_register_task(slave);

	printf("address(slave)=%p\n", slave);
	torc_init(argc, argv, MODE_MS);

	result = (double *)malloc(cnt*sizeof(double));
	ii = (double *)malloc(cnt*sizeof(double));

	double *x[cnt];

	//torc_enable_stealing();
	t0 = torc_gettime();	
	for (i=0; i<cnt; i++) {
		di = (double) (i+1);
		result[i] = 100 + i;
		if (sz == 0)
			x[i] = NULL;
		else
			x[i] = malloc(sz*sizeof(double));

		torc_create(-1, slave, 4,
			1, MPI_DOUBLE, CALL_BY_COP,
			1, MPI_DOUBLE, CALL_BY_RES,
			sz, MPI_DOUBLE, CALL_BY_RES,
			1, MPI_INT, CALL_BY_COP,
			&di, &result[i], x[i], &sz);
	}
	torc_waitall();
	t1 = torc_gettime();

	for (i=0; i<cnt; i++) {
		int j;
		for (j = 0; j < sz; j++) printf("%d %d -> %lf\n", i, j, x[i][j]);
	}

	for (i = 0; i < cnt; i++) {
		printf("Received: sqrt(%6.3f)=%6.3f\n",(double) (i+1), result[i]);
	}

	printf("Elapsed time: %.2lf seconds\n", t1-t0);
	torc_finalize();
	return 0;
}
