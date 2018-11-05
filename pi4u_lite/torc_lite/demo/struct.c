/*
 *  struct.c
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

struct tdata {
	double in;
	double out;
	char txt[256];
};

void slave(struct tdata *datain, double *out)
{
	double in;
	printf("slave2...: datain=%p, msg=%s\n", datain, datain->txt);
	sleep(1);
	in = datain->in;
	*out = sqrt(in);
	datain->out = sqrt(in);
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
	struct tdata *ttd, *td;

	torc_register_task(slave);

	printf("address(slave)=%p\n", slave);
	torc_init(argc, argv, MODE_MS);

	result = (double *)malloc(cnt*sizeof(double));
	ii = (double *)malloc(cnt*sizeof(double));
	ttd = (struct tdata *) malloc(cnt*sizeof(struct tdata));

	t0 = torc_gettime();	
	for (i=0; i<cnt; i++) {
		//td = malloc(sizeof(struct tdata));
		td = &ttd[i];
		printf("td = %p\n", td);
		td->in = (double) (i+1);
		sprintf(td->txt, "hello from task %d", i); 
		result[i] = 100 + i;
		torc_create(-1, slave, 2,
			sizeof(struct tdata)/sizeof(long), MPI_LONG, CALL_BY_REF,
			1, MPI_DOUBLE, CALL_BY_RES,
			td, &result[i]);
	}
	torc_enable_stealing();
	torc_waitall();
	torc_disable_stealing();
	t1 = torc_gettime();
#if 1
	for (i = 0; i < cnt; i++) {
		printf("Received: sqrt(%6.3f)=%6.3f - %6.3f\n",(double) (i+1), result[i], ttd[i].out);
	}
#endif
	printf("Elapsed time: %.2lf seconds\n", t1-t0);
	torc_finalize();
	return 0;
}
