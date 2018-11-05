/*
 *  broadcast.c
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

int times = 0;

int gi;
int gti[16];

double gd;
double gtd[16];

void slave()
{
	int me = torc_worker_id();
	int node = torc_node_id();

	printf("worker[%d] node[%d]: %d %f %d %f\n", me, node, gi, gd, gti[5], gtd[5]);
	sleep(1);
	return;
}

int main(int argc, char *argv[])
{
	int i;
	int ntasks;

	gi = 5;
	gd = 5.0;
	for (i = 0; i < 16; i++) {
		gti[i] = 100;
		gtd[i] = 100.0;
	}

	torc_register_task(slave);
	torc_init(argc, argv, MODE_MS);

	gi = 23;
	gd = 23.0;
	for (i = 0; i < 16; i++) {
		gti[i] = torc_worker_id() + 1000;
		gtd[i] = gti[i] + 1;
	}

	torc_broadcast(&gi, 1, MPI_INT);
	torc_broadcast(&gd, 1, MPI_DOUBLE);
	torc_broadcast(&gti, 16, MPI_INT);
	torc_broadcast(&gtd, 16, MPI_DOUBLE);

	ntasks = torc_num_workers();
	for (i=0; i<ntasks; i++) {
		torc_create(-1, slave, 0);
	}
	torc_waitall();

	torc_finalize();
	return 0;
}
