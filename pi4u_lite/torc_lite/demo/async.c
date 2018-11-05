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


int DATA_ENTRIES;	//	11
int lambda=4;

typedef struct data_s
{
	double in;
	double out;
	int node_id;
	int state;
} data_t;

data_t *result;


int find_tid()
{
	int i;
	for (i = 0; i < DATA_ENTRIES; i++)
	{
		if (result[i].state == 0)
		{
			result[i].state = 1;
			return i;
		}
	}
	return -1;
}


int taskid = 0;

pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
int completed = 0;
int exit_flag = 0;

void slave(int *ptid, double *pin);

void setexit()
{
	exit_flag = 1;
}

void callback(int *ptid, int *pnode_id, double *pout)
{
	int tid = *ptid;
	int node_id = *pnode_id;
	double out = *pout;

	pthread_mutex_lock(&m);
	completed++;
	if (completed >= 20) {
		exit_flag = 1;
	}

	printf("vvvvvvvvvvvvvvvv\n"); fflush(0);
	printf("cb: tid = %d, node_id= %d pout = %f, completed = %d\n", tid, node_id, out, completed); fflush(0);
	result[tid].out = out;
	result[tid].node_id = node_id;
	result[tid].state = 2;	// done

	if (completed % lambda == 0)
	{
		printf("!!! completed == %d\n", completed); fflush(0);
		printf("^^^^^^^^^^^^^^^^\n"); fflush(0);
	
		int i;
		for (i = 0; i < DATA_ENTRIES; i++)
		{
			if (result[i].state == 2) {
				printf("%d: %f %f %d\n", i, result[i].in, result[i].out, result[i].state);
			}
			if (result[i].state == 2) {
				if (result[i].out != sqrt(result[i].in)) {
					printf("XXX: res:%f vs ref:%f\n", result[i].out, sqrt(result[i].in));
				}
				result[i].state = 0;	// free
				result[i].in=-1;
				result[i].out=-1;
			}
		}


		if (exit_flag) {
			pthread_mutex_unlock(&m);
			return;
		}

		for (i = 0 ; i < lambda; i++)
		{
			tid = find_tid();	// result[tid].state = 1;	//pending
			double di = (double) (taskid++);
			result[tid].in = di;

			// I cannot spawn a normal task within the server thread
			printf("yy: spawning task %d for %f\n", tid, di);
			torc_create_detached(result[tid].node_id, slave, 2,
				1, MPI_INT, CALL_BY_COP,
				1, MPI_DOUBLE, CALL_BY_COP,
				&tid, &di);
		}
		pthread_mutex_unlock(&m);

	}
	else
	{
		pthread_mutex_unlock(&m);
		return;
	}

}


//void slave(int *ptid, double *pin, double *pout)
void slave(int *ptid, double *pin)
{
	int tid = *ptid;
	double in;

#if 0
	int i;
	for (i = 0; i < 10; i++)
	{
		if (exit_flag) break;
		sleep(1);
	}
	if (i != 10)
	{
		printf("Task %d was cancelled\n", tid);
	}
#else
	sleep(1);
#endif

	double out;
	double *pout = &out;
	in = *pin;
	*pout = sqrt(in);
//	printf("task %d: slave in = %f, *out = %f\n", tid, in, *pout); fflush(0);

	int node_id = torc_node_id();

	if (node_id == 0)
	{
		callback(&tid, &node_id, pout);
	}
	else
	{
		torc_create_direct(0, callback, 3,
			1, MPI_INT, CALL_BY_COP, 
			1, MPI_INT, CALL_BY_COP, 
			1, MPI_DOUBLE, CALL_BY_COP, 
			&tid, &node_id, pout);
	        torc_waitall3();
	}
}

void torc_dispatch()
{
	_torc_scheduler_loop(1);
}


int main(int argc, char *argv[])
{
	int cnt = 2*lambda-1;	// lambda + (lambda-1)
	double di;
	int i;
	double t0, t1;

	torc_register_task(slave);
	torc_register_task(callback);
	torc_register_task(setexit);

	printf("address(slave)=%p\n", slave);
	torc_init(argc, argv, MODE_MS);

	DATA_ENTRIES = 3*lambda-1;
	result = (data_t *)calloc(1, DATA_ENTRIES*sizeof(data_t));

	torc_enable_stealing();

	t0 = torc_gettime();	
	for (i=0; i<cnt; i++) {
		pthread_mutex_lock(&m);
		int tid = find_tid(); 	
		pthread_mutex_unlock(&m);

		double di = (double) (taskid++);
		result[tid].in = di;

		printf("yy: spawning task %d for %f\n", tid, di);

		torc_create_detached(tid % torc_num_workers(), slave, 2,
			1, MPI_INT, CALL_BY_COP,
			1, MPI_DOUBLE, CALL_BY_COP,
			&tid, &di);
	}

	while (!exit_flag)
	{
		torc_dispatch();
	}

	t1 = torc_gettime();

	for (i = 1; i < torc_num_nodes(); i++)
	{
		torc_create_direct(i, setexit, 0);
	        torc_waitall3();
	}

	fflush(0);

	printf("Elapsed time: %.2lf seconds\n", t1-t0);
	torc_finalize();
	return 0;
}
