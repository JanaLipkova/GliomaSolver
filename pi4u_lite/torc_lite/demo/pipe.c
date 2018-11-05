/*
 *  masterslave.c
 *  TORC_Lite
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <torc.h>
#include <unistd.h>

struct work
{
	char stage;
	int tid;
};

struct work stageA = {'A', -1};
struct work stageB = {'B', -1};
struct work stageC = {'C', -1};

void advance_and_execute(int tid)
{
	stageC.tid = stageB.tid;
	stageB.tid = stageA.tid;
	stageA.tid = tid;

	if (stageA.tid > -1) {
	} 
	if (stageB.tid > -1) {
	} 
	if (stageC.tid > -1) {
	}
	printf("[%d] task %2d: A=%2d B=%2d C=%2d\n", torc_worker_id(), tid, stageA.tid, stageB.tid, stageC.tid);
	sleep(torc_node_id()+1);
}

void task(int *p_tid)
{
	int tid = *p_tid;

	advance_and_execute(tid);

	int check = torc_fetch_work();

	if (check == 1) {
		return;
	}
	else {
		advance_and_execute(-1);
		advance_and_execute(-1);
	}
}

void spmd()
{
	int me = torc_worker_id();
	int workers = torc_num_workers();
	int cnt = 20;
	int i, iter;

	for (iter = 0; iter < 2; iter++)
	{
		MPI_Barrier(MPI_COMM_WORLD);

		if (me == 0) printf("iter = %d\n", iter);

		//if (me == 0)
		{
		for (i=me; i<cnt; i+=workers) {
			torc_task(me, (void *)task, 1, 1, MPI_INT, CALL_BY_COP, &i);
		}
		torc_i_enable_stealing();	// should be local
		}

		MPI_Barrier(MPI_COMM_WORLD);	// wait for all work to be spanwed (not actually necessary)

		torc_waitall2();	// cluster-wide 
	}
}


struct main_args
{
	int argc;
	char **argv;
};

void *posix(void *arg)
{
	struct main_args *ma = (struct main_args *) arg; 
	int i;
	double t0, t1;

	torc_register_task((void *)task);
	torc_register_task((void *)spmd);
	torc_init(ma->argc, ma->argv, MODE_MS);

	t0 = torc_gettime();	
	for (i=0; i<torc_num_workers(); i++) {
		torc_create(-1, (void *)spmd, 0);
	}
	torc_waitall();
	t1 = torc_gettime();

	printf("Elapsed time: %.2lf seconds\n", t1-t0);
	torc_finalize();
	return 0;
}

int main(int argc, char *argv[])
{
	int i;
	double t0, t1;

	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	
	struct main_args margs = {argc, argv};

	pthread_t pth;
	pthread_create(&pth, NULL, posix, &margs);
	pthread_join(pth, NULL);

	return 0;
}
