/*
 *  torc_thread.c
 *  TORC_Lite
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */
#include <torc_internal.h>
#include <torc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* External declaration from lwrte.c */
//void scheduler_loop (int);

pthread_mutex_t __m = PTHREAD_MUTEX_INITIALIZER;
int __created = 0;
//pthread_barrier_t bar;
typedef void (*sched_f)();

void *_torc_worker (void *arg)
{
	long vp_id = (long) arg;
	torc_t *rte = (torc_t *)calloc (1, sizeof(torc_t));

	rte->vp_id = vp_id;
	_lock_init (&rte->lock);
	rte->work = (sched_f) _torc_scheduler_loop;

#if DBG
	printf("[RTE %p]: NODE %d: WORKER THREAD %ld --> 0x%lx\n", rte, torc_node_id(), rte->vp_id, pthread_self()); fflush(0);
#endif
	pthread_mutex_lock(&__m);
	__created++;
	pthread_mutex_unlock(&__m);

	worker_thread[rte->vp_id] = pthread_self();

	_torc_set_vpid(vp_id);
	_torc_set_currt(rte);

	int repeat;
	while (__created < kthreads) {
		pthread_mutex_lock(&__m);
		repeat = (__created < kthreads); 
		pthread_mutex_unlock(&__m);
		if (repeat) 
			thread_sleep(10);	//sched_yield();
		else
			break;
	}

	if (vp_id == 0) {
		enter_comm_cs();
		MPI_Barrier(MPI_COMM_WORLD);
		leave_comm_cs();
	}

	if ((torc_node_id() == 0) && (vp_id == 0)) return 0;
	_torc_scheduler_loop(0);	/* never returns */

	return 0;
}

void start_worker(long id)
{
	int res;
	pthread_t pth;
	pthread_attr_t attr;
	pthread_attr_init(&attr);

	pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
/**	res = pthread_attr_setstacksize(&attr, 64*1024*1024) **/
//	printf("res = %d\n", res);

	res = pthread_create(&pth, &attr, _torc_worker, (void *)id);
	if (res == 0)
		worker_thread[id] = pth;
	else
		Error("pthread_create failed");
}

void shutdown_worker(int id)
{
	int this_node = torc_node_id();

	node_info[this_node].nworkers--;
}

pthread_mutex_t al = PTHREAD_MUTEX_INITIALIZER;
static int active_workers;


void _torc_md_init()
{
	unsigned int i;
	int this_node = torc_node_id();

	pthread_key_create(&vp_key, NULL);
	pthread_key_create(&currt_key, NULL);

	if (torc_num_nodes() > 1) {
		start_server_thread();
	}

/*	node_info[this_node].nworkers = 1;*/
	//pthread_barrier_init(&bar, NULL, kthreads);
	for (i = 1; i<kthreads; i++) {
/*		node_info[this_node].nworkers++;*/
		start_worker((long)i);
	}

	active_workers = kthreads;
}


void _torc_md_end ()
{
	unsigned int i;
	unsigned int my_vp;
	int res;

	my_vp = _torc_get_vpid();
#if DBG
	printf("worker_thread %d exits\n", my_vp); fflush(0);
#endif
	if(my_vp!=0) {
		pthread_mutex_lock(&al);
		active_workers--;
		pthread_mutex_unlock(&al);
		pthread_exit(0);
	}

	if (my_vp==0) {
		while (1)
		{
		pthread_mutex_lock(&al);
		if (active_workers == 1)
		{
			pthread_mutex_unlock(&al);
			break;
		}
		else {
			pthread_mutex_unlock(&al);
			sched_yield();
		}
		}
		/* We need a barrier here to avoid potential deadlock problems */
		enter_comm_cs();
		MPI_Barrier(MPI_COMM_WORLD);
		leave_comm_cs();

		if (torc_num_nodes() > 1) {
			shutdown_server_thread();
		}
		_torc_stats();

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		exit(0);
	}
}


void thread_sleep(int ms)
{
	struct timespec req, rem;

//	sched_yield();
//	return 0;

	req.tv_sec = ms / 1000;
	req.tv_nsec = (ms % 1000)*1E6;
	nanosleep(&req, &rem);
}

void F77_FUNC_(torc_sleep, TORC_SLEEP)(int *ms)
{
	thread_sleep(*ms);
}

void _torc_set_vpid(long vp)
{
	pthread_setspecific(vp_key, (void *)vp);
}

long _torc_get_vpid()
{
	long vp;
	vp = (long) pthread_getspecific(vp_key);

	return vp;
}

void _torc_set_currt(torc_t *task)
{
	pthread_setspecific(currt_key, (void *)task);
}

torc_t *_torc_get_currt()
{
	torc_t *task;
	task = (torc_t *) pthread_getspecific(currt_key);

	return task;
}
