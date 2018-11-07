/*
 *  torc_server.c
 *  TORC_Lite
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */
#include <torc_internal.h>
#include <torc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#define DBG	1

#define torc_desc_size sizeof(torc_t)	/* should be less */
extern MPI_Comm comm_out;

static torc_t no_work_desc;
volatile int termination_flag = 0;

int process_a_received_descriptor(torc_t *work/*, int tag1*/)
{
	int reuse = 0;
	torc_t *stolen_work, *parent;
	MPI_Status status;
	int i, istat, tag;
	char *mem;

	work->next = NULL;

#if DBG
	if (torc_node_id() == work->sourcenode) {
		printf("WHY?\n");
	}
#endif

	tag = work->sourcevpid;
	if ((tag < 0) || (tag > MAX_NVPS)) {	/* tag == MAX_NVPS occurs with asynchronous stealing */
		//Error1("Invalid message tag %d", tag);
		printf("...Invalid message tag %d from node %d [type=%d]\n", tag, work->sourcenode, work->type); fflush(0);
		MPI_Abort(MPI_COMM_WORLD, 1);
		return 1;
	}

	switch (work->type) {
	case TORC_ANSWER:
		{
#if DBG
		printf("Server %d accepted from %d, narg = %d [ANSWER]\n", torc_node_id(), work->sourcenode, work->narg);fflush(stdout);
#endif
		for (i = 0; i < work->narg; i++) { /* receive the results, if any */
			if (work->quantity[i] == 0) continue;
			if ((work->callway[i] == CALL_BY_RES) || (work->callway[i] == CALL_BY_REF)) {
				enter_comm_cs();
#if 1
				/* yyy */
				work->dtype[i] = _torc_b2mpi_type(work->btype[i]);
#endif
				istat = MPI_Recv((void *)work->localarg[i], work->quantity[i], work->dtype[i], work->sourcenode, tag, comm_out, &status);
				leave_comm_cs();
			}
#if DBG
			printf("received data from %d\n", work->sourcenode);	fflush(0);
#endif
		}
		for (i = 0; i < work->narg; i++) {
			if (work->quantity[i] == 0) continue;
			if (work->callway[i] == CALL_BY_COP2) {
				free((void *)work->localarg[i]);
				work->localarg[i] = 0;
			}
		}

		if (work->parent) {
			//printf("satisfying dependencies\n"); fflush(0);
			_torc_depsatisfy(work->parent);
			//printf("satisyedg dependencies\n"); fflush(0);
		}
		reuse = 1;
		return reuse;
		}
		break;


	case TORC_NORMAL_ENQUEUE:
		{
		if (work->homenode == torc_node_id()) {
			torc_to_i_rq(work);
			reuse = 0;
			return reuse;
		}
		else {	/* homenode != torc_node_id() */
#if DBG
			printf("Server %d accepted from %d, narg = %d [WORK]\n", torc_node_id(), work->sourcenode, work->narg);
			fflush(stdout);
#endif
			receive_arguments(work, tag);
			
#if 1
			if (work->rte_type == 20) {     /* direct execution */

				_torc_core_execution(work);
#if 0
				/* send the results back to the worker that submitted the task */
				direct_send_descriptor_answer(TORC_DIRECT_ANSWER, work->homenode, work->sourcevpid, work);
#else
				send_descriptor(work->homenode, work, TORC_ANSWER);
#endif
				reuse = 0;
				return reuse;
			}

#endif

#if DBG
			printf("SERVER: ENQ : Queue = %d Private = %d Front = %d\n", work->target_queue, work->insert_private, work->insert_in_front);
#endif

			if (work->insert_private) {	/* 1: per node - 2 : per virtual processor */
				if (work->insert_in_front)
					torc_to_i_pq(work);
				else
					torc_to_i_pq_end(work);
			}
			else {	/* 0: public queues */
				if (work->insert_in_front)
					torc_to_i_rq(work);
				else
					torc_to_i_rq_end(work);
			}

			reuse = 0;
			return reuse;
		}
	}
	break;

	case DIRECT_SYNCHRONOUS_STEALING_REQUEST:
#if DBG
		printf("Server %d received request for synchronous stealing\n", torc_node_id()); fflush(0);
#endif
		steal_attempts++;
		stolen_work = torc_i_rq_dequeue(0);
//		if (stolen_work == NULL) stolen_work = torc_i_rq_dequeue(1);
                for (i = 1; i < 10; i++) {
                        if (stolen_work == NULL) stolen_work = torc_i_rq_dequeue(i);
                        else break;
                }

		if (stolen_work != NULL) {
			direct_send_descriptor(DIRECT_SYNCHRONOUS_STEALING_REQUEST, work->sourcenode, work->sourcevpid, stolen_work);
			steal_served++;
		}
		else {
			direct_send_descriptor(DIRECT_SYNCHRONOUS_STEALING_REQUEST, work->sourcenode, work->sourcevpid, &no_work_desc);
		}
		return reuse = 0;
		break;

	case TERMINATE_LOCAL_SERVER_THREAD:
		/*thread_sleep(100);*/
		pthread_exit(0);
		reuse = 1;
		return reuse;
		break;

		
	case TERMINATE_WORKER_THREADS:
		termination_flag = 1;
		if (work->localarg[0] != torc_node_id()) {
			appl_finished++;
		}
#if DBG
		printf("Server %d will exit\n", torc_node_id()); fflush(0);
#endif
		reuse = 1;
		return reuse;
		break;

	case ENABLE_INTERNODE_STEALING:
		internode_stealing = 1;
		reuse = 1;
		break;

	case DISABLE_INTERNODE_STEALING:
		internode_stealing = 0;
		reuse = 1;
		break;
	
	case RESET_STATISTICS:
		torc_reset_statistics();
		reuse = 1;
		break;

        case TORC_BCAST:              
		{
		void *va = (void *)work->localarg[1];
		int count = work->localarg[2];
#if 1
		/* yyy */
		MPI_Datatype dtype = _torc_b2mpi_type(work->localarg[3]);
#endif
		//MPI_Datatype dtype = work->localarg[3];

		printf("TORC_BCAST: %p %d\n", va, count); fflush(0);
		enter_comm_cs();
		istat = MPI_Recv((void *)va, count, dtype, work->sourcenode, tag, comm_out, &status);
		leave_comm_cs();
		reuse = 1;
		}
		return reuse;
		break;                

	default:
		Error1("Unkown descriptor type on node %d", torc_node_id());
		break;
	}

	reuse = 1;
	return reuse;
}


/*#define DIRECT_REUSE */

void *server_loop (void *arg)
{
	torc_t *work;
	int reuse = 0;
	MPI_Status status;
	int istat;

#if DBG
	printf("Server %d begins....\n", torc_node_id()); fflush(stdout);
#endif
	memset(&no_work_desc, 0, sizeof(torc_t));
	no_work_desc.type = TORC_NO_WORK;	/* ... */

	while (1) {			
		if (!reuse)
			work = _torc_get_reused_desc();

		memset(work, 0, sizeof(torc_t));	/* lock ..?*/
#if DBG
		printf("Server %d waits for a descriptor ....\n", torc_node_id()); fflush(0);
#endif
		if (thread_safe) {
//#define _MONTE_ROSA_
#if !defined(_MONTE_ROSA_)
			istat = MPI_Recv(work, torc_desc_size, MPI_CHAR, MPI_ANY_SOURCE, MAX_NVPS, comm_out, &status);
#else
			int flag = 0;
			MPI_Request request;

			//enter_comm_cs();
			istat = MPI_Irecv(work, torc_desc_size, MPI_CHAR, MPI_ANY_SOURCE, MAX_NVPS, comm_out, &request);
			//leave_comm_cs();
			while (1) {
				//if (appl_finished == 1) pthread_exit(0);
				if (termination_flag >= 1) pthread_exit(0);
				//enter_comm_cs();
				MPI_Test(&request, &flag, &status);
				//leave_comm_cs();
				if (flag == 1) {
					break;
				} else {
					thread_sleep(yieldtime);
				}
			}
#endif
		}
		else {
			int flag = 0;
			MPI_Request request;

//			if (appl_finished == 1) {
			if (termination_flag >= 1) {
				pthread_exit(0);
			}

			enter_comm_cs();
			istat = MPI_Irecv(work, torc_desc_size, MPI_CHAR, MPI_ANY_SOURCE, MAX_NVPS, comm_out, &request);
			leave_comm_cs();
			while (1) {
				//if (appl_finished == 1) {
				if (termination_flag >=1 ) {
					printf("server threads exits!\n"); fflush(0);
					pthread_exit(0);
				}
				enter_comm_cs();
				MPI_Test(&request, &flag, &status);
				leave_comm_cs();
				if (flag == 1) {
					break;
				} else {
					thread_sleep(yieldtime);
				}
			}
		}

		reuse = process_a_received_descriptor(work /*,0 */);
		if (reuse) {
			_torc_put_reused_desc(work);
			reuse = 0;
		}

	}

	return 0;
}

/******************************************************************************/
/******************************************************************************/
static pthread_mutex_t server_thread_lock = PTHREAD_MUTEX_INITIALIZER;
static int server_thread_alive = 0;

void start_server_thread()
{
	int res;
	pthread_attr_t attr;

	pthread_attr_init(&attr);
	pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

	pthread_mutex_lock(&server_thread_lock);

	if (server_thread_alive == 1) {
		pthread_mutex_unlock(&server_thread_lock);
		return; 
	}

	res = pthread_create(&server_thread, &attr, server_loop, NULL);
	if (res != 0) {
		Error("server thread was not created!\n");
	}

	server_thread_alive = 1;

	pthread_mutex_unlock(&server_thread_lock);

}

void shutdown_server_thread()
{
	static torc_t mydata;
	
#if DBG
	printf("[%d]: Terminating local server thread....\n", torc_node_id()); fflush(0);
#endif

	pthread_mutex_lock(&server_thread_lock);
	if (server_thread_alive == 0) {
		pthread_mutex_unlock(&server_thread_lock);
		return; 
	}

	memset(&mydata, 0, sizeof(mydata));
	if (thread_safe)
		send_descriptor(torc_node_id(), &mydata, TERMINATE_LOCAL_SERVER_THREAD);
	else
		termination_flag = 1;

	pthread_join(server_thread, NULL);

	server_thread_alive = 0;
	pthread_mutex_unlock(&server_thread_lock);

}

/*************************************************************************/
/**********************     THREAD MANAGEMENT      **********************/
/*************************************************************************/

void terminate_workers()
{
	torc_t mydata;
	int node, nnodes = torc_num_nodes();
	
#if DBG
	printf("Terminating worker threads ...\n");
#endif
	memset(&mydata, 0, sizeof(torc_t));
	mydata.localarg[0] = torc_node_id();
	mydata.homenode = torc_node_id();

	for (node = 0; node < nnodes; node++) {
		if (node != torc_node_id())
			send_descriptor(node, &mydata, TERMINATE_WORKER_THREADS);
	}
}

/*************************************************************************/
/**********************     INTERNODE STEALING      **********************/
/*************************************************************************/
pthread_mutex_t sl = PTHREAD_MUTEX_INITIALIZER;

torc_t *direct_synchronous_stealing_request(int target_node)
{
	int vp = target_node;
	torc_t mydata, *work;
	
	if (termination_flag) {
		return NULL;
	}

	pthread_mutex_lock(&sl);

	work = _torc_get_reused_desc();

#if DBG
	printf("[%d] Synchronous stealing request ...\n", torc_node_id()); fflush(0);
#endif
	memset(&mydata, 0, sizeof(mydata));
	mydata.localarg[0] = torc_node_id();
	mydata.homenode = torc_node_id();

	send_descriptor(vp, &mydata, DIRECT_SYNCHRONOUS_STEALING_REQUEST);
	receive_descriptor(vp, work);	work->next = NULL;

	pthread_mutex_unlock(&sl);

	if (work->type == TORC_NO_WORK) {
		usleep(100*1000);
		/*sched_yield();*/
		/*usleep(1*1000);*/
		_torc_put_reused_desc(work);
		return NULL;
	}
	else {
		steal_hits++;
		return work;
	}
}

void torc_disable_stealing ()
{
	torc_t mydata;
	int node, nnodes = torc_num_nodes();
    
#if DBG
	printf("Disabling internode stealing ...\n"); fflush(0);
#endif
	internode_stealing = 0;
	memset(&mydata, 0, sizeof(torc_t));
	mydata.localarg[0] = torc_node_id();
	mydata.homenode = mydata.sourcenode = torc_node_id();

	for (node = 0; node < nnodes; node++) {
		if (node != torc_node_id())
			send_descriptor(node, &mydata, DISABLE_INTERNODE_STEALING);	/* OK. This descriptor is a stack variable */
	}
}



void torc_enable_stealing ()
{
	torc_t mydata;
	int node, nnodes = torc_num_nodes();
    
#if DBG
	printf("Enabling internode stealing ...\n"); fflush(0);
#endif
	internode_stealing = 1;
	memset(&mydata, 0, sizeof(torc_t));
	mydata.localarg[0] = (long) torc_node_id();
	mydata.homenode = mydata.sourcenode = torc_node_id();

	for (node = 0; node < nnodes; node++) {
		if (node != torc_node_id())
			send_descriptor(node, &mydata, ENABLE_INTERNODE_STEALING);	/* OK. This descriptor is a stack variable */
	}
}

void torc_i_enable_stealing ()
{
	internode_stealing = 1;
}

void torc_i_disable_stealing ()
{
	internode_stealing = 0;
}

/*************************************************************************/
/************************     RESET STATISTICS     ***********************/
/*************************************************************************/

void torc_reset_statistics ()
{
	torc_t mydata;
	int node, nnodes = torc_num_nodes();
    
#if DBG
	printf("Reseting statistics ...\n");
#endif
	memset(&mydata, 0, sizeof(torc_t));
	for (node = 0; node < nnodes; node++) {
		if (node != torc_node_id())
			send_descriptor(node, &mydata, RESET_STATISTICS);	/* OK. This descriptor is a stack variable */
		else
			_torc_reset_statistics();
	}
}
