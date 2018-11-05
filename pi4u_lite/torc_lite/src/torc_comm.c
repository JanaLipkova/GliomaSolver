/*
 *  torc_comm.c
 *  TORC_Lite
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */
#include <torc_internal.h>
#include <torc.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//#define MPI_Send	MPI_Ssend

pthread_mutex_t commlock = PTHREAD_MUTEX_INITIALIZER;

void enter_comm_cs()
{
        if (!thread_safe) _lock_acquire(&commlock);
//	_lock_acquire(&commlock);
}

void leave_comm_cs()
{
        if (!thread_safe) _lock_release(&commlock);
//	_lock_release(&commlock);
}

/*************************************************************************/
/**************************   DATA STRUCTURES    *************************/
/*************************************************************************/

#define torc_desc_size	sizeof(torc_t)

MPI_Comm comm_out;
/*#define comm_out MPI_COMM_WORLD*/

struct node_info *node_info;

/*************************************************************************/
/************************  INTER-NODE ROUTINES   *************************/
/*************************************************************************/

int aslr_flag = 0;
static int number_of_functions = 0;
static func_t internode_function_table[32];

void torc_register_task(void *f)
{
	internode_function_table[number_of_functions] = (func_t) f;
	number_of_functions++;
}

void F77_FUNC_(torc_register_task, TORC_REGISTER_TASK)(void *f)
{
//	printf("%s: %p\n", __FUNCTION__, f);
	torc_register_task(f);
        /* nothing to do */
}


int getfuncnum(func_t f)
{
	int i;
	for (i = 0; i < number_of_functions; i++) {
		if (f == internode_function_table[i]) return i;
	}
	return -1;
}

func_t getfuncptr(int pos)
{
/*	if (pos < 32) {*/
		return internode_function_table[pos];
/*	}*/
/*	else {*/
/*		return pos;*/	/* Not registered */
/*	}*/
}

void check_aslr()
{
	unsigned long vaddr[MAX_NODES];
	unsigned long vaddr_me = (unsigned long) check_aslr;
	int me = torc_node_id();
	int nproc = torc_num_nodes();
	
	enter_comm_cs();
	MPI_Allgather (&vaddr_me, 1, MPI_UNSIGNED_LONG,
			vaddr, 1, MPI_UNSIGNED_LONG,
			comm_out);
	leave_comm_cs();

	int i;
#if 1
	if (me == 0) {
		for (i = 0; i < torc_num_nodes(); i++) {
			printf("node %2d -> %p\n", i, (void *)vaddr[i]);
		}
	}
#endif
	for (i = 0; i < torc_num_nodes(); i++) {
		if (vaddr_me != vaddr[i]) {
			aslr_flag = 1;
			return;
		}
	}

	return;
}

int get_aslr()
{
	return aslr_flag;
}

/*************************************************************************/
/**********************   INITIALIZATION ROUTINES   **********************/
/*************************************************************************/

void _torc_comm_pre_init()
{
	static int already_called = -1;

	already_called++;
	if (already_called) {
		printf("_rte_comm_pre_init has been called already\n");
		return;
	}

	node_info =(struct node_info *) calloc(1, MAX_NODES*sizeof(struct node_info));
}


void _torc_comm_init()
{
	int i;
	int workers[MAX_NODES];
	int workers_me;
	int me = torc_node_id();
	int nproc = torc_num_nodes();
	
	workers_me = kthreads;

	check_aslr();

	enter_comm_cs();
	MPI_Allgather (&workers_me, 1, MPI_INT,
				workers, 1, MPI_INT,
				comm_out);

	MPI_Barrier(comm_out);
	leave_comm_cs();
	/* Workers have been started. The node_info array must be combined by all nodes */
	
	for (i = 0; i < torc_num_nodes(); i++) {
		/*node_info[i].nworkers = kthreads;*/	/* SMP */
		node_info[i].nworkers = workers[i];	/* SMP */
		/*printf("ni[%d].nworkes = %d\n", i, node_info[i].nworkers);*/
	}
	enter_comm_cs();
	MPI_Barrier(comm_out);
	leave_comm_cs();
	
#if DBG
	printf("[%d/%d] Node is up\n", torc_node_id(), torc_num_nodes()); fflush(0);

	/* Synchronize execution of workers */
	enter_comm_cs();
	MPI_Barrier(comm_out);
	leave_comm_cs();
#endif
}


/*************************************************************************/
/******   EXPICLIT COMMUNICATION FOR DESCRIPTORS (SEND, RECEIVE)    ******/
/*************************************************************************/

static int _torc_thread_id()
{
	if (pthread_equal(pthread_self(), server_thread))
		return MAX_NVPS;
	else
		return _torc_get_vpid();
}

void send_arguments(int node, int tag, torc_t *desc)
{
	int i;
	
	for (i = 0; i < desc->narg; i++)	{
/*		printf("ARG %d: CALL %d, QUANT %d\n", i, desc->callway[i], desc->quantity[i]); fflush(0);*/
		if (desc->quantity[i] == 0) continue;
		if ((desc->callway[i] == CALL_BY_COP)|| (desc->callway[i] == CALL_BY_VAD)) {
			if (desc->quantity[i] == 1) continue;	/* do not send anything - the value is in the descriptor */
			enter_comm_cs();
			if (desc->homenode != desc->sourcenode)
				MPI_Send(&desc->temparg[i],desc->quantity[i], desc->dtype[i],node,tag,comm_out);
			else
				MPI_Send(&desc->localarg[i],desc->quantity[i], desc->dtype[i],node,tag,comm_out);
			leave_comm_cs();
		}
		else if ((desc->callway[i] == CALL_BY_REF)||(desc->callway[i] == CALL_BY_PTR)||(desc->callway[i] == CALL_BY_COP2)) {
			enter_comm_cs();
			if (desc->homenode != desc->sourcenode) {
				MPI_Send((void *)desc->temparg[i],desc->quantity[i],desc->dtype[i],node,tag,comm_out);
			} else {
				MPI_Send((void *)desc->localarg[i],desc->quantity[i],desc->dtype[i],node,tag,comm_out);
			}
			leave_comm_cs();
		}
		else /* CALL_BY_RES */
			/*nothing*/;
	}
}

/* Send a descriptor to the target node */
void send_descriptor(int node, torc_t *desc, int type)	/* always to a server thread */
{
	int tag = _torc_thread_id();
	int i;

#if DBG
	printf("[%d] - sending to node [%d] desc -> homenode [%d], type = %d\n", torc_node_id(), node, desc->homenode, type);
#endif
	desc->sourcenode = torc_node_id();
	desc->sourcevpid = tag;	/* who sends this */
	desc->type = type;

	enter_comm_cs();
	MPI_Send(desc, torc_desc_size, MPI_CHAR,node,MAX_NVPS,comm_out);
	leave_comm_cs();
	switch (desc->type) {
		case DIRECT_SYNCHRONOUS_STEALING_REQUEST:
		case TORC_BCAST:
			return;
		case TORC_ANSWER:
			/* in case of call by reference send the data back */
			for (i = 0; i < desc->narg; i++) {
				if (desc->quantity[i] == 0) continue;
				if ((desc->callway[i] == CALL_BY_COP2) && (desc->quantity[i] > 1)) {
					free((void *)desc->temparg[i]);
				}
				if ((desc->callway[i] == CALL_BY_REF) || (desc->callway[i] == CALL_BY_RES)) {	/* send the result back */
					enter_comm_cs();
					MPI_Send((void *)desc->temparg[i],desc->quantity[i],desc->dtype[i],desc->homenode,tag,comm_out);
					leave_comm_cs();
					if (desc->quantity[i] > 1)
						free((void *)desc->temparg[i]);
				}
			}
			return;
			break;
		default:	/* TORC_NORMAL_ENQUEUE */
			if (desc->homenode == node) return;
			send_arguments(node, tag, desc);
			return;
			break;
	}
}

void direct_send_descriptor(int dummy, int sourcenode, int sourcevpid, torc_t *desc)
{
	int tag;
	desc->sourcenode = torc_node_id();
	desc->sourcevpid = MAX_NVPS;	/* the server thread responds to a stealing request from a worker*/

	if (sourcevpid == MAX_NVPS) {
		//printf("direct send to (%d, %d)\n", sourcenode, sourcevpid); fflush(0);
		sourcevpid = MAX_NVPS + 1;	/* response to server's thread direct stealing request */ 
	}

	tag = sourcevpid + 100; 

	enter_comm_cs();
	MPI_Send(desc, torc_desc_size, MPI_CHAR, sourcenode, tag, comm_out);	/* if sourcevpid == MAX_NVPS ... */
//	MPI_Ssend(desc, torc_desc_size, MPI_CHAR, sourcenode, tag, comm_out);	/* if sourcevpid == MAX_NVPS ... */
	leave_comm_cs();
	if (desc->homenode == sourcenode) {
			return;
	}

//	send_arguments(sourcenode, sourcevpid, desc);
	send_arguments(sourcenode, tag, desc);
	
/*	free(desc);*/
}

void receive_arguments(torc_t *work, int tag)
{
	int i;
	int typesize;
	char *mem;
	MPI_Status status;
	int istat;
	
	for (i = 0; i < work->narg; i++) {
#if DBG
		printf("reading arg %d (%d - %d)\n", i, work->quantity[i], work->callway[i]); fflush(0);
#endif
		if (work->quantity[i] == 0) continue;
		if ((work->quantity[i] > 1)||((work->callway[i] != CALL_BY_COP)&&(work->callway[i] != CALL_BY_VAD))) {

#if 1
			/* yyy */
			work->dtype[i] = _torc_b2mpi_type(work->btype[i]);
#endif

#if 0   /* Why? Brutus and OpenMPI */
			if (work->dtype[i] == MPI_DOUBLE_PRECISION) {
				work->dtype[i] = MPI_DOUBLE_PRECISION;
			} else if (work->dtype[i] == MPI_INTEGER) {
				work->dtype[i] = MPI_INTEGER;
			} else {
				printf("complete code for work->dtype[%d]=0x%lx\n", i, work->dtype[i]);
			}
#endif

			MPI_Type_size(work->dtype[i], &typesize);
			mem = (char *)calloc(1, work->quantity[i]*typesize);
//			printf("mem = %p\n", mem);fflush(0);
			work->temparg[i] = (INT64)mem;
			if ((work->callway[i] != CALL_BY_RES)) {	/* CALL_BY_REF etc */
//xxx				if (work->quantity[i] == 1) { printf("WHY (%d)???\n", i); fflush(0); }
				enter_comm_cs();
				istat = MPI_Recv((void *)work->temparg[i], work->quantity[i], work->dtype[i], work->sourcenode, tag, comm_out, &status);
				leave_comm_cs();
			} 
			else if (work->callway[i] == CALL_BY_RES) {
//				double setval = 0;
//				memsetvalue((void *)work->temparg[i], setval, work->quantity[i], work->dtype[i]);
			}
		}	
		else {
			work->temparg[i] = work->localarg[i];
		}
#if DBG
		printf("read+++ arg %d (%d - %d)\n", i, work->quantity[i], work->callway[i]); fflush(0);
#endif

	}
}

int receive_descriptor(int node, torc_t *rte)
{
	MPI_Status status;
	int istat;
	int tag = _torc_thread_id();

	tag = tag + 100;

	if (thread_safe) {
		enter_comm_cs();
		istat = MPI_Recv(rte, torc_desc_size, MPI_CHAR, node, tag, comm_out, &status);
		leave_comm_cs();
        } else {
		/* irecv for non thread-safe MPI libraries */
                int flag = 0;
                MPI_Request request;

                enter_comm_cs();
                istat = MPI_Irecv(rte, torc_desc_size, MPI_CHAR, node, tag, comm_out, &request);
                leave_comm_cs();
                while (1) {
			if (appl_finished == 1) { rte->type = TORC_NO_WORK; return 1;}	//pthread_exit(0);
			enter_comm_cs();
			istat = MPI_Test(&request, &flag, &status);
			leave_comm_cs();
			if (flag == 1) {        /* check of istat ? */
				//if (istat == MPI_SUCCESS)
				break;
			} else {
				thread_sleep(yieldtime);
			}
		}                
	}

	if (istat == MPI_SUCCESS)
		/*return 0*/;
	else {
		rte->type = TORC_NO_WORK;
//		printf("MPI_Recv Failed (%d)!!!\n", istat); fflush(0);
		return 1;
	}

	if (rte->type == TORC_NO_WORK) return 1;

	if (rte->homenode == torc_node_id()) {	/* the descriptor is stolen by its owner node */
	}
	else { /* homenode != torc_node_id() */
		receive_arguments(rte, tag);
	}
	return 0;
}

/*************************************************************************/
/**************************  THREAD MANAGEMENT  **************************/
/*************************************************************************/

static int _torc_num_threads(int node_id)
{
	return node_info[node_id].nworkers;	/* kthreads */
}

int _torc_total_num_threads()
{
	int i;
	int nodes = torc_num_nodes();
	int sum_vp = 0;

	for (i = 0; i < nodes; i++) {
		sum_vp += _torc_num_threads(i);
	}
	return sum_vp;
}

int global_thread_id_to_node_id(int global_thread_id)
{
	int i;
	int nodes = torc_num_nodes();
	int sum_vp = 0;

#if DBG
	printf("global_thread_id = %d\n", global_thread_id);
#endif

	for (i = 0; i < nodes; i++) {
		sum_vp += _torc_num_threads(i);
		if (global_thread_id < sum_vp) return i;
	}
	Error("target_to_node failed");
	return -1;	/* never reached */
}

int local_thread_id_to_global_thread_id(int local_thread_id)
{
	int i;
	int mynode = torc_node_id();
	int sum_vp = 0;

	for (i = 0; i < mynode; i++) {
		sum_vp += _torc_num_threads(i);
	}
	return sum_vp + local_thread_id;
}


int global_thread_id_to_local_thread_id(int global_thread_id)
{
	int i;
	int mynode = global_thread_id_to_node_id(global_thread_id);
	int sum_vp = 0;

	for (i = 0; i < mynode; i++) {
		sum_vp += _torc_num_threads(i);
	}
	return global_thread_id - sum_vp;
}

/*************************************************************************/
/**************************    BROADCASTING     **************************/
/*************************************************************************/

void torc_broadcast(void *a, long count, MPI_Datatype datatype)
{
	long mynode = torc_node_id();

	torc_t mydata;                
	int node, nnodes = torc_num_nodes();
	int tag = _torc_thread_id();

#if DBG
	printf("Broadcasting data ...\n"); fflush(0);
#endif
	memset(&mydata, 0, sizeof(torc_t));
	mydata.localarg[0] = torc_node_id();
	mydata.localarg[1] = a; 
	mydata.localarg[2] = count;
#if 1
	/* yyy */
	mydata.localarg[3] = _torc_mpi2b_type(datatype);
//	mydata.localarg[3] = datatype;                        
#endif
	mydata.homenode = mydata.sourcenode = torc_node_id();

	for (node = 0; node < nnodes; node++) {
		if (node != torc_node_id()) {
			send_descriptor(node, &mydata, TORC_BCAST);     /* OK. This descriptor is a stack variable */
			enter_comm_cs();
			MPI_Ssend(a,count,datatype,node,tag,comm_out);
			leave_comm_cs();
		}
	}
}

/* C types */
#if 0
#define T_MPI_CHAR			0
#define T_MPI_SIGNED_CHAR		1
#define T_MPI_UNSIGNED_CHAR		2
#define T_MPI_BYTE			3
#define T_MPI_WCHAR			4
#define T_MPI_SHORT			5
#define T_MPI_UNSIGNED_SHORT		6
#define T_MPI_INT			7
#define T_MPI_UNSIGNED			8
#define T_MPI_LONG			9
#define T_MPI_UNSIGNED_LONG		10
#define T_MPI_FLOAT			11
#define T_MPI_DOUBLE			12
#define T_MPI_LONG_DOUBLE		13
#define T_MPI_LONG_LONG_INT		14
#define T_MPI_UNSIGNED_LONG_LONG	15
#define T_MPI_LONG_LONG			14	//MPI_LONG_LONG_INT

/* Fortran types */
#define T_MPI_COMPLEX			15
#define T_MPI_DOUBLE_COMPLEX		16
#define T_MPI_LOGICAL			17
#define T_MPI_REAL			18
#define T_MPI_DOUBLE_PRECISION		19
#define T_MPI_INTEGER			20
#define T_MPI_2INTEGER			21

#else

enum {
/* C types */
T_MPI_CHAR = 0,
T_MPI_SIGNED_CHAR,
T_MPI_UNSIGNED_CHAR,
T_MPI_BYTE,
T_MPI_WCHAR,
T_MPI_SHORT,
T_MPI_UNSIGNED_SHORT,
T_MPI_INT,
T_MPI_UNSIGNED,
T_MPI_LONG,
T_MPI_UNSIGNED_LONG,
T_MPI_FLOAT,
T_MPI_DOUBLE,
T_MPI_LONG_DOUBLE,
T_MPI_LONG_LONG_INT,
T_MPI_UNSIGNED_LONG_LONG,

/* Fortran types */
T_MPI_COMPLEX,
T_MPI_DOUBLE_COMPLEX,
T_MPI_LOGICAL,
T_MPI_REAL,
T_MPI_DOUBLE_PRECISION,
T_MPI_INTEGER,
T_MPI_2INTEGER,

/* C types */
T_MPI_LONG_LONG = T_MPI_LONG_LONG_INT,
};

#endif

int _torc_mpi2b_type(MPI_Datatype dtype)
{
	if	(dtype == MPI_INT)		return T_MPI_INT;
	else if	(dtype == MPI_LONG)		return T_MPI_LONG;
	else if	(dtype == MPI_FLOAT)		return T_MPI_FLOAT;
	else if	(dtype == MPI_DOUBLE)		return T_MPI_DOUBLE;
	else if	(dtype == MPI_DOUBLE_PRECISION)	return T_MPI_DOUBLE_PRECISION;
	else if	(dtype == MPI_INTEGER)		return T_MPI_INTEGER;
	else 					Error("unsupported MPI data type");

        return 0;       /* never reached */
}

MPI_Datatype _torc_b2mpi_type(int btype)
{
	switch (btype)
	{
	case T_MPI_INT:			return MPI_INT; break;
	case T_MPI_LONG:		return MPI_LONG; break;
	case T_MPI_FLOAT:		return MPI_FLOAT; break;
	case T_MPI_DOUBLE:		return MPI_DOUBLE; break;
	case T_MPI_DOUBLE_PRECISION:	return MPI_DOUBLE_PRECISION; break;
	case T_MPI_INTEGER:		return MPI_INTEGER; break;
        default: 			printf("btype = %d\n", btype);
					Error("unsupported MPI datatype"); break;
        }

        return 0;       /* never reached */
}
