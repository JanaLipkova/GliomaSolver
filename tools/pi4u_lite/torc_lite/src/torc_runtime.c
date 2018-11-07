/*
 *  torc_runtime.c
 *  TORC_Lite
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */
#include <torc_internal.h>
#include <torc.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#define WAIT_COUNT 2
#define YIELD_TIME	(0)

#define A0	   
#define A1	    args[0]
#define A2	 A1,args[1]
#define A3	 A2,args[2]
#define A4	 A3,args[3]
#define A5	 A4,args[4]
#define A6	 A5,args[5]
#define A7	 A6,args[6]
#define A8	 A7,args[7]
#define A9	 A8,args[8]
#define A10	 A9,args[9]
#define A11	 A10,args[10]
#define A12	 A11,args[11]
#define A13	 A12,args[12]
#define A14	 A13,args[13]
#define A15	 A14,args[14]
#define A16	 A15,args[15]
#define A17	 A16,args[16]

#define DO_CASE(x)			\
	case x:				\
	{				\
		rte->work(A##x);	\
		break;			\
	}

/*static*/ void _torc_core_execution (torc_t *rte)
{
	VIRT_ADDR args[MAX_TORC_ARGS];
	int i;

	if (rte->work_id >= 0)
		rte->work = getfuncptr((long)rte->work_id);

	if (torc_node_id() == rte->homenode) {
		for (i=0; i<rte->narg; i++) {
			if (rte->callway[i] == CALL_BY_COP) {
				args[i] = (VIRT_ADDR) &rte->localarg[i];      /* pointer to the private copy */
			}
			else {
				args[i] = rte->localarg[i];
			}
		}
	}
	else {
		for (i=0; i<rte->narg; i++) {
			if (rte->callway[i] == CALL_BY_COP) {
				args[i] = (VIRT_ADDR) &rte->temparg[i];
			}
			else {
				args[i] = rte->temparg[i];
			}
		}
	}

        switch (rte->narg){
        DO_CASE(0);
        DO_CASE(1);
        DO_CASE(2);
        DO_CASE(3);
        DO_CASE(4);
        DO_CASE(5);
        DO_CASE(6);
        DO_CASE(7);
        DO_CASE(8);
        DO_CASE(9);
        DO_CASE(10);
        DO_CASE(11);
        DO_CASE(12);
        DO_CASE(13);
        DO_CASE(14);
        DO_CASE(15);
        DO_CASE(16);
        DO_CASE(17);
        default:
                Error("rte function with more than 17 arguments..!");
                break;
        }

}

void _torc_reset_statistics()
{
	memset(created, 0, MAX_NVPS*sizeof(unsigned long));
	memset(executed, 0, MAX_NVPS*sizeof(unsigned long));
}

void _torc_print_statistics()
{
        unsigned int i;
        unsigned long total_created = 0, total_executed = 0;

        /* Runtime statistics */
        for (i = 0; i < kthreads; i++) {
                total_created += created[i];
                total_executed += executed[i];
        }

        printf("[%2d] steals served/attempts/hits = %-3ld/%-3ld/%-3ld created = %3ld, executed = %3ld:(", torc_node_id(), 
                steal_served, steal_attempts, steal_hits, total_created, total_executed);
        for (i = 0; i < kthreads-1; i++) {
                printf("%3ld,", executed[i]); 
        }           
        printf("%3ld)\n", executed[i]); fflush(0);
}

void _torc_stats (void)
{
#if defined(TORC_STATS)
	_torc_print_statistics();
#endif
}

torc_t *_torc_self()
{
	return (torc_t *) _torc_get_currt();
}

void _torc_depadd (torc_t * rte, int ndeps)
{
	_lock_acquire(&rte->lock);
	rte->ndep += ndeps;
	_lock_release(&rte->lock);
}

static void _torc_end (void)
{
	int finalized = 0;

#if 1
	MPI_Finalized(&finalized);
	if (finalized == 1) {
		_torc_stats();
		_exit(0);
	}
#endif
	
	appl_finished=1;
	if (torc_num_nodes() > 1) {    /* notify the rest of the nodes */
		terminate_workers();
	}

	_torc_md_end ();
}

void torc_finalize() { _torc_end();}

void _torc_execute (void * arg)
{
	torc_t * me = _torc_self();
	torc_t * rte = (torc_t *) arg;
	long vp = me->vp_id;

#ifdef TORC_STATS
	if (rte->rte_type == 1) {
		executed[vp]++;
	}
#endif

	rte->vp_id = me->vp_id;
	_torc_set_currt(rte);
	_torc_core_execution(rte);
	_torc_cleanup(rte);
	_torc_set_currt(me);
}

int _torc_block (void)
{
	torc_t *rte = _torc_self();
	int self = rte->vp_id;
	int remdeps;

	_lock_acquire (&rte->lock);
	--rte->ndep;
	if (rte->ndep < 0) rte->ndep = 0;
	_lock_release (&rte->lock);
	while (1)
	{
		_lock_acquire (&rte->lock);
		remdeps = rte->ndep;
		if (remdeps<=0) {
			_lock_release (&rte->lock);
			return 1;
		}
		_lock_release (&rte->lock);
		_torc_scheduler_loop(1);
	}
	return 0;
}

/* Q: What did I do here? - A: Block until no more work exists at the cluster-layer. Useful for SPMD-like barriers */ 
int _torc_block2 (void)
{
	torc_t *rte = _torc_self();
	int self = rte->vp_id;
	int remdeps;

	_lock_acquire (&rte->lock);
	--rte->ndep;
	if (rte->ndep < 0) rte->ndep = 0;
	_lock_release (&rte->lock);
	while (1)
	{
		if (rte->ndep > 0) {
			_lock_acquire (&rte->lock);
			remdeps = rte->ndep;
			if (remdeps<=0) {
				rte->ndep = 0; 
				_lock_release (&rte->lock);
				return 1;
			}
			_lock_release (&rte->lock);
		}
		int work = _torc_scheduler_loop(1);
		if ((rte->ndep == 0) && (!work))
			return 0;
	}
	return 0;
}

void _torc_set_work_routine(torc_t *rte, void (*work)())
{
	rte->work = work;
	rte->work_id = getfuncnum(work);
	if (-1 == rte->work_id) {
#if 0
		if (get_aslr()) {
			Error1("Internode function %p not registered", work);
		}
#else
		printf("Internode function %p not registered\n", work);
#endif
	}
}

int _torc_depsatisfy (torc_t * rte)
{
	int deps;

	_lock_acquire (&rte->lock);
	deps = --rte->ndep;
	_lock_release (&rte->lock);
	return !deps;
}

extern MPI_Comm comm_out;

void _torc_opt (int argc, char *argv[])
{
	int largc = argc;
	char **largv = argv;
	char *llargv = "";
	int requested, provided;
	int initialized = 0;
	char *s;
	int val;
//	int thread_safe;

//	printf("TORC_LITE...\n"); 
#if 1
	requested = MPI_THREAD_MULTIPLE;
#else
	requested = MPI_THREAD_SINGLE;
#endif

	if (argc == 0) largv = (char **) &llargv;	/* in case argv cannot be NULL (HPMPI) */
	
	MPI_Initialized(&initialized);
	if (initialized == 0) {
		MPI_Init_thread(&largc, &largv, requested, &provided);
	}
	else {
		MPI_Query_thread(&provided);
	}

/*	printf("Thread safety: requested = %d provided = %d\n", requested, provided);*/	/*MPI_THREAD_MULTIPLE 3*/
	if (provided != MPI_THREAD_MULTIPLE)
		thread_safe = 0;
	else
		thread_safe = 1;

	kthreads = TORC_DEF_CPUS;
	s = (char *) getenv("OMP_NUM_THREADS");
	if (s != 0 && sscanf(s, "%d", &val) == 1 && val > 0)
		kthreads = val;

	s = (char *) getenv("TORC_WORKERS");
	if (s != 0 && sscanf(s, "%d", &val) == 1 && val > 0)
		kthreads = val;

	yieldtime = TORC_DEF_YIELDTIME;
	s = (char *) getenv("TORC_YIELDTIME");
	if (s != 0 && sscanf(s, "%d", &val) == 1 && val > 0)
		yieldtime = val;

	throttling_factor = -1;
	s = (char *) getenv("TORC_THROTTLING_FACTOR");
	if (s != 0 && sscanf(s, "%d", &val) == 1 && val > 0)
		throttling_factor = val;

	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_nodes);

	MPI_Barrier(MPI_COMM_WORLD);

#if 1
	int size, namelen;
	char name[MPI_MAX_PROCESSOR_NAME];
	MPI_Get_processor_name (name, &namelen);
	printf("TORC_LITE ... rank %d of %d on host %s\n", mpi_rank, mpi_nodes, name);         
#endif
	if (mpi_rank == 0) {
		/*printf("Thread safety: requested = %d provided = %d\n", requested, provided);*/	/*MPI_THREAD_MULTIPLE 3*/
		printf("The MPI implementation IS%s thread safe!\n", (thread_safe)? "": " NOT"); 
	}

	MPI_Comm_dup(MPI_COMM_WORLD, &comm_out);
	MPI_Barrier(MPI_COMM_WORLD);

	_torc_comm_pre_init();
}

/* Package initialization */
void _torc_env_init(void)
{
	rq_init ();			/* ready queues */

	_torc_md_init();	/* workers */
	_torc_comm_init();
}


torc_t *get_next_task()
{
	int self_node = torc_node_id();
	int nnodes = torc_num_nodes();
	torc_t *rte_next = NULL;
	int vp;
	int node;
	int self = torc_i_worker_id();
	
	rte_next = torc_i_pq_dequeue ();
	if (rte_next == NULL) {
		int i;
		for (i = 9; i >= 0; i--) {
			if (rte_next == NULL) rte_next = torc_i_rq_dequeue(i);
			else break;
		}

		if (internode_stealing) {
			node = (self_node + 1) % nnodes;
			while ((rte_next == NULL) && (node != self_node)) {
				rte_next = direct_synchronous_stealing_request(node);
				if (rte_next != NULL) {
					break;
				}
				node = (node + 1) % nnodes;
			}
			if (rte_next == NULL) internode_stealing = 0;
		}
	}

	return rte_next;
}


int torc_fetch_work()
{
	torc_t *task = NULL;

	task = get_next_task();
	if (task != NULL) {
		torc_to_i_pq_end(task);
		return 1;
	}
	else {
		return 0;
	}
}

void _torc_cleanup(torc_t *rte)
{
	if (rte->homenode != torc_node_id()) {
#if DBG
		printf("[%d] sending an answer to %d\n", torc_node_id(), rte->homenode);
#endif
		send_descriptor(rte->homenode, rte, TORC_ANSWER);
	}
	else {
#if DBG
		printf("[%d] satisfying deps on local inter-node desc\n", torc_node_id()); fflush(0);
#endif
		int i;
		for (i = 0; i < rte->narg; i++) {
			if ((rte->callway[i] == CALL_BY_COP2) && (rte->quantity[i]>1)) 
				if ((void *)rte->localarg[i] != NULL)
					free((void *)rte->localarg[i]);
		}

		if (rte->parent) _torc_depsatisfy(rte->parent);
	}
	_torc_put_reused_desc(rte);
}

int _torc_scheduler_loop (int once)
{
	int wait_count;
	int self = torc_i_worker_id();
	torc_t * rte_next;

	while (1) {
		rte_next = get_next_task();

		while (rte_next==NULL) {
			if (appl_finished == 1) {	/* Checking for program completion */
				_torc_md_end();
			}

			thread_sleep(yieldtime);

			rte_next = get_next_task();
			if (rte_next == NULL) {
				if (once) return 0;
				thread_sleep(yieldtime);
			}
		}

		/* Execute selected task */
		_torc_execute(rte_next);
		if (once) return 1;
	}
}


/* WTH is this? Maybe just a backup of the code? */
int _torc_scheduler_loop2 (int once)
{
	int wait_count;
	int self = torc_i_worker_id();
	torc_t * rte_next;

	while (1) {
		rte_next = get_next_task();

		while (rte_next==NULL) {
			if (appl_finished == 1) _torc_md_end();	/* Checking for program completion */

			wait_count = WAIT_COUNT;
			while (--wait_count) {
				//usleep(100*1000);
				sched_yield();
			}

			rte_next = get_next_task();
			if (rte_next == NULL)
				if (once) return 0;
		}

		/* Execute selected task */
		_torc_execute(rte_next);
		if (once) return 0;	// what is this?
	}
}


