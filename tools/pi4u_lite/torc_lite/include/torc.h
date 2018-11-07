/*
 *  torc.h
 *  TORC_Lite
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#ifndef _torc_included
#define _torc_included
#include <pthread.h>
#include <mpi.h>
#define TORC_LITE	1
/* Exported interface */

#ifdef __cplusplus
extern "C"
{
#endif

#define MODE_MW					0
#define MODE_MS					0
void torc_init (int argc, char *argv[], int ms);
void torc_reset_statistics();
typedef double torc_time_t;
torc_time_t torc_gettime();

int torc_i_worker_id(void);
int torc_i_num_workers();
int torc_worker_id();
int torc_num_workers();
int torc_getlevel();


/******  Exported Interface *******/
#define CALL_BY_COP				(int)(0x0001)	/* IN	- By copy, through pointer to private copy (C) */
#define CALL_BY_REF				(int)(0x0002)	/* INOUT- By reference */
#define CALL_BY_RES				(int)(0x0003)	/* OUT	- By result */
#define CALL_BY_PTR				(int)(0x0004)	/* IN	- By value, from address */
#define CALL_BY_VAL				(int)(0x0001)	/* IN	- By value, from address (4: C, 0: Fortran */
#define CALL_BY_COP2				(int)(0x0005)	/* IN	- By copy, through pointer to private copy (C) */
#define CALL_BY_VAD				(int)(0x0006)   /* IN   - By address - For Fortran Routines (Fortran) */

void torc_enable_stealing();
void torc_disable_stealing();
void torc_i_enable_stealing();
void torc_i_disable_stealing();

void torc_taskinit();
void torc_waitall();
void torc_waitall2();
void torc_waitall3();
void torc_tasksync();

//#ifndef __cplusplus
void torc_task (int queue, void (*f) (), int narg, ...);
void torc_task_detached (int queue, void (*f) (), int narg, ...);
void torc_task_ex (int queue, int invisible, void (*f) (), int narg, ...);
void torc_task_direct (int queue, void (*f) (), int narg, ...);
//#else
//void torc_task (int queue, void *f, int narg, ...);
//void torc_task_ex (int queue, int invisible, void *f, int narg, ...);
//void torc_task_direct (int queue, void *f, int narg, ...);
//#endif

#define torc_create		torc_task
#define torc_create_detached	torc_task_detached
#define torc_create_ex		torc_task_ex
#define torc_create_direct	torc_task_direct

int torc_node_id();
int torc_num_nodes();

void torc_broadcast(void *a, long count, MPI_Datatype dtype);
void torc_broadcast_ox(void *a, long count, int dtype);

void thread_sleep(int ms);
void torc_finalize(void);

void torc_register_task(void *f);


int torc_fetch_work();	// for dr

#ifdef __cplusplus
}
#endif

#endif

