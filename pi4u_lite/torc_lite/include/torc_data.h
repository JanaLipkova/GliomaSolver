/*
 *  torc_data.h
 *  TORC_Lite
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#ifndef _torc_data_included
#define _torc_data_included

#include "utils.h"

struct torc_data {
	/* write once - change rarely */
	unsigned int		_global_vps;
	unsigned int		_kthreads;

	unsigned int		_physcpus;
	int			_mpi_rank;
	int			_mpi_nodes;
	unsigned int		_appl_finished;

	int			_thread_safe;
	unsigned int		_internode_stealing;
	int                     _yieldtime;
	int                     _throttling_factor;

	pthread_t		_server_thread;
	pthread_t		_worker_thread[MAX_NVPS];

	/* read write */
	queue_t			_reuse_q;				
	queue_t			_private_grq;                      
	queue_t			_public_grq[10];                      

	unsigned long		_created[MAX_NVPS];
	unsigned long		_executed[MAX_NVPS];
	unsigned long		_steal_hits;
	unsigned long		_steal_served;
	unsigned long		_steal_attempts;

	pthread_key_t		_vp_key;
	pthread_key_t		_currt_key;
};

extern struct torc_data *torc_data;

#define global_vps		torc_data->_global_vps
#define kthreads		torc_data->_kthreads

#define physcpus		torc_data->_physcpus
#define mpi_rank		torc_data->_mpi_rank
#define mpi_nodes		torc_data->_mpi_nodes
#define appl_finished		torc_data->_appl_finished

#define thread_safe		torc_data->_thread_safe
#define internode_stealing	torc_data->_internode_stealing
#define yieldtime		torc_data->_yieldtime
#define throttling_factor	torc_data->_throttling_factor

#define server_thread		torc_data->_server_thread
#define worker_thread		torc_data->_worker_thread

#define reuse_q			torc_data->_reuse_q
#define private_grq		torc_data->_private_grq
#define public_grq		torc_data->_public_grq

#define created			torc_data->_created
#define executed		torc_data->_executed
#define steal_hits		torc_data->_steal_hits
#define steal_served		torc_data->_steal_served
#define steal_attempts		torc_data->_steal_attempts
#define vp_key			torc_data->_vp_key
#define currt_key		torc_data->_currt_key

#endif




