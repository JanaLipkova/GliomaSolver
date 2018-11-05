/*
 *  torc_queue.h
 *  TORC_Lite
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#ifndef _torc_queue_included
#define _torc_queue_included

/* Data structure */
QUEUE_DEFINE(torc_t,queue_t);

void rq_init ();

void torc_to_i_pq (torc_t * rte);
void torc_to_i_pq_end (torc_t * rte);
torc_t * torc_i_pq_dequeue ();

void torc_to_i_rq (torc_t * rte);
void torc_to_i_rq_end (torc_t * rte);
torc_t * torc_i_rq_dequeue (int lvl);

void torc_to_i_lrq (int which, torc_t * rte);
void torc_to_i_lrq_end (int which, torc_t * rte);
torc_t * torc_i_lrq_dequeue (int which);
torc_t * torc_i_lrq_dequeue_end (int which);
torc_t * torc_i_lrq_dequeue_inner (int which);
torc_t * torc_i_lrq_dequeue_end_inner (int which);

void _torc_put_reused_desc (torc_t * rte);
torc_t * _torc_get_reused_desc ();


void torc_to_nrq (int node, torc_t *rte);
void torc_to_nrq_end (int node, torc_t *rte);

void torc_to_rq (torc_t *rte);
void torc_to_rq_end (torc_t *rte);

void torc_to_lrq (int worker, torc_t *rte);
void torc_to_lrq_end (int worker, torc_t *rte);

void torc_to_npq (int node, torc_t *rte);
void torc_to_npq_end (int node, torc_t *rte);




#endif




