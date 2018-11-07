/*
 *  queues.h
 *  TORC_Lite
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

/*! \file queues.h
    \brief Queues management routines
 */

#ifndef __QUEUES_H__
#define __QUEUES_H__

/* locks.h should be included before queues.h */

/*
 *	Queues
 */

#define QUEUE_NODE_FIELDS(nodetype)			\
	nodetype *next, *prev			


#define QUEUE_DEFINE(nodetype,queue_type) 		\
typedef union _##queue_type {				\
	struct {					\
		_lock_t lock;				\
		nodetype *head;				\
		nodetype *tail;				\
	} q;						\
	char pad[CACHE_LINE_SIZE];			\
} queue_type


#define _queue_head(_q)		(_q)->q.head
#define _queue_tail(_q)		(_q)->q.tail

#define _queue_init(_q) {				\
	(_q)->q.head = NULL;				\
	(_q)->q.tail = NULL;				\
	_lock_init(&(_q)->q.lock);			\
}


#define _queue_destroy(_q) {				\
	_lock_destroy(&(_q)->q.lock);			\
}


#define _enqueue_head(_q, _e) {				\
	_lock_acquire(&(_q)->q.lock);			\
	(_e)->prev = NULL;				\
	(_e)->next = (_q)->q.head;			\
	if ((_q)->q.head == NULL) {			\
		(_q)->q.head = (_q)->q.tail = (_e);	\
	}						\
	else {						\
		(_q)->q.head = (_e);			\
		(_e)->next->prev = (_e);		\
	}						\
	_lock_release(&(_q)->q.lock);			\
}


#define _enqueue_tail(_q, _e) {				\
	_lock_acquire(&(_q)->q.lock);			\
	(_e)->next = NULL;				\
	(_e)->prev = (_q)->q.tail;			\
	if ((_q)->q.tail == NULL)	{		\
		(_q)->q.tail = (_q)->q.head = (_e);	\
	}						\
	else {						\
		(_q)->q.tail->next = (_e);		\
		(_q)->q.tail = (_e);			\
	}						\
	_lock_release(&(_q)->q.lock);			\
}


#define _dequeue(_q, _thread) {					\
	if ((_q)->q.head != NULL) {				\
		_lock_acquire(&(_q)->q.lock);			\
		(*(_thread)) = (_q)->q.head;			\
		if ((_q)->q.head != NULL) {			\
			(_q)->q.head = (_q)->q.head->next;	\
			if ((_q)->q.head == NULL)		\
				(_q)->q.tail = NULL;		\
			else					\
				(_q)->q.head->prev = NULL;	\
		}						\
		_lock_release(&(_q)->q.lock);			\
	}							\
}


#define _dequeue_end(_q, _thread) {				\
	_lock_acquire(&(_q)->q.lock);				\
	(*(_thread)) = (_q)->q.tail;				\
	if ((_q)->q.tail != NULL) {				\
		(_q)->q.tail = (_q)->q.tail->prev;		\
		if ((_q)->q.tail == NULL)			\
			(_q)->q.head = NULL;			\
		else						\
			(_q)->q.tail->next = NULL;		\
	}							\
	_lock_release(&(_q)->q.lock);				\
}


#define _enqueue_head_bare(_q, _e) {			\
	(_e)->prev = NULL;				\
	(_e)->next = (_q)->q.head;			\
	if ((_q)->q.head == NULL) {			\
		(_q)->q.head = (_q)->q.tail = (_e);	\
	}						\
	else {						\
		(_q)->q.head = (_e);			\
		(_e)->next->prev = (_e);		\
	}						\
}


#define _enqueue_tail_bare(_q, _e) {			\
	(_e)->next = NULL;				\
	(_e)->prev = (_q)->q.tail;			\
	if ((_q)->q.tail == NULL)	{		\
		(_q)->q.tail = (_q)->q.head = (_e);	\
	}						\
	else {						\
		(_q)->q.tail->next = (_e);		\
		(_q)->q.tail = (_e);			\
	}						\
}


#define _dequeue_bare(_q, _thread) {					\
	if ((_q)->q.head != NULL) {				\
		(*(_thread)) = (_q)->q.head;			\
		if ((_q)->q.head != NULL) {			\
			(_q)->q.head = (_q)->q.head->next;	\
			if ((_q)->q.head == NULL)		\
				(_q)->q.tail = NULL;		\
			else					\
				(_q)->q.head->prev = NULL;	\
		}						\
	}							\
}


#define _dequeue_end_bare(_q, _thread) {			\
	(*(_thread)) = (_q)->q.tail;				\
	if ((_q)->q.tail != NULL) {				\
		(_q)->q.tail = (_q)->q.tail->prev;		\
		if ((_q)->q.tail == NULL)			\
			(_q)->q.head = NULL;			\
		else						\
			(_q)->q.tail->next = NULL;		\
	}							\
}

#endif
