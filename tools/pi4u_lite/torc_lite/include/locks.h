/*
 *  locks.h
 *  TORC_Lite
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

/*! \file locks.h
    \brief Mutual exclusion routines
 */

#ifndef __LOCKS_H__
#define __LOCKS_H__

#include "ps_config.h"
/*
 *	POSIX synchronization primitves
 */

#if defined(POSIX_SPIN_LOCK) || defined(POSIX_SPIN_TRYLOCK)
#define _XOPEN_SOURCE	600
#endif

#include <sys/types.h>
#include <errno.h>
#include <pthread.h>

#if defined(POSIX_MUTEX_LOCK) || defined(POSIX_MUTEX_TRYLOCK)
typedef pthread_mutex_t _lock_t;
#define LOCK_INITIALIZER		PTHREAD_MUTEX_INITIALIZER

#define _lock_init(var)			pthread_mutex_init(var, NULL)
#define _lock_try_acquire(var)		pthread_mutex_trylock(var)
#define _lock_release(var)		pthread_mutex_unlock(var)
#define _lock_destroy(var)		pthread_mutex_destroy(var)

#if defined(POSIX_MUTEX_LOCK)
#define _lock_acquire(var)		pthread_mutex_lock(var)
#elif defined(POSIX_MUTEX_TRYLOCK)
static int _lock_acquire(_lock_t * lock)
{
	while (pthread_mutex_trylock(lock) == EBUSY) { sched_yield(); /*noop();*/ }
	return 0;
}
#endif

#elif defined(POSIX_SPIN_LOCK) || defined(POSIX_SPIN_TRYLOCK)

typedef pthread_spinlock_t _lock_t;
#define LOCK_INITIALIZER		1

#define _lock_init(var)			pthread_spin_init(var, 0)
#define _lock_try_acquire(var)		pthread_spin_trylock(var)
#define _lock_release(var)		pthread_spin_unlock(var)
#define _lock_destroy(var)		pthread_spin_destroy(var)

#if defined(POSIX_SPIN_LOCK)
#define _lock_acquire(var)		pthread_spin_lock(var)
#elif defined(POSIX_SPIN_TRYLOCK)	/* from ompi */
static int _lock_acquire(_lock_t * lock)
{
    volatile int count, delay, dummy;
    for (delay = 0; (pthread_spin_trylock(lock) == EBUSY); )
    {
      for (count = dummy = 0; count < delay; count++ )
        dummy += count;      /* To avoid compiler optimizations */
      if (delay == 0)              
        delay = 1;                            
      else
        if (delay < 10000)     /* Don't delay too much */
          delay = delay << 1;               
    }
    return 0;
}
#endif

#endif

#endif
