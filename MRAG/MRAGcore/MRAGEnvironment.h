/*
 *  MRAGEnvironment.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/21/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
namespace MRAG
{
	
//OS STUFF
#define _MRAG_OS_WINDOWS _WIN32
#define _MRAG_OS_APPLE __APPLE__

#ifdef _WIN32
#define _MRAG_OS _MRAG_OS_WINDOWS
#else 
#define _MRAG_OS _MRAG_OS_APPLE
#endif

//VIZ STUFF 
//#define _MRAG_GLUT_VIZ

	
//ARCHITECTURE STUFF
#define _MRAG_TBB

#ifdef _MRAG_TBB
		#define _MRAG_BLOCKLAB_ALLOCATOR tbb::cache_aligned_allocator
		#define _MRAG_BLOCKLAB_ALLOCATOR_HEADER "tbb/cache_aligned_allocator.h"
		
		#define _MRAG_GHOSTSCREATION_ALLOCATOR tbb::scalable_allocator
		#define _MRAG_GHOSTSCREATION_ALLOCATOR_HEADER "tbb/scalable_allocator.h"
	
		#define _MRAG_BLOCKCOLLECTION_ALLOCATOR tbb::cache_aligned_allocator
		#define _MRAG_BLOCKCOLLECTION_ALLOCATOR_HEADER "tbb/cache_aligned_allocator.h"
	
		#ifndef _MRAG_TBB_NTHREADS_HINT
			#define _MRAG_TBB_NTHREADS_HINT 6
		#endif
#endif
}


#ifdef _MRAG_TBB
#include <stdio.h>
#include <tbb/task_scheduler_init.h>
#endif
#pragma once
namespace MRAG
{
	//ENVIRONMENT: RT SETUP
	namespace Environment
	{
        /**
         * General setup of MRAG Environment. Should be called before doing stuff.
         */
		inline void setup(int threads=-1)
		{
#ifdef _MRAG_TBB
			static tbb::task_scheduler_init * init = NULL;
			
			if (init == NULL)
			{
				const int nthreads = threads==-1? _MRAG_TBB_NTHREADS_HINT : threads;
				init = new tbb::task_scheduler_init(nthreads);
				//printf("INITIALIZED THREADS=%d (_MRAG_TBB_NTHREADS_HINT is %d)\n", nthreads, _MRAG_TBB_NTHREADS_HINT);
			}
#endif

		}
	}
}
