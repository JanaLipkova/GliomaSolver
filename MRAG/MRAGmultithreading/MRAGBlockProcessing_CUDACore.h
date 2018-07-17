/*
 *  MRAGBlockProcessing_CUDACore.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/16/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
//#include <assert.h>

template <typename E>
struct BlockLab_CUDA
{
	char start[3];
	int nElements;
	short int sliceSize, rowSize;
	E * _cache;
	

	BlockLab_CUDA(const int s[3], const int e[3], E * cache):
		_cache(cache),
		nElements((e[0] - s[0])*(e[1] - s[1])*(e[2] - s[2])),
		sliceSize((e[0] - s[0])*(e[1] - s[1])),
		rowSize(e[0] - s[0])
	{
		start[0] = s[0];
		start[1] = s[1];
		start[2] = s[2];
	}

	BlockLab_CUDA(const char s[3], const char e[3], E * cache):
		_cache(cache),
		nElements((e[0] - s[0])*(e[1] - s[1])*(e[2] - s[2])),
		sliceSize((e[0] - s[0])*(e[1] - s[1])),
		rowSize(e[0] - s[0])
	{
		start[0] = s[0];
		start[1] = s[1];
		start[2] = s[2];
	}
	
	__device__
	inline E& operator()(int ix, int iy, int iz)
	{
		const int index = (ix-start[0]) + (iy-start[1])*rowSize + (iz-start[2])*sliceSize;
		
		return _cache[index];
	}

	__device__ inline int getOffset(int ix, int iy, int iz)
	{
		return (ix-start[0]) + (iy-start[1])*rowSize + (iz-start[2])*sliceSize;
	}
};
/*
struct UniformPartition
{	
	char total_work[3];
	char work_per_thread[3];
	char start[3];
	char threads[3];
	char passes;
	
	UniformPartition()
	{
	}
	
	UniformPartition(const int s[3], const int e[3], const int nThreads[3])
	{
		setup(s, e, nThreads);
	}
	
	void setup(const int s[3], const int e[3], const int nThreads[3])
	{
		threads[0] = nThreads[0];
		threads[1] = nThreads[1];
		threads[2] = nThreads[2];
		
		total_work[0] = max((int)0, e[0] - s[0]);
		total_work[1] = max((int)0, e[1] - s[1]);
		total_work[2] = max((int)0, e[2] - s[2]);
		
		work_per_thread[0] = max(1, total_work[0]/nThreads[0]);
		work_per_thread[1] = max(1, total_work[1]/nThreads[1]);
		work_per_thread[2] = max(1, total_work[2]/nThreads[2]);
		
		start[0] = s[0];
		start[1] = s[1];
		start[2] = s[2];

		const int nPassesY = 1 + total_work[1]/(work_per_thread[1]*nThreads[1]);
		const int nPassesZ = 1 + total_work[2]/(work_per_thread[2]*nThreads[2]);
		
		passes = nPassesZ > 1 ? 8 : (nPassesY > 1 ? 4 : 2);
	}
	
	__device__
	void findThreadIndex( int threadID, int output[3])
	{
		output[0] = threadID % threads[0];
		output[1] = (threadID/threads[0]) % threads[1];
		output[2] = threadID/(threads[0]*threads[1]);
	}
	
	__device__
	void getCoveredInterval(int s[3], int e[3])
	{
		s[0] = start[0];
		s[1] = start[1];
		s[2] = start[2];
		
		e[0] = start[0] + total_work[0];
		e[1] = start[1] + total_work[1];
		e[2] = start[2] + total_work[2];
	}

	void getCoveredInterval_CPU(int s[3], int e[3])
	{
		s[0] = start[0];
		s[1] = start[1];
		s[2] = start[2];
		
		e[0] = start[0] + total_work[0];
		e[1] = start[1] + total_work[1];
		e[2] = start[2] + total_work[2];
	}
	
	__device__
	bool findWorkingInterval(const int thread_index[3], const int passID, int s[3], int e[3])
	{
		if (passID >= passes) return false;
		
		const int pass[3] = {
			passID % 2, (passID/2)%2, (passID/4)%2
		};
		
		const int offset[3] = {
			work_per_thread[0]*(threads[0]*pass[0] + thread_index[0]),
			work_per_thread[1]*(threads[1]*pass[1] + thread_index[1]),
			work_per_thread[2]*(threads[2]*pass[2] + thread_index[2])
		};
		
		s[0] = start[0] + min(offset[0], total_work[0]);
		s[1] = start[1] + min(offset[1], total_work[1]);
		s[2] = start[2] + min(offset[2], total_work[2]);
		
		e[0] = start[0] + min(offset[0]+work_per_thread[0], total_work[0]);
		e[1] = start[1] + min(offset[1]+work_per_thread[1], total_work[1]);
		e[2] = start[2] + min(offset[2]+work_per_thread[2], total_work[2]);
		
		if ((e[0]-s[0])*(e[1]-s[1])*(e[2]-s[2]) > 0) 
			return true;
		else
			return false;
	}	
};
*/