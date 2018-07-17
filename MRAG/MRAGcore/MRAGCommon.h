/*
 *  MRAGCommon.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 4/23/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#pragma once

#ifndef NULL
/**
 * Define NULL if necessary.
 */
#define NULL 0
#endif

#ifndef _CUDA_SIDE
#ifndef max
/**
 * Define preprocessor max.
 */
// Gery-note: why not inline?
#define max(a,b) (((a)>=(b))? (a) : (b))
#endif
#endif

#ifndef min
/**
 * Define preprocessor min.
 */
// Gery-note: why not inline?
#define min(a,b) (((a)<(b))? (a) : (b))
#endif

/**
 * Define a default representation for real numbers.
 */

#ifndef _MRAG_REAL_
#define _MRAG_REAL_ float
#endif

typedef _MRAG_REAL_
	Real;


#pragma once

/**
 * Shall we perform debug checks?
 */
const bool cDoDebugChecks = true;

/**
 * Returns true if A => B (A implies B).
 * Useful for assertions.
 */
inline bool _implies(bool A, bool B) { return !A || B;}

#pragma mark struct I3
/**
 * An integer 3-vector.
 */
struct I3
{
	int i[3];
	
	I3(int x, int y, int z)
	{
		i[0] = x;
		i[1] = y;
		i[2] = z;
	}
	
	I3(const I3& a)
	{
		i[0] = a.i[0];
		i[1] = a.i[1];
		i[2] = a.i[2];
	}
	
	I3(){ i[0] = i[1] = i[2] = 0;}
	
	const bool operator<(const I3& a) const
	{
		return (i[0]<a.i[0] || i[0]==a.i[0] && i[1]<a.i[1] ||  i[0]==a.i[0] && i[1]==a.i[1] && i[2]<a.i[2]);
	}
};

#pragma mark struct SimpleInterval
/**
 * A simple representation of an interval with integer begin and end.
 */
struct SimpleInterval
{
	int mbegin;
	int mend;
	
	SimpleInterval(int b, int e): mbegin(b), mend(e) {}
	
	int begin() const { return mbegin;}
	int end() const { return mend;}
};

#pragma mark struct I4
/**
 * An integer 4-vector.
 */
// Gery-Note: Constructors work redundantly -> x(x_) and then x = x??
struct I4
{
	short int l,x,y,z;
	
	I4(int x_, int y_, int z_, int l_): 
		x(x_), y(y_), z(z_), l(l_)
	{
		this->x = x;
		this->y = y;
		this->z = z;
		this->l = l;
	}
	
	I4(const I4& a):
		x(a.x), y(a.y), z(a.z), l(a.l)
	{
		x = a.x;
		y = a.y;
		z = a.z;
		l = a.l;
	}
	
	I4(): x(0), y(0), z(0), l(0){}
	
	const bool operator<(const I4& i) const
	{
		return (l<i.l || l==i.l && x<i.x || l==i.l && x==i.x && y<i.y || l==i.l && x==i.x && y==i.y && z<i.z);
	}
	
	const bool operator ==(const I4& i) const
	{
		return l==i.l && x==i.x && y==i.y && z==i.z;
	}
};

namespace MRAG
{
	
#pragma mark struct WTDataType
    
struct WTDataType
{
	Real scaling;
	
	template<typename T>  WTDataType(const T& t) {scaling = (Real)t;}
	operator float() { return scaling; }
};

#pragma mark struct BlockInfo
/**
 * Helper structure to pass info of a block (like ID, position, level and a reference to the actual block).
 */
struct BlockInfo
{
	int blockID;
	void * ptrBlock;
	short int index[3];
	short int level;
    
    /** Origin of this block in all dimensions. */	
	Real origin[3];
    /** Grid spacing in all dimensions. */
	Real h[3];

	template <typename T>
#ifdef _CUDA_SIDE
	__device__
#endif
	inline void pos(T p[2], int ix, int iy) const
	{
		p[0] = origin[0] + h[0]*ix;
		p[1] = origin[1] + h[1]*iy;
	}
	
    /**
     * Returns position in physical coordinates of point in block.
     * @param T     Filled with physical coordinates.
     * @param ix    x-index (0..Block::sizeX)
     * @param iy    y-index (0..Block::sizeY)
     * @param iz    z-index (0..Block::sizeZ)
     */
	template <typename T>
#ifdef _CUDA_SIDE
	__device__
#endif
	inline void pos(T p[3], int ix, int iy, int iz) const
	{
		p[0] = origin[0] + h[0]*ix;
		p[1] = origin[1] + h[1]*iy;
		p[2] = origin[2] + h[2]*iz;
	}
	
	BlockInfo(int ID, const int idx[3], int l, const Real pos[3], const Real spacing[3]):
	blockID(ID), level(l), ptrBlock(NULL)
	{
		index[0] = idx[0];
		index[1] = idx[1];
		index[2] = idx[2];
		
		origin[0] = pos[0];
		origin[1] = pos[1];
		origin[2] = pos[2];
		
		h[0] = spacing[0];
		h[1] = spacing[1];
		h[2] = spacing[2];
	}
	
	BlockInfo(int ID, const int idx[3], int l):
	blockID(ID), level(l), ptrBlock(NULL)
	{
		index[0] = idx[0];
		index[1] = idx[1];
		index[2] = idx[2];
		
		origin[0] = 0;
		origin[1] = 0;
		origin[2] = 0;
		
		h[0] = 0;
		h[1] = 0;
		h[2] = 0;
	}
	
	BlockInfo():blockID(-1), level(-1), ptrBlock(NULL) {}
	
	BlockInfo& operator=(const BlockInfo& b)
	{
		blockID = b.blockID;
		level = b.level;
		ptrBlock = b.ptrBlock;
		
		index[0] = b.index[0];
		index[1] = b.index[1];
		index[2] = b.index[2];
		
		origin[0] = b.origin[0];
		origin[1] = b.origin[1];
		origin[2] = b.origin[2];
		
		h[0] = b.h[0];
		h[1] = b.h[1];
		h[2] = b.h[2];
		
		return *this;
	}
}; // struct BlockInfo
	
#pragma mark struct Subspace
    
struct Subspace
{	
	int start[3];
	int end[3];
	
	template<typename Integer>
	void setup(const int idx[3], const Integer stencil_start[3], const Integer stencil_end[3], const int block_size[3])
	{
		start[0] = idx[0]<0? 0 : (idx[0]==0 ? -stencil_start[0] : (block_size[0] - stencil_end[0]+1));
		start[1] = idx[1]<0? 0 : (idx[1]==0 ? -stencil_start[1] : (block_size[1] - stencil_end[1]+1));
		start[2] = idx[2]<0? 0 : (idx[2]==0 ? -stencil_start[2] : (block_size[2] - stencil_end[2]+1));
		
		end[0] = idx[0]<0? -stencil_start[0]  : (idx[0]==0 ? (block_size[0] - stencil_end[0] + 1) : block_size[0]);
		end[1] = idx[1]<0? -stencil_start[1]  : (idx[1]==0 ? (block_size[1] - stencil_end[1] + 1) : block_size[1]);
		end[2] = idx[2]<0? -stencil_start[2]  : (idx[2]==0 ? (block_size[2]- stencil_end[2] + 1) : block_size[2]);
	}
	
	template<typename Integer>
	int computeByteSize(const Integer stencil_start[3], const Integer stencil_end[3], const int sizeofElement) const 
	{
		const int nNumberOfElement = 
		(end[0] + stencil_end[0] - 1 - start[0] - stencil_start[0])*
		(end[1] + stencil_end[1] - 1 - start[1] - stencil_start[1])*
		(end[2] + stencil_end[2] - 1 - start[2] - stencil_start[2]);
	
		return nNumberOfElement*sizeofElement;
	}
}; // struct SubSpace

#pragma mark struct UniformPartition
    
struct UniformPartition
{	
	unsigned char total_work[3];
	unsigned char work_per_thread[3];
	char start[3];
	unsigned char threads[3];
	unsigned char passes_per_dim[3];
	
	UniformPartition()
	{
	}
	
	template<typename Integer>
	UniformPartition(const Integer s[3], const Integer e[3], const int nThreads[3],  const int * vMaxWorkPerThread = NULL)
	{
		setup(s, e, nThreads, vMaxWorkPerThread);
	}
	
	template<typename Integer>
	void setup(const Integer s[3], const Integer e[3], const int nThreads[3], const int * vMaxWorkPerThread = NULL)
	{
		threads[0] = nThreads[0];
		threads[1] = nThreads[1];
		threads[2] = nThreads[2];
		
		total_work[0] = max((int)0, e[0] - s[0]);
		total_work[1] = max((int)0, e[1] - s[1]);
		total_work[2] = max((int)0, e[2] - s[2]);
		
		if (vMaxWorkPerThread!= NULL)
		{
			work_per_thread[0] = min(vMaxWorkPerThread[0], max(1, total_work[0]/nThreads[0]));
			work_per_thread[1] = min(vMaxWorkPerThread[1], max(1, total_work[1]/nThreads[1]));
			work_per_thread[2] = min(vMaxWorkPerThread[2], max(1, total_work[2]/nThreads[2]));
		}
		else
		{
			work_per_thread[0] = max(1, total_work[0]/nThreads[0]);
			work_per_thread[1] = max(1, total_work[1]/nThreads[1]);
			work_per_thread[2] = max(1, total_work[2]/nThreads[2]);
		}
		
		start[0] = s[0];
		start[1] = s[1];
		start[2] = s[2];
		
		passes_per_dim[0] =  (int)ceil(total_work[0]/(double)(work_per_thread[0]*nThreads[0]));
		passes_per_dim[1] =  (int)ceil(total_work[1]/(double)(work_per_thread[1]*nThreads[1]));
		passes_per_dim[2] =  (int)ceil(total_work[2]/(double)(work_per_thread[2]*nThreads[2]));
	}
	
	const int passes_CPU()
	{
		return passes_per_dim[0]*passes_per_dim[1]*passes_per_dim[2];
	}

#ifdef _CUDA_SIDE
	__device__
#endif
	const int passes()
	{
		return passes_per_dim[0]*passes_per_dim[1]*passes_per_dim[2];
	}

#ifdef _CUDA_SIDE
	__device__
#endif
	bool findBlockWorkingArea(const int passID, int base_index[3], int end_index[3])
	{
		if (passID >= passes()) return false;

		const int pass[3] = {
			passID % passes_per_dim[0], 
			(passID/passes_per_dim[0])%passes_per_dim[1], 
			(passID/(passes_per_dim[0]*passes_per_dim[1]))
		};

		const int offset[3] = {
			work_per_thread[0]*(threads[0]*pass[0] + 0),
			work_per_thread[1]*(threads[1]*pass[1] + 0),
			work_per_thread[2]*(threads[2]*pass[2] + 0)
		};
		
		base_index[0] = start[0] + offset[0];
		base_index[1] = start[1] + offset[1];
		base_index[2] = start[2] + offset[2];
		
		end_index[0] = start[0] + min(offset[0] + threads[0]*work_per_thread[0], (int)total_work[0]);
		end_index[1] = start[1] + min(offset[1] + threads[1]*work_per_thread[1], (int)total_work[1]);
		end_index[2] = start[2] + min(offset[2] + threads[2]*work_per_thread[2], (int)total_work[2]);

		return true;
	}

#ifdef _CUDA_SIDE
	__device__
#endif
	void findThreadIndex( int threadID, int output[3])
	{
		output[0] = threadID % threads[0];
		output[1] = (threadID/threads[0]) % threads[1];
		output[2] = threadID/(threads[0]*threads[1]);
	}
	
#ifdef _CUDA_SIDE
	__device__
#endif
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
	
#ifdef _CUDA_SIDE
	__device__
#endif
	bool findWorkingInterval(const int thread_index[3], const int passID, int s[3], int e[3])
	{
		if (passID >= passes()) return false;
		
		/*const int pass[3] = {
			passID % 2, (passID/2)%2, (passID/4)%2
		};*/

		const int pass[3] = {
			passID % passes_per_dim[0], 
			(passID/passes_per_dim[0])%passes_per_dim[1], 
			(passID/(passes_per_dim[0]*passes_per_dim[1]))
		};

		const int offset[3] = {
			work_per_thread[0]*(threads[0]*pass[0] + thread_index[0]),
			work_per_thread[1]*(threads[1]*pass[1] + thread_index[1]),
			work_per_thread[2]*(threads[2]*pass[2] + thread_index[2])
		};
		
		s[0] = start[0] + min(offset[0], (int)total_work[0]);
		s[1] = start[1] + min(offset[1], (int)total_work[1]);
		s[2] = start[2] + min(offset[2], (int)total_work[2]);
		
		e[0] = start[0] + min(offset[0]+work_per_thread[0], (int)total_work[0]);
		e[1] = start[1] + min(offset[1]+work_per_thread[1], (int)total_work[1]);
		e[2] = start[2] + min(offset[2]+work_per_thread[2], (int)total_work[2]);
		
		if ((e[0]-s[0])*(e[1]-s[1])*(e[2]-s[2]) > 0) 
			return true;
		else
			return false;
	}	
}; // struct UniformPartition

	
template<int s, int e, int components=1>
struct Layer
{
	static const int start = s;
	static const int end = e;
	
	Real data[components][e-s][e-s];
	
	inline Real& operator() (int ix, int iy, int ic=0)
	{
		assert(ix >= start && ix < end);
		assert(iy >= start && iy < end);
		assert(ic >= 0 && ic<components);
		
		return data[ic][iy-s][ix-s];
	}
	
	inline Real read (int ix, int iy, int ic=0) const
	{
		assert(ix >= start && ix < end);
		assert(iy >= start && iy < end);
		assert(ic >= 0 && ic<components);
		
		return data[ic][iy-s][ix-s];
	}
	
};


}