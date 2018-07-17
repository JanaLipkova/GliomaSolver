/*
 *  MRAGBlock.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 4/23/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include "MRAGCommon.h"

namespace MRAG
{

/**
 * Represents a block of computational elements.
 * Refinement works based on blocks being split up in 4 blocks (or vice-versa).
 * cSizeX*cSizeY*cSizeZ elements of type DataType are contained in each block.
 * Data is aligned s.t. one should loop first (outermost) in z, then in y, then (innermost) in x.
 * DataType needs the following:
 *  - Default constructor (compiler-generated one may be enough)
 *  - Copy constructor (compiler-generated one may be enough)
 *  - DataType operator*(const DataType&, Real) (outside of class as we need const(!)...)
 *    (needed to refine/compress grids and to load ghosts)
 *  - void operator += (DataType)
 *    (needed to refine/compress grids and to load ghosts)
 */
template <typename DataType, int cSizeX=16, int cSizeY=1, int cSizeZ=1>
struct Block
{
    /**
     * Type of the single computational elements in the block.
     */
	typedef DataType ElementType;
	
    /**
     * The actual data of the block.
     */
	DataType data[cSizeZ][cSizeY][cSizeX];
	
    /**
     * Construct block and fill with copies of given element.
     * If no element is given the default constructor of DataType is used.
     * @param v     Element to use as template for all elements in the block.
     */
	Block(DataType v = DataType())
	{
		for(int iz=0; iz<cSizeZ; iz++)
			for(int iy=0; iy<cSizeY; iy++)
				for(int ix=0; ix<cSizeX; ix++)
					data[iz][iy][ix] = v;
	}
	
    /**
     * Get a single element from the block.
     *
     * @param ix    Index in x-direction (0 <= ix < cSizeX).
     * @param iy    Index in y-direction (0 <= iy < cSizeY).
     * @param iz    Index in z-direction (0 <= iz < cSizeZ).
     */
#ifdef _CUDA_SIDE
	__device__
#endif
	inline DataType& operator()(int ix, int iy=0, int iz=0)
	{
#ifndef _CUDA_SIDE
		if (cDoDebugChecks)
		{
			assert(ix>=0 && ix<cSizeX);
			assert(iy>=0 && iy<cSizeY);
			assert(iz>=0 && iz<cSizeZ);
			assert((DataType *)&data[iz][iy][ix] == (DataType*)data + cSizeX*cSizeY*iz + cSizeX*iy + ix);
		}
#endif
		
		return *((DataType*)(data) + ix + iy*cSizeX + iz*cSizeY*cSizeX);
		//return data[iz][iy][ix];
	}
	
    /**
     * Get a single element from the block. Does just the same as operator() but with less checks.
     *
     * @param ix    Index in x-direction (0 <= ix < cSizeX).
     * @param iy    Index in y-direction (0 <= iy < cSizeY).
     * @param iz    Index in z-direction (0 <= iz < cSizeZ).
     */
#ifdef _CUDA_SIDE
	__device__
#endif
	inline DataType& getReference(int ix, int iy=0, int iz=0)
	{
		return *((DataType *)(data + ix + iy*cSizeX + iz*cSizeY*cSizeX));
	}
	
    /**
     * Get a single element from the block based on a single index.
     * Same as block(ix,iy,iz) with i = ix + iy*cSizeX + iz*cSizeY*cSizeX)
     *
     * @param i     Single index into the block.
     */
#ifdef _CUDA_SIDE
	__device__
#endif
	inline DataType& operator[](int i)
	{
		//printf("accessing %d\n", i);
		/*if (i>=cSizeX*cSizeY*cSizeZ || i<0)
		{
			/*assert(i>=0);
			assert(i<cSizeX*cSizeY*cSizeZ);
			//abort();
		}*/
		
		return *((DataType *)data+i);
	}
	
	inline void  wtdata(Real& coeff, int code, int ix, int iy=0, int iz=0)
	{
	}
	
    /**
     * Number of elements in x-direction.
     */
	static const int sizeX = cSizeX;
    /**
     * Number of elements in y-direction.
     */
	static const int sizeY = cSizeY;
    /**
     * Number of elements in z-direction.
     */
	static const int sizeZ = cSizeZ;

    /**
     * Is x-direction relevant? (more than 1 element)
     */
	static const bool shouldProcessDirectionX = cSizeX>1;
    /**
     * Is y-direction relevant? (more than 1 element)
     */
	static const bool shouldProcessDirectionY=  cSizeY>1; 
    /**
     * Is z-direction relevant? (more than 1 element)
     */
	static const bool shouldProcessDirectionZ = cSizeZ>1; 
	
    bool shouldBeRefined() { return false; }
	bool shouldBeCompressed() { return false; }
    
    /**
     * Return size of an instance of this class.
     */
	inline static int getMemSize()  { return cSizeX*cSizeY*cSizeZ*sizeof(DataType); }
};

}
