/*
 *  MRAGBlockLab.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 6/24/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#include "MRAGrid.h"
#include "MRAGCommon.h"
#include "MRAGrid.h"
#include "MRAGBlockCollection.h"

#include "MRAGEnvironment.h"

#ifdef _MRAG_BLOCKLAB_ALLOCATOR_HEADER
#include _MRAG_BLOCKLAB_ALLOCATOR_HEADER
#endif

#include "MRAGMatrix3D.h"

#ifndef _MRAG_BLOCKLAB_ALLOCATOR
#define _MRAG_BLOCKLAB_ALLOCATOR std::allocator	
#endif

#pragma once
namespace MRAG
{
    
/**
 * Working copy of Block + Ghosts.
 * Data of original block is copied (!) here. So when changing something in
 * the lab we are not changing the original data (unless we use BlockLab::flush()).
 * @see MRAG::Block, BlockLab::prepare(), BlockLab::load()
 */
template<typename BlockType>
class BlockLab
{

protected:
	typedef typename BlockType::ElementType ElementType;
	enum eBlockLab_State {eMRAGBlockLab_Prepared, eMRAGBlockLab_Loaded, eMRAGBlockLab_Uninitialized};
	
	eBlockLab_State m_state;
	Matrix3D<ElementType, true, _MRAG_BLOCKLAB_ALLOCATOR> * m_cacheBlock;
	Real * m_weightPool;
	ElementType * m_valuePool;
	int m_weightPoolSize, m_valuePoolSize;
	int m_stencilStart[3], m_stencilEnd[3];
	const BlockCollection<BlockType>* m_refCollection;
	const BoundaryInfo* m_refBoundaryInfo;

	template<typename T> inline _MRAG_BLOCKLAB_ALLOCATOR<T> allocator() const { return _MRAG_BLOCKLAB_ALLOCATOR<T>();};
	
	template<typename T>
	void _release(T *& t) 
	{ 
		if (t != NULL)
		{
			allocator<T>().destroy(t);
			allocator<T>().deallocate(t,1); 
		}
		t = NULL; 
	}
	
	template<typename T>
	void _release_vector(T *& t, int n)
	{ 
		if (t != NULL)
			allocator<T>().deallocate(t, n); 
		
		t = NULL; 
	}
	
public:
	
	BlockLab():
	m_state(eMRAGBlockLab_Uninitialized),
	m_cacheBlock(NULL),
	m_weightPool(NULL), m_weightPoolSize(0),
	m_valuePool(NULL), m_valuePoolSize(0),
	m_refCollection(NULL), m_refBoundaryInfo(NULL)
	{
		m_stencilStart[0] = m_stencilStart[1] = m_stencilStart[2] = 0;
		m_stencilEnd[0] = m_stencilEnd[1] = m_stencilEnd[2] = 0;
	}
	
	~BlockLab()
	{
		_release(m_cacheBlock);
		_release_vector(m_weightPool, m_weightPoolSize);
		_release_vector(m_valuePool, m_valuePoolSize);
	}
	
    /**
     * Prepare the extended block.
     * @param collection    Collection of blocks in the grid (e.g. result of Grid::getBlockCollection()).
     * @param boundaryInfo  Info on the boundaries of the grid (e.g. result of Grid::getBoundaryInfo()).
     * @param stencil_start Maximal stencil used for computations at lower boundary.
     *                      Defines how many ghosts we will get in extended block.
     * @param stencil_end   Maximal stencil used for computations at lower boundary.
     *                      Defines how many ghosts we will get in extended block.
     */
    void prepare(const BlockCollection<BlockType>& collection, const BoundaryInfo& boundaryInfo, const int stencil_start[3], const int stencil_end[3])
	{
		m_refCollection = &collection;
		m_refBoundaryInfo = &boundaryInfo;
		
		assert(stencil_start[0]>= boundaryInfo.stencil_start[0]);
		assert(stencil_start[1]>= boundaryInfo.stencil_start[1]);
		assert(stencil_start[2]>= boundaryInfo.stencil_start[2]);
		assert(stencil_end[0]<= boundaryInfo.stencil_end[0]);
		assert(stencil_end[1]<= boundaryInfo.stencil_end[1]);
		assert(stencil_end[2]<= boundaryInfo.stencil_end[2]);
		
		m_stencilStart[0] = boundaryInfo.stencil_start[0];
		m_stencilStart[1] = boundaryInfo.stencil_start[1];
		m_stencilStart[2] = boundaryInfo.stencil_start[2];
		
		m_stencilEnd[0] = boundaryInfo.stencil_end[0];
		m_stencilEnd[1] = boundaryInfo.stencil_end[1];
		m_stencilEnd[2] = boundaryInfo.stencil_end[2];
		
		assert(m_stencilStart[0]<=m_stencilEnd[0]);
		assert(m_stencilStart[1]<=m_stencilEnd[1]);
		assert(m_stencilStart[2]<=m_stencilEnd[2]);
		
		if (m_cacheBlock == NULL || 
			m_cacheBlock->getSize()[0]!= BlockType::sizeX + m_stencilEnd[0] - m_stencilStart[0] -1 ||
			m_cacheBlock->getSize()[1]!= BlockType::sizeX + m_stencilEnd[1] - m_stencilStart[1] -1 ||
			m_cacheBlock->getSize()[2]!= BlockType::sizeX + m_stencilEnd[2] - m_stencilStart[2] -1 )
		{
			
			if (m_cacheBlock != NULL)
				_release(m_cacheBlock);

			m_cacheBlock = allocator< Matrix3D<ElementType,  true, _MRAG_BLOCKLAB_ALLOCATOR> >().allocate(1);
			
			allocator< Matrix3D<ElementType,  true, _MRAG_BLOCKLAB_ALLOCATOR> >().construct(m_cacheBlock, Matrix3D<ElementType,  true, _MRAG_BLOCKLAB_ALLOCATOR> ());
			
			m_cacheBlock->_Setup(BlockType::sizeX + m_stencilEnd[0] - m_stencilStart[0] -1,
								 BlockType::sizeY + m_stencilEnd[1] - m_stencilStart[1] -1,
								 BlockType::sizeZ + m_stencilEnd[2] - m_stencilStart[2] -1);

		}
		
		m_state = eMRAGBlockLab_Prepared;
	}
	
    /**
     * Load a block (incl. ghosts for it).
     * @param info  Reference to info of block to be loaded.
     */
	void load(const BlockInfo& info)
	{
		const BlockCollection<BlockType>& collection = *m_refCollection;
		const BoundaryInfo& boundaryInfo = *m_refBoundaryInfo;
		
		//0. couple of checks
		//1. load the block into the cache
		//2. create the point and the weight pools
		//3. compute the ghosts, put them into the cache
		
		//0.
		assert(m_state == eMRAGBlockLab_Prepared || m_state==eMRAGBlockLab_Loaded);
		assert(m_cacheBlock != NULL);
    
		const int nX = BlockType::sizeX;
		const int nY = BlockType::sizeY;
		const int nZ = BlockType::sizeZ;
		
		BoundaryInfoBlock& bbinfo = *boundaryInfo.boundaryInfoOfBlock.find(info.blockID)->second;
		bbinfo.lock();
		
		const vector< vector<IndexWP> >& ghosts = bbinfo.getGhosts();
		
		//1.
		{
			collection.lock(info.blockID);
			
			BlockType& block = collection[info.blockID];
			
			//*m_cacheBlock = ElementType();
			
			ElementType * ptrDestination = NULL;
			ElementType * ptrSource = &block(0);
			const size_t nCopyLength = nX * sizeof(ElementType);
			for(int iz=0; iz<nZ; iz++)
				for(int iy=0; iy<nY; iy++)
				{
					ptrDestination = &m_cacheBlock->Access(0-m_stencilStart[0], iy-m_stencilStart[1], iz-m_stencilStart[2]);
					memcpy(ptrDestination, ptrSource, nCopyLength); 
					ptrSource+=nX;
/*					for(int ix=0; ix<nX; ix++, ptrSource++, ptrDestination++)
						*ptrDestination = *ptrSource;*/
				}
			
			collection.release(info.blockID);
		}
		
		
		//2.
		{
			const int nValuePoolSize = bbinfo.getIndexPool().size();
			
			if (nValuePoolSize > m_valuePoolSize)
			{
				_release_vector(m_valuePool, m_valuePoolSize);
			
				m_valuePool = allocator<ElementType>().allocate(nValuePoolSize);
				
				m_valuePoolSize = nValuePoolSize;
			}
			
			collection.lock(bbinfo.dependentBlockIDs);

			 /* here we are going to do this:
			  for(int i=0; i<nValuePoolSize; i++)
				   m_valuePool[i] = collection[bbinfo.getIndexPool()[i].blockID][bbinfo.getIndexPool()[i].index];
				which is simply rewritten as:
			*/

			 {  
				const PointIndex  * p = &(bbinfo.getIndexPool()[0]);
				ElementType * e = &(m_valuePool[0]);
				BlockType * b;
				int oldBlockId = -1;
				
				for(int i=0;i<nValuePoolSize;++i){
				   if(oldBlockId!=(*p).blockID){
					  oldBlockId = (*p).blockID;
					  b = &(collection[oldBlockId]);
				   }
				   *e = (*b)[(*p).index];
				   ++e; ++p;
				}             
			 }
			
			collection.release(bbinfo.dependentBlockIDs);
			
			const int nWeightPoolSize = bbinfo.weightsPool.size();
			if (nWeightPoolSize > m_weightPoolSize)
			{
				_release_vector(m_weightPool, m_weightPoolSize);

				m_weightPoolSize = nWeightPoolSize;
				m_weightPool = allocator<Real>().allocate(nWeightPoolSize);
			}
				
			for(int i=0; i<nWeightPoolSize; i++)
				m_weightPool[i] = bbinfo.weightsPool[i];
		}
		
		//3.
		{
			for(int icode=0; icode<27; icode++)
			{
				if (icode == 1*1 + 3*1 + 9*1) continue;
				
				const int code[3] = { icode%3-1, (icode/3)%3-1, (icode/9)%3-1};
				
				const int s[3] = { 
					code[0]<1? (code[0]<0 ? m_stencilStart[0]:0 ) : nX, 
					code[1]<1? (code[1]<0 ? m_stencilStart[1]:0 ) : nY, 
					code[2]<1? (code[2]<0 ? m_stencilStart[2]:0 ) : nZ };
				
				const int e[3] = {
					code[0]<1? (code[0]<0 ? 0:nX ) : nX+m_stencilEnd[0]-1, 
					code[1]<1? (code[1]<0 ? 0:nY ) : nY+m_stencilEnd[1]-1, 
					code[2]<1? (code[2]<0 ? 0:nZ ) : nZ+m_stencilEnd[2]-1};
			
			
				assert(bbinfo.boundary[icode].nGhosts == (e[2]-s[2])*(e[1]-s[1])*(e[0]-s[0]));
					
				int currentghost = bbinfo.boundary[icode].start;
				
				for(int iz=s[2]; iz<e[2]; iz++)
					for(int iy=s[1]; iy<e[1]; iy++)
						for(int ix=s[0]; ix<e[0]; ix++, currentghost++)
						{
							const vector<IndexWP>& vWP= ghosts[currentghost];

							const vector<IndexWP>::const_iterator start = vWP.begin();
							const vector<IndexWP>::const_iterator end = vWP.end();
							
							ElementType ghost = ElementType();
							for(vector<IndexWP>::const_iterator it=start; it!= end; it++)
							{
								const Real w = (m_weightPool[it->weights_index[0]]*m_weightPool[it->weights_index[1]]*m_weightPool[it->weights_index[2]]);
								ghost += (m_valuePool[it->point_index]*w);
							}
							
							m_cacheBlock->Access(ix-m_stencilStart[0], iy-m_stencilStart[1], iz-m_stencilStart[2]) = ghost;
						}
			}
		}

		bbinfo.release();
		
		m_state = eMRAGBlockLab_Loaded;
	}
	
    /**
     * Write back the (changed) data in the lab overwriting the original data.
     * @param info  Reference to info of block to be flushed back.
     */
	void flush(BlockInfo& info)
	{
		const BlockCollection<BlockType>& collection = *m_refCollection;
		
		//0. couple of checks
		//1. copy back the inside data
		
		//0.
		assert(m_state == eMRAGBlockLab_Loaded);
		
		//1.
		{
			const int nX = BlockType::sizeX;
			const int nY = BlockType::sizeY;
			const int nZ = BlockType::sizeZ;
			
			BlockType& block = collection[info.blockID];
			for(int iz=0; iz<nZ; iz++)
				for(int iy=0; iy<nY; iy++)
					for(int ix=0; ix<nX; ix++)
						 block(ix, iy, iz) =  m_cacheBlock->Access(ix-m_stencilStart[0], iy-m_stencilStart[1], iz-m_stencilStart[2]);
		}
		
		m_state = eMRAGBlockLab_Prepared;
	}
	
    /**
     * Get a single element from the block.
     * stencil_start and stencil_end refer to the values passed in BlockLab::prepare().
     *
     * @param ix    Index in x-direction (stencil_start[0] <= ix < BlockType::sizeX + stencil_end[0] - 1).
     * @param iy    Index in y-direction (stencil_start[1] <= iy < BlockType::sizeY + stencil_end[1] - 1).
     * @param iz    Index in z-direction (stencil_start[2] <= iz < BlockType::sizeZ + stencil_end[2] - 1).
     */
	ElementType& operator()(int ix, int iy=0, int iz=0)
	{
		assert(m_state == eMRAGBlockLab_Loaded);
		
		const int nX = m_cacheBlock->getSize()[0];
		const int nY = m_cacheBlock->getSize()[1];
		const int nZ = m_cacheBlock->getSize()[2];
		
		assert(ix-m_stencilStart[0]>=0 && ix-m_stencilStart[0]<nX);
		assert(iy-m_stencilStart[1]>=0 && iy-m_stencilStart[1]<nY);
		assert(iz-m_stencilStart[2]>=0 && iz-m_stencilStart[2]<nZ);
		
		return m_cacheBlock->Access(ix-m_stencilStart[0], iy-m_stencilStart[1], iz-m_stencilStart[2]);
	}
	
    /** Just as BlockLab::operator() but returning a const. */
	const ElementType& read(int ix, int iy=0, int iz=0) const
	{
		assert(m_state == eMRAGBlockLab_Loaded);
		
		const int nX = m_cacheBlock->getSize()[0];
		const int nY = m_cacheBlock->getSize()[1];
		const int nZ = m_cacheBlock->getSize()[2];
		
		assert(ix-m_stencilStart[0]>=0 && ix-m_stencilStart[0]<nX);
		assert(iy-m_stencilStart[1]>=0 && iy-m_stencilStart[1]<nY);
		assert(iz-m_stencilStart[2]>=0 && iz-m_stencilStart[2]<nZ);
		
		return m_cacheBlock->Access(ix-m_stencilStart[0], iy-m_stencilStart[1], iz-m_stencilStart[2]);
	}
	
private:
	
	//forbidden
	BlockLab(const BlockLab&):
	m_state(eMRAGBlockLab_Uninitialized),
	m_cacheBlock(NULL),
	m_weightPool(NULL), m_weightPoolSize(0),
	m_valuePool(NULL), m_valuePoolSize(0) {abort();}
	
	BlockLab& operator=(const BlockLab&){abort(); return *this;}
	
};

}


