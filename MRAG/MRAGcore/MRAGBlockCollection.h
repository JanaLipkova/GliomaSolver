/*
 *  MRAGBlockCollection.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 4/28/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#include <vector>
#pragma once
#include <map>
#include <set>
using namespace std;


#ifdef _MRAG_BLOCKCOLLECTION_ALLOCATOR_HEADER
#include _MRAG_BLOCKCOLLECTION_ALLOCATOR_HEADER
#endif

#include "MRAGMatrix3D.h"

#ifndef _MRAG_BLOCKCOLLECTION_ALLOCATOR
#define _MRAG_BLOCKCOLLECTION_ALLOCATOR std::allocator	
#endif

namespace MRAG
{
    
/**
 * Collection of blocks (used by MRAG::Grid to collect blocks).
 */
template<typename BlockType_>
class BlockCollection
{
public:
	static const int nChunkSize = 1;//2*(BlockType_::shouldProcessDirectionY?2:1)*(BlockType_::shouldProcessDirectionZ?2:1) - 1;
	typedef BlockType_ BlockType;
	
    /** Default constructor creating empty collection. */
	BlockCollection();
    /** Destructor making sure, we clean up all blocks. */
	~BlockCollection();
	
    /**
     * Create or resize collection of blocks to contain certain number of blocks.
     * @param nNumberOfBlocks   Number of blocks to have in this collection.
     * @return Vector of IDs to access the blocks.
     */
	vector<int> create(const int nNumberOfBlocks=1);
    /** Remove a certain block. */
	void erase(const int ID);
    /** Resets collection to an empty one. */
	void clear();
    /** Access a certain block. */
	virtual BlockType& operator[](const int blockID) const;
	
	//REMOTE ACCESS INTERFACE
	virtual BlockType& lock(const int blockID) const { return (*this)[blockID]; }

	virtual void release(const int blockID, const bool bWriteBack=true) const {}
	
	virtual void lock(const vector<int>& blockIDs, vector<BlockType*>* output = NULL) const 
	{
		if (output != NULL)
		{
			output->resize(blockIDs.size());
			
			for(int i=0; i<blockIDs.size(); i++)
				(*output)[i] = &lock(blockIDs[i]);
		}
		else
			for(vector<int>::const_iterator it = blockIDs.begin(); it!= blockIDs.end(); it++)
				lock(*it);
	}
	
	void release(const vector<int>& blockIDs, bool bWriteBack=true) const
	{
		for(vector<int>::const_iterator it = blockIDs.begin(); it!= blockIDs.end(); it++)
			release(*it);
	}
	
	virtual float getMemorySize(bool bCountTrashAlso=false) const;
	
	const float getBlockSize() const{ return sizeof(BlockType); }
protected:

	struct Chunk
	{
		BlockType * p;
		int startID;
		int nActives;
		
		Chunk(BlockType*p_, const int n_, const int startID_): p(p_), startID(startID_), nActives(n_){}
		~Chunk() {}

		Chunk(const Chunk& a): p(a.p), startID(a.startID), nActives(a.nActives) {}
		
		Chunk& operator =(const Chunk& a)
		{
			p = (a.p);
			startID = (a.startID);
			nActives = (a.nActives);
			
			return *this;
		}
	};
	
	virtual BlockType * _allocate(int nElements) const
	{
		BlockType * ptr = allocator.allocate(nElements);
		
		for(int i=0; i<nElements; i++)
			allocator.construct(ptr+i, BlockType_());
		
		return ptr;
	}
	
	virtual void _deallocate(BlockType * ptr, int nElements) const
	{
		for(int i=0; i<nElements; i++)
			allocator.destroy(ptr+i);
		
		allocator.deallocate(ptr, nElements);
	}
	
	static int _createIDs(int n=1);
	
	vector<int> _allocateBlockInChunk(int & inoutRequestedBlocks, Chunk * chunk, int& inoutAvailableBlocksInCurrentChunk) const;
	vector<int> _allocateChunks(int & inoutRequestedBlocks, 
								set<Chunk*>& inoutChunks, map<int, BlockType*>& inoutIDToBlockPointers, 
								map<int, Chunk *>& inoutBlockIDToChunck, vector<Chunk*>& inoutRecycledChunks,  
								Chunk *& outCurrentChunk, int & outAvailableBlocksInCurrentChunk) const;
	void _emptyTrash(int nChunks=-1);
	void _trash();
	
	map<int, BlockType*> m_blockIDToBlockPointers;
	map<int, Chunk *> m_blockIDToChunck;
	set<Chunk*> m_setChunks;
	vector<Chunk*> m_trash;
	
	Chunk * m_currentChunk;
	int m_nAvailableBlocksInCurrentChunk;
	
	mutable _MRAG_BLOCKCOLLECTION_ALLOCATOR<BlockType_> allocator;
	
private:
	//forbidden
	BlockCollection(const BlockCollection&):
	m_blockIDToBlockPointers(), m_blockIDToChunck(),
	m_setChunks(), m_trash(),
	m_currentChunk(NULL),
	m_nAvailableBlocksInCurrentChunk(0){abort();}
	
	BlockCollection& operator=(BlockCollection&){abort(); return *this;}
};
	
}
#include "MRAGBlockCollection.inl"
