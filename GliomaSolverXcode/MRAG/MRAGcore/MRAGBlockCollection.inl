/*
 *  MRAGBlockCollection.inl
 *  MRAG
 *
 *  Created by Diego Rossinelli on 4/29/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#undef min
#undef max
#include <set>
#include <assert.h>
#include <math.h>

namespace MRAG
{
	template<typename BlockType_> BlockCollection<BlockType_>::BlockCollection():
		m_blockIDToBlockPointers(), m_blockIDToChunck(),
		m_setChunks(), m_trash(),
		m_currentChunk(NULL),
		m_nAvailableBlocksInCurrentChunk(0), allocator()
	{
	}
	
	template<typename BlockType_> BlockCollection<BlockType_>::~BlockCollection()
	{
		clear();
	}
	
	template<typename BlockType_>
	void BlockCollection<BlockType_>::_trash()
	{
		//old: if (m_trash.size()>(int)ceil(m_setChunks.size()*0.35))
		
		//new policy in megabytes
		const long long int maxCapacity = 50<<20;
		const long long int currentSize = m_trash.size()*sizeof(BlockType_)*nChunkSize;
		
		if (currentSize>maxCapacity)
		{
	//		printf("Too Many blocks in the trash (%.2f MB, max %.2f)\n",currentSize/1024./1024., maxCapacity/1024./1024.);
			
			int nToTrash = (int)ceil((currentSize - maxCapacity)/ 
									 (double)(m_trash.size()*sizeof(BlockType_)));
			_emptyTrash(nToTrash);
		}
	}

	template<typename BlockType_>
	void BlockCollection<BlockType_>::_emptyTrash(int nChunks)
	{
		if (nChunks == -1)
		{
			const int n = m_trash.size();
			for(int i=0; i<n; i++)
			{
				//delete [] m_trash[i]->p;
				_deallocate(m_trash[i]->p, nChunkSize);
				delete m_trash[i];
			}
			
			m_trash.clear();
			//printf("Deleted: %d\n", n);
			
		}
		else
		{
			assert(nChunks>0);
			
			const int n = std::min((int)nChunks,  (int)m_trash.size());
			for(int i=0; i<n; i++)
			{
				//delete [] m_trash.back()->p;
				_deallocate(m_trash.back()->p, nChunkSize);
				delete m_trash.back();
				
				m_trash.pop_back();
			}
		}
	}
	
	template<typename BlockType_>
	inline BlockType_& BlockCollection<BlockType_>::operator[](const int blockID) const
	{ 
		typename map<int, BlockType*>::const_iterator it = m_blockIDToBlockPointers.find(blockID);
		
		assert(it != m_blockIDToBlockPointers.end());
		
		return *it->second;
	}
	
	template<typename BlockType_>
	inline int BlockCollection<BlockType_>::_createIDs(int n)
	{
		static int scounterID = 0;
		
		int result = scounterID;
		
		scounterID+=n;
		
		return result;
	}
	
	template<typename BlockType_>
	void BlockCollection<BlockType_>::erase(const int ID)
	{
		//printf("erase %d\n", ID);
		//1. find the chunk containing this id
		//2. set the active flag as false
		//3. if the chunk is completely empty put it in the trash
		
		//1.
		typename map<int, Chunk *>::iterator itChunk = m_blockIDToChunck.find(ID);
		assert(itChunk!=m_blockIDToChunck.end());
		
		//2.
		Chunk & chunk = *itChunk->second;
		
		chunk.nActives--;

		if (&chunk == m_currentChunk)
		{
			m_currentChunk = NULL;
			m_nAvailableBlocksInCurrentChunk = 0;
			//printf("(SKIPPA IL MALEDETTO BASTARDO. CHE BASTARDO.)\n"); 
		}
		
		assert(chunk.startID<= ID);
		assert(chunk.startID + nChunkSize> ID);
		assert(chunk.nActives >= 0);
		
		//3.
		m_blockIDToChunck.erase(ID);
		m_blockIDToBlockPointers.erase(ID);
		if (chunk.nActives == 0)
		{
			m_setChunks.erase(&chunk);
			m_trash.push_back(&chunk);
			_trash();
		}
	}
	
	template<typename BlockType_>
	void BlockCollection<BlockType_>::clear()
	{
		for(typename set<Chunk *>::iterator it = m_setChunks.begin(); it!=m_setChunks.end(); it++)
			m_trash.push_back(*it);
		
		m_setChunks.clear();
		m_blockIDToBlockPointers.clear();
		m_blockIDToChunck.clear();
		
		_emptyTrash();
	}
	
	
	template<typename BlockType_>
	inline vector<int> BlockCollection<BlockType_>::_allocateBlockInChunk(int & inoutRequestedBlocks, Chunk * chunk,  int& inoutAvailableBlocksInCurrentChunk) const
	{
		std::vector<int> vIDs;
		
		assert(_implies(inoutAvailableBlocksInCurrentChunk == 0, chunk == NULL));
		assert(_implies(chunk == NULL, inoutAvailableBlocksInCurrentChunk == 0));
			   
		if (inoutAvailableBlocksInCurrentChunk == 0) return vIDs;

		const int nAllocations = std::min(inoutAvailableBlocksInCurrentChunk, inoutRequestedBlocks);
		vIDs.resize(nAllocations);
		
		const int startID = chunk->startID + nChunkSize - inoutAvailableBlocksInCurrentChunk; 
		
		
		for(int i=0; i<nAllocations; i++)
		{
			assert(startID + i >= chunk->startID );
			assert(startID + i <= chunk->startID  + nChunkSize );
			vIDs[i] = startID + i;
		}

		chunk->nActives += nAllocations;
		inoutAvailableBlocksInCurrentChunk -= nAllocations;
		inoutRequestedBlocks -= nAllocations;
		
		return vIDs;
	}
	
	template<typename BlockType_>
	inline vector<int> BlockCollection<BlockType_>::_allocateChunks(int & inoutRequestedBlocks, 
																	set<Chunk*>& inoutChunks, map<int, BlockType*>& inoutIDToBlockPointers, 
																	map<int, Chunk *>& inoutBlockIDToChunck, vector<Chunk*>& inoutRecycledChunks,  
																	Chunk *& outCurrentChunk, int & outAvailableBlocksInCurrentChunk) const
	{
		const int nChunks = (int)ceil(inoutRequestedBlocks/(double)nChunkSize);
		assert(nChunks>0);
		const int nLastBlocks = (inoutRequestedBlocks % nChunkSize);
		
		const int startID = _createIDs(nChunks*nChunkSize);
		std::vector<int> vIDs(inoutRequestedBlocks);
		
		outCurrentChunk = NULL;
		outAvailableBlocksInCurrentChunk = nChunks*nChunkSize - inoutRequestedBlocks;

		for(int iChunk=0; iChunk<nChunks; iChunk++)
		{
			Chunk * chunk = NULL;
			BlockType* blocks  = NULL;
			
			if (inoutRecycledChunks.size() == 0)
			{
				//blocks= new BlockType[nChunkSize];
				blocks = _allocate(nChunkSize);
				chunk = new Chunk(blocks, nChunkSize, startID + iChunk*nChunkSize);
			}
			else
			{
				chunk = inoutRecycledChunks.back();
				inoutRecycledChunks.pop_back();
				blocks = chunk->p;
				*chunk = Chunk(blocks, nChunkSize, startID + iChunk*nChunkSize);
				
				assert(m_setChunks.find(chunk) == m_setChunks.end());
			}
			
			if (iChunk == nChunks-1 && nLastBlocks != 0)
			{
				
				chunk->nActives = nLastBlocks;
				
				assert(chunk!=NULL);
				outCurrentChunk = chunk;
			}
				
			assert(m_setChunks.find(chunk) == m_setChunks.end());
			inoutChunks.insert(chunk);
			
			
			for(int iBlock=0; iBlock<nChunkSize; iBlock++)
			{
				const int currID = startID + iChunk*nChunkSize + iBlock;
				
				if(iChunk*nChunkSize + iBlock < inoutRequestedBlocks)
					vIDs[currID - startID] = currID;
				
				inoutIDToBlockPointers[currID] = (BlockType*)&blocks[iBlock];
				inoutBlockIDToChunck[currID] = chunk;
				assert(chunk->startID <= currID);
				assert(chunk->startID + nChunkSize > currID);
			}
		}
		
		inoutRequestedBlocks = 0;
		assert(_implies(m_nAvailableBlocksInCurrentChunk == 0, m_currentChunk == NULL));
		assert(_implies(m_currentChunk == NULL, m_nAvailableBlocksInCurrentChunk == 0));

		
		return vIDs;
	}
	
	template<typename BlockType_>
	std::vector<int> BlockCollection<BlockType_>::create(int nNumberOfBlocks)
	{
		int nRequestedBlocks = nNumberOfBlocks;
		
		vector<int> vSingleIDs = _allocateBlockInChunk(nRequestedBlocks, m_currentChunk, m_nAvailableBlocksInCurrentChunk);
		
		if (m_nAvailableBlocksInCurrentChunk==0) 
			m_currentChunk = NULL;
		
		assert(_implies(m_nAvailableBlocksInCurrentChunk == 0, m_currentChunk == NULL));
		assert(_implies(m_currentChunk == NULL, m_nAvailableBlocksInCurrentChunk == 0));
		
		if (nRequestedBlocks == 0) return vSingleIDs;
		
		vector<int> vIDs = _allocateChunks(nRequestedBlocks,
										   m_setChunks, m_blockIDToBlockPointers, m_blockIDToChunck,m_trash, 
										   m_currentChunk, m_nAvailableBlocksInCurrentChunk);
		
		vIDs.resize(vIDs.size() + vSingleIDs.size());
		
		copy(vSingleIDs.begin(), vSingleIDs.end(), vIDs.end() - vSingleIDs.size());
		
		return vIDs;
	}
	
	template<typename BlockType_>
	float BlockCollection<BlockType_>::getMemorySize(bool bCountTrashAlso) const
	{
		int memsize = 0;
				
		memsize += m_setChunks.size()*(sizeof(Chunk));
		memsize += m_blockIDToBlockPointers.size()*(sizeof(BlockType*) + sizeof(int));
		memsize += m_blockIDToChunck.size()*(sizeof(Chunk *) + sizeof(int));
		//const int r = memsize;
		for(typename set<Chunk *>::const_iterator it = m_setChunks.begin(); it!=m_setChunks.end(); it++)
			memsize += nChunkSize*sizeof(BlockType);
		
		if (bCountTrashAlso)
		{
			const int n = m_trash.size();
			for(int i=0; i<n; i++)
				memsize += nChunkSize*sizeof(BlockType);
		}
			
		return memsize/(double)(1<<20);
	}
}
