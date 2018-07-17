/*
 *  MRAGBlockCollection_CUDA.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/8/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include <map>
using namespace std;

#include <cutil.h>
#include <cuda_runtime_api.h>

#include "MRAGcore/MRAGEnvironment.h"
#include "MRAGcore/MRAGBlockCollection.h"
#include "MRAGAddressable_CUDA.h"

#ifdef _MRAG_TBB
#include <tbb/spin_mutex.h>
#endif
namespace MRAG
{
	template<typename BlockType_>
	class BlockCollection_CUDA: public BlockCollection<BlockType_>, public Addressable_CUDA
	{
	public:
		
		typedef BlockType_ BlockType;
		
		BlockCollection_CUDA(): m_mapPendingLocks(), m_mapDuplicateLocks()
		{
		}
		
		~BlockCollection_CUDA()
		{
			//no pending locks should exist
			assert(m_mapPendingLocks.size() == 0);
		}
		
		
		BlockType& operator[](const int blockID) const
		{
			//1. check that the blockID is locked
			//2. return the local pointer
			
			//1.
			typename map<int, BlockType*>::const_iterator itLocal = m_mapPendingLocks.find(blockID);
			assert(itLocal != m_mapPendingLocks.end());
			
			//2.
			return *itLocal->second;
		}
		
		//REMOTE ACCESS INTERFACE 
		//void* lock(const int blockID) const{ return (void)&lock(blockID);}

		
		BlockType& lock(const int blockID) const 
		{ 
			//1. check if the block is not locked yet, if it was return the block
			//2. locally allocate a block
			//3. perform the transfer
			//4. put the lock relation in the PendingLocks map
			//5. return the block

			try{
#ifdef _MRAG_TBB
				BlockCollection_CUDA_Mutex::scoped_lock lock(m_mutex); //scopato a porta chiusa
#endif
				//1.
				typename map<int, BlockType*>::const_iterator itLocal = m_mapPendingLocks.find(blockID);
				
				if(itLocal != m_mapPendingLocks.end())
				{
					assert(itLocal->second != NULL);
					
					m_mapDuplicateLocks[blockID]++;
					
					return *itLocal->second;
				}
				
				//2.
				BlockType * ptrLocal = this->allocator.allocate(1);
				this->allocator.construct(ptrLocal, BlockType());
				
				//3.
				typename map<int, BlockType*>::const_iterator itRemote = this->m_blockIDToBlockPointers.find(blockID);
				assert(itRemote != this->m_blockIDToBlockPointers.end());
				CUDA_SAFE_CALL(cudaMemcpy(ptrLocal, itRemote->second, sizeof(BlockType), cudaMemcpyDeviceToHost));
				
				//4.
				m_mapPendingLocks[blockID] = ptrLocal;
				
				//5.
				return *ptrLocal;
			}
			catch(...)
			{
				printf("Error in BlockCollection_CUDA::lock!\n");
				abort();
			}
		}
		
		void release(const int blockID, const bool bWriteBack=true) const 
		{
			//1. check that the block associated with the ID was locked
			//2. check that the block also exist remotely
			//3. perform the transfer
			//4. dellocate the local block and remove the lock relation from the PendingLocks map
			
#ifdef _MRAG_TBB
			BlockCollection_CUDA_Mutex::scoped_lock lock(m_mutex); //scopato a porta chiusa
#endif
			
			//1.
			typename map<int, BlockType*>::const_iterator itLocal = m_mapPendingLocks.find(blockID);
			assert(itLocal != m_mapPendingLocks.end());
			
			if(m_mapDuplicateLocks[blockID] > 0)
			{
				m_mapDuplicateLocks[blockID]--;
				
				return;
			}
			
			//2.
			typename map<int, BlockType*>::const_iterator itRemote = this->m_blockIDToBlockPointers.find(blockID);
			assert(itRemote != this->m_blockIDToBlockPointers.end());
			
			//3.
			if(bWriteBack)
			{
				CUDA_SAFE_CALL(cudaMemcpy(itRemote->second, itLocal->second, sizeof(BlockType), cudaMemcpyHostToDevice));
			}
			
			//4.
			this->allocator.destroy(itLocal->second);
			this->allocator.deallocate(itLocal->second, 1);
			m_mapPendingLocks.erase(blockID);
		}
		
		void * addressOf(const int blockID) const
		{
			typename map<int, BlockType*>::const_iterator itRemote = this->m_blockIDToBlockPointers.find(blockID);
			
			assert(itRemote != this->m_blockIDToBlockPointers.end());
			
			return (void *)itRemote->second;
		}
	
		int getNumberOfPendingLocks() const
		{
			int sum = 0;
			
			for(typename map<int, int>::const_iterator it = m_mapDuplicateLocks.begin(); it!= m_mapDuplicateLocks.end(); it++)
				sum += it->second;

			for(typename map<int, BlockType *>::const_iterator it = m_mapPendingLocks.begin(); it!= m_mapPendingLocks.end(); it++)
				sum++;
	
			return sum;
		}

	protected:
		
		BlockType * _allocate(int nElements) const
		{
			BlockType * ptr;
			
			CUDA_SAFE_CALL(cudaMalloc((void **)&ptr, sizeof(BlockType)*nElements));
			
			return ptr;
		}
		
		void _deallocate(BlockType * ptr, int nElements) const
		{
			CUDA_SAFE_CALL(cudaFree((void *)ptr));
		}
		
		// called from the CUDA blockprocessing
		void * translate(const int blockID) const
		{
			typename map<int, BlockType*>::const_iterator itRemote = this->m_blockIDToBlockPointers.find(blockID);
			assert(itRemote != this->m_blockIDToBlockPointers.end());
			
			return (void*)itRemote->second;
		}
	private:
		//forbidden
		BlockCollection_CUDA(const BlockCollection_CUDA&):BlockCollection<BlockType_>(){abort();}
		
		BlockCollection_CUDA& operator=(BlockCollection_CUDA&){abort(); return *this;}
		mutable map<int, int> m_mapDuplicateLocks;
		mutable map<int, BlockType *> m_mapPendingLocks;
		
		template<typename Btype> friend struct BlockProcessing_CUDA;
		
#ifdef _MRAG_TBB
		typedef tbb::spin_mutex BlockCollection_CUDA_Mutex;
		mutable BlockCollection_CUDA_Mutex m_mutex;
#endif
	};
}
