/*
 *  MRAGBlockProcessing_CUDA.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/10/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGBoundaryBlockPack.h"
#include "MRAGcore/MRAGBoundaryBlockInfo.h"
#include "MRAGcore/MRAGCache.h"

#include "../MRAGmultithreading/MRAGBlockCollection_CUDA.h"
#include "../MRAGmultithreading/MRAGBBPackHandler_CUDA.h"

#include <vector>
#include <map>
#include <string>
#undef min
#undef max
#include <queue>
#undef min
#undef max
#include <typeinfo>

using namespace MRAG;

void generateBoundarySubspaces(const int maxSHMEM, Subspace*& vSubspaces, int& n, const int block_size[3], const char stencil_start[3], const char stencil_end[3], const int sizeOfElement);

namespace MRAG
{
	struct IDVectorBlockInfo
	{
		priority_queue<int> blockIDs;
		
		IDVectorBlockInfo(vector<BlockInfo>& vInfo):
		blockIDs()
		{
			for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it!=vInfo.end(); it++)
				blockIDs.push(it->blockID);
		}
		
		bool operator < (const IDVectorBlockInfo& a) const 
		{ 
			priority_queue<int> q1 = blockIDs;
			priority_queue<int> q2 = a.blockIDs;
			
			const int n = min(q1.size(), q2.size());
			
			for(int i=0; i<n; i++, q1.pop(), q2.pop())
			{
				const int v1 = q1.top();
				const int v2 = q2.top();
				
				if (v1 < v2) return true;
				else if (v1 > v2) return false;
			}
			
			return (q1.size()==0 && q2.size()>0);
		}
		
		bool operator == (const IDVectorBlockInfo& a) const
		{
			priority_queue<int> q1 = blockIDs;
			priority_queue<int> q2 = a.blockIDs;
			
			if (q1.size() != q2.size()) return false;
			
			const int n = q1.size();
			
			for(int i=0; i<n; i++, q1.pop(), q2.pop())
				if (q1.top() != q2.top()) 
					return false;
			
			return true;
		}
	};
	
	struct IDBoundaryPartitions
	{
		int maxSHMEM, sizeofElement;
		int block_size[3], stencil_start[3], stencil_end[3];
		
		template<typename Integer>
		IDBoundaryPartitions(int maxSHMEM_, int sizeofElement_, const Integer stencil_start[3],  
							 const Integer stencil_end[3], const int block_size[3]):
		maxSHMEM(maxSHMEM_), sizeofElement(sizeofElement_)
		{
			this->block_size[0] = block_size[0];
			this->block_size[1] = block_size[1];
			this->block_size[2] = block_size[2];
			
			this->stencil_start[0] = stencil_start[0];
			this->stencil_start[1] = stencil_start[1];
			this->stencil_start[2] = stencil_start[2];
			
			this->stencil_end[0] = stencil_end[0];
			this->stencil_end[1] = stencil_end[1];
			this->stencil_end[2] = stencil_end[2];
		}
		
		bool operator == (const IDBoundaryPartitions& a) const
		{
			assert ( sizeof(IDBoundaryPartitions) % 4 == 0);
			
			int * ptr1 = (int*)this;
			int * ptr2 = (int*)&a;
			
			const int n = sizeof(IDBoundaryPartitions)>>2;
			
			for(int i=0; i<n; i++)
				if (ptr1[i] != ptr2[i]) return false;
			
			return true;
		}
		
		bool operator < (const IDBoundaryPartitions& a) const 
		{ 
			assert ( sizeof(IDBoundaryPartitions) % 4 == 0);
			
			int * ptr1 = (int*)this;
			int * ptr2 = (int*)&a;
			
			const int n = sizeof(IDBoundaryPartitions)>>2;
			
			for(int i=0; i<n; i++)
			{
				if (ptr1[i] < ptr2[i]) return true;
				else if (ptr1[i] > ptr2[i]) return false;
			}
			
			return false;
		}
	};
	
	struct BoundaryPartitionSet
	{
		int nSubspaces;
		int sizeSHMEM;
		UniformPartition * ptrReadLutGPU, * ptrWriteLutGPU;
		
		BoundaryPartitionSet(): ptrReadLutGPU(NULL), ptrWriteLutGPU(NULL), sizeSHMEM(0), nSubspaces(0) {}
		
		BoundaryPartitionSet(const BoundaryPartitionSet& a): 
			sizeSHMEM(a.sizeSHMEM),nSubspaces(a.nSubspaces),
			ptrReadLutGPU(a.ptrReadLutGPU), ptrWriteLutGPU(a.ptrWriteLutGPU) {}
		
		BoundaryPartitionSet& operator = (const BoundaryPartitionSet& a)
		{
			nSubspaces = a.nSubspaces;
			sizeSHMEM = a.sizeSHMEM;
			ptrReadLutGPU = a.ptrReadLutGPU;
			ptrWriteLutGPU = a.ptrWriteLutGPU;
			
			return *this;
		}
	};
	
template<typename BlockType_>
struct BlockProcessing_CUDA
{
private:
	mutable Cache<IDVectorBlockInfo, void *> m_cacheBlocksInfoGPU;
	Cache<IDBoundaryPartitions, BoundaryPartitionSet> m_cachePartitionsGPU;
	mutable map<string, map<vector<char>, void*> > m_mapProcessingGPU;
	
	BBPackHandler_CUDA packHandler;
	BlockInfo * m_ptrInfos;//, * m_ptrInfosGPU;
	
	//UniformPartition * ptrReadLutGPU, * ptrWriteLutGPU;
	
	int nAllocatedBlockInfos;

	/*void _release() // non const!
	{
		if (m_ptrInfos != NULL)
		{
			assert(m_ptrInfosGPU != NULL);
			
			CUDA_SAFE_CALL(cudaFreeHost((void *)m_ptrInfos));
			CUDA_SAFE_CALL(cudaFree((void *)m_ptrInfosGPU));

			m_ptrInfos = NULL;
			m_ptrInfosGPU = NULL;
		}
	}
	
	void _allocate(int n) // non const!
	{
		assert(m_ptrInfosGPU == NULL);
		assert(m_ptrInfos == NULL);
		
		nAllocatedBlockInfos = n;

		CUDA_SAFE_CALL(cudaMallocHost((void **)&m_ptrInfos, sizeof(BlockInfo)*nAllocatedBlockInfos));
		CUDA_SAFE_CALL(cudaMalloc((void **)&m_ptrInfosGPU, sizeof(BlockInfo)*nAllocatedBlockInfos));
	
		assert(m_ptrInfos != NULL);
		assert(m_ptrInfosGPU != NULL);
	}*/
	
	BlockInfo* _prepareBlockInfos(const BlockCollection_CUDA<BlockType_>& collection, vector<BlockInfo>& vInfo)// non const!
	{
		//1. build an id based on vInfo
		//2. check if we have worked with that id
		//3. yes: return the old pointer
		//4. no: allocate the pointer and return it
		
		//1.
		IDVectorBlockInfo _id(vInfo);
		
		//2.
		bool bCacheHit = false;
		BlockInfo *& ptrGPU = (BlockInfo *&)m_cacheBlocksInfoGPU.get(_id, bCacheHit);
		
		//3.
		//if (bCacheHit) printf("CACHE HIT in _prepareBlockInfos\n");
		if (bCacheHit) return ptrGPU;
		
		//4.
		{
			//A. if the cache slot was full, release it
			//B. re-allocate slot + workspace
			//C. fill the workspace, transmit to GPU

			const int nBlocks = vInfo.size();
			
			//A.
			if (ptrGPU != NULL) CUDA_SAFE_CALL(cudaFree((void *)ptrGPU));
			
			//B.
			CUDA_SAFE_CALL(cudaMalloc((void **)&ptrGPU, sizeof(BlockInfo)*nBlocks));
				
			if (nBlocks > nAllocatedBlockInfos)
			{
				if (m_ptrInfos != NULL)
					CUDA_SAFE_CALL(cudaFreeHost((void *)m_ptrInfos));
				
				CUDA_SAFE_CALL(cudaMallocHost((void **)&m_ptrInfos, sizeof(BlockInfo)*nBlocks));
				
				nAllocatedBlockInfos = nBlocks;
			}
			
			assert(m_ptrInfos != NULL);
			assert(ptrGPU != NULL);

			//C.
			for(int i=0; i<nBlocks; i++)
			{
				BlockInfo *info = m_ptrInfos + i;
				
				*info = vInfo[i];
				
				info->ptrBlock = collection.translate(info->blockID);
			}
			
			CUDA_SAFE_CALL(cudaMemcpy(ptrGPU, m_ptrInfos, nAllocatedBlockInfos*sizeof(BlockInfo), cudaMemcpyHostToDevice));
		
			return ptrGPU;
		}
	}
	
	template<typename Processing>
	Processing* _prepareProcessing(Processing& processing) const 
	{
		//1. create the id
		//2. look for a cached value
		//3. not found: store the new value
		//4. return the pointer
		
		string sProcessingName = typeid(Processing).name();
		map<vector<char>, void *>& mapID2ptr = m_mapProcessingGPU[sProcessingName];
		
		const int byteSize = sizeof(Processing);
		
		//1.
		vector<char> footprint(byteSize);
		
		{
			char * ptr = (char*)&processing;
			for(int i=0; i<byteSize; i++, ptr++)
				footprint[i] = *ptr;
			//printf("Type: %s\n", sProcessingName.data());
			//printf("Footprint: ");
			//for(int i=0; i<byteSize; i++)
			//		printf("%x",footprint[i]);
			//printf("\n");
			int ii;
			//cin>>ii;
		}			
		
		Processing * processingGPU = (Processing *)mapID2ptr[footprint];
		
		//2.
		if (processingGPU == NULL) // cache miss
		{
			//3.
			CUDA_SAFE_CALL(cudaMalloc((void **)&processingGPU, sizeof(Processing)));
			CUDA_SAFE_CALL(cudaMemcpy((void*)processingGPU, (void*)&processing, sizeof(Processing), cudaMemcpyHostToDevice));
			
			mapID2ptr[footprint] = processingGPU;
		}
		else
			;//printf("CACHE HIT in _prepareProcessing\n");
		
		//4.
		return processingGPU;
	}
	
	template<typename Processing>
	void _prepareBoundaryPartitions(const int nThreads[3], Processing& processing, Subspace * vSubspaces, const int nSubspaces, const int maxSizeSHMEM,
									UniformPartition *& ptrReadLutGPU, UniformPartition *& ptrWriteLutGPU, int& sizeSHMEM) const  
	{
		static int oldNumberOfSubspaces = 0;
		static UniformPartition * ptrReadLut = NULL;
		static UniformPartition * ptrWriteLut = NULL;
		
			if (ptrReadLutGPU != NULL)
			{
				assert(ptrWriteLutGPU != NULL);
				CUDA_SAFE_CALL( cudaFree((void *)ptrReadLutGPU));
				CUDA_SAFE_CALL( cudaFree((void *)ptrWriteLutGPU));
				
				ptrReadLutGPU = NULL;
				ptrWriteLutGPU = NULL;
			}
			else
				assert(ptrWriteLutGPU == NULL);
			
			CUDA_SAFE_CALL( cudaMalloc((void **)&ptrReadLutGPU, sizeof(UniformPartition)*nSubspaces));
			CUDA_SAFE_CALL( cudaMalloc((void **)&ptrWriteLutGPU, sizeof(UniformPartition)*nSubspaces));

		if (ptrReadLut == NULL || nSubspaces != oldNumberOfSubspaces)
		{
			if (ptrReadLut != NULL)
			{
				assert(ptrWriteLut != NULL);
				CUDA_SAFE_CALL( cudaFreeHost((void *)ptrReadLut));
				CUDA_SAFE_CALL( cudaFreeHost((void *)ptrWriteLut));
				
				ptrReadLut = NULL;
				ptrWriteLut = NULL;
			}
			else
				assert(ptrWriteLut == NULL);
			
			CUDA_SAFE_CALL( cudaMallocHost((void **)&ptrReadLut, sizeof(UniformPartition)*nSubspaces));
			CUDA_SAFE_CALL( cudaMallocHost((void **)&ptrWriteLut, sizeof(UniformPartition)*nSubspaces));
			
			oldNumberOfSubspaces = nSubspaces;
		}
			assert(ptrReadLut != NULL);
			assert(ptrWriteLut != NULL);
			assert(ptrReadLutGPU != NULL);
			assert(ptrWriteLutGPU != NULL);

		for(int i=0; i<nSubspaces; i++)
		{
			Subspace& subspace = vSubspaces[i];
			
			ptrWriteLut[i].setup(subspace.start, subspace.end, nThreads);
			
			const int read_start[3] = {
				subspace.start[0] + processing.stencil_start[0],
				subspace.start[1] + processing.stencil_start[1],
				subspace.start[2] + processing.stencil_start[2]
			};
			
			const int read_end[3] = {
				subspace.end[0] + processing.stencil_end[0] - 1,
				subspace.end[1] + processing.stencil_end[1] - 1,
				subspace.end[2] + processing.stencil_end[2] - 1
			};
			
			ptrReadLut[i].setup(read_start, read_end, nThreads);
		}
		
		CUDA_SAFE_CALL( cudaMemcpy(ptrReadLutGPU, ptrReadLut, sizeof(UniformPartition)*nSubspaces, cudaMemcpyHostToDevice));
		CUDA_SAFE_CALL( cudaMemcpy(ptrWriteLutGPU, ptrWriteLut, sizeof(UniformPartition)*nSubspaces, cudaMemcpyHostToDevice));
		
		//find minimum shared mem space
		{
			int maxCoverage = 0;
			for(int iSubspace=0;iSubspace<nSubspaces; iSubspace++)
			{
				int s[3] = {0,0,0};
				int e[3] = {0,0,0};
				
				ptrReadLut[iSubspace].getCoveredInterval_CPU(s,e);
				
				maxCoverage = max(maxCoverage, (e[0]-s[0])*(e[1]-s[1])*(e[2]-s[2]));
			}
			
			sizeSHMEM = maxCoverage*sizeof(typename BlockType_::ElementType);
		}
			
		//printf("SHARED MEMORY PER BLOCK: %f KB (%f MB)\n", sizeSHMEM/1024., sizeSHMEM/1024./1024.);
		
		assert(sizeSHMEM <= maxSizeSHMEM);
	}

	
public:
	BlockProcessing_CUDA(): m_ptrInfos(NULL), /*m_ptrInfosGPU(NULL),*/ nAllocatedBlockInfos(0) /*,
	ptrReadLutGPU(NULL), ptrWriteLutGPU(NULL) */ , m_cachePartitionsGPU(), m_cacheBlocksInfoGPU()
	{
		m_cachePartitionsGPU.setup(10);
		m_cacheBlocksInfoGPU.setup(10);
	}
	
	~BlockProcessing_CUDA(){ /**_release(); */}
	
	template <typename Collection , typename Processing>
	void process_mmt(vector<BlockInfo>& vInfo, const Collection& c, Processing& p)
	{
		try
		{
			//1.trasmute man in animal
			//2.invoke the beast
			
			//1.
			const int nThreads[3] = {2,2,1}; //1,32,1
			
			const int nBlocks = vInfo.size();
			if (nBlocks==0) return;
			
			const BlockCollection_CUDA<BlockType_>& collection = dynamic_cast<const BlockCollection_CUDA<BlockType_>&>(c);
			assert(collection.getNumberOfPendingLocks() == 0);

			BlockInfo * ptrInfosGPU = _prepareBlockInfos(collection, vInfo);
			
			Processing * processingGPU = _prepareProcessing(p);
			
			//2
			//printf("process_mmt_CUDA\n");
		 	process_mmt_CUDA(ptrInfosGPU, nBlocks, p, processingGPU, nThreads);
		}
		catch(...)
		{
			printf("Exception in BlockProcessing_CUDA::process_mmt!\n");
			abort();
		}
	}

	
	template <typename Processing, typename Collection>
	void process_mmt(vector<BlockInfo>& vInfo, Collection& c, const BoundaryInfo& b, Processing& p)
	{
		try
		{
			//1. translate man in animal BlocksInfo
			//2. translate man in animal BBInfo
			//3. prepare the lut for the patches
			//4. invoke the beast
			
			//1.
			const int nThreads[3] = {4,4,1};//16,16,1
			
			const int nBlocks = vInfo.size();
			if (nBlocks==0) return;
			
			const BlockCollection_CUDA<BlockType_>& collection = dynamic_cast<const BlockCollection_CUDA<BlockType_>&>(c);
			assert(collection.getNumberOfPendingLocks() == 0);

			BlockInfo * ptrInfosGPU = _prepareBlockInfos(collection, vInfo);
			
			Processing * processingGPU = _prepareProcessing(p);
			
			//2.
			BBPackHandler_CUDA::ArrayOfBoundaryPack& resultBBPacks = packHandler.translate(vInfo, b, collection);
			
			//3
			UniformPartition * ptrReadLutGPU = NULL;
			UniformPartition* ptrWriteLutGPU = NULL;
			int nSubspaces = 0;
			int sizeSHMEM = 0;
			
			{
				typedef BlockType_ B;
				typedef typename B::ElementType E;
				
				const int block_size[3] = {
					B::sizeX, B::sizeY, B::sizeZ
				};
				
				const int KB = 1024;
				const int maxSHMEM = 1*KB;

				IDBoundaryPartitions _id(maxSHMEM, sizeof(E), p.stencil_start, p.stencil_end, block_size);
				
				bool bCacheHit = false;
				BoundaryPartitionSet& partition_set = m_cachePartitionsGPU.get(_id, bCacheHit);
				
				if (!bCacheHit)
				{
					if (partition_set.ptrReadLutGPU != NULL)
					{
						assert(partition_set.ptrWriteLutGPU);
						
						CUDA_SAFE_CALL( cudaFree((void *)partition_set.ptrReadLutGPU));
						CUDA_SAFE_CALL( cudaFree((void *)partition_set.ptrWriteLutGPU));
						
						partition_set.ptrReadLutGPU = NULL;
						partition_set.ptrWriteLutGPU = NULL;
					}
					
					Subspace * vSubspaces = NULL;
					
					generateBoundarySubspaces(maxSHMEM, vSubspaces, partition_set.nSubspaces , block_size, p.stencil_start, p.stencil_end, sizeof(E));
					
					_prepareBoundaryPartitions(nThreads, p, vSubspaces, partition_set.nSubspaces, maxSHMEM,  partition_set.ptrReadLutGPU, partition_set.ptrWriteLutGPU, partition_set.sizeSHMEM );

					delete [] vSubspaces;
				}
				else
					;//printf("CACHE HIT in _prepareBoundaryPartitions\n");
				
				nSubspaces = partition_set.nSubspaces;
				sizeSHMEM = partition_set.sizeSHMEM;
				ptrReadLutGPU = partition_set.ptrReadLutGPU;
				ptrWriteLutGPU = partition_set.ptrWriteLutGPU;
			}
			
			//4.
			//printf("process_mmtgathering_CUDA\n");
			process_mmtgathering_CUDA(ptrInfosGPU, resultBBPacks.ptrPacks, nBlocks, 
									  ptrReadLutGPU, ptrWriteLutGPU, nSubspaces, p, processingGPU, nThreads, sizeSHMEM);
		}
		catch(...)
		{
			printf("Exception in BlockProcessing_CUDA::process_mmt(vector<BlockInfo>& vInfo, Collection& c, BoundaryInfo& b, Processing& p)!\n");
			abort();
		}
	}
};
}



