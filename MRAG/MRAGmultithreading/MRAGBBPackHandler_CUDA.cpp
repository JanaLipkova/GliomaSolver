/*
 *  MRAGBBPackHandler_CUDA.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/14/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#include "MRAGBBPackHandler_CUDA.h"


namespace MRAG
{
	BBPackHandler_CUDA::ArrayOfBoundaryPack::ArrayOfBoundaryPack(const int n): 
	nPacks(n), bFinalized(false), ptrBoundaryPacksCPU(NULL), ptrBoundaryPacksGPU(NULL), ptrPacks(NULL)
	{
		CUDA_SAFE_CALL(cudaMallocHost((void **)&ptrBoundaryPacksCPU, sizeof(BoundaryBlockPack*)*nPacks));

		assert(ptrBoundaryPacksCPU != NULL);
	}
	
	BBPackHandler_CUDA::ArrayOfBoundaryPack::~ArrayOfBoundaryPack()
	{
		CUDA_SAFE_CALL(cudaFreeHost((void *)ptrBoundaryPacksCPU));
		
		if(bFinalized && ptrBoundaryPacksGPU != NULL)
			CUDA_SAFE_CALL(cudaFree((void *)ptrBoundaryPacksGPU));
	}
	
	BoundaryBlockPack *& BBPackHandler_CUDA::ArrayOfBoundaryPack::operator[](int i)
	{
		assert(!bFinalized);
		assert(i>=0);
		assert(i<nPacks);
		
		return ptrBoundaryPacksCPU[i];
	}
	
	void BBPackHandler_CUDA::ArrayOfBoundaryPack::finalize()
	{
		assert(!bFinalized);
		
		CUDA_SAFE_CALL(cudaMalloc((void **)&ptrBoundaryPacksGPU, sizeof(BoundaryBlockPack*)*nPacks));
		CUDA_SAFE_CALL(cudaMemcpy(ptrBoundaryPacksGPU, ptrBoundaryPacksCPU, sizeof(BoundaryBlockPack*)*nPacks, cudaMemcpyHostToDevice));
		ptrPacks = ptrBoundaryPacksGPU;
		
		bFinalized = true;
	}

	
	void BBPackHandler_CUDA::_erase(int blockID, BoundaryBlockPack * ptrGPU)
	{
		//1. couple of checks
		//2. call free
		//3. remove the pointers from the datastructures
		
		//1.
		map<int, BoundaryBlockPack *>::iterator itPack = m_mapID2Pack.find(blockID);
		assert(itPack != m_mapID2Pack.end());
		
		map<BoundaryBlockPack*, BoundaryBlockPack *>::iterator itPackCPU = m_mapGPU2CPU.find(ptrGPU);
		assert(itPackCPU != m_mapGPU2CPU.end());
		
		map<int, FootPrint>::iterator itFootPrint = m_mapID2FootPrint.find(blockID);
		assert(itFootPrint != m_mapID2FootPrint.end());
		
		//2.
		CUDA_SAFE_CALL(cudaFree((void *)ptrGPU));
		CUDA_SAFE_CALL(cudaFreeHost((void *)itPackCPU->second));
		
		//3.
		m_mapID2Pack.erase(itPack);
		m_mapGPU2CPU.erase(itPackCPU);
		m_mapID2FootPrint.erase(itFootPrint);
	}

	BoundaryBlockPack* BBPackHandler_CUDA::_translate(int blockID, BoundaryInfoBlock& in_BBI, const char stencil_start[3], const char stencil_end[3])
	{
		//1. store the footprint
		//2. create the BBPack on the CPU
		//3. replicate it to the GPU, also by translating the right pointers
		//4. return the GPU pointer
		
		//1.
		m_mapID2FootPrint[blockID] = in_BBI.dependentBlockIDs;
		
		//2. BUTCHER BUTCHER
		BoundaryBlockPack * pack = NULL;

		{
			typedef Real WeightType;
			typedef short int GhostOffset;
			
			const int nWeights = in_BBI.weightsPool.size();
			const int nPoints = in_BBI.indexPool.size();
			const int nGhosts = in_BBI.ghosts.size();
			
			vector<short int> vGhostOffsets(nGhosts+1);
			int nInstructions = 0;
			
			typedef BoundaryInfoBlock::ReconstructionInfo ReconstructionInfo;
			int counter=0;
			for(vector<ReconstructionInfo>::const_iterator itGhost = in_BBI.ghosts.begin(); itGhost!= in_BBI.ghosts.end(); counter++, itGhost++)
			{
				vGhostOffsets[counter] = nInstructions;
				nInstructions += itGhost->size();
			}
			vGhostOffsets[counter] = nInstructions;
			assert(counter == nGhosts);
			
			const unsigned int byteSizes[5] = {
				sizeof(BoundaryBlockPack), 
				nWeights*sizeof(WeightType), 
				nPoints*sizeof(PointerToPoint),
				(nGhosts+1)*sizeof(GhostOffset),
				nInstructions*sizeof(InstructionItem)
			};
			
			const int nSizeInBytes = byteSizes[0] + byteSizes[1] + byteSizes[2] + byteSizes[3] + byteSizes[4];
			
			//char * ptr = (char*)malloc(nSizeInBytes);
			char * ptr = NULL;
			CUDA_SAFE_CALL(cudaMallocHost((void **)&ptr, nSizeInBytes));

			assert(ptr != NULL);
			assert((int)(ptr) % 4 == 0);
			
			pack = (BoundaryBlockPack *)ptr;
			
			pack->byteSize = nSizeInBytes;
			pack->nWeights = nWeights;
			pack->nPoints = nPoints;
			
			for(int i=0; i<27; i++) pack->boundary_start[i] = in_BBI.boundary[i].start;
			for(int i=0; i<3; i++) pack->stencil_start[i] = stencil_start[i];
			for(int i=0; i<3; i++) pack->stencil_end[i] = stencil_end[i];
			
			ptr += byteSizes[0];
			pack->ptrWeights = (Real*)ptr;
			for(int i=0; i<nWeights; i++)
				pack->ptrWeights[i] = in_BBI.weightsPool[i];
			
			ptr += byteSizes[1];
			pack->ptrPoints = (PointerToPoint*)ptr;
			for(int i=0; i<nPoints; i++)
			{
				pack->ptrPoints[i].ptrBlock = remoteCollection->addressOf(in_BBI.indexPool[i].blockID);
				pack->ptrPoints[i].index = in_BBI.indexPool[i].index;
			}
			
			ptr += byteSizes[2];
			pack->ptrGhostOffsets = (GhostOffset *)ptr;
			std::copy(vGhostOffsets.begin(), vGhostOffsets.end(), pack->ptrGhostOffsets);
			
			ptr += byteSizes[3];
			pack->ptrInstructions = (InstructionItem*)ptr;
			counter = 0;
			for(vector<ReconstructionInfo>::const_iterator itGhost = in_BBI.ghosts.begin(); itGhost!= in_BBI.ghosts.end(); itGhost++)
			{
				for(vector<IndexWP>::const_iterator it = itGhost->begin(); it!= itGhost->end(); it++, counter++)
				{
					InstructionItem& dest = pack->ptrInstructions[counter];
					
					dest.idxWeight[0] = it->weights_index[0];
					dest.idxWeight[1] = it->weights_index[1];
					dest.idxWeight[2] = it->weights_index[2];
					dest.idxPoint = it->point_index;
				}
			}
			assert(vGhostOffsets.back() == counter);
		}

		//3.
		BoundaryBlockPack * packGPU = NULL;
		{
			CUDA_SAFE_CALL(cudaMalloc((void **)&packGPU, pack->byteSize));
			
			const int drift = (char *)packGPU - (char *)pack;
			pack->ptrWeights = (Real *)((char *)pack->ptrWeights + drift);
			pack->ptrPoints = (PointerToPoint *)((char *)pack->ptrPoints + drift);
			pack->ptrGhostOffsets = (short int *)((char *)pack->ptrGhostOffsets + drift);
			pack->ptrInstructions = (InstructionItem *)((char *)pack->ptrInstructions + drift);

			CUDA_SAFE_CALL(cudaMemcpy(packGPU, pack, pack->byteSize, cudaMemcpyHostToDevice));
			
#if __DEVICE_EMULATION__
				memset(pack, 0, pack->byteSize);
#endif
			
			m_mapID2Pack[blockID] = packGPU;
			m_mapGPU2CPU [packGPU] = pack; // GPU Pack -> CPU pack
		}
		
		//4.
		return packGPU;
	}

}