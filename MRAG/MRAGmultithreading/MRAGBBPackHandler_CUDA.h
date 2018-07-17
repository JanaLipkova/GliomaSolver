/*
 *  MRAGBBPackHandler.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/14/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include <vector>
#include <map>

#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGBoundaryBlockPack.h"
#include "MRAGcore/MRAGBoundaryBlockInfo.h"
#include "MRAGcore/MRAGCache.h"
#include "MRAGAddressable_CUDA.h"

//#define __DEVICE_EMULATION__ 1
#include <cutil.h>
#include <cuda_runtime_api.h>

using namespace std;

namespace MRAG
{
	
	class BBPackHandler_CUDA
	{
		typedef vector<int> FootPrint;
		
		void _erase(int blockID, BoundaryBlockPack * ptrGPU);
		BoundaryBlockPack* _translate(int blockID, BoundaryInfoBlock& in_BBI, const char stencil_start[3], const char stencil_end[3]);
		
	public:
		
		struct ArrayOfBoundaryPack
		{
		public:
			ArrayOfBoundaryPack(const int n);
			
			~ArrayOfBoundaryPack();
			
			BoundaryBlockPack *& operator[](int i);
			
			void finalize();

			int nPacks;
			BoundaryBlockPack ** ptrPacks;
			
		
			
			BoundaryBlockPack ** ptrBoundaryPacksGPU;
			BoundaryBlockPack ** ptrBoundaryPacksCPU;
			bool bFinalized;
			
		private:
			//forbidden
			ArrayOfBoundaryPack(const ArrayOfBoundaryPack& b):
			bFinalized(false), ptrPacks(NULL), ptrBoundaryPacksGPU(NULL), ptrBoundaryPacksCPU(NULL), nPacks(0) { abort(); }
			ArrayOfBoundaryPack& operator=(const ArrayOfBoundaryPack& b){abort(); return *this;}
			
		};
		
		template <typename Collection>
		ArrayOfBoundaryPack& translate(const vector<BlockInfo>& vInfo, const BoundaryInfo& boundaryInfo, const Collection& collection)
		{
			//1. setup
			//2. check wheter there are obsolete packs, remove them
			//3. iterating over vInfo, translate every BBI into a BBPack
			//4. return the collected result
			
			//1.
			if (m_result != NULL)
				delete m_result;
			
			m_result = new ArrayOfBoundaryPack(vInfo.size());
			
			ArrayOfBoundaryPack& result = *m_result;
			
			remoteCollection = &collection;
			
			//2.
			{
				vector<BoundaryBlockPack *> vToErase;
				vector<int> vToEraseID;
				for(map<int, BoundaryBlockPack*>::const_iterator it = m_mapID2Pack.begin(); it != m_mapID2Pack.end(); it++)
				{
					map<int, BoundaryInfoBlock*>::const_iterator itBBI =  boundaryInfo.boundaryInfoOfBlock.find(it->first);
					
					const bool bFound = itBBI != boundaryInfo.boundaryInfoOfBlock.end();
					
					if (!bFound)
					{
						vToErase.push_back(it->second);
						vToEraseID.push_back(it->first);
					}
				}
				
				for(int i=0; i<vToEraseID.size(); i++)
					_erase(vToEraseID[i], vToErase[i]);
			}
			
			//3.
			{
				for(int i=0; i<vInfo.size(); i++)
				{
					const int blockID = vInfo[i].blockID;
					
					map<int, BoundaryInfoBlock*>::const_iterator itBBI = boundaryInfo.boundaryInfoOfBlock.find(blockID);
					assert(itBBI != boundaryInfo.boundaryInfoOfBlock.end());
					
					map<int, BoundaryBlockPack *>::const_iterator itBBPack = m_mapID2Pack.find(blockID);
					
					const bool bFound = itBBPack != m_mapID2Pack.end();
					
					if (bFound)
					{
						map<int, FootPrint>::const_iterator itFootPrint= m_mapID2FootPrint.find(blockID);
						assert(itFootPrint != m_mapID2FootPrint.end());
						
						vector<int>& dependentBlockIDs = itBBI->second->dependentBlockIDs;
						const FootPrint& f = itFootPrint->second;
						
						bool bSameFootPrint = true;
						for(int j=0; j<dependentBlockIDs.size() && bSameFootPrint; j++)
							bSameFootPrint &= (dependentBlockIDs[j] == f[j]);
						
						if (bSameFootPrint)
							result[i]  = itBBPack->second;
						else
						{
							_erase(blockID, itBBPack->second);
							result[i]  = _translate(blockID, *itBBI->second, boundaryInfo.stencil_start, boundaryInfo.stencil_end);
						}
					}
					else
						result[i]  = _translate(blockID, *itBBI->second, boundaryInfo.stencil_start, boundaryInfo.stencil_end);
				}
			}
			
			
			//4.
			//look for cached values
			{
				const int n = result.nPacks;
				
				vector<void *> _id(n);
				
				for(int i=0; i<n; i++)
					_id[i] = result[i];
				
				bool bCacheHit = false;
				BoundaryBlockPack **& ptrGPU = m_cacheBoundaryPacksGPU.get(_id, bCacheHit);
				
				if (bCacheHit)
				{
					//printf("CACHE HIT in BBPackHandler_CUDA::translate()\n");
					result.ptrPacks = ptrGPU;
					result.ptrBoundaryPacksGPU = NULL;
					result.bFinalized = true;
				}
				else
				{
					if (ptrGPU != NULL) CUDA_SAFE_CALL(cudaFree((void*)ptrGPU));
					
					result.finalize();
					
					ptrGPU = result.ptrPacks;
					result.ptrBoundaryPacksGPU = NULL; 
				}
			}
					
			
			{
				int s=0;
				for(map<int, BoundaryBlockPack *>::const_iterator it =m_mapID2Pack.begin(); it!=m_mapID2Pack.end(); it++ )
				{
				//	printf("size: %d\n", m_mapGPU2CPU[it->second]->byteSize);
					s+= m_mapGPU2CPU[it->second]->byteSize;
				}
				
				printf("stored packs total size: %f KB\n", s/1024.);
			
			}
			return result;
		}
		
		BBPackHandler_CUDA(): remoteCollection(NULL), m_result(NULL),
			m_mapID2FootPrint(), m_mapID2Pack(), m_mapGPU2CPU(), m_cacheBoundaryPacksGPU()
		{
			m_cacheBoundaryPacksGPU.setup(10);
		}
		
	private:
		
		Cache< vector<void*>, BoundaryBlockPack **> m_cacheBoundaryPacksGPU;
		
		map<int, FootPrint> m_mapID2FootPrint;
		map<int, BoundaryBlockPack *> m_mapID2Pack; // Block ID -> Pack in the GPU
		map<BoundaryBlockPack *, BoundaryBlockPack *> m_mapGPU2CPU; // GPU Pack -> CPU pack
		
		ArrayOfBoundaryPack * m_result;
		
		const Addressable_CUDA * remoteCollection;
		
		//forbidden
		BBPackHandler_CUDA(const BBPackHandler_CUDA& b):
		m_mapID2FootPrint(), m_mapID2Pack(), m_mapGPU2CPU(),
		m_result(NULL), remoteCollection(NULL) { abort(); }
		BBPackHandler_CUDA& operator=(const BBPackHandler_CUDA& b){abort(); return *this;}
	};
}
