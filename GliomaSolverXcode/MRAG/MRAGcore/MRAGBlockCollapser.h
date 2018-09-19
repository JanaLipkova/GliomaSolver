/*
 *  BlockCollapser.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 5/21/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */  
#include "MRAGBlockLab.h"
#include "MRAGCompressionPlan.h"
#pragma once
namespace MRAG
{
	template <typename WaveletsType, typename BlockCollectionType>
	class BlockCollapser
	{
		typedef typename BlockCollectionType::BlockType BlockType;
		typedef typename BlockType::ElementType ElementType;
	protected:
		BlockLab<BlockType> m_blockLab;
		
		
		template <typename BlockLabType>
		inline ElementType _filter(BlockLabType& lab, const int dx, const int dy, const int dz) const
		{
			const int s[3] = {
				!BlockType::shouldProcessDirectionY?0 : WaveletsType::HaSupport[0],
				!BlockType::shouldProcessDirectionY?0 : WaveletsType::HaSupport[0],
				!BlockType::shouldProcessDirectionZ?0 : WaveletsType::HaSupport[0]
			};
			
			const int e[3] = {
				!BlockType::shouldProcessDirectionX?1 : WaveletsType::HaSupport[1],
				!BlockType::shouldProcessDirectionY?1 : WaveletsType::HaSupport[1],
				!BlockType::shouldProcessDirectionZ?1 : WaveletsType::HaSupport[1]
			};
			
			const int o[3] = {
				!BlockType::shouldProcessDirectionX?0:dx*2,
				!BlockType::shouldProcessDirectionY?0:dy*2,
				!BlockType::shouldProcessDirectionZ?0:dz*2
			};
			
			ElementType result = ElementType();
			for(int iz=s[2]; iz<e[2]; iz++)
			{
				const double wZ = (!BlockType::shouldProcessDirectionZ)?1: WaveletsType::getHa(iz);
				for(int iy=s[1]; iy<e[1]; iy++)
				{
					const double wY = (!BlockType::shouldProcessDirectionY)?1: WaveletsType::getHa(iy);
					const double wZwY = wZ*wY;
					
					for(int ix=s[0]; ix<e[0]; ix++)//double)
						result += lab.read(o[0] - ix, o[1] - iy, o[2] - iz)* ( wZwY*WaveletsType::getHa(ix) );
				}
			}
			
			return result;
		}
		
		template <typename BlockLabType>
		void _collapseBlocks(BlockLabType& lab, CompressionPlan::Collapse& collapseInfo, vector<BlockType *>& source, BlockType & dest)
		{
			const int nSourceX = BlockType::sizeX/2;
			const int nSourceY = !BlockType::shouldProcessDirectionY? 1 : BlockType::sizeY/2;
			const int nSourceZ = !BlockType::shouldProcessDirectionZ? 1 : BlockType::sizeZ/2;
			
			for(int iSource =0; iSource<collapseInfo.nSources; iSource++)
			{
				//const I3 idx = collapseInfo.
				const int *abs_idx = collapseInfo.source_index[iSource].i;
				const int dest_idx[3] = {abs_idx[0]/2, abs_idx[1]/2, abs_idx[2]/2};
				const int idx[3] = { abs_idx[0] - 2*dest_idx[0], abs_idx[1] - 2*dest_idx[1], abs_idx[2] - 2*dest_idx[2]};
				const int s[3] = {idx[0]*nSourceX, idx[1]*nSourceY, idx[2]*nSourceZ};
				const int e[3] = {s[0] + nSourceX, s[1] + nSourceY , s[2] + nSourceZ};
				
				BlockInfo info(collapseInfo.source_blockIDs[iSource], abs_idx, collapseInfo.source_level );
				lab.load(info);
	
				for(int iz=s[2]; iz<e[2]; iz++)
					for(int iy=s[1]; iy<e[1]; iy++)
						for(int ix=s[0]; ix<e[0]; ix++)
							dest(ix, iy, iz) = _filter(lab, ix-s[0], iy-s[1], iz-s[2]);
			}
		}
		
	public:
		BlockCollapser(): m_blockLab() {}
	
		virtual CompressionResult collapse(BlockCollection<BlockType>& collection, BoundaryInfo& boundaryInfo, const CompressionPlan& compressionPlan, vector<BlockCollapseInfo>& vNewIDs)
		{
			return _collapse(m_blockLab, collection, boundaryInfo, compressionPlan, vNewIDs);
		}
		
		template<typename BlockLabType>
		CompressionResult _collapse(BlockLabType& lab, BlockCollection<BlockType>& collection, BoundaryInfo& boundaryInfo, const CompressionPlan& compressionPlan, vector<BlockCollapseInfo>& vNewIDs)
		{
			//0. setup	
			//1. allocate the new blocks
			//2. fill them
			//3. remove the old blocks
			//4. output the result
			
			//0.
			{
				const int stencilStart[3] = {
					!BlockType::shouldProcessDirectionX?0: 1 - WaveletsType::HaSupport[1],
					!BlockType::shouldProcessDirectionY?0: 1 - WaveletsType::HaSupport[1],
					!BlockType::shouldProcessDirectionZ?0: 1 - WaveletsType::HaSupport[1]
				};
				
				const int stencilEnd[3] = {
					!BlockType::shouldProcessDirectionX?1: - WaveletsType::HaSupport[0],
					!BlockType::shouldProcessDirectionY?1: - WaveletsType::HaSupport[0],
					!BlockType::shouldProcessDirectionZ?1: - WaveletsType::HaSupport[0]
				};	
				
				lab.prepare(collection, boundaryInfo, stencilStart, stencilEnd);
			}
			
			//1.
			vector<int> vIDs = collection.create(compressionPlan.nCollapses);
			
			//2.
			int nCollapsed = 0;
			for(int c=0; c<compressionPlan.nCollapses; c++)
			{
				CompressionPlan::Collapse& collapseInfo = compressionPlan.vCollapseArray[c];
				
				vector<BlockType *>source;
				for(int s=0; s<collapseInfo.nSources; s++)
					source.push_back(&collection.lock(collapseInfo.source_blockIDs[s]));
				
				_collapseBlocks(lab, collapseInfo, source, collection.lock(vIDs[c]));
				collection.release(vIDs[c]);

				for(int s=0; s<collapseInfo.nSources; s++)
					collection.release(collapseInfo.source_blockIDs[s]);

				nCollapsed += source.size();
			}
			
			//3.
			for(int c=0; c<compressionPlan.nCollapses; c++)
			{
				CompressionPlan::Collapse& collapseInfo = compressionPlan.vCollapseArray[c];
				
				for(int s=0; s<collapseInfo.nSources; s++)
					collection.erase(collapseInfo.source_blockIDs[s]);
			}
			
			//4.
			for(int c=0; c<compressionPlan.nCollapses; c++)
			{
				BlockCollapseInfo blockinfo = {c, vIDs[c]};
				vNewIDs.push_back(blockinfo);
			}
			
			CompressionResult result = {nCollapsed, vIDs.size()};
			return result;
		}
	};
	
	template <typename WaveletsType, typename BlockCollectionType, typename BlockLabType>
	class BlockCollapserGenericLab : public BlockCollapser<WaveletsType, BlockCollectionType>
	{
		typedef typename BlockCollectionType::BlockType BlockType;
		typedef typename BlockType::ElementType ElementType;
	protected:
		
		BlockLabType m_myGenericLab;
	public:
		BlockCollapserGenericLab(): m_myGenericLab(), BlockCollapser<WaveletsType, BlockCollectionType>() {} 
		
		CompressionResult collapse(BlockCollectionType& collection, BoundaryInfo& boundaryInfo, const CompressionPlan& compressionPlan, vector<BlockCollapseInfo>& vNewIDs)
		{
			return _collapse(m_myGenericLab, collection, boundaryInfo, compressionPlan, vNewIDs);
		}
	};
}
