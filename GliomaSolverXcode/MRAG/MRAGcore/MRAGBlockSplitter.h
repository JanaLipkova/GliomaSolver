/*
 *  MRAGBlockSplitter.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 5/27/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include "MRAGRefinementPlan.h"

namespace MRAG
{	
	template <typename WaveletsType, typename BlockCollectionType>
	class BlockSplitter
		{
			typedef typename BlockCollectionType::BlockType BlockType;
			typedef typename BlockType::ElementType ElementType;
		protected:
			BlockLab<BlockType> m_blockLab;
			
			template<typename BlockLabType> 
			void _splitBlocks(BlockLabType& lab, BlockInfo& info, const BoundaryInfoBlock& bbinfo,  const BlockType &source,  vector<BlockType *>&dest, vector<RefinementPlanNode*> refinementinfo)
			{
				lab.load(info);
				
				const int nX = BlockType::sizeX;
				const int nY = BlockType::sizeY;
				const int nZ = BlockType::sizeZ;
				
				assert(refinementinfo.size()<=8);
				BlockType* destMat[2][2][2] = {0,0,0,0, 0,0,0,0};
				
				for(int i=0; i<dest.size(); i++)
				{
					//ElementType e;
					//e = 1.0;
					*dest[i] = ElementType();
					
					const int  * idx = refinementinfo[i]->relative_index;
					destMat[idx[0]][idx[1]][idx[2]] = dest[i];
					
					assert(idx[0] == 0 || idx[0] == 1);
					assert(idx[1] == 0 || idx[1] == 1);
					assert(idx[2] == 0 || idx[2] == 1);
				}
				
				const int s[3] = {
					(int)ceil(WaveletsType::HsSupport[0]*0.5), 
					BlockType::shouldProcessDirectionY? (int)ceil(WaveletsType::HsSupport[0]*0.5) : 0,
					BlockType::shouldProcessDirectionZ? (int)ceil(WaveletsType::HsSupport[0]*0.5) : 0
				};
				
				const int e[3] = {
					(1 + (int)floor((WaveletsType::HsSupport[1]-1)*0.5 + nX)),
					BlockType::shouldProcessDirectionY? (1 + (int)floor((WaveletsType::HsSupport[1]-1)*0.5 + nY)) : 1,
					BlockType::shouldProcessDirectionZ? (1 + (int)floor((WaveletsType::HsSupport[1]-1)*0.5 + nZ)) : 1
				};
				
				for(int sz=s[2]; sz<e[2]; sz++)
					for(int sy=s[1]; sy<e[1]; sy++)
						for(int sx=s[0]; sx<e[0]; sx++)
						{
							const int influence_start[3] = {
								!BlockType::shouldProcessDirectionX?0 : std::max(0,(2*sx + 1 - WaveletsType::HsSupport[1])),
								!BlockType::shouldProcessDirectionY?0 : std::max(0,(2*sy + 1 - WaveletsType::HsSupport[1])),
								!BlockType::shouldProcessDirectionZ?0 : std::max(0,(2*sz + 1 - WaveletsType::HsSupport[1])),
							};
							
							const int influence_end[3] = {
								!BlockType::shouldProcessDirectionX?1 : std::min(2*nX,(2*sx + 1 - WaveletsType::HsSupport[0])),
								!BlockType::shouldProcessDirectionY?1 : std::min(2*nY,(2*sy + 1 - WaveletsType::HsSupport[0])),
								!BlockType::shouldProcessDirectionZ?1 : std::min(2*nZ,(2*sz + 1 - WaveletsType::HsSupport[0])),
							};
														
							for(int dz=influence_start[2]; dz<influence_end[2];  dz++)
							{
								const double wZ = (!BlockType::shouldProcessDirectionZ)? 1: WaveletsType::getHs(2*sz - dz);
								
								for(int dy=influence_start[1]; dy<influence_end[1];  dy++)
								{
									const double wY = (!BlockType::shouldProcessDirectionY)?1: WaveletsType::getHs(2*sy - dy);
									const double wZwY = wZ*wY;
									
									for(int dx=influence_start[0]; dx<influence_end[0];  dx++)
									{
										const int idx[3] = {
											(!BlockType::shouldProcessDirectionX)?0 : dx/BlockType::sizeX,
											(!BlockType::shouldProcessDirectionY)?0 : dy/BlockType::sizeY,
											(!BlockType::shouldProcessDirectionZ)?0 : dz/BlockType::sizeZ
										};
										
										(*destMat[idx[0]][idx[1]][idx[2]])(dx-idx[0]*nX,dy-idx[1]*nY,dz-idx[2]*nZ) += lab.read(sx,sy,sz)*( wZwY*WaveletsType::getHs(2*sx - dx) );
									}
								}
							}
						}
			}
			
		public:
			
			BlockSplitter(): m_blockLab() {} 
			
			virtual RefinementResult split(BlockCollectionType& collection, BoundaryInfo& boundaryInfo, const RefinementPlan& refinementPlan, vector<RefinementReport>& vNewIDs)
			{
				return _split(m_blockLab, collection, boundaryInfo, refinementPlan, vNewIDs);
			}
			
			template<typename BlockLabType>
			RefinementResult _split(BlockLabType& lab, BlockCollectionType& collection, BoundaryInfo& boundaryInfo, const RefinementPlan& refinementPlan, vector<RefinementReport>& vNewIDs)
			{
				//0. setup
				//1. allocate the new blocks
				//2. fill them
				//3. remove the old blocks
				//4. output the result
				
				//0.
				{
					const int stencilStart[3] = {
						(int)ceil(WaveletsType::HsSupport[0]*0.5), 
						BlockType::shouldProcessDirectionY? (int)ceil(WaveletsType::HsSupport[0]*0.5) : 0,
						BlockType::shouldProcessDirectionZ? (int)ceil(WaveletsType::HsSupport[0]*0.5) : 0
					};
					
					const int stencilEnd[3] = {
						(2 + (int)floor((WaveletsType::HsSupport[1]-1)*0.5 )),
						BlockType::shouldProcessDirectionY? (2 + (int)floor((WaveletsType::HsSupport[1]-1)*0.5 )) : 1,
						BlockType::shouldProcessDirectionZ? (2 + (int)floor((WaveletsType::HsSupport[1]-1)*0.5 )) : 1					
					};
					
					lab.prepare(collection, boundaryInfo, stencilStart, stencilEnd);
				}
				
				//1.
				//vector<int> vIDs;  
				//vIDs.reserve(refinementPlan.nNewBlocks);
				vector<int> vIDs = collection.create(refinementPlan.nNewBlocks);
				
				//2.
				int iCurrentBlock = 0;
				const int nRefinements = refinementPlan.refinements.size();
				for(int r=0; r<nRefinements; r++)
				{
					SingleRefinement& refinement =  *refinementPlan.refinements[r];
					
					const BlockType& source = collection.lock(refinement.block_info_source.blockID);//collection[refinement.blockID];
					
					//{
					//	vector<int> vLocalIDs = collection.create(refinement.children.size());
					//	for(int i=0; i<vLocalIDs.size(); i++)
					//		vIDs.push_back(vLocalIDs[i]);
					//}
					
					vector<BlockType *> destinations;
					vector<int> lockedChildren;
					lockedChildren.reserve(refinement.children.size());
					
					for(int d=0; d<refinement.children.size(); d++,iCurrentBlock++)  
					{
						//printf("******ID : %d\n", vIDs[iCurrentBlock]);
						collection.lock(vIDs[iCurrentBlock]);
						lockedChildren.push_back(vIDs[iCurrentBlock]);
						destinations.push_back(&collection[vIDs[iCurrentBlock]]);
					}
					
					const BoundaryInfoBlock & bbinfo = *boundaryInfo.boundaryInfoOfBlock[refinement.block_info_source.blockID];
					
					const int idx_fake[3] = {0,0,0};
					BlockInfo info(refinement.block_info_source.blockID, idx_fake, refinement.children[0]->level-1);
					
					_splitBlocks(lab, refinement.block_info_source, bbinfo, source, destinations, refinement.children);
					 
					collection.release(refinement.block_info_source.blockID);
					for(vector<int>::const_iterator it = lockedChildren.begin(); it!=lockedChildren.end(); it++)
						collection.release(*it);
				}
				
				//3.
				for(int r=0; r<nRefinements; r++)
				{
					SingleRefinement& refinement =  *refinementPlan.refinements[r];
					collection.erase(refinement.block_info_source.blockID);
				}
				
				//4.
				iCurrentBlock = 0;
				for(int r=0; r<nRefinements; r++)
				{
					RefinementReport refinementInfo;
					
					refinementInfo.blockID = refinementPlan.refinements[r]->block_info_source.blockID;
					
					for(int d=0; d<refinementPlan.refinements[r]->children.size(); d++)
						refinementInfo.childrenIDs.push_back(vIDs[iCurrentBlock++]);
					
					vNewIDs.push_back(refinementInfo);
				}
				
				RefinementResult result;
				result.nChildrenBlocks = vIDs.size();
				result.nCollapsedParentBlocks = nRefinements;
				
				return result;
			}
		};
	
	template <typename WaveletsType, typename BlockCollectionType, typename BlockLabType>
	class BlockSplitterGenericLab : public BlockSplitter<WaveletsType, BlockCollectionType>
	{
		typedef typename BlockCollectionType::BlockType BlockType;
		typedef typename BlockType::ElementType ElementType;
	protected:
		
		BlockLabType m_myGenericLab;
		
	public:
		BlockSplitterGenericLab(): m_myGenericLab(), BlockSplitter<WaveletsType, BlockCollectionType>(){}
		
		RefinementResult split(BlockCollectionType& collection, BoundaryInfo& boundaryInfo, const RefinementPlan& refinementPlan, vector<RefinementReport>& vNewIDs)
		{
			return _split(m_myGenericLab, collection, boundaryInfo, refinementPlan, vNewIDs);
		}
	};
}
