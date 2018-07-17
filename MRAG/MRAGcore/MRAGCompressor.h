/*
 *  MRAGCompressor.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 4/24/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include "MRAGGridNode.h"
#include "MRAGCompressionPlan.h"

namespace MRAG
{
	struct GridNodeCollapseInfo
	{
		int collapseID;
		GridNode * node;
	};
	
	class Compressor
	{
		const int m_nMaxLevelJump;
		
		template<typename HierarchyType, typename NeighborhoodType>
		GridNode* _findRealNode(const HierarchyType& hierarchy, NeighborhoodType& neighborhood, GridNode* node)
		{
			//this function finds the "true" grid node when the input grid node is a ghost
			
			typename NeighborhoodType::const_iterator itNode = neighborhood.find(const_cast<GridNode*>(node));
			const bool bNotGhost = (itNode != neighborhood.end());
			
			if (bNotGhost) 
				return node;
			else 
			{
				//look for the parent's children. the one with the same blockID is the one.				
				typename HierarchyType::const_iterator itParent = hierarchy.find(node->parent);
				
				assert(itParent != hierarchy.end());
				
				vector<GridNode *> children = itParent->second;
				
				int iChild;
				for(iChild =0; iChild<children.size(); iChild++)
					if (children[iChild]->blockID == node->blockID) break;
				
				assert(children[iChild]->blockID == node->blockID);
				
				return children[iChild];
				
			}
		}
		
	public:
		const int getMaxLevelJump() const
		{
			return m_nMaxLevelJump;
		}
		
		Compressor(int nMaxLevelJump = 1):
		m_nMaxLevelJump(nMaxLevelJump)
		{
		}
		
		template<typename HierarchyType, typename NeighborhoodType>
		CompressionPlan * createPlan(const HierarchyType& hierarchy, NeighborhoodType& neighborhood, vector<GridNodeCollapseInfo>& vToCollapse)
		{
			//1. check the nodes whose children are full and ready to be compressed
			//2. generate a compression plan
			
			//1.
			for(typename HierarchyType::const_iterator it=hierarchy.begin(); it!= hierarchy.end(); it++)
			{
				GridNode * node = it->first;
				
				if (node == NULL || node->isEmpty==false || it->second[0]->isEmpty == true) continue;
				
				const vector<GridNode *>& children = it->second;
				bool bEveryChildrenIsFull = true;
				bool bShouldCompress = true;
				
				
				for(int i=0; i<children.size(); i++)
				{
					bEveryChildrenIsFull &= children[i]->isEmpty==false;
					bShouldCompress &= children[i]->shouldBeCompressed==true;
					
					if  (!bEveryChildrenIsFull || !bShouldCompress) break;
				}
				
				if (bEveryChildrenIsFull && bShouldCompress)
				{
					bool bAccettableLevelJump = true;
					
					const int iCandidateLevel = it->first->level;
					
					for(int i=0; i<children.size(); i++)
					{
						const vector<GridNode*>& neighbors =  neighborhood[children[i]];
						
						for(int j = 0; j<neighbors.size(); j++)
							bAccettableLevelJump &= (abs(neighbors[j]->level - iCandidateLevel)<=m_nMaxLevelJump);// && _findRealNode(hierarchy, neighborhood, neighbors[j])->shouldBeCompressed;
						
						if (!bAccettableLevelJump) break;
					}
						
					if (bAccettableLevelJump)
					{
						GridNodeCollapseInfo collapseinfo = {vToCollapse.size(), node};
						vToCollapse.push_back(collapseinfo);
					}
				}
			}
			
			//2.
			CompressionPlan * plan = new CompressionPlan(vToCollapse.size());
			const int nCollapses = plan->nCollapses;
	
			for(int i=0; i<nCollapses; i++)
			{
				CompressionPlan::Collapse& collapse = plan->vCollapseArray[i];
				typename HierarchyType::const_iterator itChildren = hierarchy.find(vToCollapse[i].node);
				assert(itChildren!=hierarchy.end());
				
				const vector<GridNode *>& children = itChildren->second;//hierarchy[vToCollapse[i].node];
				
				collapse.collapseID = i;
				collapse.nSources = children.size();
				collapse.source_level = children[0]->level;
				for(int s=0;s<collapse.nSources; s++)
				{
					const GridNode * node = children[s];
					
					assert(node->isEmpty == false);
					assert(node->level == collapse.source_level);
					collapse.source_blockIDs[s] = node->blockID;
					collapse.source_index[s].i[0] = node->index[0];
					collapse.source_index[s].i[1] = node->index[1];
					collapse.source_index[s].i[2] = node->index[2];
				}
			}
			
			return plan;
		}
	};
}
