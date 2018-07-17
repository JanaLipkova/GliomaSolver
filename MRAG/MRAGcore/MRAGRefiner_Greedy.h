/*
 *  MRAGRefiner_Greedy.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/21/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "MRAGRefiner.h"

namespace MRAG
{

	class Refiner_Greedy: public Refiner
	{
		
		
	public:
		Refiner_Greedy(int nMaxLevelJump=1): Refiner(nMaxLevelJump)
		{
		}
		
		virtual RefinementPlan* createPlan(const HierarchyType& hierarchy, const NeighborhoodType& neighborhood, const bool vProcessingDirections[3], vector<NodeToRefine>& vRefinements)
		{
			//1. create a fresh plan
			//2. iterate over the leafs: if a leaf should be refined (and the jump resolution is not too high), create new children leaves
			//3. return it
			
			//1.
			vRefinements.clear();
			RefinementPlan * plan = new RefinementPlan;
			
			//2.
			int nNewChildren = 0;
			for(HierarchyType::const_iterator it=hierarchy.begin(); it!=hierarchy.end(); it++)
			{
				GridNode * node = it->first;
				
				if(node == NULL || node->isEmpty) continue;
				
				bool bJumpInLevelTooHigh = false;
				{
					const int candidate_level = node->level + 1;
					NeighborhoodType::const_iterator itNode = neighborhood.find(node);
					const vector<GridNode*>& neighbors = itNode->second;
					
					for (vector<GridNode*>::const_iterator n=neighbors.begin(); n!=neighbors.end() && bJumpInLevelTooHigh==false; n++)
						bJumpInLevelTooHigh |= (abs((*n)->level - candidate_level) > m_nMaxLevelJump);
				}
				
				if(node->shouldBeRefined && !bJumpInLevelTooHigh)
				{
					SingleRefinement singleRefinement;
					singleRefinement.blockID = node->blockID;
					vector<RefinementPlanNode *>& vChildren = singleRefinement.children;
					
					for(int i=0; i<8; i++)
					{
						const int offset[3] = {i&1, (i>>1)&1, (i>>2)&1};
						
						if ((offset[0]==0 || vProcessingDirections[0]) && 
							(offset[1]==0 || vProcessingDirections[1]) &&
							(offset[2]==0 || vProcessingDirections[2]))
						{
							RefinementPlanNode * child = new RefinementPlanNode;
							child->level = node->level + 1;
							
							child->index[0] = node->index[0]*2 + offset[0];
							child->index[1] = node->index[1]*2 + offset[1];
							child->index[2] = node->index[2]*2 + offset[2];
							
							child->relative_index[0] = offset[0];
							child->relative_index[1] = offset[1];
							child->relative_index[2] = offset[2];
							
							vChildren.push_back(child);
							
							nNewChildren++;
						}
					}
					
					plan->refinements.push_back(singleRefinement);
					
					vRefinements.push_back(NodeToRefine(node, vRefinements.size()));
				}
			}
			plan->nNewBlocks = nNewChildren;
			
			//3.
			return plan;
		}
	};
}
