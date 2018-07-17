/*
 *  MRAGRefiner_SpaceExtension.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 9/9/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */


#pragma once

#include "MRAGcore/MRAGRefiner.h"

namespace MRAG
{
	class Refiner_SpaceExtension: public Refiner
	{
		int m_nMaxLevel;
	public:
		Refiner_SpaceExtension(int nMaxLevelJump=1, int nMaxLevel=-1): Refiner(nMaxLevelJump), m_nMaxLevel(nMaxLevel)
		{
		}
		
		virtual RefinementPlan* createPlan(const HierarchyType& hierarchy, const NeighborhoodType& neighborhood, const bool vProcessingDirections[3], vector<NodeToRefine>& vRefinements)
		{
			//this isessentially a wrapper
			//1. iterate over the blocks
			//2. if a block should be refined, check its neighbor blocks
			//3. if a neighbor block is coarser, refine that too
			//4. prevent that a block on the finest level is refined any further
			//5. ready to generate a refinement plan
			
			//1.
			for(NeighborhoodType::const_iterator itS= neighborhood.begin(); itS!=neighborhood.end(); itS++)
			{
				//2.
				
				const vector<GridNode *>& n = itS->second;
				
				for(vector<GridNode *>::const_iterator itD = n.begin(); itD != n.end(); itD++)
				{
					//3.
					if (neighborhood.find(*itD) != neighborhood.end()) // neighbor block not a ghost
						if ((*itD)->level > itS->first->level) // neighbor block coarser
						{
							itS->first->shouldBeRefined = true;
							break;
						}
				}
			}
			
			//4.
			if (m_nMaxLevel>=0)
			{
				for(NeighborhoodType::const_iterator itS= neighborhood.begin(); itS!=neighborhood.end(); itS++)
				{
					assert(itS->first->level <= m_nMaxLevel);
					
					if(itS->first->level == m_nMaxLevel)
						itS->first->shouldBeRefined = false;
				}
			}
			
			//5.
			return Refiner::createPlan(hierarchy, neighborhood, vProcessingDirections, vRefinements);
		}
	};
}
