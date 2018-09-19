/*
 *  MRAGRefiner.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 4/24/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <map>
#include <vector>
#include <set>
using namespace std;

#include "MRAGRefinementPlan.h"
#include "MRAGGridNode.h"

namespace MRAG
{
	class Refiner
	{
	protected:
		const int m_nMaxLevelJump;
		
		template <typename C, typename K> const bool _isFound(const C& c, const K& k) const
		{
			typename C::const_iterator it = c.find(k);
			return (it != c.end() );
		}
		
	public:
		
		const int getMaxLevelJump() const
		{
			return m_nMaxLevelJump;
		}
		
		Refiner(int nMaxLevelJump=1): m_nMaxLevelJump(nMaxLevelJump) 
		{
		}
		
		virtual RefinementPlan* createPlan(const HierarchyType& hierarchy, const NeighborhoodType& neighborhood, const bool vProcessingDirections[3], vector<NodeToRefine>& vRefinements)
		{
			//1. create a fresh plan
			//2. iterate over the leafs: if a leaf should be refined put it in a victim list.
			//3. chain reaction on the victim list to guarantee the max resolution jump
			//4. create a plan from the victim list
			
			//1.
			typedef set<const GridNode *> VictimSet;
			vRefinements.clear();
			RefinementPlan * plan = new RefinementPlan;
			
			//2.
			VictimSet victims;
			
			for(HierarchyType::const_iterator it=hierarchy.begin(); it!=hierarchy.end(); it++)
			{
				const GridNode * node = it->first;
				
				if (node != NULL && !node->isEmpty && node->shouldBeRefined) victims.insert(node);
			}
			
			//3.
			{
				//3.a allocation
				//3.b spot new victims from the old ones
				//3.c if there are no new victim anymore, this is the end of pt 3.
				//3.d put the new victims into the 'victims' set
				
				//3.a
				VictimSet tmp1, tmp2;
				VictimSet& new_victims = tmp1;
				VictimSet& old_victims = tmp2;
				vector<const GridNode *> ghost_nodes;
		
				for(VictimSet::iterator victim = victims.begin();  victim!=victims.end(); victim++)
						old_victims.insert(*victim);
				
				do
				{
					//3.b
					{
						new_victims.clear();
						for(VictimSet::iterator victim = old_victims.begin();  victim!=old_victims.end(); victim++)
						{
							assert((*victim)->isEmpty == false);
							
							const int new_victim_level = (*victim)->level + 1;
							
							NeighborhoodType::const_iterator itNode = neighborhood.find(const_cast<GridNode*>(*victim));
							const bool bFound = (itNode != neighborhood.end());
							
							const vector<GridNode*> * pneighbors = NULL;
							
							if (!bFound) //it's a ghost
							{
								//look for the parent's children. the one with the same blockID is the one.
								//remove then the ghosts from the vicims
								
								HierarchyType::const_iterator itParent = hierarchy.find((*victim)->parent);
								
								assert(itParent != hierarchy.end());
								
								vector<GridNode *> children = itParent->second;
						
								int iChild;
								for(iChild =0; iChild<children.size(); iChild++)
									if (children[iChild]->blockID == (*victim)->blockID) break;
								
								assert(children[iChild]->blockID == (*victim)->blockID);
								
								NeighborhoodType::const_iterator itRealNode = neighborhood.find(children[iChild]);
								
								ghost_nodes.push_back(*victim);
								victims.insert(children[iChild]);
								
								pneighbors = &itRealNode->second;
							}
							else
								pneighbors = &itNode->second;
							 
							const vector<GridNode*>& neighbors = *pneighbors;
							
							for (vector<GridNode*>::const_iterator n=neighbors.begin(); n!=neighbors.end(); n++)
							{
								const GridNode& node = **n;
								const bool bJumpInLevelTooHigh = (abs(node.level - new_victim_level) > m_nMaxLevelJump);
								const bool bNewVictim = !_isFound(victims, &node);
								
								if (bJumpInLevelTooHigh && bNewVictim)
								{
									new_victims.insert(&node);
									victims.insert(&node);
								}
							}
						}
					}
					
					//3.c
					if (new_victims.size() == 0) break;
					
					//3.d
					std::swap(new_victims, old_victims);
				
				} while(true);
				
				
				for(int i=0; i<ghost_nodes.size(); i++)
					victims.erase(ghost_nodes[i]);
			}
			
			//4.
			{
				int nNewChildren = 0;
				for(VictimSet::const_iterator victim = victims.begin();  victim!=victims.end(); victim++)
				{
					const GridNode * node = *victim;
					
					SingleRefinement& singleRefinement = *(plan->createEntry());
					singleRefinement.block_info_source.blockID = node->blockID;
					singleRefinement.block_info_source.index[0] = node->index[0];
					singleRefinement.block_info_source.index[1] = node->index[1];
					singleRefinement.block_info_source.index[2] = node->index[2];
					singleRefinement.block_info_source.level = node->level;

					for(int i=0; i<8; i++)
					{
						const int offset[3] = {i&1, (i>>1)&1, (i>>2)&1};
						
						if ((offset[0]==0 || vProcessingDirections[0]) && 
							(offset[1]==0 || vProcessingDirections[1]) &&
							(offset[2]==0 || vProcessingDirections[2]))
						{
							RefinementPlanNode& child = singleRefinement.createEntry();
							child.level = node->level + 1;
							
							child.index[0] = node->index[0]*2 + offset[0];
							child.index[1] = node->index[1]*2 + offset[1];
							child.index[2] = node->index[2]*2 + offset[2];
							
							child.relative_index[0] = offset[0];
							child.relative_index[1] = offset[1];
							child.relative_index[2] = offset[2];
							
							nNewChildren++;
						}
					}
					
					vRefinements.push_back(NodeToRefine(node, vRefinements.size()));
				}
				
				plan->nNewBlocks = nNewChildren;
			}
			
		return plan;
	 }
	};
}
