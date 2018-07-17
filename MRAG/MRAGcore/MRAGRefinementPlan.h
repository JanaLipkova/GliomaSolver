/*
 *  MRAGRefinementPlan.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 4/24/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "MRAGGridNode.h"

#include <list>
#include <map>
using namespace std;

namespace MRAG
{
	struct RefinementPlanNode
	{
		int index[3];
		int relative_index[3];
		int level;
	};
	
	struct RefinementReport
	{
		RefinementReport(): blockID(0), childrenIDs() {} 
		int blockID;
		vector<int> childrenIDs;
	};
	
	struct SingleRefinement
	{
		//int blockID;
		BlockInfo block_info_source;
		vector<RefinementPlanNode*> children;
		//vector<int> vBlockIDsToRemove;
		
		SingleRefinement(): block_info_source(), children() {}
		
		SingleRefinement( SingleRefinement& c):block_info_source(), children() {abort();}
		
		SingleRefinement(const SingleRefinement& c):block_info_source(), children() { abort(); }
		
		RefinementPlanNode& createEntry()
		{
			RefinementPlanNode * result = new RefinementPlanNode;
			
			children.push_back(result);

			return *result;
		}
		
		~SingleRefinement()
		{
			for(int i=0; i<children.size(); i++)
				delete children[i];
		}
	};
	
	struct RefinementPlan
	{
		vector<SingleRefinement*> refinements;

		int nNewBlocks;
		
		SingleRefinement * createEntry()
		{
			SingleRefinement* c = new SingleRefinement;
			
			refinements.push_back(c);
			
			return c;
		}
		
		RefinementPlan():refinements(), nNewBlocks(0){}
		
		~RefinementPlan()
		{
			for(int i=0; i<refinements.size(); i++)
				delete refinements[i];
		}
	};
	
	struct RefinementResult
	{
	private:
		bool bFailed;
		
	public:
		int nChildrenBlocks;
		int nCollapsedParentBlocks;
		
		const bool hasFailed() const { return bFailed; }
		
		RefinementResult(const bool bFailed_=false):nChildrenBlocks(0), nCollapsedParentBlocks(0), bFailed(bFailed_){}
		
		RefinementResult operator += (const RefinementResult& r)
		{
			nChildrenBlocks += r.nChildrenBlocks;
			nCollapsedParentBlocks += r.nCollapsedParentBlocks;
			bFailed = bFailed || r.bFailed;
			
			return *this;
		}
	};
	
	struct NodeToRefine
	{
		const GridNode * node;
		int refinementID;
		NodeToRefine(const GridNode * n,int ID): node(n), refinementID(ID){}
	};
}
