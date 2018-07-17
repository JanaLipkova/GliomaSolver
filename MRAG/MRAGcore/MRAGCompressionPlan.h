/*
 *  MRAGCompressionPlan.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 5/20/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#include "MRAGCommon.h"

#pragma once
namespace MRAG
{
	struct CompressionPlan
	{
		struct Collapse
		{
			int collapseID;
			int source_blockIDs[8];
			I3 source_index[8];
			int source_level;
			int nSources;
			
			Collapse(): collapseID(0), source_level(-1), nSources(0) {} 
		};
		
		int nCollapses;
		Collapse* vCollapseArray;
		
		CompressionPlan(int n): nCollapses(n), vCollapseArray(NULL)
		{
			vCollapseArray = new Collapse[nCollapses];
		}
		
		~CompressionPlan()
		{
			delete [] vCollapseArray;
		}
		
	private:
		CompressionPlan(const CompressionPlan&):nCollapses(0), vCollapseArray(NULL)  { abort(); }
		CompressionPlan& operator=(const CompressionPlan&) { abort(); return *this;}
	};
	
	struct CompressionResult
	{
		int nCollapsedBlocks;
		int nNewBlocks;
	};

	struct BlockCollapseInfo
	{
		int collapseID;
		int newBlockID;
	};

}

