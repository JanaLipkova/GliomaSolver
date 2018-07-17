/*
 *  MRAGSpaceTimeSorter.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 5/28/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once
#undef max
#undef min
#include <vector>

#include <stack>
using namespace std;

#include <string>
#include <fstream>
#include "MRAGcore/MRAGrid.h"
#include "MRAGcore/MRAGBlock.h"
#include "MRAGcore/MRAGGridNode.h"
#include "MRAGcore/MRAGBlockCollection.h"
#include "MRAGcore/MRAGBlockLab.h"
#include "MRAGcore/MRAGHelperMacros.h"

namespace MRAG
{
	template<typename TWavelets, typename TBlock, typename TProjector = dummy_projector, int nChannels=1 >
	class IO_BaseClass
	{
	public:
		// Constructor/destructor
		IO_BaseClass();
		~IO_BaseClass();
		
		// Typedefs
		typedef MRAG::Grid<TWavelets, TBlock> GridType;
		typedef map<GridNode *, vector<GridNode *> > HierarchyType;
		typedef typename TBlock::ElementType ElementType;
			
		// Virtual methods
		virtual void Write( GridType & inputGrid, string fileName ) = 0;
		virtual void Read( GridType & inputGrid, string fileName ) = 0;
	};
	
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	IO_BaseClass< TWavelets, TBlock, TProjector, nChannels >::
	IO_BaseClass()
	{
	}
	
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	IO_BaseClass< TWavelets, TBlock, TProjector, nChannels >::
	~IO_BaseClass()
	{
	}
}
