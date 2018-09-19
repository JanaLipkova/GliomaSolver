/*
 *  MRAG_STDTestL1.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/16/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "MRAGCommon.h"
#include "MRAGrid.h"
   
namespace MRAG
{
	
template <typename Wavelets, typename Block>
class MRAG_STDTestL1 //LEVEL 1: crashes, memory leaks 
{
protected:
	virtual bool run()
	{
		// OVERLOAD THIS
		
		try
		{
			const int nBlocks = 20;
			MRAG::Grid<Wavelets, Block> * grid = new MRAG::Grid<Wavelets, Block>(nBlocks,nBlocks,1);
			
			delete grid;
		}
		catch(...)
		{
			return false;
		}
		
		return true;
	}
	
public:
	static void runTests()
	{
		MRAG_STDTestL1<Wavelets, Block> test;
		
		test.run();
	}
};


}