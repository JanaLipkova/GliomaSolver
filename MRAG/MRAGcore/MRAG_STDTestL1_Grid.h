/*
 *  MRAG_STDTestGrid.h
 *  MRAG
 *
 *  Created by basil bayati on 7/21/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "MRAG_STDTestL1.h"

namespace MRAG
{
	
template <typename Wavelets, typename Block>
class MRAG_STDTestGrid : public MRAG_STDTestL1< Wavelets, Block >  //LEVEL 1: crashes, memory leaks
{
protected:
	bool run(int nX, int nY, int nZ)
	{		
		try
		{
			cout << "Trying nX: " << nX << " nY: " << nY << " nZ: " << nZ << endl;
			
			MRAG::Grid<Wavelets, Block> * grid = new MRAG::Grid<Wavelets, Block>(nX,nY,nZ);
			delete grid;
			
			cout << endl;
		}
		catch(...)
		{
			cout << "exception here at: nX: " << nX << " nY: " << nY << " nZ: " << nZ << endl;
			abort();
		}
		
		return true;
	}
public:
	static void runTests()
	{
		MRAG_STDTestGrid<Wavelets, Block> test;
		
		int maxGridSize = 6;
		
		for (int i = 0; i < maxGridSize; ++i)
			for (int j = 0; j < maxGridSize; ++j)
				for (int k = 0; k < maxGridSize; ++k)
				{
					test.run( pow(2, i), pow(2, j), pow(2, k) );
				}
		
	}
};


}