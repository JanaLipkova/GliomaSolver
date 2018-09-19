/*
 *  MRAG_STDTestL2.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/25/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

/*
 *  MRAG_STDTestGrid.h
 *  MRAG
 *
 *  Created by basil bayati on 7/21/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
# pragma once

#define _BINARYTEST

#include "../MRAGio/MRAG_IO_BaseClass.h"
#ifndef _BINARYTEST
 #include "../MRAGio/MRAG_IO_Native.h"
#else
 #include "../MRAGio/MRAG_IO_Binary.h"
#endif
#include "../MRAGio/MRAG_IO_VTK.h"
#include "../MRAGcore/MRAGRefiner.h"

namespace MRAG
{
	
template <typename TWavelets, typename TBlock>
class MRAG_STDTestL3_IO   //LEVEL 3
{
protected:
	MRAG_STDTestL3_IO(){};
	~MRAG_STDTestL3_IO(){};
	
	static const int sX=4;
	static const int sY=4;
	static const int sZ=4;
	
	void _ic_func(Real x[3], typename TBlock::ElementType & outE)
	{
		outE=4.0*(0.5-x[0])*(0.5-x[0]);
		//outE=x[0]+10.0*x[1]+100.0*x[2];
	}
	
	
	bool run(void)
	{		
		try
		{
			// Test 1:	I am printing out a grid (culo) with its structure, then I read it, then I print it again (culoCopy) to see if
			//			culo and culoCopy are the same
			
			// Test 2: I refine one block of the grid and I print it out as a VTK format.
			
			// Test 1
			// -----------------------------------------------
			
			// Instatiate
#ifndef _BINARYTEST
			MRAG::IO_Native< TWavelets, TBlock > IO;
#else
			MRAG::IO_Binary< TWavelets, TBlock > IO;
#endif

			// Create a grid which is not meant to be used, to mess up blockIDs
			MRAG::Grid<TWavelets, TBlock> uselessGrid(sX,sY,sZ);
		
			// Create grid to be written
			MRAG::Grid<TWavelets, TBlock> grid(sX,sY,sZ);
			
			// Initialize grid
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				TBlock& block = grid.getBlockCollection()[info.blockID];
				for(int iz=0; iz<TBlock::sizeZ; iz++)
					for(int iy=0; iy<TBlock::sizeY; iy++)
						for(int ix=0; ix<TBlock::sizeX; ix++)
						{
							Real x[3];
							info.pos(x,ix,iy,iz); 
							_ic_func(x,block(ix,iy,iz));
						}
			}
			
			// Write it
			//what is this? bah unsigned  int channelToBeUsed = 1;
			IO.Write( grid, "culo" );
			
			// Initialize new grid to store read information and read printed file
			MRAG::Grid<TWavelets, TBlock> loadedGrid(sX,sY,sZ);			
			IO.Read( loadedGrid, "culo" );
			
			// Re-print it to check that it is equal to the origial one
			IO.Write( loadedGrid, "culoCopy" );
			
			//Verification:
			typename TBlock::ElementType sumerr(0);
			//compare loaded data with Initial Data
			{
			vector<BlockInfo> vInfo = loadedGrid.getBlocksInfo();
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				TBlock& block = loadedGrid.getBlockCollection()[info.blockID];
				for(int iz=0; iz<TBlock::sizeZ; iz++)
					for(int iy=0; iy<TBlock::sizeY; iy++)
						for(int ix=0; ix<TBlock::sizeX; ix++)
						{
							Real x[3];
							info.pos(x,ix,iy,iz);
							typename TBlock::ElementType InitialData;
							_ic_func(x,InitialData);
							if( block(ix,iy,iz)!=InitialData )
							{
								printf("Something went wrong: at x=[ %g,%g,%g], initial: %g, readin: %g \n",x[0],x[1],x[2],InitialData,block(ix,iy,iz));
							}
							sumerr+=fabs(InitialData-block(ix,iy,iz));


						}
				
			}
			}
			
			
			
			if(sumerr>std::numeric_limits<float>::epsilon())
			{
				std::cout  << "Test 1 FAILED! Data is not the same" <<std::endl;
				std::cout << "Error is: " << sumerr <<std::endl;

			}
			else {
				std::cout  << "Test 1 succesfull! Data is the same" <<std::endl;
				std::cout << "Error is: " <<sumerr <<std::endl;

			}
			
		
			
			
			
					
			
			// Test 2
			// -----------------------------------------------
			
			// Refine grid
			Refiner refiner;
			grid.setRefiner(&refiner);
			set<int> blocksToRefine;
			blocksToRefine.insert(grid.getBlocksInfo()[0].blockID);
			grid.refine(blocksToRefine);
			
			// Initialize
			vInfo = grid.getBlocksInfo();
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				TBlock& block = grid.getBlockCollection()[info.blockID];
				for(int iz=0; iz<TBlock::sizeZ; iz++)
					for(int iy=0; iy<TBlock::sizeY; iy++)
						for(int ix=0; ix<TBlock::sizeX; ix++)
						{
							double x[3];
							info.pos(x,ix,iy,iz); 
							block(ix,iy,iz) = 4.0*(0.5-x[0])*(0.5-x[0]);
						}
			}
						
			// Instatiate VTK writer and print the new file
			MRAG::IO_VTK< TWavelets, TBlock > vtk;
			vtk.Write( grid, "culo" );
		}
		catch(...)
		{
			cout << "exception here" << endl;
			abort();
			return false;
		}
		
		return true;
	}
public:
	static void runTests()
	{
		MRAG_STDTestL3_IO<TWavelets, TBlock> test;
		test.run();
	}
};


}