/*
 *  MRAGSpaceTimeSorter.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 5/28/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include <vector>
#include <stack>
using namespace std;

#include <string>
#include <fstream>
#include "MRAGcore/MRAGrid.h"
#include "MRAGcore/MRAGBlock.h"
#include "MRAGcore/MRAGGridNode.h"
#include "MRAGcore/MRAGBlockCollection.h"
#include "MRAGio/MRAG_IO_BaseClass.h"


namespace MRAG
{
	
	template<typename TWavelets, typename TBlock, typename TProjector = dummy_projector, int nChannels=1>
	class IO_VTK : IO_BaseClass<TWavelets, TBlock, TProjector, nChannels>
	{
	public:
		// Constructor/destructor
		IO_VTK();
		~IO_VTK();
		
		// Typedefs
		typedef MRAG::IO_BaseClass<TWavelets, TBlock, TProjector, nChannels> SuperClass;
		typedef typename SuperClass::ElementType ElementType;
		typedef MRAG::Grid<TWavelets, TBlock> GridType;
		typedef map<GridNode *, vector<GridNode *> > HierarchyType;
			
		// Virtual methods
		virtual void Write( GridType & inputGrid, string fileName );
		virtual void Read( GridType & inputGrid, string fileName );
	};
	
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	IO_VTK< TWavelets, TBlock, TProjector, nChannels >::
	IO_VTK()
	{
	}
	
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	IO_VTK< TWavelets, TBlock, TProjector, nChannels >::
	~IO_VTK()
	{
	}
	
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	void
	IO_VTK< TWavelets, TBlock, TProjector, nChannels >::
	Write( GridType & inputGrid, string fileName )
	{
		// Open output file
		fileName += ".vtk";
		ofstream output( fileName.c_str() );
		
		// Calculate total number of points and cells
		vector<BlockInfo> vInfo = inputGrid.getBlocksInfo();
		unsigned int totalNumberOfPoints = 0;
		unsigned int verticesPerCell = 0;
		unsigned int totalNumberOfCells = 0;
		if( TBlock::sizeZ == 1 )
		{
			totalNumberOfPoints = vInfo.size() * (TBlock::sizeX+1) * (TBlock::sizeY+1);
			totalNumberOfCells = vInfo.size() * (TBlock::sizeX) * (TBlock::sizeY);
			verticesPerCell = 4;
		}
		
		if( TBlock::sizeZ > 1 )
		{
			totalNumberOfPoints = vInfo.size() * (TBlock::sizeX+1) * (TBlock::sizeY+1) * (TBlock::sizeZ+1);
			totalNumberOfCells = vInfo.size() * (TBlock::sizeX) * (TBlock::sizeY) * (TBlock::sizeZ);
			verticesPerCell = 8;
		}
		
		// Preapare header
		output << "# vtk DataFile Version 2.0\n";
		output << "LaMinkia\n";
		output << "ASCII\n";
		output << "DATASET UNSTRUCTURED_GRID\n";
				
		// Print out points coordinates
		output << "POINTS " << totalNumberOfPoints << " float\n";
		for(int i=0; i<vInfo.size(); i++)
		{
			BlockInfo& info = vInfo[i];
			//don't need that for locations of points: TBlock& block = inputGrid.getBlockCollection()[info.blockID];
			if( TBlock::sizeZ > 1 )
			{
				for(int iz=0; iz<TBlock::sizeZ+1; iz++)
					for(int iy=0; iy<TBlock::sizeY+1; iy++)
						for(int ix=0; ix<TBlock::sizeX+1; ix++)
						{
							double x[3];
							info.pos(x, ix, iy, iz);
							for(int j = 0; j < 3; j++){ output << x[j] << "  "; }
							output << endl;
						}
			}
			if( TBlock::sizeZ == 1 )
			{
				for(int iy=0; iy<TBlock::sizeY+1; iy++)
					for(int ix=0; ix<TBlock::sizeX+1; ix++)
					{
						double x[2];
						info.pos(x, ix, iy);
						for(int j = 0; j < 2; j++){ output << x[j] << "  "; }
						output << 0.0;
						output << endl;
					}
			}
		}

		// Print out connectivity
		output << "CELLS " << totalNumberOfCells << " " << totalNumberOfCells * (verticesPerCell+1) << std::endl;
		unsigned int counter = 0;
		for(int i=0; i<vInfo.size(); i++)
		{
			if( TBlock::sizeZ > 1 )
			{
				for(int iz=0; iz<TBlock::sizeZ; iz++)
					for(int iy=0; iy<TBlock::sizeY; iy++)
						for(int ix=0; ix<TBlock::sizeX; ix++)
						{
							output << verticesPerCell << "  ";
							for(int j = 0; j < verticesPerCell; j++)
							{
								int shift[3] = { j&1, (j>>1)&1, (j>>2)&1 };
								int ixx = ix + shift[0];
								int iyy = iy + shift[1];
								int izz = iz + shift[2];
								int pointIndex = counter*(TBlock::sizeX+1)*(TBlock::sizeY+1)*(TBlock::sizeZ+1) + izz*(TBlock::sizeX+1)*(TBlock::sizeY+1) + iyy*(TBlock::sizeX+1) + ixx;
								output << pointIndex << "  ";
							}
							output << std::endl;
						}
				counter += 1;
			}
			if( TBlock::sizeZ == 1 )
			{
				for(int iy=0; iy<TBlock::sizeY; iy++)
					for(int ix=0; ix<TBlock::sizeX; ix++)
						{
							output << verticesPerCell << "  ";
							for(int j = 0; j < verticesPerCell; j++)
							{
								int shift[3] = { j&1, (j>>1)&1, (j>>2)&1 };
								int ixx = ix + shift[0];
								int iyy = iy + shift[1];
								int izz = 0;
								int pointIndex = counter*(TBlock::sizeX+1)*(TBlock::sizeY+1)*(TBlock::sizeZ) + izz*(TBlock::sizeX+1)*(TBlock::sizeY+1) + iyy*(TBlock::sizeX+1) + ixx;
								output << pointIndex << "  ";
							}
							output << std::endl;
						}
				counter += 1;
			}
		}
		
		// Print out cell type
		output << "CELL_TYPES " << totalNumberOfCells << std::endl;
		if( TBlock::sizeZ == 1 ){ for(int i=0; i<totalNumberOfCells; i++){ output << 8 << std::endl; } }
		if( TBlock::sizeZ > 1 ){ for(int i=0; i<totalNumberOfCells; i++){ output << 11 << std::endl; } }
		
		// Print out look up table
		output << "POINT_DATA " << totalNumberOfPoints << std::endl;
		output << "SCALARS scalars float 1"  << std::endl;
		output << "LOOKUP_TABLE default"  << std::endl;
		BlockLab<TBlock> lab;
		int steStart[3] ={0,0,0};
		int steEnd[3] = {2,2,2};
		//added this for 2D support:
		if (TBlock::sizeZ==1)
		{
			steEnd[2]=1;
		}
		

		lab.prepare(inputGrid.getBlockCollection(), inputGrid.getBoundaryInfo(),steStart,steEnd);
		for(int i=0; i<vInfo.size(); i++)
		{
			BlockInfo& info = vInfo[i];
			lab.load(info);
			//using the blocklab, don't need this: \TBlock& block = inputGrid.getBlockCollection()[info.blockID];
			if( TBlock::sizeZ > 1 )
			{
				for(int iz=0; iz<TBlock::sizeZ+1; iz++)
					for(int iy=0; iy<TBlock::sizeY+1; iy++)
						for(int ix=0; ix<TBlock::sizeX+1; ix++)
						{
							output << TProjector::template Project<ElementType, nChannels>( lab(ix,iy,iz) ) << "  ";
						}
				output << std::endl;
			}
			if( TBlock::sizeZ == 1 )
			{
				//not used: unsigned int iz = 0;
				for(int iy=0; iy<TBlock::sizeY+1; iy++)
					for(int ix=0; ix<TBlock::sizeX+1; ix++)
					{
						output << TProjector::template Project<ElementType, nChannels>( lab(ix,iy) ) << "  ";
					}
				output << std::endl;
			}
		}
			
		// Close output file
		output.close();
		cout << "Output file " << fileName << " printed!" << endl;
	}
	
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	void
	IO_VTK< TWavelets, TBlock, TProjector, nChannels >::
	Read( GridType & inputGrid, string fileName )
	{
		cout << "I dont read vtk shit!" << endl;
	}
	
}
