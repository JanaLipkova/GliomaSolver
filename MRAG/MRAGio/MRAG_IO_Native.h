/*
 *  IO_Native.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 5/28/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 *	Changed by Manfred Quack on 09/09/09: _WriteBlock and _ReadBlock functions (to facilitate _Binary version)
 *
 */

#pragma once
#include <vector>
#include <stack>
using namespace std;

#include <string>
#include <fstream>
#include "../MRAGcore/MRAGrid.h"
#include "../MRAGcore/MRAGBlock.h"
#include "../MRAGcore/MRAGGridNode.h"
#include "../MRAGcore/MRAGBlockCollection.h"
#include "MRAG_IO_BaseClass.h"


namespace MRAG
{
	
	template<typename TWavelets, typename TBlock, typename TProjector = dummy_projector, int nChannels=1>
	class IO_Native : IO_BaseClass<TWavelets, TBlock, TProjector, nChannels>
	{
	public:
		// Constructor/destructor
		IO_Native();
		~IO_Native();
		
		// Typedefs
		typedef MRAG::IO_BaseClass<TWavelets, TBlock, TProjector, nChannels> SuperClass;
		typedef typename SuperClass::ElementType ElementType;
		typedef MRAG::Grid<TWavelets, TBlock> GridType;
		typedef map<GridNode *, vector<GridNode *> > HierarchyType;
			
		// Virtual methods
		virtual void Write( GridType & inputGrid, string fileName );
		virtual void Read( GridType & inputGrid, string fileName );
	
	protected:
	
		void _Write( GridType & inputGrid, string fileName, FILE* binaryOutputFile=NULL );
		void _Read( GridType & inputGrid, string fileName, FILE* binaryInputFile=NULL );

		
		//virtual:
		virtual void _WriteBlock(TBlock& block, std::ofstream & outputstream,FILE * outputBinary);
		virtual void _ReadBlock(TBlock& block, std::ifstream & inputstream,FILE * inputBinary);

	
	};
	
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	IO_Native< TWavelets, TBlock, TProjector, nChannels >::
	IO_Native()
	{
	}
	
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	IO_Native< TWavelets, TBlock, TProjector, nChannels >::
	~IO_Native()
	{
	}

    
    //Virtual Write function.
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	void
	IO_Native< TWavelets, TBlock, TProjector, nChannels >::
	Write( GridType & inputGrid, string fileName )
	{
		_Write( inputGrid, fileName );
	}

	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	void
	IO_Native< TWavelets, TBlock, TProjector, nChannels >::
	Read( GridType & inputGrid, string fileName )
	{
		_Read( inputGrid, fileName );
	}

	//Virtual _WriteBlock Method:
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	void IO_Native< TWavelets, TBlock, TProjector, nChannels >::
	_WriteBlock(TBlock& block, std::ofstream & outputstream,FILE* outputBinary)
	{
		outputstream.precision(16);
		for(int iz=0; iz<TBlock::sizeZ; iz++)
			for(int iy=0; iy<TBlock::sizeY; iy++)
				for(int ix=0; ix<TBlock::sizeX; ix++)
				{
					outputstream << TProjector::template Project<ElementType, nChannels>( block(ix,iy,iz) ) << "  ";
				}
    }

	//Virtual _ReadBlock Method:
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	void IO_Native< TWavelets, TBlock, TProjector, nChannels >::
	_ReadBlock(TBlock& block, std::ifstream & inputstream,FILE* inputBinary)
	{
		for(int iz=0; iz<TBlock::sizeZ; iz++)
			for(int iy=0; iy<TBlock::sizeY; iy++)
				for(int ix=0; ix<TBlock::sizeX; ix++)
				{
					ElementType value;
					inputstream >> value;
					//printf("value=%f\n", value);
					block(ix,iy,iz) = value;
				}
	}


	
    //Protected internal _Write function
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	void
	IO_Native< TWavelets, TBlock, TProjector, nChannels >::
	_Write( GridType & inputGrid, string fileName, FILE* binaryOutputFile )
	{
		// Open output file
		fileName += ".txt";
		ofstream output( fileName.data() );
		
		if(!output)
		{
			printf("FILE NOT FOUND!\n");
			abort();
		}
		
		// Map pointers to blocks
		map<GridNode*,int> mappingNodes;
		unsigned int counter = 0;
		for(HierarchyType::const_iterator it=inputGrid.m_hierarchy.begin(); it!=inputGrid.m_hierarchy.end(); it++)
		{
			mappingNodes[it->first] = counter;
			counter += 1;
		}
		
		// Map blockIDs
		vector<BlockInfo> vInfo = inputGrid.getBlocksInfo();
		map<int,int> mappingIDs;
		for(int i=0; i<vInfo.size(); i++)
		{
			BlockInfo& info = vInfo[i];
			mappingIDs[info.blockID] = i;
		}

		// Loop over hierarchy and print out information to reconstruct grid topology
		output << inputGrid.m_hierarchy.size() << endl << endl;
		output << mappingNodes[NULL] << endl << endl;
		for(HierarchyType::const_iterator it=inputGrid.m_hierarchy.begin(); it!=inputGrid.m_hierarchy.end(); it++)
		{
			if( it->first == NULL )
			{	
				// If the node is the root I print out fake values for level, isEmpty, blockID, etc. They will not be used
				// is just to mantain the format regular
				output << mappingNodes[it->first] << endl;
				output << "0" << endl;
				output << "1" << endl;
				output << "-1" << endl;
				output << "0" << endl;
				for( unsigned int i = 0; i < 3; i++ ){ output << "0" << "  "; }
				output << endl;
				output << it->second.size() << "  ";
				for( unsigned int i = 0; i < it->second.size(); i++ ){ output << mappingNodes[ it->second[i] ] << "  "; }
				output << endl;
				output << endl;
			}
			else
			{
				output << mappingNodes[it->first] << endl;
				output << it->first->level << endl;
				output << it->first->isEmpty << endl;
				output << mappingIDs[it->first->blockID] << endl;
				output << mappingNodes[it->first->parent] << endl;
				for( unsigned int i = 0; i < 3; i++ ){ output << it->first->index[i] << "  "; }
				output << endl;
				output << it->second.size() << "  ";
				for( unsigned int i = 0; i < it->second.size(); i++ ){ output << mappingNodes[ it->second[i] ] << "  "; }
				output << endl;
				output << endl;
			}
		}
		

		// Print out data contained in blocks
		for(int i=0; i<vInfo.size(); i++)
		{
			const BlockInfo& info = vInfo[i];
			TBlock& block = inputGrid.getBlockCollection().lock(info.blockID);
			output << endl;
			output << mappingIDs[info.blockID] << endl;

			_WriteBlock(block,output,binaryOutputFile);
			
			output << endl;		
			inputGrid.getBlockCollection().release(info.blockID);
		}
		
		
		// Close output file
		output.close();
		cout << "Output file " << fileName << " printed!" << endl;
	}
	
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	void
	IO_Native< TWavelets, TBlock, TProjector, nChannels >::
	_Read( GridType & inputGrid, string fileName, FILE* binaryInputFile )
	{
		inputGrid._dispose();
		
		// Open input file
		fileName += ".txt";
		ifstream input( fileName.c_str() );
		
		if(!input)
		{
			printf("FILE NOT FOUND!\n");
			abort();
		}
		
		//The following code is new, it has been rewritten (by diego)
		
		// Map ID to nodes, load hierarchy in the form of parentID -> childrenIDs
		map<int, GridNode*> mapNodeIDToNode;
		map<int, int> mapChildIDToParentID;
		map<int, int> mapCounterToBlockID;
		
		// Read topology structure
		int nLeaves = 0;
		int numberOfNodes = 0;
		int rootPosition = 0;
		int nodeID = 0;
		int level = 0;
		bool isEmpty = 0;
		int blockID = 0;
		int parentID = 0;
		int index[3] = {0,0,0};
		int childrenNum = 0;
		
		input >> numberOfNodes;
		
		if(numberOfNodes==0)
		{
			printf("numberOfNodes == 0, this is a failure.\nDid you compile this program in Debug Mode with GCC <4.0? STD C++ DEBUG NOT SUPPORTED with that compiler\n");
			abort();
		}
		input >> rootPosition;
		for( unsigned int i = 0; i < numberOfNodes; i++ )
		{	
			input >> nodeID;
			input >> level;
			input >> isEmpty;
			input >> blockID;
			input >> parentID;
			input >> index[0]; 
			input >> index[1]; 
			input >> index[2]; 
			input >> childrenNum;
			
			vector<int> childrenID(childrenNum);
			for( unsigned int j = 0; j < childrenNum; j++ )
			{
				unsigned int value = 0;
				input >> value;
				childrenID[j] = value;
			}
			
			nLeaves += (int)(!isEmpty);
			
			mapNodeIDToNode[nodeID] = new GridNode(isEmpty, NULL, blockID, index[0], index[1], index[2], level);
			mapChildIDToParentID[nodeID] = parentID;
		}
		
		//attach parent pointers, create hierarchy
		{
			delete mapNodeIDToNode[rootPosition];
			
			mapNodeIDToNode[rootPosition] = NULL;
			
			for( unsigned int i = 0; i < numberOfNodes; i++ )
			{
				if (i == rootPosition) continue;
				
				GridNode * child = mapNodeIDToNode[i];
				GridNode * parent = mapNodeIDToNode[mapChildIDToParentID[i]];
				
				child->parent = parent;
				
				inputGrid.m_hierarchy[child];
				
				inputGrid.m_hierarchy[parent].push_back(child);
			}
		}
		
		//allocate blocks, attach them to the leaves
		{
			const int nBlocks = nLeaves;
			
			vector<int> vBlockIDs = inputGrid.m_blockCollection.create(nBlocks);
			
			for( unsigned int i = 0; i < nBlocks; i++ )
				mapCounterToBlockID[i] = vBlockIDs[i];
			
			for(map<int, GridNode*>::const_iterator it = mapNodeIDToNode.begin(); it != mapNodeIDToNode.end(); it++)
				if (it->second != NULL && !it->second->isEmpty)
					it->second->blockID = mapCounterToBlockID[it->second->blockID];
		}
		
		// Read block data
		{
			const int nBlocks = nLeaves;
			for(int i=0; i<nBlocks; i++)
			{
				input >> blockID;

				TBlock& block = inputGrid.getBlockCollection().lock(mapCounterToBlockID[blockID]);
				
				_ReadBlock(block,input,binaryInputFile);
				
				inputGrid.getBlockCollection().release(mapCounterToBlockID[blockID]);
			}
		}
		
		inputGrid._refresh(true);

		// Close input file
		input.close();
		cout << "Input file " << fileName << " read!" << endl;
	}
	
}
