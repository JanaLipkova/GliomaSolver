/*
 *  IO_Binary.h
 *  MRAG
 *
 *  Created by Manfred Quack on 09/09/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *  Changed 13.September: Now using inheritance to prevent code-replication.
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
#include "MRAG_IO_Native.h"


namespace MRAG
{
	
	template<typename TWavelets, typename TBlock, typename TProjector = dummy_projector, int nChannels=1>
	class IO_Binary : IO_Native<TWavelets, TBlock, TProjector, nChannels>
	{
	public:
		// Constructor/destructor
		IO_Binary();
		~IO_Binary();
		
		// Typedefs
		typedef MRAG::IO_BaseClass<TWavelets, TBlock, TProjector, nChannels> SuperClass;
		typedef typename SuperClass::ElementType ElementType;
		typedef MRAG::Grid<TWavelets, TBlock> GridType;
		typedef map<GridNode *, vector<GridNode *> > HierarchyType;
			
		// Virtual methods
		virtual void Write( GridType & inputGrid, string fileName );
		virtual void Read( GridType & inputGrid, string fileName );
	
	protected:
		
		//virtual:
		virtual void _WriteBlock(TBlock& block, std::ofstream & outputstream, FILE * outputBinary );
		virtual void _ReadBlock(TBlock& block, std::ifstream & inputstream, FILE * inputBinary );

	
	};
	
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	IO_Binary< TWavelets, TBlock, TProjector, nChannels >::
    IO_Binary()
	{
	}
	
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	IO_Binary< TWavelets, TBlock, TProjector, nChannels >::
	~IO_Binary()
	{
	}

    
    //Virtual Write function.
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	void
	IO_Binary< TWavelets, TBlock, TProjector, nChannels >::
	Write( GridType & inputGrid, string fileName )
	{
		int status=-1;
		// Open output file
		string fileNameBinary = fileName+".mrg";
		const char fmode='w';
		FILE* binaryout=fopen(fileNameBinary.c_str(),&fmode);
		assert (binaryout!=NULL);
		
		//Calling Base-Class _Write function, but now with a pointer to the binary File.
		_Write( inputGrid, fileName,binaryout);
		status=fclose(binaryout);
		if(status==0)
		{
			cout << "Sucessfully wrote to binary file " << fileNameBinary << endl;
		}
		else {
			cout << "Something went wrong writing to binary file " << fileNameBinary << endl;
			abort();
		}
	}


	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	void
	IO_Binary< TWavelets, TBlock, TProjector, nChannels >::
	Read( GridType & inputGrid, string fileName )
	{
		int status=-1;
		string fileNameBinary = fileName+".mrg";
		const char fmode='r';
		FILE* binaryin=fopen(fileNameBinary.c_str(),&fmode);
		assert(binaryin!=NULL);
		//Calling Base-Class _Read function, but now with a pointer to the binary File.
		_Read( inputGrid, fileName,binaryin);
		status=fclose(binaryin);
		if(status==0)
		{
			cout << "Sucessfully read from binary file " << fileNameBinary << endl;
		}
		else {
			cout << "Something went wrong reading from binary file " << fileNameBinary << endl;
			abort();
		}

	}

	//Virtual _WriteBlock Method:
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	void IO_Binary< TWavelets, TBlock, TProjector, nChannels >::
	_WriteBlock(TBlock& block, std::ofstream & outputstream, FILE * outputBinary )
	{

		assert(outputBinary!=NULL);
		block.serialize(outputBinary);
    }

	//Virtual _WriteBlock Method:
	template<typename TWavelets, typename TBlock, typename TProjector, int nChannels>
	void IO_Binary< TWavelets, TBlock, TProjector, nChannels >::
	_ReadBlock(TBlock& block, std::ifstream & inputstream,FILE * inputBinary)
	{
		
		assert(inputBinary!=NULL);
		block.deserialize(inputBinary);
	}


	
}//namespace
