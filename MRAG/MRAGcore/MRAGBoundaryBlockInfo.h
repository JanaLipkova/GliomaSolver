 /*
 *  BoundaryBlockInfo.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 4/24/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once
#undef min
#undef max
#include <vector>
#undef min
#undef max
#include <map>
#undef min
#undef max
#include <set>
#undef min
#undef max

using namespace std;

#include "MRAGCommon.h"
#include "MRAGHuffmanEncoder.h"

namespace MRAG
{
	struct PointIndex
	{
		int blockID;
		/*unsigned short*/  int index;
	};
	
	struct IndexWP
	{
		/*unsigned short*/ int point_index;
		unsigned char weights_index[3];
	};
	
	struct BlockOfGhosts
	{
		unsigned int start;
		unsigned int nGhosts;
	};
	
	struct FaceInfoHeader
	{
		int start[3], end[3];
		unsigned int start_data;
	};
	
	struct BoundaryBlockHeader
	{
		int memsize;
		FaceInfoHeader faces[6];
		int points;
		int weights;
		unsigned int start_weights;
	};
	
	struct BoundaryInfoHeader
	{
		int nBlocks;
		int memsize;
	};
	

	struct BoundaryInfoBlock
	{	
		enum BBIState {BBIState_Initialized, BBIState_Locked, BBIState_Unlocked};
	private:
		BBIState state;
		int nLocks;
		int block_size[3];
		vector<PointIndex> indexPool;
		
		typedef vector<IndexWP> ReconstructionInfo;
		vector< ReconstructionInfo > ghosts;
	
		HuffmanEncoder<unsigned short> encodedInstructionSizes;
		Encoder<unsigned char> encodedInstructionItemsWs;
		Encoder<unsigned short> encodedInstructionItemsPts;
		Encoder<unsigned char> vBlockID_encodedPointIndices3D;
		vector< pair<int, int> > vBlockID_Points;
		bool bCompressed;
		
		void _compress();
		void _decompress();
		void _discard_decompression();
	public:
		BlockOfGhosts boundary[27];
		vector<int> dependentBlockIDs;
		vector<double> weightsPool;
		
		void lock()
		{
			if (state == BBIState_Unlocked)
				
			if(bCompressed)
				_decompress();
			
			nLocks++;
			
			state = BBIState_Locked;
		}
		
		inline const vector<PointIndex>& getIndexPool() const
		{
			assert(state == BBIState_Initialized || state == BBIState_Locked);
			return indexPool;
		}
		
		const vector< ReconstructionInfo >& getGhosts() const
		{
			assert(state == BBIState_Initialized || state == BBIState_Locked);
			return ghosts;
		}
		
		void release()
		{
			assert(nLocks>0);
			assert(state == BBIState_Initialized || state == BBIState_Locked);
			
			nLocks--;
			
			if(nLocks==0)
			{
				const double oldMB = getMemorySize();
				if (state == BBIState_Initialized  && oldMB>128. && !bCompressed)
				{
					_compress();
					_discard_decompression();
					
					bCompressed = true;
					
					const double newMB = getMemorySize();
					printf("Compression factor %f, (%f MB -> %f MB)\n", oldMB/newMB, oldMB, newMB);
				}
				else if (bCompressed)
					_discard_decompression();
				
				
			}
			
			state = BBIState_Unlocked;
		}
		
		void * createBBPack();
		
		float getMemorySize() const;
		
		BoundaryInfoBlock(const int block_size_[3]): 
			indexPool(), weightsPool(), ghosts(), 
			dependentBlockIDs(), state(BBIState_Initialized), nLocks(1),
			encodedInstructionSizes(), encodedInstructionItemsWs(), encodedInstructionItemsPts(),
			vBlockID_encodedPointIndices3D(), vBlockID_Points(), bCompressed(false)
		{
			block_size[0] = block_size_[0];
			block_size[1] = block_size_[1];
			block_size[2] = block_size_[2];
		}
		
		~BoundaryInfoBlock()
		{
			indexPool.clear();
			weightsPool.clear();
			dependentBlockIDs.clear();
			for(int i=0; i<ghosts.size(); i++) 
				ghosts[i].clear();
			
			ghosts.clear();
		}
		
		template<typename WaveletsType, typename BlockType> friend	class MRAG_BBInfoCreator;
	};
		
	struct BoundaryInfo
	{
		char stencil_start[3], stencil_end[3];
		map<int, BoundaryInfoBlock*> boundaryInfoOfBlock;
		
		//set<int> invalidBlocks;
		
		float getMemorySize() const;
		void clear();
		void erase(int blockID, bool bCritical=true);
		BoundaryInfo():boundaryInfoOfBlock(){}
		~BoundaryInfo() { clear(); }
	};
}
