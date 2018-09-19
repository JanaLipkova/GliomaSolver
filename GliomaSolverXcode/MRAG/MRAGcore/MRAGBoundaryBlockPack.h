/*
 *  MRAGBoundaryBlockPack.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/13/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#include <stdio.h>
#pragma once
#include "MRAGCommon.h"

namespace MRAG
{
	struct PointerToPoint
	{
		void * ptrBlock;
		int index;
	};
	
	struct InstructionItem
	{
		short int idxPoint;
		short int idxWeight[3];
	};
	
	struct BoundaryBlockPack
	{
		char stencil_start[3], stencil_end[3];
		short int boundary_start[27];
		int byteSize;
		int nWeights, nPoints;
		
		Real *ptrWeights; 
		PointerToPoint *ptrPoints; 
		short int *ptrGhostOffsets; 
		InstructionItem *ptrInstructions;
		
		template <typename B>
#ifdef _CUDA_SIDE
		__device__
#endif
		void findCodeAndOffset(int input_pos[3], int& offset, int& code)
		{
			int dir[3] = {
				B::sizeX>1? floor(input_pos[0]/(float)B::sizeX) : 0,
				B::sizeY>1? floor(input_pos[1]/(float)B::sizeY) : 0,
				B::sizeZ>1? floor(input_pos[2]/(float)B::sizeZ) : 0
			 };
			
			code = (dir[2]+1)*9 + (dir[1]+1)*3 + (dir[0]+1);
			
			const int size[3] = {
				dir[0]<0? -stencil_start[0] : (dir[0]==0 ? B::sizeX : stencil_end[0]-1),
				dir[1]<0? -stencil_start[1] : (dir[1]==0 ? B::sizeY : stencil_end[1]-1),
				dir[2]<0? -stencil_start[2] : (dir[2]==0 ? B::sizeZ : stencil_end[2]-1)				
			};
			
			const int s[3] = {
				dir[0]<0? stencil_start[0] : (dir[0]==0 ? 0 : B::sizeX),
				dir[1]<0? stencil_start[1] : (dir[1]==0 ? 0 : B::sizeY),
				dir[2]<0? stencil_start[2] : (dir[2]==0 ? 0 : B::sizeZ)
			};
			
			offset =	(input_pos[0] - s[0]) + 
						(input_pos[1] - s[1])*size[0] + 
						(input_pos[2] - s[2])*size[0]*size[1]; 
		}
	};
	
	template <typename BlockType, typename ElementType>
#ifdef _CUDA_SIDE
	__device__
#endif
	inline void constructGhosts(BoundaryBlockPack& pack, int ipos[3], ElementType& poolElement)
	{	
		int code, ig;
		pack.findCodeAndOffset<BlockType>(ipos, ig, code);
		
		short int * ghost_offsets = pack.ptrGhostOffsets + pack.boundary_start[code];
			
		poolElement = ElementType();
		
		const int ghostStart = ghost_offsets[ig];
		const int nInstructions = ghost_offsets[ig+1] - ghostStart;

		InstructionItem * itemArray = &pack.ptrInstructions[ghostStart];
		
		for(int instr=0; instr<nInstructions; instr++)
		{
			InstructionItem item = itemArray[instr];

			Real w = (	pack.ptrWeights[item.idxWeight[0]]*
						pack.ptrWeights[item.idxWeight[1]]*
						pack.ptrWeights[item.idxWeight[2]]);
			
			PointerToPoint point = pack.ptrPoints[item.idxPoint];

			ElementType contrib = (*((BlockType *)point.ptrBlock))[point.index];
			
			poolElement += (contrib*w);
		}
	}
}
