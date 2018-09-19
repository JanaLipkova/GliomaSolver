/*
 *  MRAGBlockProcessing_CUDA.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/13/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include <vector>
#include <set>

#include "MRAGBlockProcessing_CUDA.h"


void split(const Subspace& subspace, const int idx[3], vector<Subspace*>& vOutput)  
{
	assert(idx[0]!=0 || idx[1]!=0 || idx[2]!=0);
	
	//printf("split\n");
	
	const bool shouldSplit[3] = {
		idx[0] == 0 && subspace.end[0]-subspace.start[0]>1,
		idx[1] == 0 && subspace.end[1]-subspace.start[1]>1,
		idx[2] == 0 && subspace.end[2]-subspace.start[2]>1
	};
	
	const int split_size[3] = {
		1+(int)(idx[0] == 0), 
		1+(int)(idx[1] == 0), 
		1+(int)(idx[2] == 0)
	};
	
	const int m[3] = {
		subspace.end[0] - (subspace.end[0] - subspace.start[0])/2,
		subspace.end[1] - (subspace.end[1] - subspace.start[1])/2,
		subspace.end[2] - (subspace.end[2] - subspace.start[2])/2
	};
	
	const int lut_start[3][2] = {
		subspace.start[0], shouldSplit[0]? m[0] : subspace.start[0],
		subspace.start[1], shouldSplit[1]? m[1] : subspace.start[1],
		subspace.start[2], shouldSplit[2]? m[2] : subspace.start[2]
	};
	
	const int lut_end[3][2] = {
		shouldSplit[0]? m[0] : subspace.end[0], subspace.end[0],
		shouldSplit[1]? m[1] : subspace.end[1], subspace.end[1],
		shouldSplit[2]? m[2] : subspace.end[2], subspace.end[2]
	};
	
	int i[3];
	for(i[2]=0; i[2]<split_size[2]; i[2]++)
		for(i[1]=0; i[1]<split_size[1]; i[1]++)
			for(i[0]=0; i[0]<split_size[0]; i[0]++)
			{
				Subspace * child = new Subspace;
				
				child->start[0] = lut_start[0][i[0]];
				child->start[1] = lut_start[1][i[1]];
				child->start[2] = lut_start[2][i[2]];
				
				child->end[0] = lut_end[0][i[0]];
				child->end[1] = lut_end[1][i[1]];
				child->end[2] = lut_end[2][i[2]];
				
				assert((child->end[0]-child->start[0])*(child->end[1]-child->start[1])*(child->end[2]-child->start[2]) > 0);
				
				vOutput.push_back(child);
			}
}

void generateBoundarySubspaces(const int maxSHMEM, Subspace*& vSubspaces, int& n, const int block_size[3], const char stencil_start[3], const char stencil_end[3], const int sizeOfElement)
{
	//A. setup
	//B. create suitable subspaces that fits in the shared memory
	//C. convert the subspaces into uniform partitions

	//A.
	const int nByteSizeBound = maxSHMEM;

	//B.
	set<Subspace *> boundaries[27];
	for(int code=0; code<27; code++)
	{
		if (code == 1 + 3 + 9) continue;
		
		const int idx[3] = { code%3 - 1, (code/3)%3 - 1 , (code/9)%3 - 1}; 
		
		Subspace * initial_subspace = new Subspace();
		
		initial_subspace->setup(idx, stencil_start, stencil_end, block_size);
		
		{
			vector<Subspace *> tmp1, tmp2;
			vector<Subspace *>& vSource = tmp1;
			vector<Subspace *>& vDest = tmp2;
			
			vSource.push_back(initial_subspace);
			
			do
			{
				vDest.clear();
				
				for(vector<Subspace *>::iterator it = vSource.begin(); it != vSource.end(); it++)
				{
					Subspace * subspace = *it;
					
					const int size = subspace->computeByteSize(stencil_start, stencil_end, sizeOfElement);

					if ( size > nByteSizeBound)
					{
						split(*subspace, idx, vDest);
						
						delete subspace;
						
						subspace = NULL;
					}
					else
						vDest.push_back(subspace);
					
					*it = NULL;
				}
				
				std::swap(vSource, vDest);
			}
			while (vDest.size() != vSource.size());
			
			for( vector<Subspace *>::iterator it = vSource.begin(); it != vSource.end(); it++)
			{
				const Subspace& s = **it;
				
				if ((s.end[0]-s.start[0])*(s.end[1]-s.start[1])*(s.end[2]-s.start[2]) == 0)
				{
					delete *it;
					
					continue;
				}

				boundaries[code].insert(*it);
			}
		}
	}
	
	//C.
	{
		int nSubspaces = 0;
		for(int i=0; i<27; i++)
			nSubspaces += boundaries[i].size();
		
		n = nSubspaces;
		vSubspaces = new Subspace[nSubspaces];
		
		Subspace* ptrDest = vSubspaces;
		
		for(int i=0; i<27; i++)
			for(set<Subspace *>::iterator it = boundaries[i].begin(); it != boundaries[i].end(); it++, ptrDest++)
			{
				memcpy(ptrDest, *it, sizeof(Subspace));
				
				delete *it;
			}
	}
}
