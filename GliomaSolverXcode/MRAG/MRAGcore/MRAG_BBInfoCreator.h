/*
 *  MRAG_BBInfoCreator.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/25/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "MRAGGridNode.h"
#include "MRAGCommon.h"

#include <vector>
#include <map>
#include <math.h>
#include <assert.h>

#ifdef _MRAG_GHOSTSCREATION_ALLOCATOR_HEADER
#include _MRAG_GHOSTSCREATION_ALLOCATOR_HEADER
#endif

#ifndef _MRAG_GHOSTSCREATION_ALLOCATOR
#define _MRAG_GHOSTSCREATION_ALLOCATOR std::allocator	
#endif

#include "MRAG_SmartBlockFinder.h"
#include "MRAGMatrix3D.h"


using namespace std;
namespace MRAG
{

template<typename WaveletsType, typename BlockType>
class MRAG_BBInfoCreator
{
	typedef BlockType B;
	typedef WaveletsType W;
	
	template<typename T> inline _MRAG_GHOSTSCREATION_ALLOCATOR<T> allocator() const { return _MRAG_GHOSTSCREATION_ALLOCATOR<T>();}
	
	static const bool bVerbose = false;
	
	struct BastardGhost
	{
		template<typename T> inline _MRAG_GHOSTSCREATION_ALLOCATOR<T> allocator() const { return _MRAG_GHOSTSCREATION_ALLOCATOR<T>();}
		
		typedef Matrix3D<BastardGhost *, true, _MRAG_GHOSTSCREATION_ALLOCATOR> MatrixOfRequests;
		
		enum BastardGhost_Status {
			BastardGhost_Unresolved=0,
			BastardGhost_Unresolved_PendingRequests=1, 
			BastardGhost_ResolvedKnown=2, 
			BastardGhost_ResolvedUnknown=3
		};
		
		BastardGhost_Status status, candidate_status;
		
		vector<BastardGhost* , _MRAG_GHOSTSCREATION_ALLOCATOR<BastardGhost*> > customers;
		MatrixOfRequests * matRequests;
		I3 index;
		short int block_level;
		double * vWeights[3];
		int nWeightSize[3], nMissingReferences, iPoolIndex;
		
		void resolvedLeaf(bool bKnown, int iPoolIndex)
		{
			assert(status == BastardGhost_Unresolved);

			this->iPoolIndex = iPoolIndex;
			
			status = bKnown? BastardGhost_ResolvedKnown : BastardGhost_ResolvedUnknown;
			
			for(int c=0; c<customers.size(); c++)
				customers[c]->resolved(this, status==BastardGhost_ResolvedKnown? true : false);
		}
		
		void resolved(BastardGhost* retainer, bool bKnown)
		{
			bool bFound = false;
			for(int i =0; i<matRequests->getNumberOfElements();i++)
				bFound |= matRequests->LinAccess(i) == retainer;
			
			assert(bFound);
			assert(status == BastardGhost_Unresolved_PendingRequests);
			assert(nMissingReferences > 0);
				
			nMissingReferences--;
			
			candidate_status = bKnown? BastardGhost_ResolvedKnown : BastardGhost_ResolvedUnknown;
			
			if (nMissingReferences == 0)
			{
				status = candidate_status;
				
				for(int c=0; c<customers.size(); c++)
					customers[c]->resolved(this, status==BastardGhost_ResolvedKnown? true : false);
			}
		}
		
		template <typename FrontType ,typename BufferLayer, typename LayerAllocator, template <typename V, typename A> class BufferType>
		void generateRequests(SmartBlockFinder& smartFinder, BufferType<BufferLayer,LayerAllocator>& buffer, FrontType& newFront, int minLevel, int lying_level)
		{
			assert(status == BastardGhost_Unresolved);
			
			const bool bAnalysis = lying_level>block_level;
			assert(lying_level != block_level);
			
			int s[3] = {
				B::shouldProcessDirectionX? (bAnalysis? index.i[0]*2 - W::HaSupport[1] + 1 : (int)ceil((index.i[0] + W::HsSupport[0])*0.5) ):0,
				B::shouldProcessDirectionY? (bAnalysis? index.i[1]*2 - W::HaSupport[1] + 1 : (int)ceil((index.i[1] + W::HsSupport[0])*0.5) ):0,
				B::shouldProcessDirectionZ? (bAnalysis? index.i[2]*2 - W::HaSupport[1] + 1 : (int)ceil((index.i[2] + W::HsSupport[0])*0.5) ):0
			};
			
			int e[3] = {
				B::shouldProcessDirectionX? (bAnalysis? index.i[0]*2 - W::HaSupport[0] + 1 : 1+(int)floor((index.i[0] + W::HsSupport[1]-1)*0.5) ):1,
				B::shouldProcessDirectionY? (bAnalysis? index.i[1]*2 - W::HaSupport[0] + 1 : 1+(int)floor((index.i[1] + W::HsSupport[1]-1)*0.5) ):1,
				B::shouldProcessDirectionZ? (bAnalysis? index.i[2]*2 - W::HaSupport[0] + 1 : 1+(int)floor((index.i[2] + W::HsSupport[1]-1)*0.5) ):1
			};
			
			const int request_block_level = block_level + (bAnalysis?+1:-1);
			
			matRequests = allocator<MatrixOfRequests>().allocate(1);
			allocator<MatrixOfRequests>().construct(matRequests,  MatrixOfRequests());
			matRequests->_Setup(e[0]-s[0], e[1]-s[1], e[2]-s[2]);
			
			*matRequests = NULL;
			
			for(int c=0; c<3; c++)
			{
				nWeightSize[c] = e[c]-s[c];
			
				vWeights[c] = allocator<double>().allocate(nWeightSize[c]);
				
				if (nWeightSize[c] == 1)
				{
					vWeights[c][0] = 1;
					continue;
				}
				
				if (bAnalysis)
					for(int i=s[c]; i<e[c]; i++)
						vWeights[c][i-s[c]] = W::getHa(index.i[c]*2-i);
				else
					for(int i=s[c]; i<e[c]; i++)
						vWeights[c][i-s[c]] = W::getHs(2*i-index.i[c]);
			}
			
			nMissingReferences = 0;
			
			int i[3];
			for(i[2]=s[2]; i[2]<e[2]; i[2]++)
			for(i[1]=s[1]; i[1]<e[1]; i[1]++)
			for(i[0]=s[0]; i[0]<e[0]; i[0]++)
			{
				const I3 key(i[0], i[1], i[2]);
				/*i[0]*((bAnalysis || !B::shouldProcessDirectionX)?1:2),
							 i[1]*((bAnalysis || !B::shouldProcessDirectionY)?1:2),
							 i[2]*((bAnalysis || !B::shouldProcessDirectionZ)?1:2)* /);//i[1], i[2]); */
				
				const double w = vWeights[0][i[0]-s[0]]*vWeights[1][i[1]-s[1]]*vWeights[2][i[2]-s[2]];
				
				if (w == 0.0) continue;
				
				typename BufferLayer::const_iterator it = buffer[request_block_level- minLevel].find(key);
				
				BastardGhost * request = NULL;
				
				if (it != buffer[request_block_level- minLevel].end()) 
					request = it->second;
				else
				{
					request = allocator<BastardGhost>().allocate(1);
					allocator<BastardGhost>().construct(request, BastardGhost(request_block_level, i));
					
					newFront.push_back(request);
					buffer[request_block_level - minLevel][key] = request;
				}
				
				matRequests->Access(i[0]-s[0], i[1]-s[1], i[2]-s[2]) = request;
				
				if (request->status == BastardGhost::BastardGhost_Unresolved || request->status ==BastardGhost_Unresolved_PendingRequests)
				{
					request->customers.push_back(this);
					nMissingReferences++;
				}
				else abort();
			}
			
			status = BastardGhost::BastardGhost_Unresolved_PendingRequests;
		}
		
		template<int iPass, typename WeightsSet, typename MappingW2I>
		void collect(WeightsSet& wSet, MappingW2I& mapW2I, vector<IndexWP>& info, vector<double>& weightsPool, double wX=1, double wY=1, double wZ=1, bool bVerbose=false)
		{ 
			if (matRequests == NULL) // leaf
			{
				if (iPass==1)
				{
					const int nInfo = info.size();
					info.resize(nInfo+1);
					IndexWP& iwp = info[nInfo];
					iwp.point_index = iPoolIndex;
					
					if ( wX==1 && wY==1 && wZ==1)
					{
						iwp.weights_index[0] = 0;
						iwp.weights_index[1] = 0;
						iwp.weights_index[2] = 0;
					}
					else
					{
						iwp.weights_index[0] = mapW2I[wX];
						iwp.weights_index[1] = mapW2I[wY];
						iwp.weights_index[2] = mapW2I[wZ];
					}
					
				}
				else if (iPass==0)
				{
					wSet.insert(wX);
					wSet.insert(wY);
					wSet.insert(wZ);
				}
			}
			else
			{
				const int n[3] = {
					matRequests->getSize()[0],
					matRequests->getSize()[1],
					matRequests->getSize()[2]
				};
				
				int i[3];
				for(i[2]=0; i[2]<n[2]; i[2]++)
				for(i[1]=0; i[1]<n[1]; i[1]++)
				for(i[0]=0; i[0]<n[0]; i[0]++)
				{
					BastardGhost* request =  matRequests->Access(i[0], i[1], i[2]);
					if (request!=NULL) 
						request->template collect<iPass>(wSet, mapW2I, info, weightsPool, wX*vWeights[0][i[0]], wY*vWeights[1][i[1]], wZ*vWeights[2][i[2]], bVerbose);
				}
			}
		}
		

		BastardGhost(const int block_level_, const int point_index[3]):
			status(BastardGhost_Unresolved), candidate_status(BastardGhost_ResolvedKnown),
			block_level(block_level_), nMissingReferences(0), index(), iPoolIndex(0),
			customers(), matRequests(NULL)
		{
			index.i[0] = point_index[0];
			index.i[1] = point_index[1];
			index.i[2] = point_index[2];
			
			vWeights[0] = NULL;
			vWeights[1] = NULL;
			vWeights[2] = NULL;
		}
		
		//used by the allocator policy
		BastardGhost(const BastardGhost& b):
		status(BastardGhost_Unresolved), candidate_status(BastardGhost_ResolvedKnown),
		block_level(b.block_level), nMissingReferences(0), index(), iPoolIndex(0),
		customers(), matRequests(NULL)
		{
			index.i[0] = b.index.i[0];
			index.i[1] = b.index.i[1];
			index.i[2] = b.index.i[2];
			
			vWeights[0] = NULL;
			vWeights[1] = NULL;
			vWeights[2] = NULL;
		}
		
		~BastardGhost()
		{	
			if (matRequests!=NULL)
			{
				allocator<MatrixOfRequests>().destroy(matRequests);
				allocator<MatrixOfRequests>().deallocate(matRequests,1);
			}
			
			if (vWeights[0] != NULL) allocator<double>().deallocate(vWeights[0], nWeightSize[0]);
			if (vWeights[1] != NULL) allocator<double>().deallocate(vWeights[1], nWeightSize[1]);
			if (vWeights[2] != NULL) allocator<double>().deallocate(vWeights[2], nWeightSize[2]);
			
			vWeights[0] = NULL;
			vWeights[1] = NULL;
			vWeights[2] = NULL;
		}
				
	private:
		
		//forbidden
		BastardGhost& operator=(const BastardGhost&){ abort(); return *this;}
	};
	

	
private:
	
	int stencil_start[3];
	int stencil_end[3];
	int block_size[3];

	
	void _computeEasyGhosts_Analysis(BoundaryInfoBlock& bb, const GridNode& b, const GridNode& n, int code,
									 vector<BastardGhost *>& bastards, int & nFoundEasyGhosts) const;
	void _computeEasyGhosts_Synthesis(BoundaryInfoBlock& bb, const GridNode& b, const GridNode& n, int code,
									 vector<BastardGhost *>& bastards, int & nFoundEasyGhosts) const;
	void _computeEasyGhosts_SameLevel(BoundaryInfoBlock& bb, const GridNode& b, const GridNode& n, int code,
									 vector<BastardGhost *>& bastards, int & nFoundEasyGhosts) const;
	
	void _resolveBastardGhosts(SmartBlockFinder& smartGuy, BoundaryInfoBlock& bb, const GridNode& b, const vector<GridNode*>& neighbors,
									  const vector<BastardGhost *>& bastards) const;
	
public:
	MRAG_BBInfoCreator(const int stencil_start[3], const int  stencil_end[3], const int block_size[3]);
	
	BoundaryInfoBlock* createBoundaryInfoBlock(const GridNode& b, const vector<GridNode*>& neighbors) const;
};

}

#include "MRAG_BBInfoCreator.inl"
