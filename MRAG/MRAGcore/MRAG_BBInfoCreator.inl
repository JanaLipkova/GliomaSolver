/*
 *  MRAG_BBInfoCreator.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/25/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#include "MRAGGridNode.h"
#include "MRAGBoundaryBlockInfo.h"
#include "MRAG_BBInfoCreator.h"
#include <set>

//HELPER
template <typename W, typename Real, bool bAnalysisFilter>
inline  const Real  * cookFilter(const int level_difference, const int  * filter_support);

namespace MRAG
{
	inline int _computeNeighborCode(const GridNode& b, const GridNode& n) 
	{
		const double dilate = pow(2.,-(n.level - b.level));
		
		const double p[3] = {
			dilate*(n.index[0]+0.5) - (b.index[0]+0.5),
			dilate*(n.index[1]+0.5) - (b.index[1]+0.5),
			dilate*(n.index[2]+0.5) - (b.index[2]+0.5)
		};
		
		const double pmax = std::max(fabs(p[0]), std::max(fabs(p[1]), fabs(p[2])));
		
		const int d[3] = {
			(int)(p[0]/pmax), 
			(int)(p[1]/pmax), 
			(int)(p[2]/pmax)
		};
		
		return d[0]+1 + (d[1]+1)*3 + (d[2]+1)*9;
	}
	

	template<typename WaveletsType, typename BlockType>
	void MRAG_BBInfoCreator<WaveletsType, BlockType>::_computeEasyGhosts_SameLevel(BoundaryInfoBlock& bb, const GridNode& b, const GridNode& n, int code, 
																				   vector<BastardGhost *>& bastards, int & nFoundEasyGhosts) const 
	{
		//1. find the intersection
		//2. put the real points in the index pool
		//3. copy the info in the right ghosts
		
		//1.
		const int iS_bspace[3] = {
			std::max((n.index[0] - b.index[0])*B::sizeX, stencil_start[0]),
			std::max((n.index[1] - b.index[1])*B::sizeY, stencil_start[1]),
			std::max((n.index[2] - b.index[2])*B::sizeZ, stencil_start[2])
		};
		
		const int iE_bspace[3] = {
			std::min((n.index[0]+1 - b.index[0])*B::sizeX, B::sizeX + stencil_end[0] - 1),
			std::min((n.index[1]+1 - b.index[1])*B::sizeY, B::sizeY + stencil_end[1] - 1),
			std::min((n.index[2]+1 - b.index[2])*B::sizeZ, B::sizeZ + stencil_end[2] - 1)
		};
		
		const int iS_nspace[3] = {
			iS_bspace[0] + (b.index[0] - n.index[0])*B::sizeX,
			iS_bspace[1] + (b.index[1] - n.index[1])*B::sizeY,
			iS_bspace[2] + (b.index[2] - n.index[2])*B::sizeZ
		};
		
		const int iE_nspace[3] = {
			iE_bspace[0] + (b.index[0] - n.index[0])*B::sizeX,
			iE_bspace[1] + (b.index[1] - n.index[1])*B::sizeY,
			iE_bspace[2] + (b.index[2] - n.index[2])*B::sizeZ
		};

		
		nFoundEasyGhosts = 
			std::max(0,(iE_bspace[0]-iS_bspace[0]))*
			std::max(0,(iE_bspace[1]-iS_bspace[1]))*
			std::max(0,(iE_bspace[2]-iS_bspace[2]));
		
		if (nFoundEasyGhosts == 0) return;
		
		vector<PointIndex>& poolIndex = bb.indexPool;
		const int iIndexStart = poolIndex.size();
		const int nNewPoints = nFoundEasyGhosts;
		
		//2.
		{
			poolIndex.resize(iIndexStart + nNewPoints);
			
			int src[3];
			int counter = 0;
			for(src[2]=iS_nspace[2]; src[2]<iE_nspace[2]; src[2]++)
			for(src[1]=iS_nspace[1]; src[1]<iE_nspace[1]; src[1]++)
			for(src[0]=iS_nspace[0]; src[0]<iE_nspace[0]; src[0]++)
			{
				PointIndex& pti = poolIndex[iIndexStart + counter];
				pti.blockID = n.blockID;
				pti.index = src[0] + src[1]*B::sizeX + src[2]*B::sizeY*B::sizeX;
				counter++;
			}

			assert(counter == nNewPoints);
		}
		
		//3.
		{
			const int iGhostStart = bb.boundary[code].start;
			
			const int dir[3] = {code%3 - 1, (code/3) %3 -1, (code/9) %3 -1};
			
			const int ghosts_start[3] = {
				dir[0]<0? stencil_start[0] : (dir[0]==0 ? 0 : B::sizeX),
				dir[1]<0? stencil_start[1] : (dir[1]==0 ? 0 : B::sizeY),
				dir[2]<0? stencil_start[2] : (dir[2]==0 ? 0 : B::sizeZ)
			};
			
			const int ghosts_end[3] = {
				dir[0]<0? 0 : (dir[0]==0 ? B::sizeX : (B::sizeX + stencil_end[0] - 1)),
				dir[1]<0? 0 : (dir[1]==0 ? B::sizeY : (B::sizeY + stencil_end[1] - 1)),
				dir[2]<0? 0 : (dir[2]==0 ? B::sizeZ : (B::sizeZ + stencil_end[2] - 1))
			};
			
			const int ghosts_size[3] = {
				ghosts_end[0] - ghosts_start[0],
				ghosts_end[1] - ghosts_start[1],
				ghosts_end[2] - ghosts_start[2]
			};

			int counter = 0;
			int dst[3];
			for(dst[2]=iS_bspace[2]; dst[2]<iE_bspace[2]; dst[2]++)
			for(dst[1]=iS_bspace[1]; dst[1]<iE_bspace[1]; dst[1]++)
			for(dst[0]=iS_bspace[0]; dst[0]<iE_bspace[0]; dst[0]++)
			{
				vector<IndexWP>& ghost = bb.ghosts[iGhostStart+
												   dst[0]-ghosts_start[0] +
												   (dst[1]-ghosts_start[1])*ghosts_size[0]+
												   (dst[2]-ghosts_start[2])*ghosts_size[0]*ghosts_size[1]];
							
				const int index = iIndexStart + counter++;
				
				ghost.resize(1);
				
				IndexWP& iwp = ghost[0];
				iwp.weights_index[0] = 0;
				iwp.weights_index[1] = 0;
				iwp.weights_index[2] = 0;
				iwp.point_index = index;
			}
		}
	}
	
	template<typename WaveletsType, typename BlockType>
	void MRAG_BBInfoCreator<WaveletsType, BlockType>::_computeEasyGhosts_Synthesis(BoundaryInfoBlock& bb, const GridNode& b, const GridNode& n, int code, 
																				  vector<BastardGhost *>& bastards, int & nFoundEasyGhosts) const 
	{
		//1. find the "reconstructable" area in n-space into the b-space
		//2. find the intersection between 2. and the ghosts
		//3. allocate the bastards (complement of 3.)
		//4. if 3. is void return
		//5. cook the filter + put it in the wpool
		//6. recognize the needed point in n + put them in the ipool
		//7. scatter to know how to resize the ghosts
		//8. scatter again to fill the ghosts
				
		const int iGhostStart = bb.boundary[code].start;
		
		const int dir[3] = {code%3 - 1, (code/3) %3 -1, (code/9) %3 -1};
		
		const int ghosts_start[3] = {
			dir[0]<0? stencil_start[0] : (dir[0]==0 ? 0 : B::sizeX),
			dir[1]<0? stencil_start[1] : (dir[1]==0 ? 0 : B::sizeY),
			dir[2]<0? stencil_start[2] : (dir[2]==0 ? 0 : B::sizeZ)
		};
		
		const int ghosts_end[3] = {
			dir[0]<0? 0 : (dir[0]==0 ? B::sizeX : (B::sizeX + stencil_end[0] - 1)),
			dir[1]<0? 0 : (dir[1]==0 ? B::sizeY : (B::sizeY + stencil_end[1] - 1)),
			dir[2]<0? 0 : (dir[2]==0 ? B::sizeZ : (B::sizeZ + stencil_end[2] - 1))
		};
		
		const int ghosts_size[3] = {
			ghosts_end[0] - ghosts_start[0],
			ghosts_end[1] - ghosts_start[1],
			ghosts_end[2] - ghosts_start[2]
		};
		
		//1.
		const int inv_level_difference = -(n.level - b.level);
		
		assert(inv_level_difference>0);

		const int sReconstructable_bspace[3] = {
			B::shouldProcessDirectionX ? (((n.index[0]*B::sizeX -1+1 - W::HsSupport[0] - 1 << inv_level_difference) + W::HsSupport[0] + 1) - b.index[0]*B::sizeX) : 0,
			B::shouldProcessDirectionY ? (((n.index[1]*B::sizeY -1+1 - W::HsSupport[0] - 1 << inv_level_difference) + W::HsSupport[0] + 1) - b.index[1]*B::sizeY) : 0,
			B::shouldProcessDirectionZ ? (((n.index[2]*B::sizeZ -1+1 - W::HsSupport[0] - 1 << inv_level_difference) + W::HsSupport[0] + 1) - b.index[2]*B::sizeZ) : 0
		};
		
		const int eReconstructable_bspace[3] = {
			B::shouldProcessDirectionX ? (((1+n.index[0])*B::sizeX- W::HsSupport[1] + 1 << inv_level_difference) + W::HsSupport[1] - 1 - b.index[0]*B::sizeX) : 1,
			B::shouldProcessDirectionY ? (((1+n.index[1])*B::sizeY- W::HsSupport[1] + 1 << inv_level_difference) + W::HsSupport[1] - 1 - b.index[1]*B::sizeY) : 1,
			B::shouldProcessDirectionZ ? (((1+n.index[2])*B::sizeZ- W::HsSupport[1] + 1 << inv_level_difference) + W::HsSupport[1] - 1 - b.index[2]*B::sizeZ) : 1,
		};
		
		//2.
		const int sEasy_bspace[3] = {
			std::max(sReconstructable_bspace[0], ghosts_start[0]),
			std::max(sReconstructable_bspace[1], ghosts_start[1]),
			std::max(sReconstructable_bspace[2], ghosts_start[2])
		};
		
		const int eEasy_bspace[3] = {
			std::min(eReconstructable_bspace[0], ghosts_end[0]),
			std::min(eReconstructable_bspace[1], ghosts_end[1]),
			std::min(eReconstructable_bspace[2], ghosts_end[2])
		};
		
		nFoundEasyGhosts = 
			std::max(0,(eEasy_bspace[0]-sEasy_bspace[0]))*
			std::max(0,(eEasy_bspace[1]-sEasy_bspace[1]))*
			std::max(0,(eEasy_bspace[2]-sEasy_bspace[2]));
		
		//3.
		const int ncoverage_start[3] = {
			(0 + n.index[0]*B::sizeX << inv_level_difference) - b.index[0]*B::sizeX,
			(0 + n.index[1]*B::sizeY << inv_level_difference) - b.index[1]*B::sizeY,
			(0 + n.index[2]*B::sizeZ << inv_level_difference) - b.index[2]*B::sizeZ
		};
		
		const int ncoverage_end[3] = {
			((1+n.index[0])*B::sizeX << inv_level_difference) - b.index[0]*B::sizeX,
			((1+n.index[1])*B::sizeY << inv_level_difference) - b.index[1]*B::sizeY,
			((1+n.index[2])*B::sizeZ << inv_level_difference) - b.index[2]*B::sizeZ		
		};
		
		const int nb_intersec_start[3] = {
			std::max(ncoverage_start[0], ghosts_start[0]),
			std::max(ncoverage_start[1], ghosts_start[1]),
			std::max(ncoverage_start[2], ghosts_start[2])		
		};
		
		const int nb_intersec_end[3] = {
			std::min(ncoverage_end[0], ghosts_end[0]),
			std::min(ncoverage_end[1], ghosts_end[1]),
			std::min(ncoverage_end[2], ghosts_end[2])		
		};
		
		int i[3];
		int nBastards = 0;
		for(i[2]= nb_intersec_start[2]; i[2]<nb_intersec_end[2]; i[2]++)
		for(i[1]= nb_intersec_start[1]; i[1]<nb_intersec_end[1]; i[1]++)
		for(i[0]= nb_intersec_start[0]; i[0]<nb_intersec_end[0]; i[0]++)
		{
			const bool bInsideEasyRange = 
				i[0] >= sEasy_bspace[0] && i[0]<eEasy_bspace[0] &&
				i[1] >= sEasy_bspace[1] && i[1]<eEasy_bspace[1] &&
				i[2] >= sEasy_bspace[2] && i[2]<eEasy_bspace[2];
			
			if (bInsideEasyRange)
				i[0] = eEasy_bspace[0];
			
			if (i[0]>= nb_intersec_end[0]) continue;
			
			nBastards++;
		
			BastardGhost * bastard = allocator<BastardGhost>().allocate(1);
			allocator<BastardGhost>().construct(bastard, BastardGhost( b.level, i));
			
			bastards.push_back(bastard);
		}
		
		//4.
		if (nFoundEasyGhosts == 0) return;
		
		//5.
		const int filter_support[2] =  {
			(- W::HsSupport[1] + 1 << inv_level_difference) + W::HsSupport[1] - 1,
			(1 - W::HsSupport[0] - 1 << inv_level_difference) + W::HsSupport[0] + 1,
		};
		
		const Real * weights = cookFilter< W, Real, false >(inv_level_difference, (const int *)filter_support);
		
		vector<double>& wPool = bb.weightsPool;
		const int iWStart = bb.weightsPool.size();
		const int iWeightsToAdd = filter_support[1] - filter_support[0];
		
		wPool.resize(iWStart + iWeightsToAdd);
		
		for(int i=0; i<iWeightsToAdd; i++)
			wPool[iWStart + i] = weights[i];
		
		//6.
		const int sEasy_nspace[3] = {
			B::shouldProcessDirectionX ? ((b.index[0]*B::sizeX + sEasy_bspace[0] - filter_support[1]>>inv_level_difference)+1 - n.index[0]*B::sizeX) : 0,
			B::shouldProcessDirectionY ? ((b.index[1]*B::sizeY + sEasy_bspace[1] - filter_support[1]>>inv_level_difference)+1 - n.index[1]*B::sizeY) : 0,
			B::shouldProcessDirectionZ ? ((b.index[2]*B::sizeZ + sEasy_bspace[2] - filter_support[1]>>inv_level_difference)+1 - n.index[2]*B::sizeZ) : 0
		};
		
		const int eEasy_nspace[3] = {
			B::shouldProcessDirectionX ? ((b.index[0]*B::sizeX + eEasy_bspace[0] - filter_support[0]-1>>inv_level_difference) - n.index[0]*B::sizeX +1) : 1,
			B::shouldProcessDirectionY ? ((b.index[1]*B::sizeX + eEasy_bspace[1] - filter_support[0]-1>>inv_level_difference) - n.index[1]*B::sizeX +1) : 1,
			B::shouldProcessDirectionZ ? ((b.index[2]*B::sizeX + eEasy_bspace[2] - filter_support[0]-1>>inv_level_difference) - n.index[2]*B::sizeX +1) : 1
		};
		
		const int sIteration_nspace[3] = {
			std::max(0, sEasy_nspace[0]), 
			std::max(0, sEasy_nspace[1]), 
			std::max(0, sEasy_nspace[2]), 
		};
		
		const int eIteration_nspace[3] = {
			std::min((int)B::sizeX, eEasy_nspace[0]), 
			std::min((int)B::sizeY, eEasy_nspace[1]), 
			std::min((int)B::sizeZ, eEasy_nspace[2]), 
		};
		
		const int size_iteration_nspace[3] = {
			eIteration_nspace[0] - sIteration_nspace[0],
			eIteration_nspace[1] - sIteration_nspace[1],
			eIteration_nspace[2] - sIteration_nspace[2]
		};
			
		Matrix3D<char, true, _MRAG_GHOSTSCREATION_ALLOCATOR> matNeeded(size_iteration_nspace[0], size_iteration_nspace[1], size_iteration_nspace[2]);
		Matrix3D<int, true, _MRAG_GHOSTSCREATION_ALLOCATOR> matKey(size_iteration_nspace[0], size_iteration_nspace[1], size_iteration_nspace[2]);
		Matrix3D<int, true, _MRAG_GHOSTSCREATION_ALLOCATOR> matDest(eEasy_bspace[0] - sEasy_bspace[0], eEasy_bspace[1] - sEasy_bspace[1], eEasy_bspace[2] - sEasy_bspace[2]);
		matNeeded = 0;
		
		int currPointsToAdd = 0;
		
		{
			int s[3], d[3];
			for(s[2]=sIteration_nspace[2]; s[2]<eIteration_nspace[2]; s[2]++)
			for(s[1]=sIteration_nspace[1]; s[1]<eIteration_nspace[1]; s[1]++)
			for(s[0]=sIteration_nspace[0]; s[0]<eIteration_nspace[0]; s[0]++)
			{
				const int s_dest[3] = {
					std::max(sEasy_bspace[0], (n.index[0]*B::sizeX + s[0]<<inv_level_difference) + filter_support[0] - b.index[0]*B::sizeX),
					std::max(sEasy_bspace[1], (n.index[1]*B::sizeY + s[1]<<inv_level_difference) + filter_support[0] - b.index[1]*B::sizeY),
					std::max(sEasy_bspace[2], (n.index[2]*B::sizeZ + s[2]<<inv_level_difference) + filter_support[0] - b.index[2]*B::sizeZ)
				};
				
				const int e_dest[3] = {
					std::min(eEasy_bspace[0], (n.index[0]*B::sizeX + s[0]<<inv_level_difference) + filter_support[1] - b.index[0]*B::sizeX),
					std::min(eEasy_bspace[1], (n.index[1]*B::sizeY + s[1]<<inv_level_difference) + filter_support[1] - b.index[1]*B::sizeY),
					std::min(eEasy_bspace[2], (n.index[2]*B::sizeZ + s[2]<<inv_level_difference) + filter_support[1] - b.index[2]*B::sizeZ)
				};
				
				const int s_bspace[3] = {
					(n.index[0]*B::sizeX + s[0]<<inv_level_difference) - b.index[0]*B::sizeX,
					(n.index[1]*B::sizeY + s[1]<<inv_level_difference) - b.index[1]*B::sizeY,
					(n.index[2]*B::sizeZ + s[2]<<inv_level_difference) - b.index[2]*B::sizeZ
				};
			
				for(d[2]=s_dest[2]; d[2]<e_dest[2]; d[2]++)
				for(d[1]=s_dest[1]; d[1]<e_dest[1]; d[1]++)
				for(d[0]=s_dest[0]; d[0]<e_dest[0]; d[0]++)
				{
					const double w = 
						(B::shouldProcessDirectionX ? weights[d[0]-s_bspace[0]-filter_support[0]] : 1)*
						(B::shouldProcessDirectionY ? weights[d[1]-s_bspace[1]-filter_support[0]] : 1)*
						(B::shouldProcessDirectionZ ? weights[d[2]-s_bspace[2]-filter_support[0]] : 1);
					
					if (fabs(w)>0)
					{
						
						matNeeded.Access(s[0]-sIteration_nspace[0], s[1]-sIteration_nspace[1], s[2]-sIteration_nspace[2]) = 1;
						currPointsToAdd++;
						goto _exit_dest;
					}
				}
				
			_exit_dest:		
				int iDummy;
				iDummy = 1;
			}
		}
		
		vector<PointIndex>& indexPool = bb.indexPool;
		
		const int nPointsToAdd = currPointsToAdd;
		const int iStartI = indexPool.size();
		
		indexPool.resize(iStartI + nPointsToAdd);
		
		//A. fill the information of the index pool
		{
			int counter = 0;
			const int nElements = matNeeded.getNumberOfElements();
			for(int i=0; i<nElements; i++)
				if(matNeeded.LinAccess(i)==1) matKey.LinAccess(i) = iStartI + counter++;
			
			counter = 0;
			for(int i=0; i<nElements; i++)
				if(matNeeded.LinAccess(i) == 1)
				{
					PointIndex& pti = indexPool[iStartI + counter++];
					
					const int n_index[3] = {
						(i % size_iteration_nspace[0]) + sIteration_nspace[0], 
						(i/size_iteration_nspace[0] % size_iteration_nspace[1]) + sIteration_nspace[1], 
						(i/size_iteration_nspace[0]/size_iteration_nspace[1] % size_iteration_nspace[2]) + sIteration_nspace[2]};
					
					const int index = n_index[0] + n_index[1]*B::sizeX +  n_index[2]*B::sizeX*B::sizeY;
					
					pti.index = index;
					pti.blockID = n.blockID;
				}
		}
		
		//B. compute the degree of the scattering graph for the destination points
		{
			matDest = 0;
			
			int s[3], d[3];
			for(s[2]=sIteration_nspace[2]; s[2]<eIteration_nspace[2]; s[2]++)
			for(s[1]=sIteration_nspace[1]; s[1]<eIteration_nspace[1]; s[1]++)
			for(s[0]=sIteration_nspace[0]; s[0]<eIteration_nspace[0]; s[0]++)
			{
				const int s_dest[3] = {
					std::max(sEasy_bspace[0], (n.index[0]*B::sizeX + s[0]<<inv_level_difference) + filter_support[0] - b.index[0]*B::sizeX),
					std::max(sEasy_bspace[1], (n.index[1]*B::sizeY + s[1]<<inv_level_difference) + filter_support[0] - b.index[1]*B::sizeY),
					std::max(sEasy_bspace[2], (n.index[2]*B::sizeZ + s[2]<<inv_level_difference) + filter_support[0] - b.index[2]*B::sizeZ)
				};
				
				const int e_dest[3] = {
					std::min(eEasy_bspace[0], (n.index[0]*B::sizeX + s[0]<<inv_level_difference) + filter_support[1] - b.index[0]*B::sizeX),
					std::min(eEasy_bspace[1], (n.index[1]*B::sizeY + s[1]<<inv_level_difference) + filter_support[1] - b.index[1]*B::sizeY),
					std::min(eEasy_bspace[2], (n.index[2]*B::sizeZ + s[2]<<inv_level_difference) + filter_support[1] - b.index[2]*B::sizeZ)
				};
				
				const int s_bspace[3] = {
					(n.index[0]*B::sizeX + s[0]<<inv_level_difference) - b.index[0]*B::sizeX,
					(n.index[1]*B::sizeY + s[1]<<inv_level_difference) - b.index[1]*B::sizeY,
					(n.index[2]*B::sizeZ + s[2]<<inv_level_difference) - b.index[2]*B::sizeZ
				};
				
				for(d[2]=s_dest[2]; d[2]<e_dest[2]; d[2]++)
				for(d[1]=s_dest[1]; d[1]<e_dest[1]; d[1]++)
				for(d[0]=s_dest[0]; d[0]<e_dest[0]; d[0]++)
				{
					const double w = 
					(B::shouldProcessDirectionX ? weights[d[0]-s_bspace[0]-filter_support[0]] : 1)*
					(B::shouldProcessDirectionY ? weights[d[1]-s_bspace[1]-filter_support[0]] : 1)*
					(B::shouldProcessDirectionZ ? weights[d[2]-s_bspace[2]-filter_support[0]] : 1);

					if (fabs(w)>0)
						matDest.Access(d[0] - sEasy_bspace[0], d[1] - sEasy_bspace[1], d[2] - sEasy_bspace[2])++;
				}
			}
			
			for(d[2]=sEasy_bspace[2]; d[2]<eEasy_bspace[2]; d[2]++)
			for(d[1]=sEasy_bspace[1]; d[1]<eEasy_bspace[1]; d[1]++)
			for(d[0]=sEasy_bspace[0]; d[0]<eEasy_bspace[0]; d[0]++)
			{
				int val = matDest.Access(d[0] - sEasy_bspace[0], d[1] - sEasy_bspace[1], d[2] - sEasy_bspace[2]);
				
				if (val>0)
				{
					const int ghost_index = iGhostStart+
							d[0]-ghosts_start[0] +
							(d[1]-ghosts_start[1])*ghosts_size[0]+
							(d[2]-ghosts_start[2])*ghosts_size[0]*ghosts_size[1];
				
					bb.ghosts[ghost_index].resize(val);
				}
			}
			
			//C. now fill the ghosts
			matDest = 0;
			for(s[2]=sIteration_nspace[2]; s[2]<eIteration_nspace[2]; s[2]++)
			for(s[1]=sIteration_nspace[1]; s[1]<eIteration_nspace[1]; s[1]++)
			for(s[0]=sIteration_nspace[0]; s[0]<eIteration_nspace[0]; s[0]++)
			{
				const int s_dest[3] = {
					std::max(sEasy_bspace[0], (n.index[0]*B::sizeX + s[0]<<inv_level_difference) + filter_support[0] - b.index[0]*B::sizeX),
					std::max(sEasy_bspace[1], (n.index[1]*B::sizeY + s[1]<<inv_level_difference) + filter_support[0] - b.index[1]*B::sizeY),
					std::max(sEasy_bspace[2], (n.index[2]*B::sizeZ + s[2]<<inv_level_difference) + filter_support[0] - b.index[2]*B::sizeZ)
				};
				
				const int e_dest[3] = {
					std::min(eEasy_bspace[0], (n.index[0]*B::sizeX + s[0]<<inv_level_difference) + filter_support[1] - b.index[0]*B::sizeX),
					std::min(eEasy_bspace[1], (n.index[1]*B::sizeY + s[1]<<inv_level_difference) + filter_support[1] - b.index[1]*B::sizeY),
					std::min(eEasy_bspace[2], (n.index[2]*B::sizeZ + s[2]<<inv_level_difference) + filter_support[1] - b.index[2]*B::sizeZ)
				};
				
				const int s_bspace[3] = {
					(n.index[0]*B::sizeX + s[0]<<inv_level_difference) - b.index[0]*B::sizeX,
					(n.index[1]*B::sizeY + s[1]<<inv_level_difference) - b.index[1]*B::sizeY,
					(n.index[2]*B::sizeZ + s[2]<<inv_level_difference) - b.index[2]*B::sizeZ
				};
				
				for(d[2]=s_dest[2]; d[2]<e_dest[2]; d[2]++)
				for(d[1]=s_dest[1]; d[1]<e_dest[1]; d[1]++)
				for(d[0]=s_dest[0]; d[0]<e_dest[0]; d[0]++)
				{
					const double w = 
					(B::shouldProcessDirectionX ? weights[d[0]-s_bspace[0]-filter_support[0]] : 1)*
					(B::shouldProcessDirectionY ? weights[d[1]-s_bspace[1]-filter_support[0]] : 1)*
					(B::shouldProcessDirectionZ ? weights[d[2]-s_bspace[2]-filter_support[0]] : 1);
					
					if (fabs(w)>0)
					{
						int& destHit = matDest.Access(d[0] - sEasy_bspace[0], d[1] - sEasy_bspace[1], d[2] - sEasy_bspace[2]);

						const int ghost_index = iGhostStart+
						d[0]-ghosts_start[0] +
						(d[1]-ghosts_start[1])*ghosts_size[0]+
						(d[2]-ghosts_start[2])*ghosts_size[0]*ghosts_size[1];
						
						IndexWP& iwp = bb.ghosts[ghost_index][destHit];
						
						iwp.point_index = matKey.Access(s[0]-sIteration_nspace[0], s[1]-sIteration_nspace[1], s[2]-sIteration_nspace[2]);
						iwp.weights_index[0] = B::shouldProcessDirectionX ? iWStart + d[0]-s_bspace[0]-filter_support[0] : 0;
						iwp.weights_index[1] = B::shouldProcessDirectionY ? iWStart + d[1]-s_bspace[1]-filter_support[0] : 0;
						iwp.weights_index[2] = B::shouldProcessDirectionZ ? iWStart + d[2]-s_bspace[2]-filter_support[0] : 0;
						
						destHit++;
					}
				}
			}
		}

		allocator<Real>().deallocate(const_cast<Real*>(weights), filter_support[1] - filter_support[0]);
	}
	
template<typename WaveletsType, typename BlockType>
void MRAG_BBInfoCreator<WaveletsType, BlockType>::_computeEasyGhosts_Analysis(BoundaryInfoBlock& bb, const GridNode& b, const GridNode& n, int code, 
																			  vector<BastardGhost *>& bastards, int & nFoundEasyGhosts) const 
{
	//1. compute the filter width to used to do a level_difference-analysis filter 
	//2. compute the range of the ghosts that can be easily computed
	//3. put the rest as bastards
	//4. if the range is non zero compute the weight of each point (it does not change!)
	//5. transfer this info in the bb

	//1.
	const int dir[3] = {code%3 - 1, (code/3) %3 -1, (code/9) %3 -1};
	
	const int ghosts_start[3] = {
		dir[0]<0? stencil_start[0] : (dir[0]==0 ? 0 : B::sizeX),
		dir[1]<0? stencil_start[1] : (dir[1]==0 ? 0 : B::sizeY),
		dir[2]<0? stencil_start[2] : (dir[2]==0 ? 0 : B::sizeZ)
	};
	
	const int ghosts_end[3] = {
		dir[0]<0? 0 : (dir[0]==0 ? B::sizeX : (B::sizeX + stencil_end[0] - 1)),
		dir[1]<0? 0 : (dir[1]==0 ? B::sizeY : (B::sizeY + stencil_end[1] - 1)),
		dir[2]<0? 0 : (dir[2]==0 ? B::sizeZ : (B::sizeZ + stencil_end[2] - 1))
	};
	
	const int level_difference = n.level - b.level;
	
	assert(level_difference>0);
	
	const int filter_support[2] =  { //( in the n- index space)
		(- W::HaSupport[1] + 1 << level_difference) + W::HaSupport[1] - 1,
		(1 - W::HaSupport[0] - 1 << level_difference) + W::HaSupport[0] + 1,
	}; 

	const int easy_bspace_start[3] = {
		B::shouldProcessDirectionX ? (int)ceil(pow(2.0, -level_difference)*((0 - filter_support[0]) + n.index[0]*B::sizeX) - b.index[0]*B::sizeX) : 0,
		B::shouldProcessDirectionY ? (int)ceil(pow(2.0, -level_difference)*((0 - filter_support[0]) + n.index[1]*B::sizeY) - b.index[1]*B::sizeY) : 0,
		B::shouldProcessDirectionZ ? (int)ceil(pow(2.0, -level_difference)*((0 - filter_support[0]) + n.index[2]*B::sizeZ) - b.index[2]*B::sizeZ) : 0
	};
	
	const int easy_bspace_end[3] = {
		B::shouldProcessDirectionX ? (1+(int)floor(pow(2.0, -level_difference)*((B::sizeX - filter_support[1]) + n.index[0]*B::sizeX)) - b.index[0]*B::sizeX) : 1,
		B::shouldProcessDirectionY ? (1+(int)floor(pow(2.0, -level_difference)*((B::sizeY - filter_support[1]) + n.index[1]*B::sizeY)) - b.index[1]*B::sizeY) : 1,
		B::shouldProcessDirectionZ ? (1+(int)floor(pow(2.0, -level_difference)*((B::sizeZ - filter_support[1]) + n.index[2]*B::sizeZ)) - b.index[2]*B::sizeZ) : 1	
	};

	//2.
	const int easy_range_start[3] = {
		std::max(easy_bspace_start[0], ghosts_start[0]),
		std::max(easy_bspace_start[1], ghosts_start[1]),
		std::max(easy_bspace_start[2], ghosts_start[2])
	};
	
	const int easy_range_end[3] = {
		std::min(easy_bspace_end[0], ghosts_end[0]),
		std::min(easy_bspace_end[1], ghosts_end[1]),
		std::min(easy_bspace_end[2], ghosts_end[2])
	};
	
	nFoundEasyGhosts = 
		std::max(0,(easy_range_end[0]-easy_range_start[0]))*
		std::max(0,(easy_range_end[1]-easy_range_start[1]))*
		std::max(0,(easy_range_end[2]-easy_range_start[2]));
	
	//3.
	{
		const int entire_bspace_start[3] = {
			B::shouldProcessDirectionX ? ((n.index[0]*B::sizeX-1>>level_difference)+1 - b.index[0]*B::sizeX) : 0,
			B::shouldProcessDirectionY ? ((n.index[1]*B::sizeY-1>>level_difference)+1 - b.index[1]*B::sizeY) : 0,
			B::shouldProcessDirectionZ ? ((n.index[2]*B::sizeZ-1>>level_difference)+1 - b.index[2]*B::sizeZ) : 0,
		};
		
		const int entire_bspace_end[3] = {
			B::shouldProcessDirectionX ? (((n.index[0]+1)*B::sizeX-1>>level_difference)+1 - b.index[0]*B::sizeX) : 1,
			B::shouldProcessDirectionY ? (((n.index[1]+1)*B::sizeY-1>>level_difference)+1 - b.index[1]*B::sizeY) : 1,
			B::shouldProcessDirectionZ ? (((n.index[2]+1)*B::sizeZ-1>>level_difference)+1 - b.index[2]*B::sizeZ) : 1,	
		};
		

		const int s[3] = {
			std::max(entire_bspace_start[0], ghosts_start[0]),
			std::max(entire_bspace_start[1], ghosts_start[1]),
			std::max(entire_bspace_start[2], ghosts_start[2])			
		};
		
		const int e[3] = {
			std::min(entire_bspace_end[0], ghosts_end[0]),
			std::min(entire_bspace_end[1], ghosts_end[1]),
			std::min(entire_bspace_end[2], ghosts_end[2])			
		};
		
		int nBastards = 0;
		int i[3];
		for(i[2]=s[2]; i[2]<e[2]; i[2]++)
		for(i[1]=s[1]; i[1]<e[1]; i[1]++)
		for(i[0]=s[0]; i[0]<e[0]; i[0]++)
		{
			const bool bInsideEasyRange = 
			i[0] >= easy_range_start[0] && i[0]<easy_range_end[0] &&
			i[1] >= easy_range_start[1] && i[1]<easy_range_end[1] &&
			i[2] >= easy_range_start[2] && i[2]<easy_range_end[2];
			
			if (bInsideEasyRange) i[0] = easy_range_end[0];
			
			if (i[0]>= e[0]) continue;
			
			nBastards++;
			
			BastardGhost * bastard = allocator<BastardGhost>().allocate(1);
			allocator<BastardGhost>().construct(bastard, BastardGhost( b.level, i));
		
			bastards.push_back(bastard);
		}
				
		assert(nBastards+  nFoundEasyGhosts == (e[0]-s[0])*(e[1]-s[1])*(e[2]-s[2]));
	}
	

	
	//4.
	if (nFoundEasyGhosts == 0) return;
	
	{
		//A. cook the filter 
		//B. put the filter in the bbinfo
		//C. check which neighbor points are affected 
		
		//A.
		const Real * weights = cookFilter<W,Real, true>(level_difference, (const int *)filter_support);
		
		//B.
		vector<double> & wPool = bb.weightsPool;
		
		const int nFilterSize = filter_support[1]-filter_support[0];
		const int iStartW = wPool.size();
		
		{
			wPool.resize(iStartW + nFilterSize + 1);
			
			for(int i=0; i<nFilterSize; i++)
				wPool[i+iStartW] = weights[i];
		}
			
		//C.
		{
			const int tmp_start[3] = {
				((b.index[0]*B::sizeX + easy_range_start[0]) << level_difference) - n.index[0]*B::sizeX,
				((b.index[1]*B::sizeY + easy_range_start[1]) << level_difference) - n.index[1]*B::sizeY,
				((b.index[2]*B::sizeZ + easy_range_start[2]) << level_difference) - n.index[2]*B::sizeZ			
			};
			
			const int tmp_end[3] = {
				((b.index[0]*B::sizeX + easy_range_end[0]-1) << level_difference) - n.index[0]*B::sizeX+1,
				((b.index[1]*B::sizeY + easy_range_end[1]-1) << level_difference) - n.index[1]*B::sizeY+1,
				((b.index[2]*B::sizeZ + easy_range_end[2]-1) << level_difference) - n.index[2]*B::sizeZ+1			
			};
			
			const int easy_nspace_start[3] = {
				B::shouldProcessDirectionX ? (tmp_start[0] + filter_support[0])  : 0,
				B::shouldProcessDirectionY ? (tmp_start[1] + filter_support[0])  : 0,
				B::shouldProcessDirectionZ ? (tmp_start[2] + filter_support[0])  : 0
			};
			
			const int easy_nspace_end[3] = {
				B::shouldProcessDirectionX ? (-1 +  tmp_end[0] + filter_support[1]) : 1,
				B::shouldProcessDirectionY ? (-1 +  tmp_end[1] + filter_support[1]) : 1,
				B::shouldProcessDirectionZ ? (-1 +  tmp_end[2] + filter_support[1]) : 1	
			};
			
			const int easy_size[3] = {
				easy_nspace_end[0] - easy_nspace_start[0],
				easy_nspace_end[1] - easy_nspace_start[1],
				easy_nspace_end[2] - easy_nspace_start[2],
			};
			
			Matrix3D<int> matPoints(easy_range_end[0] - easy_range_start[0],
									easy_range_end[1] - easy_range_start[1],
									easy_range_end[2] - easy_range_start[2]);
			
			Matrix3D<char> matNeeded(easy_size[0], easy_size[1], easy_size[2]);
			Matrix3D<int> matKey(easy_size[0], easy_size[1], easy_size[2]);
			
			matPoints = 0;
			matNeeded = 0;
			
			int d[3],s[3];
			
			for(d[2]=easy_range_start[2]; d[2]<easy_range_end[2]; d[2]++)
			for(d[1]=easy_range_start[1]; d[1]<easy_range_end[1]; d[1]++)
			for(d[0]=easy_range_start[0]; d[0]<easy_range_end[0]; d[0]++)
			{
				const int d_nspace[3] = {
					((b.index[0]*B::sizeX + d[0]) << level_difference) - n.index[0]*B::sizeX,
					((b.index[1]*B::sizeY + d[1]) << level_difference) - n.index[1]*B::sizeY,
					((b.index[2]*B::sizeZ + d[2]) << level_difference) - n.index[2]*B::sizeZ
				};
				
				const int src_start[3] = {
					B::shouldProcessDirectionX ? (d_nspace[0] + filter_support[0]) : 0,
					B::shouldProcessDirectionY ? (d_nspace[1] + filter_support[0]) : 0,
					B::shouldProcessDirectionZ ? (d_nspace[2] + filter_support[0]) : 0,
				};
				
				const int src_end[3] = {
					B::shouldProcessDirectionX ? (d_nspace[0] + filter_support[1]) : 1,
					B::shouldProcessDirectionY ? (d_nspace[1] + filter_support[1]) : 1,
					B::shouldProcessDirectionZ ? (d_nspace[2] + filter_support[1]) : 1,
				};
				
				int counter = 0;
				for(s[2]=src_start[2]; s[2]<src_end[2]; s[2]++)
				for(s[1]=src_start[1]; s[1]<src_end[1]; s[1]++)
				for(s[0]=src_start[0]; s[0]<src_end[0]; s[0]++)
				{
					const double w = 
						(B::shouldProcessDirectionX ? weights[s[0]-src_start[0]-filter_support[0]] : 1)*
						(B::shouldProcessDirectionY ? weights[s[1]-src_start[1]-filter_support[0]] : 1)*
						(B::shouldProcessDirectionZ ? weights[s[2]-src_start[2]-filter_support[0]] : 1);
						
					if ( w > 0 ) 
					{
						matNeeded.Access(s[0]-easy_nspace_start[0], 
										 s[1]-easy_nspace_start[1], 
										 s[2]-easy_nspace_start[2]) = 1;
						counter++;
					}
				}
				
				matPoints.Access(d[0]-easy_range_start[0],
								 d[1]-easy_range_start[1],
								 d[2]-easy_range_start[2]) = counter;
				
				assert(counter>0);
			}
			
			vector<PointIndex>& poolIndex = bb.indexPool;
			const int iIndexStart = poolIndex.size();
			
			//C-2. compute how many source points are needed
			int sum = 0;
			const int nElements = matNeeded.getNumberOfElements();
			for(int i=0; i<nElements; i++)
				sum += matNeeded.LinAccess(i);
		
			const int nPoints = sum;
			
			poolIndex.resize(iIndexStart + nPoints);
			
			//C-3. fill the index pool
			{
				int counter = 0;
				const int nElements = matNeeded.getNumberOfElements();
				for(int i=0; i<nElements; i++)
					if (matNeeded.LinAccess(i) == 1)
					{
						const int n_index[3] = {
							(i % easy_size[0]) + easy_nspace_start[0], 
							(i/easy_size[0] % easy_size[1]) + easy_nspace_start[1], 
							(i/easy_size[0]/easy_size[1] % easy_size[2]) + easy_nspace_start[2]};
						
						const int index = n_index[0] + n_index[1]*B::sizeX +  n_index[2]*B::sizeX*B::sizeY;
						
						const int key = iIndexStart + counter++;
						
						matKey.LinAccess(i) = key;
		
						PointIndex& pti = poolIndex[key];
						pti.blockID = n.blockID;
						pti.index = index;
					}
			}
		
		
			const int iGhostStart = bb.boundary[code].start;
			
			const int dir[3] = {code%3 - 1, (code/3) %3 -1, (code/9) %3 -1};
			
			const int ghosts_start[3] = {
				dir[0]<0? stencil_start[0] : (dir[0]==0 ? 0 : B::sizeX),
				dir[1]<0? stencil_start[1] : (dir[1]==0 ? 0 : B::sizeY),
				dir[2]<0? stencil_start[2] : (dir[2]==0 ? 0 : B::sizeZ)
			};
			
			const int ghosts_end[3] = {
				dir[0]<0? 0 : (dir[0]==0 ? B::sizeX : (B::sizeX + stencil_end[0] - 1)),
				dir[1]<0? 0 : (dir[1]==0 ? B::sizeY : (B::sizeY + stencil_end[1] - 1)),
				dir[2]<0? 0 : (dir[2]==0 ? B::sizeZ : (B::sizeZ + stencil_end[2] - 1))
			};

			const int ghosts_size[3] = {
				ghosts_end[0] - ghosts_start[0],
				ghosts_end[1] - ghosts_start[1],
				ghosts_end[2] - ghosts_start[2]
			};
			
			//C-4. now fill the ghosts
			for(d[2]=easy_range_start[2]; d[2]<easy_range_end[2]; d[2]++)
			for(d[1]=easy_range_start[1]; d[1]<easy_range_end[1]; d[1]++)
			for(d[0]=easy_range_start[0]; d[0]<easy_range_end[0]; d[0]++)
			{
				vector<IndexWP>& ghost = bb.ghosts[iGhostStart+
							d[0]-ghosts_start[0] +
							(d[1]-ghosts_start[1])*ghosts_size[0]+
							(d[2]-ghosts_start[2])*ghosts_size[0]*ghosts_size[1]];
				
				const int d_nspace[3] = {
					((b.index[0]*B::sizeX + d[0]) << level_difference) - n.index[0]*B::sizeX,
					((b.index[1]*B::sizeY + d[1]) << level_difference) - n.index[1]*B::sizeY,
					((b.index[2]*B::sizeZ + d[2]) << level_difference) - n.index[2]*B::sizeZ
				};
				
				const int src_start[3] = {
					B::shouldProcessDirectionX ? (d_nspace[0] + filter_support[0]) : 0,
					B::shouldProcessDirectionY ? (d_nspace[1] + filter_support[0]) : 0,
					B::shouldProcessDirectionZ ? (d_nspace[2] + filter_support[0]) : 0,
				};
				
				const int src_end[3] = {
					B::shouldProcessDirectionX ? (d_nspace[0] + filter_support[1]) : 1,
					B::shouldProcessDirectionY ? (d_nspace[1] + filter_support[1]) : 1,
					B::shouldProcessDirectionZ ? (d_nspace[2] + filter_support[1]) : 1,
				};
				
				ghost.resize( matPoints.Access(d[0]-easy_range_start[0],
											   d[1]-easy_range_start[1],
											   d[2]-easy_range_start[2]) );
				
				int counter = 0;
				for(s[2]=src_start[2]; s[2]<src_end[2]; s[2]++)
				for(s[1]=src_start[1]; s[1]<src_end[1]; s[1]++)
				for(s[0]=src_start[0]; s[0]<src_end[0]; s[0]++)
				{
					const double w[3] = {
						(B::shouldProcessDirectionX ? weights[s[0]-src_start[0]-filter_support[0]] : 1),
						(B::shouldProcessDirectionY ? weights[s[1]-src_start[1]-filter_support[0]] : 1),
						(B::shouldProcessDirectionZ ? weights[s[2]-src_start[2]-filter_support[0]] : 1)						
					};
					
					if ( w[0]*w[1]*w[2] > 0 ) 
					{
						IndexWP& iwp = ghost[counter];
						iwp.point_index = matKey.Access(s[0]-easy_nspace_start[0], 
														s[1]-easy_nspace_start[1], 
														s[2]-easy_nspace_start[2]);

						iwp.weights_index[0] = B::shouldProcessDirectionX ? (iStartW + s[0]-src_start[0]-filter_support[0]) : 0;
						iwp.weights_index[1] = B::shouldProcessDirectionY ? (iStartW + s[1]-src_start[1]-filter_support[0]) : 0;
						iwp.weights_index[2] = B::shouldProcessDirectionZ ? (iStartW + s[2]-src_start[2]-filter_support[0]) : 0;
						counter++;
					}
				}
			}
		}
		
		allocator<Real>().deallocate(const_cast<Real*>(weights), filter_support[1] - filter_support[0]);
	}
}
	
template<typename WaveletsType, typename BlockType>
MRAG_BBInfoCreator<WaveletsType, BlockType>::MRAG_BBInfoCreator(const int stencil_start[3], const int  stencil_end[3], const int block_size[3])
{
	this->stencil_start[0] = stencil_start[0];
	this->stencil_start[1] = stencil_start[1];
	this->stencil_start[2] = stencil_start[2];
	
	this->stencil_end[0] = stencil_end[0];
	this->stencil_end[1] = stencil_end[1];
	this->stencil_end[2] = stencil_end[2];
	
	this->block_size[0] = block_size[0];
	this->block_size[1] = block_size[1];
	this->block_size[2] = block_size[2];
}

template<typename WaveletsType, typename BlockType>
BoundaryInfoBlock* MRAG_BBInfoCreator<WaveletsType, BlockType>::createBoundaryInfoBlock(const GridNode& b, const vector<GridNode*>& neighbors) const 
{
	//1. allocate BoundaryInfoBlock, fill dependentBlockIDs
	//2. allocate the ghosts needed
	//3. put 1.0 in the wpool
	//4. look for easy ghosts
	//5. solve for the bastards
	
	//1.
	BoundaryInfoBlock * bbinfo = new BoundaryInfoBlock(block_size);
	
	{
		set<int> s;
		for(int i=0; i<neighbors.size(); i++)
			s.insert(neighbors[i]->blockID);
		
		bbinfo->dependentBlockIDs.resize(s.size());
		
		copy(s.begin(), s.end(), bbinfo->dependentBlockIDs.begin());
	}

	//2.
	{
		int nTotalGhosts = 0;
		for(int code=0; code<27; code++)
		{
			if (code == 1 + 3 + 9) continue;
			
			const int d[3] = {code%3 - 1, (code/3) %3 -1, (code/9) %3 -1};
			
			const int s[3] = {
				d[0]<0? stencil_start[0] : (d[0]==0 ? 0 : B::sizeX),
				d[1]<0? stencil_start[1] : (d[1]==0 ? 0 : B::sizeY),
				d[2]<0? stencil_start[2] : (d[2]==0 ? 0 : B::sizeZ)
			};
			
			const int e[3] = {
				d[0]<0? 0 : (d[0]==0 ? B::sizeX : (B::sizeX + stencil_end[0] - 1)),
				d[1]<0? 0 : (d[1]==0 ? B::sizeY : (B::sizeY + stencil_end[1] - 1)),
				d[2]<0? 0 : (d[2]==0 ? B::sizeZ : (B::sizeZ + stencil_end[2] - 1))
			};
			
			const int nBunchOfGhosts =  (e[0]-s[0])*(e[1]-s[1])*(e[2]-s[2]);
			assert(nBunchOfGhosts>= 0);
			bbinfo->boundary[code].start = nTotalGhosts;
			bbinfo->boundary[code].nGhosts = nBunchOfGhosts;
			
			nTotalGhosts += nBunchOfGhosts;
		}
		
		bbinfo->ghosts.resize(nTotalGhosts);
	}

	//3.
	bbinfo->weightsPool.push_back(1.0);
	
	//4.
	vector<BastardGhost *> bastards;
	
	set<int> coveredCodes;
	for(int i=0; i<27; i++)
		coveredCodes.insert(i);
	
	int nEasyOnes = 0;

	for(vector<GridNode*>::const_iterator itNeighbor = neighbors.begin(); itNeighbor != neighbors.end(); itNeighbor++)
	{
		const GridNode& n = **itNeighbor;
		
		int code = _computeNeighborCode(b, n);
		
		if (n.level > b.level)
			_computeEasyGhosts_Analysis(*bbinfo, b, n, code, bastards, nEasyOnes);
		else if (n.level < b.level)
			_computeEasyGhosts_Synthesis(*bbinfo, b, n, code, bastards, nEasyOnes);
		else if (n.level == b.level)
			_computeEasyGhosts_SameLevel(*bbinfo, b, n, code, bastards, nEasyOnes);
		
		coveredCodes.erase(code);
	}
	
	for(set<int>::iterator itCode = coveredCodes.begin(); itCode!=coveredCodes.end(); itCode++ )
	{
		const int code = *itCode;
		
		if (code == 1 + 3 + 9) continue;
		
		const int d[3] = {code%3 - 1, (code/3) %3 -1, (code/9) %3 -1};
		
		const int s[3] = {
			d[0]<0? stencil_start[0] : (d[0]==0 ? 0 : B::sizeX),
			d[1]<0? stencil_start[1] : (d[1]==0 ? 0 : B::sizeY),
			d[2]<0? stencil_start[2] : (d[2]==0 ? 0 : B::sizeZ)
		};
		
		const int e[3] = {
			d[0]<0? 0 : (d[0]==0 ? B::sizeX : (B::sizeX + stencil_end[0] - 1)),
			d[1]<0? 0 : (d[1]==0 ? B::sizeY : (B::sizeY + stencil_end[1] - 1)),
			d[2]<0? 0 : (d[2]==0 ? B::sizeZ : (B::sizeZ + stencil_end[2] - 1))
		};
		
		
		int i[3];
		for(i[2]=s[2]; i[2]<e[2]; i[2]++)
		for(i[1]=s[1]; i[1]<e[1]; i[1]++)
		for(i[0]=s[0]; i[0]<e[0]; i[0]++)
		{
			BastardGhost * bastard = allocator<BastardGhost>().allocate(1);
			allocator<BastardGhost>().construct(bastard, BastardGhost( b.level, i));
			bastards.push_back(bastard);
		}
	}
		
	if (bastards.size()>0)
	{
		SmartBlockFinder smartGuy(b, neighbors, stencil_start, stencil_end, block_size);
		
		_resolveBastardGhosts(smartGuy, *bbinfo, b, neighbors, bastards);
	}
	
	if (bVerbose) printf("Weights: %d, Points:%d\n", (int)bbinfo->weightsPool.size(), (int)bbinfo->indexPool.size());
	
	for(typename vector<BastardGhost *>::iterator it= bastards.begin(); it!= bastards.end(); it++)
	{
		allocator<BastardGhost>().destroy(*it);
		allocator<BastardGhost>().deallocate(*it,1);
		*it = NULL;
	}
	
	{
		PointIndex pointindex;
		IndexWP indexwp;
		if (powf(2,sizeof(pointindex.index)*8-1) <  B::sizeX* B::sizeY* B::sizeZ)
		{
			printf("OUCH: (1<<(sizeof(pointindex.index)*8-1) <  B::sizeX* B::sizeY* B::sizeZ)\n");
			printf("-> OUCH: %f < %d\n", powf(2,sizeof(pointindex.index)*8-1), B::sizeX* B::sizeY* B::sizeZ);
			abort();
		}
		
		if (powf(2, sizeof(indexwp.weights_index[0])*8) <= bbinfo->weightsPool.size())
		{
			printf("OUCH: (1<<(sizeof(indexwp.weights_index[0])*8) <= bbinfo->weightsPool.size()\n");
			printf("-> OUCH: %f < %d\n", powf(2, sizeof(indexwp.weights_index[0])*8) , (int)bbinfo->weightsPool.size());
			abort();
		}
	}
	
	bbinfo->release();
	
	return bbinfo;
}
	
template<typename WaveletsType, typename BlockType>
void MRAG_BBInfoCreator<WaveletsType, BlockType>::_resolveBastardGhosts(
			SmartBlockFinder& smartFinder, BoundaryInfoBlock& bb, const GridNode& b, 
			const vector<GridNode*>& neighbors, const vector<BastardGhost *>& bastards) const
{
	//1. setup
	//2. if the new front is empty exit
	//3.	pop a ghost from the old front:
	//4-a.		If the ghost is resolved using the smartFinder, back-track the solution to the parents in the buffer, go to 3.
	//4-b.		The ghost is not found: find the level where the ghost is lying.
	//5.		if ghost.level < b.level generate synthesis ghosts otherwise analysis (abort if lying level = block level)
	//6.		put the generated ghosts in the new front if they not exist already
	//7.	swap the old front with the new front, goto 2.
	//8.	collect the shit
	//9.	free the shit
	
	//1.
	typedef _MRAG_GHOSTSCREATION_ALLOCATOR<BastardGhost *> ABG;
	typedef vector<BastardGhost *, ABG> GhostVector;
	typedef _MRAG_GHOSTSCREATION_ALLOCATOR< std::pair<const I3, BastardGhost*> > ABGMap;
	typedef map<I3, BastardGhost*, std::less<I3>, ABGMap > GhostMap;
	typedef _MRAG_GHOSTSCREATION_ALLOCATOR<GhostMap> ABGBuffer;
	typedef vector< GhostMap, ABGBuffer > GhostBuffer;
	typedef set<double, std::less<double>, _MRAG_GHOSTSCREATION_ALLOCATOR<double> > WeightSet;
	typedef map<double, int, std::less<double>, _MRAG_GHOSTSCREATION_ALLOCATOR<std::pair< const double, int > > > WeightToIndexMap;

    const int start_level = smartFinder.minLevel;
	
	GhostBuffer buffer(smartFinder.maxLevel - start_level+1);

	GhostVector pezzodimerda(bastards.size());
	copy(bastards.begin(), bastards.end(), pezzodimerda.begin());
	
	GhostVector tmp;
	tmp.reserve(bastards.size());
	
	GhostVector& newFront = tmp ;
	GhostVector& oldFront = pezzodimerda;
	
	GhostVector resolvedBastards;
	tmp.reserve(4*bastards.size());
	
	vector<PointIndex>& indexPool = bb.indexPool;
	
	const int iBastardStart = indexPool.size();
	
	//2.
	while (oldFront.size()>0)
	{
		newFront.clear();
	
		//3.
		while(oldFront.size()>0)
		{
			BastardGhost* bastard = oldFront.back();
			oldFront.pop_back();
			
			//4.
			const GridNode * node = NULL;
			PointIndex * pti = smartFinder.findData(bastard->block_level, bastard->index, &node);
			
			if (pti!=NULL)
			{
				//4-a
				indexPool.push_back(*pti);
				resolvedBastards.push_back(bastard);
			}
			else //4-b, 5., 6.
				bastard->generateRequests(smartFinder, buffer, newFront, start_level, node->level);
		}
		
		//7.
		std::swap(oldFront, newFront);
	}
	
	//7-b. trigger the solution computation from the leaf
	{
		const int nBastards = resolvedBastards.size();
		for(int i=0; i<nBastards; i++)
			resolvedBastards[i]->resolvedLeaf(true, iBastardStart+i);
	}
	
	//7-c. check that every ghost is resolved
	bool bAbort = false;
	for(int l=0; l<buffer.size(); l++)
	{
		GhostMap& m = buffer[l];
		for(typename GhostMap::const_iterator it = m.begin(); it!=m.end(); it++)
		{
			if (it->second->status != BastardGhost::BastardGhost_ResolvedKnown)
			{
				printf("ERROR!! (level %d) key: %d %d %d value: ghost level %d index %d %d %d\n", 
					   l, it->first.i[0], it->first.i[1], it->first.i[2],
					   l, it->second->index.i[0], it->second->index.i[1], it->second->index.i[2]);
				if (it->second->status != BastardGhost::BastardGhost_Unresolved) printf("no pending requests!! status: %d\n", it->second->status);
				bAbort = true;
			}
		}
	}
	if (bAbort)
	{
		printf("Aborting in MRAG_BBInfoCreator<WaveletsType, BlockType>::_resolveBastardGhosts/\n");
		printf("const GridNode& b = (index=%d %d %d, level = %d)\n", b.index[0], b.index[1], b.index[2], b.level);
		printf("Neighbors are: \n");
		for(int i=0; i<neighbors.size(); i++)
			printf("neighbor %d: (index=%d %d %d, level = %d)\n", i, neighbors[i]->index[0], neighbors[i]->index[1], neighbors[i]->index[2], neighbors[i]->level);
		
		printf("end of neighbors.\n");
		abort();
	}

	//7-d. emotional shit
	if (bVerbose) printf("EVERYTHING WENT OK\n");
	
	WeightSet wSet;
	WeightToIndexMap mapW2I;
	vector<double>& weightsPool = bb.weightsPool;
	
	//8.
	for(int iPass=0; iPass<2; iPass++)
	{
		for(int i=0; i<bastards.size(); i++)
		{
			BastardGhost * ghost = bastards[i];
			const int ghost_index[3] =  {ghost->index.i[0], ghost->index.i[1], ghost->index.i[2]};
			
			const int d[3] = {
				(int)floor(ghost_index[0]/(double)block_size[0]),
				(int)floor(ghost_index[1]/(double)block_size[1]),
				(int)floor(ghost_index[2]/(double)block_size[2])
			};
			
			const int code = d[0]+1 + (d[1]+1)*3 + (d[2]+1)*9;
			
			const int s[3] = {
				d[0]<0? stencil_start[0] : (d[0]==0 ? 0 : B::sizeX),
				d[1]<0? stencil_start[1] : (d[1]==0 ? 0 : B::sizeY),
				d[2]<0? stencil_start[2] : (d[2]==0 ? 0 : B::sizeZ)
			};
			
			const int e[3] = {
				d[0]<0? 0 : (d[0]==0 ? B::sizeX : (B::sizeX + stencil_end[0] - 1)),
				d[1]<0? 0 : (d[1]==0 ? B::sizeY : (B::sizeY + stencil_end[1] - 1)),
				d[2]<0? 0 : (d[2]==0 ? B::sizeZ : (B::sizeZ + stencil_end[2] - 1))
			};
			
			const int size_chunk[3] = {e[0] - s[0], e[1] - s[1], e[2] - s[2]};
			const int rel_pos[3] = {ghost_index[0] - s[0], ghost_index[1] - s[1], ghost_index[2] - s[2]};
			const int index = bb.boundary[code].start + rel_pos[0] + rel_pos[1]*size_chunk[0] + rel_pos[2]*size_chunk[0]*size_chunk[1];
			
			assert(index - bb.boundary[code].start < bb.boundary[code].nGhosts);
			
			if (iPass==0)
				ghost->template collect<0>(wSet, mapW2I, bb.ghosts[index], weightsPool, 1,1,1,false);
			else
				ghost->template collect<1>(wSet, mapW2I, bb.ghosts[index], weightsPool, 1,1,1,false);
		}
		
		if (iPass==0)
		{
			//1. put wSet in the w pool
			//2. generate the mapping w2i
			
			const int iStart = weightsPool.size();
			const int nNewWeights = wSet.size();
			int c;
			
			//1.
			weightsPool.resize(iStart + nNewWeights);
			
			c = iStart;
			for(WeightSet::iterator it = wSet.begin(); it!=wSet.end(); it++)
				weightsPool[c++] = *it;
		
			//2.
			c = iStart + nNewWeights - 1;
			int counter = 0;
			for(vector<double>::const_reverse_iterator it = weightsPool.rbegin(); counter<nNewWeights; counter++, it++)
				mapW2I[*it] = c--;
		}
	}
	
	//9.
	for(typename GhostBuffer::iterator itV = buffer.begin(); itV!=buffer.end(); itV++)
	{
		for(typename GhostMap::iterator itM = itV->begin(); itM!=itV->end(); itM++)
		{
			allocator<BastardGhost>().destroy(itM->second);
			allocator<BastardGhost>().deallocate(itM->second,1);
			
			itM->second = NULL;
		}
		
		itV->clear();
	}
	buffer.clear();
}		
	
}

template <typename W, typename Real, bool bAnalysisFilter>
inline  const Real  * cookFilter(const int level_difference, const int  * filter_support)
{
	typedef _MRAG_GHOSTSCREATION_ALLOCATOR<Real> RealAllocator;
	
	const int offset = 0 - filter_support[0];
	const int start = filter_support[0];
	const int end = filter_support[1];
	
	const int sH = bAnalysisFilter? W::HaSupport[0] : W::HsSupport[0];
	const int eH = bAnalysisFilter? W::HaSupport[1] : W::HsSupport[1];
	
	RealAllocator alloc;
	
	Real * weights = alloc.allocate(filter_support[1] - filter_support[0]);
	Real * tmp = alloc.allocate(filter_support[1] - filter_support[0]);
	Real * H = alloc.allocate(eH-sH);

	int i, s, d, dl, eD, sD, eS,sS;
	Real * dest, *filt, * src; 
	
	//alone in the dark
	weights += offset;
	tmp += offset;
	H -= sH;
	
	for(i=sH; i<eH; i++) H[i] = bAnalysisFilter? W::getHa(i) : W::getHs(i);
	for(i=start; i<end; i++) weights[i] = 0.0;
	for(i=sH; i<eH; i++) weights[-i] = H[i];

	sS = -eH + 1;
	eS = -sH + 1;
	sD = 2*sS - eH + 1;
	eD = 2*eS - 1 - sH;

	for(dl = level_difference; dl>1; dl--)
	{
		for(d=sD, dest=tmp+sD; d<eD; d++) *dest++ = 0.0;
		
		for(s=sS, src=weights + sS; s<eS; s++, src++)
			for(d=sH, dest=tmp+2*s-eH+1, filt=H+eH-1; d<eH; d++)
					*dest++ += *src * *filt--;
		
		sS += sS - eH + 1;
		eS += eS - sH - 1;
		sD += sD - eH + 1;
		eD += eD - sH - 1;
		
		dest = tmp;
		tmp = weights;
		weights = dest;
	}
	
	weights -= offset;
	tmp -= offset;
	H += sH;
	
	alloc.deallocate(tmp, filter_support[1] - filter_support[0]);
	alloc.deallocate(H, eH-sH);
	
	return weights;
}

