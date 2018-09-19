/*
 *  MRAG_SmartBlockFinder.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 8/19/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once

namespace MRAG
{
	struct SmartBlockFinder
	{
		template<typename T> inline  static _MRAG_GHOSTSCREATION_ALLOCATOR<T> allocator()  { return _MRAG_GHOSTSCREATION_ALLOCATOR<T>();}
		
	public:
		static const bool bVerbose = false;
		
		int stencil_start[3], stencil_end[3], block_size[3];
		int minLevel, maxLevel, referenceLevel;
		int nNeighborhoodSize;
		int nSubStart[27];
		
		GridNode ** vNeighborhood;
		const GridNode * refCenterNode;
		PointIndex * result;
		
		SmartBlockFinder(const GridNode& b,  const vector<GridNode*>& neighbors, 
						 const int stencil_start[3], const int  stencil_end[3], const int block_size[3]);
		
		~SmartBlockFinder()
		{
			allocator<PointIndex>().deallocate(result,1);
			allocator<GridNode *>().deallocate(vNeighborhood, nNeighborhoodSize);
		}
		
		const GridNode * findLyingNode(const int level, const int index[3]) const;
		
		const int findLyingLevel(const int level, const int index[3]) const;
		
		PointIndex * findData(const int level, I3& i, const GridNode ** node=NULL) const;
		
	private:
		
		//forbidden
		SmartBlockFinder(const SmartBlockFinder&): 
		vNeighborhood(0), nNeighborhoodSize(0), minLevel(-1), maxLevel(-1), 
		referenceLevel(-1),  refCenterNode(NULL), result(NULL){abort();}
		
		SmartBlockFinder& operator=(const SmartBlockFinder&){abort(); return *this;}
	};

	inline SmartBlockFinder::SmartBlockFinder(const GridNode& b,  const vector<GridNode*>& neighbors,
									   const int stencil_start[3], const int  stencil_end[3], const int block_size[3]):
		vNeighborhood(NULL), nNeighborhoodSize(0), minLevel(b.level), maxLevel(b.level), referenceLevel(b.level), refCenterNode(&b), result(NULL)
	{
		result = allocator<PointIndex>().allocate(1);
		
		for(int i=0; i<neighbors.size(); i++)
		{
			const GridNode& n = *neighbors[i];
			
			minLevel = std::min(minLevel, n.level);
			maxLevel = std::max(maxLevel, n.level);
		}
		
		this->stencil_start[0] = stencil_start[0];
		this->stencil_start[1] = stencil_start[1];
		this->stencil_start[2] = stencil_start[2];
		
		this->stencil_end[0] = stencil_end[0];
		this->stencil_end[1] = stencil_end[1];
		this->stencil_end[2] = stencil_end[2];
		
		this->block_size[0] = block_size[0];
		this->block_size[1] = block_size[1];
		this->block_size[2] = block_size[2];
		
		//compute  substart
		int nTotalElements = 0;
		{
			const int length = (1<< maxLevel - b.level);
			
			for(int code=0; code<27; code++)
			{
				if (code == 1 + 3 + 9) continue;
				
				const int d[3] = {code%3 - 1, (code/3) %3 -1, (code/9) %3 -1};
				
				const int size_chunk[3] = {
					d[0]==0? (block_size[0]>1 ? length : 1) : (block_size[0]>1 ? 1 : 0),
					d[1]==0? (block_size[1]>1 ? length : 1) : (block_size[1]>1 ? 1 : 0),
					d[2]==0? (block_size[2]>1 ? length : 1) : (block_size[2]>1 ? 1 : 0)
				};
				
				const int nBunchOfElements =  size_chunk[0]*size_chunk[1]*size_chunk[2];

				nSubStart[code] = nTotalElements;

				nTotalElements += nBunchOfElements;
			}
		}

		
		//creation of the neighborhood
		{
			const int length = (1<< maxLevel - b.level);
			
			nNeighborhoodSize = 
					(block_size[0]>1 ? length+2 : 1)*(block_size[1]>1 ? length+2 : 1)*(block_size[2]>1 ? length+2 : 1) - 
					(block_size[0]>1 ? length : 1)*(block_size[1]>1 ? length : 1)*(block_size[2]>1 ? length : 1);

			assert(nNeighborhoodSize == nTotalElements);
			
			vNeighborhood = allocator<GridNode*>().allocate(nNeighborhoodSize);
			
			for(int i=0; i<nNeighborhoodSize; i++) 
				vNeighborhood[i] = NULL;

			const int s_b_finest[3] = {
				block_size[0]>1? (b.index[0] << maxLevel - b.level) : 0,
				block_size[1]>1? (b.index[1] << maxLevel - b.level) : 0,
				block_size[2]>1? (b.index[2] << maxLevel - b.level) : 0,
			};			
						
			for(int i=0; i<neighbors.size(); i++)
			{
				//1. map n o the finest level (b-space centered)
				//2. compute the intersection between n_finestlevel and the box (-1, length+1)
				//3. iterate over the intersection and put there the n-pointer
				
				GridNode& n = *neighbors[i];
				
				const int s_n_finest[3] = {
					block_size[0]>1? (n.index[0] << maxLevel - n.level) - s_b_finest[0] : 0,
					block_size[1]>1? (n.index[1] << maxLevel - n.level) - s_b_finest[1] : 0,
					block_size[2]>1? (n.index[2] << maxLevel - n.level) - s_b_finest[2] : 0
				};
				
				const int e_n_finest[3] = {
					block_size[0]>1? s_n_finest[0] + (1 << maxLevel - n.level) : 1,
					block_size[1]>1? s_n_finest[1] + (1 << maxLevel - n.level) : 1,
					block_size[2]>1? s_n_finest[2] + (1 << maxLevel - n.level) : 1
				};
				
				const int s_intersection[3] = {
					std::max(s_n_finest[0], -1),
					std::max(s_n_finest[1], -1),
					std::max(s_n_finest[2], -1),
				};
				
				const int e_intersection[3] = {
					std::min(e_n_finest[0], length + 1),
					std::min(e_n_finest[1], length + 1),
					std::min(e_n_finest[2], length + 1),
				};
				
				int d[3];
				for(d[2]=s_intersection[2]; d[2]<e_intersection[2]; d[2]++)
				for(d[1]=s_intersection[1]; d[1]<e_intersection[1]; d[1]++)
				for(d[0]=s_intersection[0]; d[0]<e_intersection[0]; d[0]++)
				{
					const int dir[3] = {
						(int)floor(d[0]/(double)length),
						(int)floor(d[1]/(double)length),
						(int)floor(d[2]/(double)length)
					};
					
					const int code = (dir[2]+1)*9 + (dir[1]+1)*3 + (dir[0]+1);
					
					const int size_chunk[3] = {
						dir[0]==0? (block_size[0]>1 ? length : 1) : (block_size[0]>1 ? 1 : 0),
						dir[1]==0? (block_size[1]>1 ? length : 1) : (block_size[1]>1 ? 1 : 0),
						dir[2]==0? (block_size[2]>1 ? length : 1) : (block_size[2]>1 ? 1 : 0)
					};
					
					const int start_chunk[3] = {
						dir[0]<0? -1 : (dir[0]==0 ? 0 : length),
						dir[1]<0? -1 : (dir[1]==0 ? 0 : length),
						dir[2]<0? -1 : (dir[2]==0 ? 0 : length)
					};
					
					const int dest_pos[3] = {
						d[0] - start_chunk[0],
						d[1] - start_chunk[1],
						d[2] - start_chunk[2],
					};
					
					const int row_size = size_chunk[0];
					const int slice_size = size_chunk[0]*size_chunk[1];
					const int dest_index =  nSubStart[code] + dest_pos[0] + row_size*dest_pos[1] + slice_size*dest_pos[2];
					
					assert(nSubStart[code]<= dest_index && dest_index <  nSubStart[code] + size_chunk[0]*size_chunk[1]*size_chunk[2]);
			
					vNeighborhood[ dest_index ] = &n;
				}
			}
		}
		
		for(int c=0; c<27; c++)
		{
			const int length = (1<< maxLevel - b.level);
			const int dir[3] = {
				(c%3) - 1, (c/3 %3) - 1, (c/9 %3) - 1
			};
			
			const int size_chunk[3] = {
				dir[0]==0? (block_size[0]>1 ? length : 1) : (block_size[0]>1 ? 1 : 0),
				dir[1]==0? (block_size[1]>1 ? length : 1) : (block_size[1]>1 ? 1 : 0),
				dir[2]==0? (block_size[2]>1 ? length : 1) : (block_size[2]>1 ? 1 : 0)
			};
			
			const int n= size_chunk[0]*size_chunk[1]*size_chunk[2];
	
			for(int i=nSubStart[c]; i<n+nSubStart[c]; i++)
			{
				if (c == 1+3+9) continue;
				assert(vNeighborhood[i]!= NULL);
				if (vNeighborhood[i]!= NULL && c!=13) if (bVerbose) printf("vNeighborhood[%d]->level = %d\n", i, vNeighborhood[i]->level);
			}
		}
	}
	

	
	inline const GridNode * SmartBlockFinder::findLyingNode(const int level, const int index[3]) const
	{
		//1. map the point to the finest level
		//2. assert that it falls inside the feasible area
		//3. if the code is 9+3+1 return b.l
		//4. otherwise look in the lut neighborhood
		
		const int level_difference = maxLevel - level;
		const int length = 1<< maxLevel - refCenterNode->level;
		
		//1.
		const int candidate_i[3] = {
			block_size[0]>1? (int)floor((index[0]<<level_difference)/(double)block_size[0]) : 0,
			block_size[1]>1? (int)floor((index[1]<<level_difference)/(double)block_size[1]) : 0,
			block_size[2]>1? (int)floor((index[2]<<level_difference)/(double)block_size[2]) : 0
		};
		

		//2.
		const int i[3] = {
			candidate_i[0]<-1?-1 : (candidate_i[0]>length? length: candidate_i[0]),
			candidate_i[1]<-1?-1 : (candidate_i[1]>length? length: candidate_i[1]),
			candidate_i[2]<-1?-1 : (candidate_i[2]>length? length: candidate_i[2])
		};
		
		assert(i[0]>= -1 && i[0]<length+1);
		assert(i[1]>= -1 && i[1]<length+1);
		assert(i[2]>= -1 && i[2]<length+1);
		
		//3.
		const int dir[3] = {
			floor(i[0]/(double)length),
			floor(i[1]/(double)length),
			floor(i[2]/(double)length)
		};

		const int code = (dir[2]+1)*9 + (dir[1]+1)*3 + (dir[0]+1);
		
		if (code == 1 + 3 + 9) return refCenterNode;
		
		//4.
		const int size_chunk[3] = {
			dir[0]==0? (block_size[0]>1 ? length : 1) : (block_size[0]>1 ? 1 : 0),
			dir[1]==0? (block_size[1]>1 ? length : 1) : (block_size[1]>1 ? 1 : 0),
			dir[2]==0? (block_size[2]>1 ? length : 1) : (block_size[2]>1 ? 1 : 0)
		};
		
		const int s[3] = {
			dir[0]<0? -1 : (dir[0]==0 ? 0 : length),
			dir[1]<0? -1 : (dir[1]==0 ? 0 : length),
			dir[2]<0? -1 : (dir[2]==0 ? 0 : length)
		};
		
		const int row_size = size_chunk[0];
		const int slice_size = size_chunk[0]*size_chunk[1];
		
		GridNode ** destination = vNeighborhood + nSubStart[code];
		
		const GridNode * result = destination[ (i[0]-s[0]) + row_size*(i[1]-s[1]) + slice_size*(i[2]-s[2]) ];
		
		//assert(result->level == level);
		const int level_diff_gridnode = maxLevel - result->level;
		const int level_diff_centernode = maxLevel - refCenterNode->level;
		for(int i=0; i<3; i++)
		{
			assert((result->index[i] << level_diff_gridnode)*block_size[i] <= (index[i]<<level_difference)  + (refCenterNode->index[i]<<level_diff_centernode)*block_size[i] || block_size[i]==1);
			assert((result->index[i]+1 << level_diff_gridnode)*block_size[i] > (index[i]<<level_difference) + (refCenterNode->index[i]<<level_diff_centernode)*block_size[i] || block_size[i]==1);
		}
		
		return result;
	}
	
	inline const int SmartBlockFinder::findLyingLevel(const int level, const int index[3]) const
	{
		return findLyingNode(level, index)->level;
	}
	
	inline PointIndex * SmartBlockFinder::findData(const int level, I3& ghost, const GridNode ** node_output) const
	{
		//1. ask for the lying level
		//2. if it is not the same, return null
		//3. compute the point index
		
		const GridNode * node = findLyingNode(level, ghost.i);
		const int lying_level = node->level;
		
		if (node_output != NULL) *node_output = node;
		
		const double dilate = pow(0.5, refCenterNode->level- node->level);
		
		if (lying_level != level) return NULL;
		
		const int d[3] = {
			ghost.i[0] - block_size[0]*node->index[0] + (int)floor(dilate*block_size[0]*refCenterNode->index[0]),
			ghost.i[1] - block_size[1]*node->index[1] + (int)floor(dilate*block_size[1]*refCenterNode->index[1]),
			ghost.i[2] - block_size[2]*node->index[2] + (int)floor(dilate*block_size[2]*refCenterNode->index[2])
		};
				
		assert(d[0]<block_size[0]);
		assert(d[1]<block_size[1]);
		assert(d[2]<block_size[2]);
		
		assert(d[0]>=0);
		assert(d[1]>=0);
		assert(d[2]>=0);
		
		result->index = d[0] + d[1]*block_size[0] +  d[2]*block_size[0]*block_size[1];
		assert(result->index<block_size[0]*block_size[1]*block_size[2]);
		result->blockID = node->blockID;
		
		return result;
	}
}
