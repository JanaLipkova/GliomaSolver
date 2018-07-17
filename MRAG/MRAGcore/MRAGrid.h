/*
 *  MRAGrid.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 4/23/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once
#undef min
#undef max
#include <vector>
#include <map>
#undef min
#undef max
#include <set>
using namespace std;

#include "MRAGEnvironment.h"
#include "MRAGGridNode.h"
#include "MRAGCompressor.h"
#include "MRAGRefinementPlan.h"
#include "MRAGCommon.h"
#include "MRAGBoundaryBlockInfo.h"
#include "MRAGBlockCollection.h"
#include "MRAGMatrix2D.h"
#include "MRAGBlockCollapser.h"
#include "MRAGBlockSplitter.h"

namespace MRAG
{
	class Refiner;
	typedef map<GridNode *, vector<GridNode *> > HierarchyType;
	typedef map<GridNode *, vector<GridNode *> > NeighborhoodType;
	
    /**
     * The grid.
     * Central structure and first that should be created in a simulation.
     * Manages blocks (using MRAG::BlockCollection) and all info needed to work with them.
     * @param WaveletType   Kind of wavelet used in this grid (used for compression, refinement...)
     * @param BlockType     Kind of block (incl. size of it) used in this grid.
     * @see MRAG::Block
     */
	template <typename WaveletType, typename BlockType>
	class Grid
	{	
	public:
		typedef BlockType GridBlockType;
		typedef typename BlockType::ElementType ElementType;
		
		//CLASS MANAGEMENT
        /**
         * Default constructor creating a grid with a certain number of blocks.
         * Corresponds to Grid(nBlocksX, nBlocksY, nBlocksZ, maxStencil[2][3], collection, bVerbose)
         * with maxStencil = {{-1,-1,-1}, {2,2,2}}.
         */
		Grid(int nBlocksX, int nBlocksY=1, int nBlocksZ=1, BlockCollection<BlockType>* collection = NULL, bool bVerbose=true);
        /**
         * Create a grid with a certain number of blocks and a maximal stencil for ghost-computations.
         * @param nBlocksX      Number of blocks in x-dimension.
         * @param nBlocksY      Number of blocks in y-dimension.
         * @param nBlocksZ      Number of blocks in z-dimension.
         * @param maxStencil    Maximal stencil used for computations (wavelet-stencil taken care of automatically).
         *                      Defines how many ghosts we will need for computations.
         * @param collection    (optional, default=NULL) MRAG::BlockCollection to be used in this grid (if NULL a fresh one is created).
         *                      The given collection (containing all blocks) will be resized (not reset!) if needed.
         * @param bVerbose      (optional, default=true) Enables verbose output (currently not used).
         */
		Grid(int nBlocksX, int nBlocksY, int nBlocksZ, const int maxStencil[2][3], BlockCollection<BlockType>* collection = NULL, bool bVerbose=true);
        /** Default destructor cleaning up allocated stuff (if a (non-NULL) collection was passed in the constructor, that one is not deleted). */
		~Grid();
		float getMemorySize() const;
		
		//TASK - ACCESS
		int getCurrentMaxLevel() const { return _computeMaxLevel(m_hierarchy); }
		int getCurrentMinLevel() const { return _computeMinLevel(m_blockAtLevel); }
		/** Return info on all blocks in the grid. @see MRAG::BlockInfo */
        vector<BlockInfo> getBlocksInfo() const;
        /** Return the collection of block defining all blocks in the grid. @see MRAG::BlockCollection */
		const BlockCollection<BlockType>& getBlockCollection() { return m_blockCollection; }
		BoundaryInfo& getBoundaryInfo(int * output_stencil_start=NULL, int * output_stencil_end=NULL);
		BoundaryInfo* createBoundaryInfo(const int requested_stencil_start[3], const int requested_stencil_end[3]) const;
		
		CompressionResult compress(const set<int>& blocksToCompress, int& nCollapsedBlocks);
		void setCompressor(Compressor * compressor);
		
		RefinementResult refine(const set<int>& blocksToRefine);
		void setRefiner(Refiner * refiner);
		
		void setBlockCollapser( BlockCollapser< WaveletType, BlockCollection<BlockType> > * collapser );
		void setBlockSplitter( BlockSplitter<  WaveletType, BlockCollection<BlockType> >* splitter  );
	
	protected:
		enum eGridStatus {eGridStatus_Locked, eGridStatus_Initialized};
		
		//fine tasks
		void _collapse(HierarchyType& hierarchy, GridNode * parent, int newBlockID) const;
		
		//TASK - const
		int _computeMaxLevel(const HierarchyType& hierarchy) const;
		int _computeMinLevel(const vector<vector<BlockInfo> >& blockAtLevel) const;
		void _computeNeighborhood(const HierarchyType& hierarchy, NeighborhoodType& neighborhood, map<GridNode *, map<int, GridNode *> > & ghostNodes) const;
		void _computeBlockAtLevel(const HierarchyType& hierarchy, vector<vector<BlockInfo> >& blockAtLevel) const;
		void _computeBoundaryInfo(BoundaryInfo& binfo, const int requested_stencil_start[3], const int requested_stencil_end[3], vector<GridNode*>& vNodesToCompute) const;
		void _computeMaxStencilUsed(int start[3], int end[3]) const;
		void _checkResolutionJumpCondition(int maxJump, const int * stencil_start = NULL, const int * stencil_end = NULL) const;
		int	 _computeMaxLevelJump() const;
		
		//TASK - non const
		inline void _refresh(bool bUpdateLazyData = false)
		{
			_computeBlockAtLevel(m_hierarchy, m_blockAtLevel);
			_computeNeighborhood(m_hierarchy, m_neighborhood, m_ghostNodes);
				
			m_mapGhost2Node.clear();
			for(map< GridNode *, map<int, GridNode *> >::const_iterator itGridNode = m_ghostNodes.begin(); itGridNode!=m_ghostNodes.end(); itGridNode++)
			{
				const map<int, GridNode *> & currentGhosts = itGridNode->second;
				for(map<int, GridNode *>::const_iterator itGhostNode = currentGhosts.begin(); itGhostNode!=currentGhosts.end(); itGhostNode++)
					m_mapGhost2Node[itGhostNode->second]=itGridNode->first;
			}
			
			if (bUpdateLazyData)
			{
				m_boundaryInfo.clear();

				m_setInvalidBBInfo.clear();
				
				for(HierarchyType::iterator it=m_neighborhood.begin(); it!=m_neighborhood.end(); it++)
					m_setInvalidBBInfo.insert(it->first);
			}
		}
		
		virtual void _setup(int nBlocksX, int nBlocksY, int nBlocksZ);
		
		void _dispose();
		
		//blocks info
		vector<vector<BlockInfo> > m_blockAtLevel;
		BlockCollection<BlockType>& m_blockCollection;
		bool m_bCollectionOwner;
		
		//node info
		HierarchyType m_hierarchy;
		NeighborhoodType m_neighborhood;
		map<GridNode *, map<int, GridNode *> > m_ghostNodes;
		map<GridNode *, GridNode *> m_mapGhost2Node;
		
		//boundary info
		int m_maxStencil[2][3];
		bool m_vPeriodicDirection[3];
		set<GridNode *> m_setInvalidBBInfo;
		BoundaryInfo m_boundaryInfo;
		
		//class info
		Real m_vPosition[3], m_vSize[3];
		bool m_vProcessingDirections[3];
		eGridStatus m_status;
		bool m_bVerbose;
		
		Refiner * m_refRefiner;
		Compressor *  m_refCompressor;
		BlockCollapser< WaveletType, BlockCollection<BlockType> >* m_blockCollapser;
		BlockSplitter<  WaveletType, BlockCollection<BlockType> >* m_blockSplitter;
		
		//friend class GridViewer;
		friend class SpaceTimeSorter;
		template <typename W, typename B>  friend class MRAG_STDTestL1;
		template <typename W, typename B>  friend class MRAG_STDTestL1_BoundaryInfo;
		template <typename W, typename B, typename P, int C> friend class IO_BaseClass;
		template <typename W, typename B, typename P, int C > friend class IO_Native;
		template <typename W, typename B, typename P, int C > friend class IO_Binary;
		template <typename W, typename B, typename P, int C > friend class IO_VTK;
	private:
		
		//forbidden
		Grid(const Grid&):m_blockAtLevel(), 
		m_blockCollection(*(new BlockCollection<BlockType>())), 
		m_bCollectionOwner(true),
		m_blockCollapser(),//blocks
		m_hierarchy(), m_neighborhood(), m_setInvalidBBInfo(), m_boundaryInfo(), m_mapGhost2Node(), m_ghostNodes(), //nodes
		m_status(eGridStatus_Initialized), m_bVerbose(true),
		m_refRefiner(NULL), m_refCompressor(NULL){abort();}
		
		Grid& operator =(const Grid&){abort(); return *this;}
	};
}

#include "MRAGrid.inl"

