
#include "MRAGrid.h"
#include "MRAGMatrix3D.h"
#include "MRAGridHelpers.h"
//#include "MRAGBoundaryBlockHelpers.h"
#include "MRAG_BBInfoCreator.h"
#include "MRAGRefiner.h"

#ifdef _MRAG_TBB
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#endif

namespace MRAG
{
	
#pragma mark -
#pragma mark Class Management
	
	template <typename WaveletType, typename BlockType>
	Grid<WaveletType, BlockType>::Grid(int nBlocksX, int nBlocksY, int nBlocksZ, BlockCollection<BlockType>* collection, bool bVerbose):
		m_blockAtLevel(), 
		m_blockCollection((collection != NULL)? *collection : *(new BlockCollection<BlockType>())), 
		m_bCollectionOwner(collection == NULL),
		m_blockCollapser(), m_blockSplitter(),
		m_hierarchy(), m_neighborhood(), m_setInvalidBBInfo(), m_boundaryInfo(), m_mapGhost2Node(), m_ghostNodes(), //nodes
		m_status(eGridStatus_Initialized), m_bVerbose(bVerbose),
		m_refRefiner(NULL), m_refCompressor(NULL)
	{
		m_blockCollapser = new BlockCollapser< WaveletType, BlockCollection<BlockType> >;
		m_blockSplitter = new BlockSplitter<  WaveletType, BlockCollection<BlockType> >;
		
		m_vProcessingDirections[0] = BlockType::shouldProcessDirectionX;
		m_vProcessingDirections[1] = BlockType::shouldProcessDirectionY;
		m_vProcessingDirections[2] = BlockType::shouldProcessDirectionZ;
		
		m_vPeriodicDirection[0] = true;
		m_vPeriodicDirection[1] = true;
		m_vPeriodicDirection[2] = true;
		
		m_vPosition[0] = 0;
		m_vPosition[1] = 0;
		m_vPosition[2] = 0;
		
		m_vSize[0] = 1;
		m_vSize[1] = 1;
		m_vSize[2] = 1;
		
		m_maxStencil[0][0] = m_maxStencil[0][1] = m_maxStencil[0][2] = -1;
		m_maxStencil[1][0] = m_maxStencil[1][1] = m_maxStencil[1][2] = 2;

		
		_setup(nBlocksX, nBlocksY, nBlocksZ);
	}

	template <typename WaveletType, typename BlockType>
	Grid<WaveletType, BlockType>::Grid(int nBlocksX, int nBlocksY, int nBlocksZ, const int maxStencil[2][3], BlockCollection<BlockType>* collection, bool bVerbose):
		m_blockAtLevel(), 
		m_blockCollection((collection != NULL)? *collection : *(new BlockCollection<BlockType>())), 
		m_bCollectionOwner(collection == NULL),
		m_blockCollapser(NULL), m_blockSplitter(NULL),
		m_hierarchy(), m_neighborhood(), m_setInvalidBBInfo(), m_boundaryInfo(), m_mapGhost2Node(), m_ghostNodes(), //nodes
		m_status(eGridStatus_Initialized), m_bVerbose(bVerbose),
		m_refRefiner(NULL), m_refCompressor(NULL)
	{
		m_blockCollapser = new BlockCollapser< WaveletType, BlockCollection<BlockType> >;
		m_blockSplitter = new BlockSplitter<  WaveletType, BlockCollection<BlockType> >;
		
		m_vProcessingDirections[0] = BlockType::shouldProcessDirectionX;
		m_vProcessingDirections[1] = BlockType::shouldProcessDirectionY;
		m_vProcessingDirections[2] = BlockType::shouldProcessDirectionZ;
		
		m_vPeriodicDirection[0] = true;
		m_vPeriodicDirection[1] = true;
		m_vPeriodicDirection[2] = true;
		
		m_vPosition[0] = 0;
		m_vPosition[1] = 0;
		m_vPosition[2] = 0;
		
		m_vSize[0] = 1;
		m_vSize[1] = 1;
		m_vSize[2] = 1;
		
		for(int i=0; i<2; i++)
			for(int j=0; j<3; j++)
				m_maxStencil[i][j] = maxStencil[i][j];

		_setup(nBlocksX, nBlocksY, nBlocksZ);
	}
		
	template <typename WaveletType, typename BlockType>
	void Grid<WaveletType, BlockType>::_setup(int nBlocksX, int nBlocksY, int nBlocksZ)
	{	
		int nBlocks[3] = {nBlocksX, nBlocksY, nBlocksZ};

		//find l
		int iMaxLevel = MRAG::MRAGridHelpers::_computeDesiredLevel(nBlocks);
		
		//create blocks
		std::vector<int> vIDs = m_blockCollection.create(nBlocksX*nBlocksY*nBlocksZ);
		assert(vIDs.size() > 0);
		
		//create hierarchy
		m_hierarchy.clear();
		MRAG::MRAGridHelpers::_createChildren(vIDs, m_hierarchy, NULL, 0, nBlocks, 0, iMaxLevel);
		assert(vIDs.size() == 0);
		
		//create blockAtLevel, neighborhood
		_refresh(true);
	}
		
	template <typename WaveletType, typename BlockType>
	Grid<WaveletType, BlockType>::~Grid()
	{
		//1. dispose data
		//2. delete block collection
		//3. set pointers to NULL

		//1.
		_dispose();

		//2.
		if (m_bCollectionOwner)
			delete &m_blockCollection;

		//3.
		m_refRefiner = NULL;
		m_refCompressor = NULL;
	}

	template <typename WaveletType, typename BlockType>
	void Grid<WaveletType, BlockType>::_dispose()
	{
		//1. release blocks
		//2. deallocate GridNodes in hierarchy
		//3. deallocate Ghost GridNodes in ghosts
		//4. deallocate boundary info
		//5. clear the data structures

		//1.
		for(NeighborhoodType::const_iterator it = m_neighborhood.begin(); it != m_neighborhood.end(); it++)
			m_blockCollection.erase(it->first->blockID);

		//2.
		for(HierarchyType::iterator it=m_hierarchy.begin(); it!=m_hierarchy.end(); it++)
		{
			const GridNode * node = it->first;
			if(node!=NULL) delete node;
			
			it->second.clear();
		}
		
		//3.
		for(map<GridNode *, map<int, GridNode *> >::iterator it=m_ghostNodes.begin(); it!=m_ghostNodes.end(); it++)
			for(map<int, GridNode *>::iterator itGhosts=it->second.begin(); itGhosts!=it->second.end(); itGhosts++)
				delete itGhosts->second;

		//4.
		m_setInvalidBBInfo.clear();
		m_boundaryInfo.clear();
		
		//5.
		m_hierarchy.clear();
		m_neighborhood.clear();
		m_ghostNodes.clear();
		m_blockAtLevel.clear();
		
	}

	template <typename WaveletType, typename BlockType>
	float Grid<WaveletType, BlockType>::getMemorySize() const
	{
		int memsize = sizeof(this);
		
		for(int i=0;i<m_blockAtLevel.size(); i++)
			memsize += sizeof(int)*m_blockAtLevel[i].size();
		
		memsize += 1024*1024*m_blockCollection.getMemorySize();
		
		memsize+= m_hierarchy.size()*(sizeof(HierarchyType::key_type) + sizeof(HierarchyType::mapped_type));
		
		memsize+= m_neighborhood.size()*(sizeof(NeighborhoodType::key_type) + sizeof(NeighborhoodType::mapped_type));
		
		memsize+= 1024*1024*m_boundaryInfo.getMemorySize();
		
		return memsize/(double)(1<<20);
	}
		
	template <typename WaveletType, typename BlockType>
	vector<BlockInfo> Grid<WaveletType, BlockType>::getBlocksInfo() const 
	{
		vector<BlockInfo> result;
		
		for(typename HierarchyType::const_iterator it=m_hierarchy.begin(); it!=m_hierarchy.end(); it++)
		{
			const GridNode * node = it->first;
			if(node!=NULL && !node->isEmpty)
			{
				const double dilate = pow(2.0, -node->level);
				
				const Real h[3] = {
					dilate*m_vSize[0]/BlockType::sizeX, 
					m_vProcessingDirections[1]? dilate*m_vSize[1]/BlockType::sizeY  : 1,
					m_vProcessingDirections[2]? (dilate*m_vSize[2]/BlockType::sizeZ) : 1
				};

				const Real p[3] = {
					m_vPosition[0] + node->index[0]*dilate*m_vSize[0] + (WaveletType::bIsCellCentered?0.5*h[0] : 0),
					m_vPosition[1] + node->index[1]*dilate*m_vSize[1] + (WaveletType::bIsCellCentered?0.5*h[1] : 0),
					m_vPosition[2] + node->index[2]*dilate*m_vSize[2] + (WaveletType::bIsCellCentered?0.5*h[2] : 0)
				};
				
				result.push_back(BlockInfo(node->blockID, node->index, node->level, p, h));
			}
		}
		
		return result;
	}
	
#pragma mark -
#pragma mark Refinement

	template <typename WaveletType, typename BlockType>
	void Grid<WaveletType, BlockType>::setBlockSplitter( BlockSplitter<  WaveletType, BlockCollection<BlockType> >* splitter )
	{
		if (splitter == NULL)
		{
			printf("void Grid<WaveletType, BlockType>::setBlockSplitter: NULL as arguments!\n");
			abort();
		}
		
		assert(m_blockSplitter != NULL);
		delete m_blockSplitter;
		
		m_blockSplitter = splitter;
	}
	
	template <typename WaveletType, typename BlockType>
	void Grid<WaveletType, BlockType>::setRefiner(Refiner * refiner)
	{
		_checkResolutionJumpCondition(refiner->getMaxLevelJump());
		
		m_refRefiner=refiner;
	}
	
	template <typename WaveletType, typename BlockType>
	RefinementResult Grid<WaveletType, BlockType>::refine(const set<int>& blocksToRefine) 
	{
		if (blocksToRefine.size() == 0) return RefinementResult();
		
		//1. iterate over the nodes, look which one has an associated block to be refined
		//2. create a refinement plan
		//3. realize the plan
		//4. update the hieararchy, update boundary info
		//5. update the invalid node set
		//6. return how many blocks we created
		
		//1.
		for(typename HierarchyType::const_iterator it=m_hierarchy.begin(); it!=m_hierarchy.end(); it++)
		{
			GridNode * node = it->first;
			
			if (node==NULL || node->isEmpty ) continue;
			
			const bool bFound = blocksToRefine.find(node->blockID) != blocksToRefine.end();
			node->shouldBeRefined = (node->isEmpty)? false: bFound;
		}
		
		//2.
		vector<NodeToRefine> vToRefine;
		RefinementPlan * plan = m_refRefiner->createPlan(m_hierarchy, m_neighborhood, m_vProcessingDirections, vToRefine);
		if (plan == NULL)
			return RefinementResult(true);
		
		//3.
		vector<RefinementReport> vRefinementReport;
		RefinementResult result = m_blockSplitter->split(m_blockCollection, getBoundaryInfo(), *plan, vRefinementReport);
		//printf("vRefinementReport size: %d\n", vRefinementReport.size());
		
		//4.
		vector<GridNode *> newNodes;
		{
			const int nRefinements = plan->refinements.size();
			for(int i=0; i<nRefinements; i++)
			{
				GridNode * parent =  const_cast<GridNode *>(vToRefine[i].node);
				
				vector<RefinementPlanNode *>& vNodes = plan->refinements[i]->children;
				vector<int>& blockIDs = vRefinementReport[i].childrenIDs;
				
				int nChildren = vNodes.size();
				
				for(int c=0; c<nChildren; c++)
				{
					GridNode * node = new GridNode(false, parent, blockIDs[c], 
												   vNodes[c]->index[0], vNodes[c]->index[1], vNodes[c]->index[2], vNodes[c]->level);
					
					m_hierarchy[parent].push_back(node);
					m_hierarchy[node]=vector<GridNode*>();
					
					newNodes.push_back(node);
				}
				
				
				parent->isEmpty = true;
				parent->shouldBeRefined = false;
				
				m_boundaryInfo.erase(parent->blockID);
				parent->blockID = -1;
			}
		}
		delete plan;
		
		_refresh();
		
		//5.
		set<int> blockToErase;
		for(vector<GridNode *>::const_iterator itNewNode= newNodes.begin(); itNewNode!=newNodes.end(); itNewNode++)
		{
			m_setInvalidBBInfo.insert(*itNewNode);
			
			NeighborhoodType::const_iterator itNeighbors = m_neighborhood.find(*itNewNode);
			
			assert(itNeighbors != m_neighborhood.end());
			
			const vector<GridNode*>& neighbors = itNeighbors->second;
			
			for(vector<GridNode*>::const_iterator it = neighbors.begin(); it!=neighbors.end(); it++)
			{
				GridNode * node = NULL;
				
				NeighborhoodType::const_iterator itN =  m_neighborhood.find(*it);
				
				if (itN != m_neighborhood.end())
					node = itN->first;
				else
				{
					map<GridNode*, GridNode*>::const_iterator itG = m_mapGhost2Node.find(*it);
					
					assert(itG != m_mapGhost2Node.end());
					
					node = itG->second;
				}
				
				blockToErase.insert(node->blockID);
				m_setInvalidBBInfo.insert(node);
			}
		}
		
		for(set<int>::const_iterator it=blockToErase.begin(); it!=blockToErase.end(); it++)
			m_boundaryInfo.erase(*it, false);
		//6.
		return result;
	}
	
#pragma mark -
#pragma mark Compression
	
	template <typename WaveletType, typename BlockType>
	void Grid<WaveletType, BlockType>::setBlockCollapser( BlockCollapser< WaveletType, BlockCollection<BlockType> > * collapser )
	{
		if (collapser == NULL)
		{
			printf("void Grid<WaveletType, BlockType>::setBlockCollapser: NULL as arguments!\n");
			abort();
		}
		
		assert(m_blockCollapser != NULL);
		delete m_blockCollapser;
		
		m_blockCollapser = collapser;
	}
	
	template <typename WaveletType, typename BlockType>
	void Grid<WaveletType, BlockType>::setCompressor(Compressor * compressor)
	{
		_checkResolutionJumpCondition(compressor->getMaxLevelJump());
		
		m_refCompressor=compressor;
	}
		
	template <typename WaveletType, typename BlockType>
	CompressionResult Grid<WaveletType, BlockType>::compress(const set<int>& blocksToCompress, int& nCollapsedBlocks)
	{
		//1. iterate over full nodes, fill the compress flag
		//2. create a compression plan
		//3. collapse the blocks
		//4. generate map CollapseID -> BlockID
		//5. update boundaryinfo and the invalid set of BBInfo
		//6. update hieararchy
		//7. refresh the other data structures depending on hierarchy
		
		//1.
		for(typename HierarchyType::const_iterator it=m_hierarchy.begin(); it!=m_hierarchy.end(); it++)
		{
			GridNode * node = it->first;
			if(node!=NULL)
					node->shouldBeCompressed = (node->isEmpty)? false: (blocksToCompress.find(node->blockID) != blocksToCompress.end());
		}
		
		//2.
		assert(m_refCompressor != NULL);
		vector<GridNodeCollapseInfo> vToCollapse;
		CompressionPlan * plan = m_refCompressor->createPlan(m_hierarchy, m_neighborhood, vToCollapse);
		
		//3.
		vector<BlockCollapseInfo> vBlockCollapseInfo;
		CompressionResult result = m_blockCollapser->collapse(m_blockCollection, getBoundaryInfo(),  *plan, vBlockCollapseInfo);
		
		//4
		map<int, int> mapCollapseIDBlockID;
		for(int i=0; i<vBlockCollapseInfo.size(); i++)
			mapCollapseIDBlockID[vBlockCollapseInfo[i].collapseID] = vBlockCollapseInfo[i].newBlockID;
		
		//5.
		vector<GridNode *> newNodes;
		newNodes.reserve(vToCollapse.size());
		for(int i=0; i<vToCollapse.size(); i++)
		{
			vector<GridNode *>& v = m_hierarchy[vToCollapse[i].node];
			
			for(vector<GridNode *>::iterator it=v.begin(); it!=v.end(); it++)
			{
				m_boundaryInfo.erase((*it)->blockID);
				(*it)->blockID = -1;
			}
			
			newNodes.push_back(vToCollapse[i].node);
		}
		
		//6.
		for(int i=0; i<vToCollapse.size(); i++)
			_collapse(m_hierarchy, vToCollapse[i].node, mapCollapseIDBlockID[vToCollapse[i].collapseID]);
		delete plan;
		
		//7.
		_refresh();
		
		/*for(vector<GridNode *>::const_iterator itNewNode= newNodes.begin(); itNewNode!=newNodes.end(); itNewNode++)
		{
			vector<GridNode*>& neighbors = m_neighborhood[*itNewNode];
			
			m_setInvalidBBInfo.insert(*itNewNode);
			m_setInvalidBBInfo.insert(neighbors.begin(), neighbors.end());
		}*/
		set<int> blockToErase;
		for(vector<GridNode *>::const_iterator itNewNode= newNodes.begin(); itNewNode!=newNodes.end(); itNewNode++)
		{
			m_setInvalidBBInfo.insert(*itNewNode);
			
			NeighborhoodType::const_iterator itNeighbors = m_neighborhood.find(*itNewNode);
			
			assert(itNeighbors != m_neighborhood.end());
			
			const vector<GridNode*>& neighbors = itNeighbors->second;
			
			for(vector<GridNode*>::const_iterator it = neighbors.begin(); it!=neighbors.end(); it++)
			{
				GridNode * node = NULL;
				
				NeighborhoodType::const_iterator itN =  m_neighborhood.find(*it);
				
				if (itN != m_neighborhood.end())
					node = itN->first;
				else
				{
					map<GridNode*, GridNode*>::const_iterator itG = m_mapGhost2Node.find(*it);
					
					assert(itG != m_mapGhost2Node.end());
					
					node = itG->second;
				}
				
				blockToErase.insert(node->blockID);
				m_setInvalidBBInfo.insert(node);
			}
		}
		
		for(set<int>::const_iterator it=blockToErase.begin(); it!=blockToErase.end(); it++)
			m_boundaryInfo.erase(*it, false);
		
		nCollapsedBlocks = vBlockCollapseInfo.size();
		
		return result;
	}

	template <typename WaveletType, typename BlockType>
	void Grid<WaveletType, BlockType>::_collapse(HierarchyType& hierarchy, GridNode * parent, int newBlockID) const
	{
		vector<GridNode*>& children = hierarchy[parent];
		
		for(int i=0; i<children.size(); i++)
		{
			assert(children[i]->shouldBeCompressed == true && children[i]->isEmpty == false);
		}
		
		for(int i=0; i<children.size(); i++)
		{
			GridNode * ptr = children[i];
			hierarchy.erase(children[i]);
			delete ptr;
		}
		
		hierarchy[parent].clear();
		
		parent->blockID = newBlockID;
		parent->isEmpty = false;
	}

#pragma mark -
#pragma mark Boundary Info Computation
	
	template <typename WaveletType, typename BlockType>
	BoundaryInfo& Grid<WaveletType, BlockType>::getBoundaryInfo(int * output_stencil_start, int * output_stencil_end)
	{
		//0. couple of checks
		//1. compute the stencil used in the grid
		//2. fill the output
		//3. if m_bAlreadyComputedBoundaryInfo is true, it means that the work has already been done
		//4. otherwise go and do it
		
		//0.
		if (m_refRefiner == NULL)
			_checkResolutionJumpCondition(_computeMaxLevelJump());
		
		//1.
		int stencil_start[3] = {0,0,0};
		int stencil_end[3] = {0,0,0};
		_computeMaxStencilUsed(stencil_start, stencil_end);
		
		//2.
		if (output_stencil_start != NULL)
		{
			output_stencil_start[0] = stencil_start[0];
			output_stencil_start[1] = stencil_start[1];
			output_stencil_start[2] = stencil_start[2];
		}
		
		if (output_stencil_end != NULL)
		{
			output_stencil_end[0] = stencil_start[0];
			output_stencil_end[1] = stencil_start[1];
			output_stencil_end[2] = stencil_start[2];
		}
		
		//3.
		if (m_setInvalidBBInfo.size() == 0) return m_boundaryInfo;
		//printf("Found %d invalid blocks\n", m_setInvalidBBInfo.size() );
		
		//4.
		vector<GridNode*> vNodes(m_setInvalidBBInfo.size());
		copy(m_setInvalidBBInfo.begin(), m_setInvalidBBInfo.end(), vNodes.begin());
		
		_computeBoundaryInfo(m_boundaryInfo, stencil_start, stencil_end, vNodes);
		
		m_setInvalidBBInfo.clear();
		
		return m_boundaryInfo;
	}
	
	template <typename WaveletType, typename BlockType>
	BoundaryInfo* Grid<WaveletType, BlockType>::createBoundaryInfo(const int requested_stencil_start[3], const int requested_stencil_end[3]) const
	{
		_checkResolutionJumpCondition(_computeMaxLevelJump(), requested_stencil_start, requested_stencil_end);
		
		BoundaryInfo * binfo = new BoundaryInfo;
		
		vector<GridNode*> vNodes(m_neighborhood.size());
		vector<GridNode*>::iterator itVector = vNodes.begin();
		
		for(HierarchyType::const_iterator it=m_hierarchy.begin(); it!=m_hierarchy.end(); it++)
		{
			if (it->first == NULL || it->first->isEmpty) continue;
			
			*itVector = it->first;
			itVector++;
		}
		
		_computeBoundaryInfo(*binfo, requested_stencil_start, requested_stencil_end, vNodes);
		
		return binfo;
	}
	
	template <typename WaveletType, typename BlockType>
	struct Body_CreateBoundaryInfo
	{
		struct ParallelItem
		{
			const GridNode* node;
			const vector<GridNode*>& neighbors;
			BoundaryInfoBlock * bbi;
			
			ParallelItem(GridNode* node_,  const vector<GridNode*>& neighbors_): node(node_), neighbors(neighbors_), bbi(NULL) {}
		};
		
		vector<ParallelItem*> vWorkingList;
		BoundaryInfo& m_binfo; 
		const int (&requested_stencil_start)[3];
		const int (&requested_stencil_end)[3];
		
		Body_CreateBoundaryInfo(const NeighborhoodType& neighborhood, BoundaryInfo& binfo, 
								const int requested_stencil_start_[3], const int requested_stencil_end_[3], vector<GridNode*>& vNodesToCompute ):
		m_binfo(binfo), vWorkingList(), 
		requested_stencil_start((const int (&)[3])requested_stencil_start_),
		requested_stencil_end((const int (&)[3])requested_stencil_end_)
		{
			
			vWorkingList.reserve(vNodesToCompute.size());
			
			for(vector<GridNode*>::const_iterator it = vNodesToCompute.begin(); it!=vNodesToCompute.end(); it++)
			{
				GridNode * node = *it;
				
				NeighborhoodType::const_iterator itNeighbors = neighborhood.find(node);
				
				assert(itNeighbors != neighborhood.end());
				
				map<int, BoundaryInfoBlock*>::iterator itBBI = binfo.boundaryInfoOfBlock.find(node->blockID);
				
				if (itBBI!= binfo.boundaryInfoOfBlock.end()) delete itBBI->second;
				
				vWorkingList.push_back(new ParallelItem(node, itNeighbors->second));
			}
		}
		
		template<typename BlockedRange>
		void operator()(const BlockedRange& r) const
		{
			const int block_size[3] = {
				BlockType::sizeX,
				BlockType::sizeY,
				BlockType::sizeZ
			};
			
			const int requested_stencil_start[3] = {
				m_binfo.stencil_start[0] ,
				m_binfo.stencil_start[1] ,
				m_binfo.stencil_start[2]
			};
			
			const int requested_stencil_end[3] = {
				m_binfo.stencil_end[0] ,
				m_binfo.stencil_end[1] ,
				m_binfo.stencil_end[2]
			};
			
			
			MRAG_BBInfoCreator<WaveletType, BlockType> bbcreator(requested_stencil_start, requested_stencil_end, block_size);
			//printf("=============:: %d %d\n", r.begin(), r.end());
			for(int i=r.begin(); i!=r.end(); i++)
				vWorkingList[i]->bbi = bbcreator.createBoundaryInfoBlock(*vWorkingList[i]->node, vWorkingList[i]->neighbors);
		}	
		
		void collect()
		{		
			for(typename vector<ParallelItem*>::const_iterator it = vWorkingList.begin(); it != vWorkingList.end(); it++)
			{
				const int blockID = (*it)->node->blockID;
				
				map<int, BoundaryInfoBlock*>::iterator itBBI = m_binfo.boundaryInfoOfBlock.find(blockID);
				
				if (itBBI== m_binfo.boundaryInfoOfBlock.end())
					m_binfo.boundaryInfoOfBlock[blockID] = (*it)->bbi;
				else
					itBBI->second = (*it)->bbi;
			}
			
			for(typename vector<ParallelItem*>::iterator it = vWorkingList.begin(); it != vWorkingList.end(); it++)
			{
				delete *it;
				*it = NULL;
			}
		}
	};
	
	template <typename WaveletType, typename BlockType>
	void Grid<WaveletType, BlockType>::_computeBoundaryInfo(BoundaryInfo& binfo, const int requested_stencil_start[3], const int requested_stencil_end[3], vector<GridNode*>& vNodesToCompute ) const
	{
		binfo.stencil_start[0] = requested_stencil_start[0];
		binfo.stencil_start[1] = requested_stencil_start[1];
		binfo.stencil_start[2] = requested_stencil_start[2];
		
		binfo.stencil_end[0] = requested_stencil_end[0];
		binfo.stencil_end[1] = requested_stencil_end[1];
		binfo.stencil_end[2] = requested_stencil_end[2];
		
		const int n = vNodesToCompute.size();
		Body_CreateBoundaryInfo<WaveletType, BlockType> body(m_neighborhood, binfo, requested_stencil_start, requested_stencil_end, vNodesToCompute);

#ifdef _MRAG_TBB
		//printf("(Grid::_computeBoundaryInfo: TBB)\n");
		const int nThreads = _MRAG_TBB_NTHREADS_HINT;
		tbb::parallel_for(tbb::blocked_range<size_t>(0, n,std::max(1,n/nThreads)), body);
#else
		body(SimpleInterval(0, n));
#endif
		
		body.collect();
	}

#pragma mark -
#pragma mark Helpers const

	template <typename WaveletType, typename BlockType>
	int Grid<WaveletType, BlockType>::_computeMaxLevel(const HierarchyType& hierarchy) const
	{
		int lmax = 0;
		
		for(typename HierarchyType::const_iterator it=hierarchy.begin(); it!=hierarchy.end(); it++)
			if (it->first != NULL)
				lmax = std::max(lmax, it->first->level);
		
		return lmax;
	}
	
	template <typename WaveletType, typename BlockType>
	int Grid<WaveletType, BlockType>::_computeMinLevel(const vector<vector<BlockInfo> >& blockAtLevel) const
	{		
		int l;
		
		for(l=0; l<blockAtLevel.size(); l++)
			if (blockAtLevel[l].size()>0) break;
		
		assert(l<blockAtLevel.size());
		
		return l;
	}
	
	
	template <typename WaveletType, typename BlockType>
	int Grid<WaveletType, BlockType>::_computeMaxLevelJump() const
	{
		int jump = 0;
		for(NeighborhoodType::const_iterator it = m_neighborhood.begin(); it != m_neighborhood.end(); it++)
		{
			const vector<GridNode*>& neighbors = it->second;
			
			for (vector<GridNode*>::const_iterator n=neighbors.begin(); n!=neighbors.end(); n++)
				jump = std::max( jump, abs(it->first->level - (*n)->level) );
		}
		
		return jump;
	}
	
	template <typename WaveletType, typename BlockType>
	void Grid<WaveletType, BlockType>::_computeBlockAtLevel(const HierarchyType& hierarchy, vector<vector<BlockInfo> >& blockAtLevel) const
	{
		vector<BlockInfo> blocks = getBlocksInfo();
		
		blockAtLevel.clear();
		
		const int levels = 1+_computeMaxLevel(hierarchy);
		
		blockAtLevel.resize(levels);
		
		for(vector<BlockInfo>::const_iterator it=blocks.begin(); it!=blocks.end(); it++)
			blockAtLevel[it->level].push_back(*it);
	}
	
	template <typename WaveletType, typename BlockType>
	void Grid<WaveletType, BlockType>::_computeMaxStencilUsed(int stencil_start[3], int stencil_end[3]) const 
	{
		const int support_needed_fwt[2] = {
			1+min(-WaveletType::GaSupport[1], -WaveletType::HaSupport[1]),
			max(-WaveletType::GaSupport[0], -WaveletType::HaSupport[0])
		};
		
		const int support_needed_split[2] = {
			(int)ceil(WaveletType::HsSupport[0]*0.5),
			(2 + (int)floor((WaveletType::HsSupport[1]-1)*0.5 ))
		};
		
		const int support_needed_collapse[2] = {
			1 - WaveletType::HaSupport[1],
			- WaveletType::HaSupport[0]
		};
		
		const int wavelet_left_supp = std::min(-1, std::min(support_needed_fwt[0], std::min(support_needed_split[0], support_needed_collapse[0])));
		const int wavelet_right_supp = std::max(2, std::max(support_needed_fwt[1], std::max(support_needed_split[1], support_needed_collapse[1])));
		
		stencil_start[0] = m_vProcessingDirections[0]? std::min(wavelet_left_supp, m_maxStencil[0][0]) : 0;
		stencil_start[1] = m_vProcessingDirections[1]? std::min(wavelet_left_supp, m_maxStencil[0][1]) : 0;
		stencil_start[2] = m_vProcessingDirections[2]? std::min(wavelet_left_supp, m_maxStencil[0][2]) : 0;
		
		stencil_end[0] = m_vProcessingDirections[0]? std::max(wavelet_right_supp, m_maxStencil[1][0]) : 1;
		stencil_end[1] = m_vProcessingDirections[1]? std::max(wavelet_right_supp, m_maxStencil[1][1]) : 1;
		stencil_end[2] = m_vProcessingDirections[2]? std::max(wavelet_right_supp, m_maxStencil[1][2]) : 1;
	}
	
	template <typename WaveletType, typename BlockType>
	void Grid<WaveletType, BlockType>::_checkResolutionJumpCondition(int maxJump, const int * client_stencil_start, const int * client_stencil_end) const
	{
		int stencil_start[3] = {0,0,0};
		int stencil_end[3] = {0,0,0};
		
		_computeMaxStencilUsed(stencil_start, stencil_end);
		
		if (client_stencil_start != NULL)
		{
			assert(client_stencil_end != NULL);
			
			stencil_start[0] = std::min(stencil_start[0], client_stencil_start[0]);
			stencil_start[1] = std::min(stencil_start[1], client_stencil_start[1]);
			stencil_start[2] = std::min(stencil_start[2], client_stencil_start[2]);
			
			stencil_end[0] = std::max(stencil_end[0], client_stencil_end[0]);
			stencil_end[1] = std::max(stencil_end[1], client_stencil_end[1]);
			stencil_end[2] = std::max(stencil_end[2], client_stencil_end[2]);
		}
		
		const int w[3] = {
			std::max((int)abs(stencil_start[0]), -1 + stencil_end[0]),
			std::max((int)abs(stencil_start[1]), -1 + stencil_end[1]),
			std::max((int)abs(stencil_start[2]), -1 + stencil_end[2])
		};
		
		const int minimalBlockSize[3] = {
			std::max(1, w[0]<<maxJump),
			std::max(1, w[1]<<maxJump),
			std::max(1, w[2]<<maxJump)
		};
		
		const bool bAbort = 
		BlockType::sizeX<minimalBlockSize[0] || 
		BlockType::sizeY<minimalBlockSize[1] || 
		BlockType::sizeZ<minimalBlockSize[2] ;
		
		if (bAbort)
		{
			printf("Grid::_checkResolutionJumpCondition: ABORT! Blocksize + resolution jump condition not satisfied for this stencil width!\n");
			printf("stencil width = %d %d %d, minimalBlockSize = %d %d %d, jump-res = %d, used block size = %d %d %d\n",
				   w[0]+1, w[1]+1, w[2]+1, minimalBlockSize[0], minimalBlockSize[1], minimalBlockSize[2], maxJump, 
				   BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ);
			
			abort();
		}
	}
	
	template <typename WaveletType, typename BlockType>
	void Grid<WaveletType, BlockType>::_computeNeighborhood(const HierarchyType& hierarchy, NeighborhoodType& neighborhood, map<GridNode *, map<int, GridNode *> > & ghostNodes) const
	{
		//0. create a suitable data-structure for the bucket list
		//1. generate bucket list of nodes depending on their locations
		//2. use it to generate the neighborhood
		//3. clear it cause i am a good boy
		
		//0.
		{
			for(map<GridNode *, map<int, GridNode *> >::iterator it1 = ghostNodes.begin(); it1!= ghostNodes.end(); it1++)
			{
				map<int, GridNode *> & m = (it1->second);
				for(map<int, GridNode *>::iterator it2 = m.begin(); it2!= m.end(); it2++)
					delete it2->second;
			}
			
			ghostNodes.clear();
			
			
			for(NeighborhoodType::iterator it1 = neighborhood.begin(); it1!= neighborhood.end(); it1++)
				it1->second.clear();
			
			neighborhood.clear();
		}
		
		
		typedef MRAG::MRAGridHelpers::FlattenedBlock FlattenedBlock;
		typedef vector< const FlattenedBlock* > B;
		
		//1.	
		const int max_level = _computeMaxLevel(hierarchy);
		const int reference_level = max_level;
		const int reference_size = 1<<reference_level;
		
		const int nPeriodicBoundaries = (int)(m_vPeriodicDirection[0]) + (int)(m_vPeriodicDirection[1]) + (int)(m_vPeriodicDirection[2]);
		
		MRAG::Matrix3D< B > & bucket = *(new MRAG::Matrix3D< B >(
																 m_vProcessingDirections[0]?reference_size:1, 
																 m_vProcessingDirections[1]?reference_size:1, 
																 m_vProcessingDirections[2]?reference_size:1));
		
		B pool;
		
		{
			int nAdded = 0;
			for(HierarchyType::const_iterator it=hierarchy.begin(); it!=hierarchy.end(); it++)
			{
				if (it->first == NULL || it->first->isEmpty) continue;
				
				const GridNode & node = *(it->first);
				
				const int level_difference = max_level - node.level;
				const int level_offset = max_level - reference_level;
				
				const int finest_start[3] = {
					m_vProcessingDirections[0] ? ((node.index[0]<<level_difference)) : 0, 
					m_vProcessingDirections[1] ? ((node.index[1]<<level_difference)) : 0, 
					m_vProcessingDirections[2] ? ((node.index[2]<<level_difference)) : 0
				};
				
				const int finest_end[3] = {
					m_vProcessingDirections[0] ? ((node.index[0]+1<<level_difference)) : 1, 
					m_vProcessingDirections[1] ? ((node.index[1]+1<<level_difference)) : 1, 
					m_vProcessingDirections[2] ? ((node.index[2]+1<<level_difference)) : 1
				};
				
				const FlattenedBlock * block = new FlattenedBlock(finest_start, finest_end, &node);
				pool.push_back(block);
				
				const int reference_start[3] = {
					finest_start[0]>>level_offset, 
					finest_start[1]>>level_offset, 
					finest_start[2]>>level_offset
				};
				
				const int reference_end[3] = {
					reference_start[0] + std::max((int)1, (finest_end[0]>>level_offset) - reference_start[0]), 
					reference_start[1] + std::max((int)1, (finest_end[1]>>level_offset) - reference_start[1]),
					reference_start[2] + std::max((int)1, (finest_end[2]>>level_offset) - reference_start[2])
				};
				
				nAdded  += (reference_end[2]-reference_start[2])*(reference_end[1]-reference_start[1])*(reference_end[0]-reference_start[0]);
				
				
				for(int iz = reference_start[2]; iz<reference_end[2]; iz++)
					for(int iy = reference_start[1]; iy<reference_end[1]; iy++)
						for(int ix = reference_start[0]; ix<reference_end[0]; ix++)
							bucket.Access(ix,iy,iz).push_back(block);
				
				int ghosts = 0;
				if (nPeriodicBoundaries>0)
				{
					const bool vProcessAndPeriodic[3] = { 
						m_vProcessingDirections[0] && m_vPeriodicDirection[0], 
						m_vProcessingDirections[1] && m_vPeriodicDirection[1], 
						m_vProcessingDirections[2] && m_vPeriodicDirection[2]
					};
					
					const int n = (1 << node.level);
					
					if (node.index[0] !=0 && node.index[1] !=0 && node.index[2] !=0 &&
						node.index[0]!=n-1 && node.index[1]!=n-1 && node.index[2]!=n-1) continue;
					
					for(int code=0; code<27; code++)
					{
						if (code == 1 + 3 + 9) continue;
						
						const int d[3] = {-(code%3-1), -((code/3)%3-1), -((code/9)%3-1)};
						int idx[3] = {node.index[0] + d[0], node.index[1] + d[1], node.index[2] + d[2]};
						
						if (!(idx[0]<0 || idx[0]>= n ||
							  idx[1]<0 || idx[1]>= n ||
							  idx[2]<0 || idx[2]>= n )) continue;
						
						if (vProcessAndPeriodic[0]) idx[0] = (idx[0]+n)%n;
						if (vProcessAndPeriodic[1]) idx[1] = (idx[1]+n)%n;
						if (vProcessAndPeriodic[2]) idx[2] = (idx[2]+n)%n;
						
						if (idx[0]<0 || idx[0]>= n ||
							idx[1]<0 || idx[1]>= n ||
							idx[2]<0 || idx[2]>= n ) continue;
						ghosts++;
						GridNode * ghostNode = new GridNode(node.isEmpty, node.parent, node.blockID, node.index[0] - d[0]*n, node.index[1] - d[1]*n, node.index[2] - d[2]*n, node.level);
						ghostNodes[(GridNode*)&node][code] = ghostNode;
					}
					
				}
			}
			//printf("Total Added: %d\n", nAdded);
			
		}
		
		if (nPeriodicBoundaries>0)
		{
			const int nX = bucket.getSize()[0];
			const int nY = bucket.getSize()[1];
			const int nZ = bucket.getSize()[2];
			
			for(int iz = 0; iz<nZ; iz++)
				for(int iy = 0; iy<nY; iy++)
					for(int ix = 0; ix<nX; ix++)
					{
						//current block-list
						if (iz !=0 && iy !=0 && ix !=0 &&
							iz!=nZ-1 && iy!=nY-1 && ix!=nX-1)
							ix = nX-1;
						
						
						assert(iz==nZ-1 || iy==nY-1 || ix==nX-1 || ix == 0 || iy == 0 || iz == 0);
						
						for(int i=0; i<27; i++)
						{
							const int dx = i%3-1;
							const int dy = (i/3)%3-1;
							const int dz = (i/9)%3-1;
							
							//if (dz*4 + dy*2 + dx<=0) continue;
							if (i == 1 + 3 + 9) continue;
							
							//const int original_idx[3] = {ix + dx, iy + dy, iz + dz};
							int idx[3] = {ix + dx, iy + dy, iz + dz};
							
							if (m_vProcessingDirections[0] && m_vPeriodicDirection[0]) idx[0] = (idx[0]+nX)%nX;
							if (m_vProcessingDirections[1] && m_vPeriodicDirection[1]) idx[1] = (idx[1]+nY)%nY;
							if (m_vProcessingDirections[2] && m_vPeriodicDirection[2]) idx[2] = (idx[2]+nZ)%nZ;
							
							if (idx[0]<0 || idx[0]>= nX ||
								idx[1]<0 || idx[1]>= nY ||
								idx[2]<0 || idx[2]>= nZ ) continue;
							
							vector< const FlattenedBlock* >& v = bucket.Access(idx[0], idx[1], idx[2]);
							
							const int nElements = v.size();
							for(int iElement=0; iElement<nElements; iElement++)
							{
								if(ghostNodes[(GridNode*)v[iElement]->target][i] == NULL) continue;
								
								
								GridNode& g = *ghostNodes[(GridNode*)v[iElement]->target][i];
								
								const int level_difference = max_level - g.level;
								//const int level_offset = max_level - reference_level;
								
								
								/*	int sign[3] = {
								 original_idx[0]-idx[0]<0?-1:1,
								 original_idx[1]-idx[1]<0?-1:1,
								 original_idx[2]-idx[2]<0?-1:1 };
								 
								 if (original_idx[0]==idx[0])sign[0]=0;
								 if (original_idx[1]==idx[1])sign[1]=0;
								 if (original_idx[2]==idx[2])sign[2]=0;
								 
								 const int displacement[3] = {
								 sign[0]*(1<<max_level),
								 sign[1]*(1<<max_level),
								 sign[2]*(1<<max_level)
								 };
								 
								 const int s[3] = {
								 v[iElement]->index[0]+ displacement[0], 
								 v[iElement]->index[1]+ displacement[1], 
								 v[iElement]->index[2]+ displacement[2]
								 };
								 
								 const int e[3] = {
								 s[0] + v[iElement]->size[0], 
								 s[1] + v[iElement]->size[1], 
								 s[2] + v[iElement]->size[2]
								 };*/
								
								//const int level_offset = max_level - reference_level;
								
								const int finest_start[3] = {
									m_vProcessingDirections[0] ? ((g.index[0]<<level_difference)) : 0, 
									m_vProcessingDirections[1] ? ((g.index[1]<<level_difference)) : 0, 
									m_vProcessingDirections[2] ? ((g.index[2]<<level_difference)) : 0
								};
								
								const int finest_end[3] = {
									m_vProcessingDirections[0] ? ((g.index[0]+1<<level_difference)) : 1, 
									m_vProcessingDirections[1] ? ((g.index[1]+1<<level_difference)) : 1, 
									m_vProcessingDirections[2] ? ((g.index[2]+1<<level_difference)) : 1
								};
								
								
								//printf("TOOO FASST\n");
								//FlattenedBlock * block = new FlattenedBlock(s, e, v[iElement]->target);
								FlattenedBlock * block = new FlattenedBlock(finest_start, finest_end, ghostNodes[(GridNode*)v[iElement]->target][i]);
								pool.push_back(block);
								bucket.Access(ix,iy,iz).push_back(block);
							}
						}
						
					}
		}
		
		//2.
		{
			int nComparisons = 0;
			int nNeighborsConnections = 0;
			
			
			
			const int nX = bucket.getSize()[0];
			const int nY = bucket.getSize()[1];
			const int nZ = bucket.getSize()[2];
			
			for(int iz = 0; iz<nZ; iz++)
				for(int iy = 0; iy<nY; iy++)
					for(int ix = 0; ix<nX; ix++)
					{
						//current block-list
						const B & currentlist = bucket.Access(ix,iy,iz);
						
						for(int i=0; i<27; i++)
						{
							const int dx = i%3-1;
							const int dy = (i/3)%3-1;
							const int dz = (i/9)%3-1;
							
							if (dz*4 + dy*2 + dx<0) continue;
							
							const int idx[3] = {ix + dx, iy + dy, iz + dz};
							
							
							if (idx[0]<0 || idx[0]>= nX ||
								idx[1]<0 || idx[1]>= nY ||
								idx[2]<0 || idx[2]>= nZ ) continue;
							
							
							const B & neighbors = bucket.Access(idx[0],idx[1],idx[2]);
							if (neighbors.size() == 0) continue;
							
							const int n1 = currentlist.size();
							const int n2 = neighbors.size();
							
							if (&currentlist != &neighbors)
								for (int b2=0; b2<n2; b2++)
									for (int b1=0; b1<n1; b1++)
										FlattenedBlock::compareAndAdd(*currentlist[b1], *neighbors[b2], neighborhood,nComparisons, nNeighborsConnections);
							else
								for (int b2=0; b2<n2; b2++)
									for (int b1=0; b1<=b2; b1++)
										FlattenedBlock::compareAndAdd(*currentlist[b1], *currentlist[b2], neighborhood,nComparisons, nNeighborsConnections);
						}
					}
			
			//printf("Done N=%dx%dx%d Compared:%d Found:%d\n", nX, nY, nZ, nComparisons, nNeighborsConnections);
		}
		
		//3.
		{
			for(map< GridNode *, map<int, GridNode *> >::const_iterator itGridNode = ghostNodes.begin(); itGridNode!=ghostNodes.end(); itGridNode++)
			{
				const map<int, GridNode *> & currentGhosts = itGridNode->second;
				for(map<int, GridNode *>::const_iterator itGhostNode = currentGhosts.begin(); itGhostNode!=currentGhosts.end(); itGhostNode++)
					neighborhood.erase(itGhostNode->second);
			}
			
			for(typename B::iterator e = pool.begin(); e!= pool.end(); e++)
			{
				delete *e;
				*e = NULL;
			}
			
			delete &bucket;
		}
	}
	
}

