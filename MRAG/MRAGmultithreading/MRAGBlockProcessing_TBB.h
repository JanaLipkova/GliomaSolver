/*
 *  MRAGBlockProcessing_TBB.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 9/13/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include <vector>
#include <map>

//#define TBB_DO_THREADING_TOOLS 1
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/pipeline.h"
#include "tbb/concurrent_queue.h"

#include "MRAGcore/MRAGEnvironment.h"
#pragma once
#ifdef _MRAG_TBB
#include "MRAGBlockProcessing_SingleCPU.h"
#undef max

using namespace std;

using namespace tbb;

namespace MRAG
{
	namespace Multithreading
	{
        /**
         * Functor to actually perform the operations on the blocks.
         * See MRAG::Multithreading::DummyBlockFunctor for a sample ProcessingMT type.
         */
		template <typename BlockType, template <typename BB> class Lab, typename Collection, typename ProcessingMT, int nSlots>
		class BlockProcessingMT_TBB
		{
			Collection& collection;
			BoundaryInfo& boundaryInfo;
			ProcessingMT& processing;
		
			const BlockInfo * ptrInfos;
			
			concurrent_bounded_queue<Lab<BlockType> *>& m_availableLabs;
			
		public:
			BlockProcessingMT_TBB(concurrent_bounded_queue<Lab<BlockType> *>& availableLabs, const BlockInfo * ptrInfos_, Collection& collection_, 
							  BoundaryInfo& boundaryInfo_, ProcessingMT& processing_):
				collection(collection_), boundaryInfo(boundaryInfo_), 
				processing(processing_), 
				ptrInfos(ptrInfos_), 
				m_availableLabs(availableLabs)
			{
				for(int i=0; i<nSlots; i++)
				{
					Lab<BlockType> * lab = NULL;
					availableLabs.pop(lab);
					assert(lab!=NULL);
					
                    // stencil_start and stencil_end required for ProcessingMT
					lab->prepare(collection_, boundaryInfo_, processing_.stencil_start, processing_.stencil_end);
					
					availableLabs.push(lab);
				}
			}
		
			template <typename BlockedRange>
			void operator()(const BlockedRange& r) const
			{
				Lab<BlockType>* lab = NULL;
				m_availableLabs.pop(lab);
				assert(lab != NULL);
				
				const int nBlocks = r.end() - r.begin();
				const BlockInfo* v = ptrInfos + r.begin();
				
				for(int iB=0; iB<nBlocks; iB++)
				{
					const BlockInfo& info = v[iB];
					BlockType& block = *(BlockType*)info.ptrBlock;
					
					lab->load(info);
                    
					// operator()(LabType&, const BlockInfo&, BlockType&) required for ProcessingMT
					processing(*lab, info, block);
				}
				
				m_availableLabs.push(lab);	
			}	
			
			BlockProcessingMT_TBB(const BlockProcessingMT_TBB& p):
			collection(p.collection), boundaryInfo(p.boundaryInfo), processing(p.processing), 
			ptrInfos(p.ptrInfos), m_availableLabs(p.m_availableLabs){}
		
		private:
			//forbidden
			BlockProcessingMT_TBB& operator=(const BlockProcessingMT_TBB& p){abort(); return *this;}
		}; /* BlockProcessingMT_TBB */
		
        /**
         * Functor to actually perform the operations on the blocks.
         * See MRAG::Multithreading::DummySimpleBlockFunctor for a sample ProcessingMT type.
         */
		template <typename BlockType, typename ProcessingMT>
		class BlockProcessingMT_Simple_TBB
		{
			ProcessingMT& processing;
			const BlockInfo * ptrInfos;
			
		public:
			
			BlockProcessingMT_Simple_TBB(const BlockInfo * ptrInfos_, ProcessingMT& processing_):
			processing(processing_), ptrInfos(ptrInfos_){}
			
			BlockProcessingMT_Simple_TBB(const BlockProcessingMT_Simple_TBB& p):
			processing(p.processing), ptrInfos(p.ptrInfos){}
			
			//forbidden
			BlockProcessingMT_Simple_TBB& operator=(const BlockProcessingMT_Simple_TBB& p){abort(); return *this;}

			
			template <typename BlockedRange>
			void operator()(const BlockedRange& r) const
			{
				const int nBlocks = r.end() - r.begin();
				const BlockInfo* v = ptrInfos + r.begin();
				
                // copy constructor required for ProcessingMT
				ProcessingMT p = processing;
				for(int iB=0; iB<nBlocks; iB++)
				{
					const BlockInfo& info = v[iB];
					BlockType& block = *(BlockType*)(info.ptrBlock);
                    // operator()(const BlockInfo&, BlockType&) required for ProcessingMT
					p(info, block);
				}
			}	
		}; /* BlockProcessingMT_Simple_TBB */

        /**
         * Process blocks with tbb (using threads).
         * @see BlockProcessing_TBB::process
         */
		template <typename BlockType>
		class BlockProcessing_TBB
		{
			static map<string, vector<void *> > s_cachedResources;
			static BlockInfo * s_ptrInfos;
			static int s_nBlocks;
	
		protected:
			template <typename Lab>
			static void _getResources(concurrent_bounded_queue<Lab*>& resources, const int nSlots)
			{
				//0. checks
				//1. find if the resources are already in the cache
				//2. if not allocate them, store in the cache
				//3. fill the queue
				
				//0.
				//assert(resources.size() == 0);
				assert(resources.empty());
				
				//1.
				string key = typeid(Lab).name();
				map<string, vector<void *> >::const_iterator it = s_cachedResources.find(key);
				const bool bCacheMiss = (it == s_cachedResources.end());
				
				if (bCacheMiss)
				{
					//2.
					typedef cache_aligned_allocator<Lab> lab_allocator;
					
					vector<void*>& new_resources = s_cachedResources[key];
					new_resources.reserve(nSlots);
					
					for(int i=0; i<nSlots; i++)
					{
						Lab* t = lab_allocator().allocate(1);
						Lab* lab = new((void*)(t)) Lab();
						
						//3.
						new_resources.push_back(lab);
						resources.push(lab);
					}
				}
				else //3.
					for(vector<void *>::const_iterator itElem = it->second.begin(); itElem!=it->second.end(); itElem++)
						resources.push((Lab*)*itElem);
			}
			
			template<typename Collection>
			static const BlockInfo * _prepareBlockInfos(vector<BlockInfo>& vInfo, Collection& collection)
			{
				typedef cache_aligned_allocator<BlockInfo> info_allocator;
				
				const int nCurrBlocks = vInfo.size();
				
				if (s_nBlocks < nCurrBlocks)
				{
					if (s_ptrInfos != NULL)
						info_allocator().deallocate(s_ptrInfos, s_nBlocks);
					
					s_ptrInfos = info_allocator().allocate(nCurrBlocks);
					s_nBlocks = nCurrBlocks;
				}
				
				BlockInfo * currInfo = s_ptrInfos;
				for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++, currInfo++)
					*currInfo = *it;

				currInfo = s_ptrInfos;
				for(int i=0; i<nCurrBlocks; i++, currInfo++)
					currInfo->ptrBlock = &collection.lock(currInfo->blockID);
				
				return s_ptrInfos; 
			}
			
			template<typename Collection>
			static void _releaseBlockPointers(vector<BlockInfo>& vInfo , Collection& collection)
			{
				for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++)
					collection.release(it->blockID);
			}
			
		public:
			BlockProcessing_TBB()
			{
			}
			
			virtual ~BlockProcessing_TBB()
			{
			}
			
            /**
             * Process blocks in parallel using parallel_for (see tbb-doc) to split up the work.
             * @param vInfo         Info of all the blocks (to be processed) in the grid
             *                      (e.g. result of Grid::getBlocksInfo()).
             * @param c             Collection of all the blocks (to be processed) in the grid
             *                      (e.g. result of Grid::getBlockCollection()).
             * @param p             Functor processing the block.
             *                      See MRAG::Multithreading::DummySimpleBlockFunctor for details.
             * @param nGranularity  Granularity for the parallel_for (see tbb-doc).
             *                      Optional: if not set, auto_partitioner (see tbb-doc) will be used.
             */
			template <typename Processing, typename Collection>
			static void process(vector<BlockInfo>& vInfo, Collection& c, Processing& p, 
								int nGranularity = -1)
			{								
				const bool bAutomatic = nGranularity<0;
				
				const BlockInfo* infos = _prepareBlockInfos(vInfo, c);
				
				BlockProcessingMT_Simple_TBB<BlockType,Processing> body(infos, p);
				
				if (bAutomatic)
					parallel_for(blocked_range<size_t>(0,vInfo.size()), body,  auto_partitioner());
				else
					parallel_for(blocked_range<size_t>(0,vInfo.size(), nGranularity), body);
				
				_releaseBlockPointers(vInfo, c);
			}
			
            /**
             * Process blocks in parallel using parallel_for (see tbb-doc) to split up the work.
             * @param vInfo         Info of all the blocks (to be processed) in the grid
             *                      (e.g. result of Grid::getBlocksInfo()).
             * @param c             Collection of all the blocks (to be processed) in the grid
             *                      (e.g. result of Grid::getBlockCollection()).
             * @param b             Info on the boundaries of the grid (e.g. result of Grid::getBoundaryInfo())
             * @param p             Functor processing the block.
             *                      See MRAG::Multithreading::DummyBlockFunctor for details.
             * @param nGranularity  Granularity for the parallel_for (see tbb-doc).
             *                      Optional: if not set, auto_partitioner (see tbb-doc) will be used.
             */
			template <template <typename Btype> class Lab, typename Processing, typename Collection>
			static void process(vector<BlockInfo>& vInfo, Collection& c, BoundaryInfo& b, Processing& p, 
								int nGranularity = -1)
			{
				const int nSlots= (int)(_MRAG_TBB_NTHREADS_HINT);
				
				concurrent_bounded_queue<Lab<BlockType> *> resources;
				_getResources(resources, nSlots);
				
				const BlockInfo* infos = _prepareBlockInfos(vInfo, c);
				
				BlockProcessingMT_TBB<BlockType, Lab, Collection, Processing, nSlots> body(resources, infos, c, b, p) ;
				
				const bool bAutomatic = nGranularity<0;
				if (bAutomatic)
					parallel_for(blocked_range<size_t>(0,vInfo.size()), body,  auto_partitioner());
				else
					parallel_for(blocked_range<size_t>(0,vInfo.size(), nGranularity), body);
				
				_releaseBlockPointers(vInfo, c);
			}
		}; /* BlockProcessing_TBB */
		
		template <typename BlockType>
		map<string, vector<void *> > BlockProcessing_TBB<BlockType>::s_cachedResources;

		template <typename BlockType>
		BlockInfo * BlockProcessing_TBB<BlockType>::s_ptrInfos = NULL;
		
		template <typename BlockType>
		int BlockProcessing_TBB<BlockType>::s_nBlocks = 0;

#pragma mark class BlockProcessing_Pipeline_TBB
		/**
         * Enhanced blockprocessing designed for efficient gathering.
         */
		template <typename BlockType, template <typename Btype> class LabType, int nSlots>
		class BlockProcessing_Pipeline_TBB: public BlockProcessing_TBB<BlockType>
		{
			LabType<BlockType> m_LabSlots[nSlots];
			
			template<typename Processing>
			struct Token
			{
				BlockInfo* info;
				Processing* processing;
				LabType<BlockType>* lab;
				BlockType* block;
				int slotID;
				
				Token(): info(NULL), lab(NULL), processing(NULL), block(NULL), slotID(-1) {}
				
				void process() 
				{
					lab->load(*info);
					const Processing& p = *processing;
                    // operator()(LabType&, const BlockInfo&, BlockType&) required for Processing
					p(*lab, *info, *block);
				}
			};
			
			template<typename Processing>
			class Filter_PickABlock: public tbb::filter
			{				
				Token<Processing> tokens[nSlots];
				const BlockInfo * ptrCurrentInfo;
				
				
				concurrent_bounded_queue<int> availableBlocks;
			public:
				concurrent_bounded_queue<int> availableSlots;
				
				Filter_PickABlock(LabType<BlockType>* labSlots,  const BlockInfo* vInfo, const int nBlocks, Processing& p):
					filter(false),	
					ptrCurrentInfo(vInfo),
					availableSlots(), availableBlocks()
				{
					availableSlots.set_capacity(nSlots);
					for(int i=0; i<nSlots; i++)
					{
						tokens[i].processing = &p;
						tokens[i].lab = &labSlots[i];
						availableSlots.push(i);
					}
					
					availableBlocks.set_capacity(nBlocks);
					for(int i=0; i<nBlocks; i++)
						availableBlocks.push(i);
				}
				
				~Filter_PickABlock()
				{
					for(int i=0; i<nSlots; i++)
					{
						tokens[i].processing = NULL;
						tokens[i].lab = NULL;
					}
				}
				
				void* operator()(void *)
				{
					int iBlock = -1;
					//if( !availableBlocks.pop_if_present(iBlock) ) return NULL;
					if( !availableBlocks.try_pop(iBlock) ) return NULL; 
					assert(iBlock>=0);
					
					int slotID = -1;
					availableSlots.pop(slotID);
					assert(slotID>=0 && slotID<nSlots);
					
					Token<Processing>& currenToken = tokens[slotID];
					
					currenToken.info = const_cast<BlockInfo*>(&ptrCurrentInfo[iBlock]);
					currenToken.block = (BlockType*)currenToken.info->ptrBlock;
					currenToken.slotID = slotID;

					return &currenToken;
				}
			private:
				//forbidden
				Filter_PickABlock(const Filter_PickABlock& p):filter(false),ptrCurrentInfo(NULL),availableSlots(), availableBlocks()
				{abort();}
				
				Filter_PickABlock& operator=(const Filter_PickABlock& p){abort(); return *this;}
			};
			
			template<typename Processing>
			class Filter_Process: public tbb::filter
			{
			public:
				Filter_Process(): filter(false){}
				
				void * operator()(void * ptr_token)
				{
					Token<Processing>& token = *(Token<Processing> *)ptr_token;
					token.process();
					return &token;
				}
			};
			
			template<typename Processing>
			class Filter_Release: public tbb::filter
			{
				concurrent_bounded_queue<int>& availableSlots;
			public:
				Filter_Release(Filter_PickABlock<Processing>& pickAblock): 
				filter(false),
				availableSlots(pickAblock.availableSlots)
				{
				}
				
				void * operator()(void *ptr_token) 
				{
					Token<Processing>& token = *(Token<Processing> *)ptr_token;
					
					availableSlots.push(token.slotID);
					
					return NULL; 
				}
			};
			
		public:
			BlockProcessing_Pipeline_TBB(): BlockProcessing_TBB<BlockType>(){}
			~BlockProcessing_Pipeline_TBB(){}
			
            /**
             * Alias for process<LabType>(vInfo, c, b, p, nGranularity) of BlockProcessing_TBB.
             * @see BlockProcessing_TBB::process()
             */
			template <typename Processing, typename Collection>
			static void process(vector<BlockInfo>& vInfo, Collection& c, BoundaryInfo& b, Processing& p, int nGranularity = -1)
			{
				BlockProcessing_TBB<BlockType>::template process<LabType>(vInfo, c, b, p, nGranularity);
			}
			
            /**
             * Alias for process<LabType>(vInfo, c, p, nGranularity) of BlockProcessing_TBB.
             * @see BlockProcessing_TBB::process()
             */
			template <typename Processing, typename Collection>
			static void process(vector<BlockInfo>& vInfo, Collection& c, Processing& p, int nGranularity = -1)
			{
				BlockProcessing_TBB<BlockType>::process(vInfo, c, p, nGranularity);
			}
			
            /**
             * Process blocks in parallel using a pipeline (see tbb-doc) to split up the work.
             * @param vInfo         Info of all the blocks (to be processed) in the grid
             *                      (e.g. result of Grid::getBlocksInfo()).
             * @param c             Collection of all the blocks (to be processed) in the grid
             *                      (e.g. result of Grid::getBlockCollection()).
             * @param b             Info on the boundaries of the grid (e.g. result of Grid::getBoundaryInfo()).
             * @param p             Functor processing the block.
             *                      See MRAG::Multithreading::DummyBlockFunctor for details.
             * @param nGranularity  Granularity for the parallel_for (see tbb-doc).
             *                      Optional: if not set, auto_partitioner (see tbb-doc) will be used.
             */
			template <typename Processing, typename Collection>
			void pipeline_process(vector<BlockInfo>& vInfo, Collection& c, BoundaryInfo& b, Processing& p)
			{
				//1. prepare resources
				//2. allocate filters
				//3. compose pipeline
				//4. run
				//5. dispose
				
				//1.
				if (vInfo.size() == 0) return;
				
				for(int i=0; i<nSlots; i++) {
                    // stencil_start and stencil_end required for Processing
					m_LabSlots[i].prepare(c, b, p.stencil_start, p.stencil_end);
                }
				
				const int nBlocks = vInfo.size();
				const BlockInfo * infos = BlockProcessing_TBB<BlockType>::_prepareBlockInfos(vInfo, c);//new BlockInfo[nBlocks];
				
				//2.
				Filter_PickABlock<Processing> filter1(m_LabSlots, infos, nBlocks, p);
				Filter_Process<Processing> filter2;
				Filter_Release<Processing> filter3(filter1);
				
				//3.
				pipeline mypipeline;
				mypipeline.add_filter(filter1);
				mypipeline.add_filter(filter2);
				mypipeline.add_filter(filter3);
				
				//4.
				mypipeline.run(nSlots);
				
				//5.
				mypipeline.clear();
				this->_releaseBlockPointers(vInfo, c);
			}
		}; /* BlockProcessing_Pipeline_TBB */
	} /* namespace Multithreading */
} /* namespace MRAG */

#endif
