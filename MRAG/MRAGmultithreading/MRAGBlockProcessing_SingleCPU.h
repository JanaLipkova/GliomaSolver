/*
 *  MRAGBlockProcessing_SingleCPU.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 9/13/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */


#pragma once
#include "MRAGcore/MRAGCommon.h"

namespace MRAG
{
	namespace Multithreading
	{
        /**
         * Sample Functor acting on simple blocks without ghosts (in this sample we do nothing).
         * This is just to show what functions should be provided by Processing functors
         * passed to process(vector< BlockInfo > &, Collection &, Processing&, ...) 
         * as in MRAG::Multithreading::BlockProcessing_SingleCPU, 
         *       MRAG::Multithreading::BlockProcessing_TBB and 
         *       MRAG::Multithreading::BlockProcessing_Pipeline_TBB
         */
        class DummySimpleBlockFunctor {
        public:
            // Note for programmers: functor is used in the end as ProcessingMT type in BlockProcessingMT_Simple...
            /**
             * Copy constructor (may not be needed as compiler-generated one may be enough).
             */
            DummySimpleBlockFunctor(const DummySimpleBlockFunctor& c) {
                // copy content
                memcpy(this, &c, sizeof(DummySimpleBlockFunctor) );
                // debug
                printf("DummySimpleBlockFunctor - Copy constructor.\n");
            }
            
            /**
             * Default constructor (may not be needed as compiler-generated one may be enough).
             */
            DummySimpleBlockFunctor() {
                // debug
                printf("DummySimpleBlockFunctor - Default constructor.\n");
            }
            
            /**
             * Perform some operation on a simple block (in this sample we do nothing).
             * @param info  Reference to info of current block.
             * @param b     Reference to actual block to be used as input and output.
             */
            template<typename BlockType>
            inline void operator()(const BlockInfo& info, BlockType& b) const {
                // define some shorthands for types
                typedef BlockType B;
                typedef typename BlockType::ElementType E;
                
                // looping in this order for efficiency reasons (see MRAG::Block)
                for(int iz=0; iz<B::sizeZ; iz++) {
                    for(int iy=0; iy<B::sizeY; iy++) {
                        for(int ix=0; ix<B::sizeX; ix++) {
                            // get element in block (indices for b may not be < 0 or >= size)
                            E element = b(ix, iy, iz);
                            // do something with it (here we do...nothing)
                            b(ix, iy, iz) = b(ix, iy, iz);
                        }
                    }
                }
                
                // alternative way of looping
                const int n = B::sizeZ*B::sizeY*B::sizeX;
                E* ptrE = &(b[0]);
                for(int iE=0; iE<n; iE++, ptrE++) {
                    // ptrE points to current element
                }
                
                // debug
                printf("DummySimpleBlockFunctor - Operation on block %d @ %d, %d, %d\n", 
                       info.blockID, info.index[0], info.index[1], info.index[2]);
            }
        };
        
        /**
         * Sample Functor acting on extended blocks with ghosts (in this sample we do nothing).
         * This is just to show what functions should be provided by Processing functors
         * passed to process(vector< BlockInfo > &, Collection &, BoundaryInfo &, Processing&, ...) 
         * and to pipeline_process(vector< BlockInfo > &, Collection &, BoundaryInfo &, Processing &)
         * as in MRAG::Multithreading::BlockProcessing_SingleCPU, 
         *       MRAG::Multithreading::BlockProcessing_TBB and 
         *       MRAG::Multithreading::BlockProcessing_Pipeline_TBB
         */
        class DummyBlockFunctor {
        public:
            // Note for programmers: functor is used in the end as ProcessingMT type in BlockProcessingMT...
            /**
             * Maximal stencil used for computations at lower boundary.
             * Defines how many ghosts we will need for computations.
             */
            int stencil_start[3];
            /**
             * Maximal stencil used for computations at upper boundary.
             * Defines how many ghosts we will need for computations.
             */
            int stencil_end[3];
            
            /**
             * Copy constructor (may not be needed as compiler-generated one may be enough).
             */
            DummyBlockFunctor(const DummyBlockFunctor& c) {
                // copy content
                memcpy(this, &c, sizeof(DummyBlockFunctor) );
                // debug
                printf("DummyBlockFunctor - Copy constructor.\n");
            }
            
            /**
             * Default constructor (needed to set at least stencils).
             */
            DummyBlockFunctor() {
                // set stencil
                stencil_start[0] = stencil_start[1] = -1; stencil_start[2] = 0;
                stencil_end[0] = stencil_end[1] = 2; stencil_end[2] = 1;
                // debug
                printf("DummyBlockFunctor - Default constructor.\n");
            }
            
            /**
             * Perform some operation on an extended block (in this sample we do nothing).
             * @param i     Reference to extended block (including ghosts as in MRAG::BlockLab) to be used as input.
             * @param info  Reference to info of current block.
             * @param o     Reference to actual block to be used as output.
             */
            template<typename LabType, typename BlockType>
            inline void operator()(LabType& i, const BlockInfo& info, BlockType& o) const {
                // define some shorthands for types
                typedef BlockType B;
                typedef typename BlockType::ElementType E;
                
                // looping in this order for efficiency reasons (see MRAG::Block)
                for(int iz=0; iz<B::sizeZ; iz++) {
                    for(int iy=0; iy<B::sizeY; iy++) {
                        for(int ix=0; ix<B::sizeX; ix++) {
                            // get element in block (indices may be < 0 or >= size)
                            E element = i(ix, iy, iz);
                            // do something with it (here we do...nothing)
                            // (indices for o may not be < 0 or >= size)
                            o(ix, iy, iz) = i(ix, iy, iz);
                        }
                    }
                }
                // debug
                printf("DummyBlockFunctor - Operation on block %d @ %d, %d, %d\n", 
                       info.blockID, info.index[0], info.index[1], info.index[2]);
            }
        };
        
        /**
         * Functor to actually perform the operations on the blocks.
         * See MRAG::Multithreading::DummySimpleBlockFunctor for a sample ProcessingMT type.
         */
		template <typename BlockType, typename ProcessingMT>
		class BlockProcessingMT_Simple
		{
			ProcessingMT& processing;
			vector<BlockInfo>& vInfo;
			BlockType ** ptrBlocks;
			
		public:
			BlockProcessingMT_Simple(vector<BlockInfo>& vInfo_, ProcessingMT& processing_, BlockType ** ptrs):
			processing(processing_), vInfo(vInfo_), ptrBlocks(ptrs) {}
			
			template <typename BlockedRange>
			void operator()(const BlockedRange& r) const
			{
				typedef BlockType B;
				typedef typename B::ElementType E;
				
				const int nBlocks = r.end() - r.begin();
				BlockInfo* v = new BlockInfo[nBlocks];
				
				for (int i=0; i<nBlocks; i++)
					v[i] = vInfo[i+r.begin()];
				
                // copy constructor required for ProcessingMT
				ProcessingMT p = processing;
				
				for (int iB=0; iB<nBlocks; iB++) {
                    // operator()(const BlockInfo&, BlockType&) required for ProcessingMT
					p(v[iB], *ptrBlocks[iB+r.begin()]);
                }
				
				delete [] v;
			}		
		};
		
        /**
         * Functor to actually perform the operations on the blocks.
         * See MRAG::Multithreading::DummyBlockFunctor for a sample ProcessingMT type.
         */
		template <typename BlockType, template <typename BB> class Lab, typename Collection, typename ProcessingMT>
		class BlockProcessingMT
		{
			vector<BlockInfo>& vInfo;
			Collection& collection;
			BoundaryInfo& boundaryInfo;
			ProcessingMT& processing;
			BlockType ** ptrBlocks;
			
		public:
			BlockProcessingMT(vector<BlockInfo>& vInfo_, Collection& collection_, 
							  BoundaryInfo& boundaryInfo_, ProcessingMT& processing_, BlockType ** ptrs):
			vInfo(vInfo_), collection(collection_), boundaryInfo(boundaryInfo_), processing(processing_), ptrBlocks(ptrs){}
			
			template <typename BlockedRange>
			void operator()(const BlockedRange& r) const
			{
				typedef BlockType B;
				typedef typename B::ElementType E;
				
				const int nBlocks = r.end() - r.begin();
				BlockInfo* v = new BlockInfo[nBlocks];
				for(int i=0; i<nBlocks; i++) v[i] = vInfo[i+r.begin()];
				
				Lab<B> lab;
                // copy constructor required for ProcessingMT
				ProcessingMT p = processing;
                // stencil_start and stencil_end required for ProcessingMT
				lab.prepare(collection, boundaryInfo, processing.stencil_start, processing.stencil_end);
				
				for(int iB=0; iB<nBlocks; iB++)
				{
					BlockInfo& info = v[iB];
					B& b = * ptrBlocks[r.begin() + iB];
					lab.load(info);
					// operator()(LabType&, const BlockInfo&, BlockType&) required for ProcessingMT
					p(lab, info, b);
				}
				
				delete [] v;
			}		
		};
		
        /**
         * Process blocks on a single core.
         * @see BlockProcessing_SingleCPU::process
         */
		template <typename BlockType>
		class BlockProcessing_SingleCPU
		{
			template<typename Collection>
			static BlockType** createBlockPointers(vector<BlockInfo>& vInfo, Collection& collection)
			{
				BlockType** result = new BlockType*[vInfo.size()];
				
				for(int i=0; i<vInfo.size(); i++)
					result[i] = &collection.lock(vInfo[i].blockID);
				
				return result;
			}
			
			template<typename Collection>
			static void destroyBlockPointers(BlockType**& ptrs, vector<BlockInfo>& vInfo , Collection& collection)
			{
				for(int i=0; i<vInfo.size(); i++)
					collection.release(vInfo[i].blockID);
				
				delete [] ptrs;
				ptrs = NULL;
			}
			
		public:
			
            /**
             * Process blocks one after the other.
             * @param vInfo         Info of all the blocks (to be processed) in the grid
             *                      (e.g. result of Grid::getBlocksInfo()).
             * @param c             Collection of all the blocks (to be processed) in the grid
             *                      (e.g. result of Grid::getBlockCollection()).
             * @param p             Functor processing the block.
             *                      See MRAG::Multithreading::DummySimpleBlockFunctor for details.
             * @param dummy         Optional and never used (just to match signature of tbb-versions).
             */
			template <typename Processing, typename Collection>
			static void process(vector<BlockInfo>& vInfo, Collection& c, Processing& p, int dummy = -1)
			{
				BlockType** ptrs =  createBlockPointers(vInfo, c);
				
				BlockProcessingMT_Simple<BlockType,Processing> bps(vInfo, p, ptrs);
				bps(SimpleInterval(0, vInfo.size()));
				
				destroyBlockPointers(ptrs, vInfo, c);
			}
			
            /**
             * Process blocks one after the other.
             * @param vInfo         Info of all the blocks (to be processed) in the grid
             *                      (e.g. result of Grid::getBlocksInfo()).
             * @param c             Collection of all the blocks (to be processed) in the grid
             *                      (e.g. result of Grid::getBlockCollection()).
             * @param b             Info on the boundaries of the grid (e.g. result of Grid::getBoundaryInfo())
             * @param p             Functor processing the block.
             *                      See MRAG::Multithreading::DummyBlockFunctor for details.
             * @param dummy         Optional and never used (just to match signature of tbb-versions).
             */
			template <template <typename Btype> class Lab, typename Processing, typename Collection>
			static void process(vector<BlockInfo>& vInfo, Collection& c, BoundaryInfo& b, Processing& p, int dummy = -1)
			{
				BlockType** ptrs =  createBlockPointers(vInfo, c);
				
				BlockProcessingMT<BlockType, Lab, Collection, Processing> bpg(vInfo, c, b, p, ptrs);
				bpg(SimpleInterval(0, vInfo.size()));
				
				destroyBlockPointers(ptrs, vInfo, c);
			}
			
            /**
             * Alias for process<BlockLab>(vInfo, c, b, p).
             * @see BlockProcessing_SingleCPU::process()
             */
			template <typename Processing, typename Collection>
			void pipeline_process(vector<BlockInfo>& vInfo, Collection& c, BoundaryInfo& b, Processing& p)
			{
				process<BlockLab>(vInfo, c, b, p);
			}
			
			template <typename B> friend class BlockProcessing_TBB; 
		};
	} /* namespace Multithreading */
} /* namespace MRAG */

