/*
 *  MRAGWrapperCUDA.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/10/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

//shitty macro for the stupid nvcc compiler, I hate this

#define _DECLARATION_CudaWrapper_Gathering(BLOCKTYPE, PROCESSINGTYPE) \
\
extern "C" void process_mmtgathering_CUDA_##BLOCKTYPE##PROCESSINGTYPE(BlockInfo * vInfo, BoundaryBlockPack** vBBPacks, const int nBlocks, UniformPartition * vReadLUT,  UniformPartition * vWriteLUT, const int nSubspaces, PROCESSINGTYPE<BLOCKTYPE>& p, PROCESSINGTYPE<BLOCKTYPE>* pGPU, const int nThreads[3], const int shmemsize); \
\
void process_mmtgathering_CUDA(BlockInfo * vInfo, BoundaryBlockPack** vBBPacks, const int nBlocks, UniformPartition * vReadLUT,  UniformPartition * vWriteLUT, const int nSubspaces, PROCESSINGTYPE<BLOCKTYPE>& p, PROCESSINGTYPE<BLOCKTYPE>* pGPU, const int nThreads[3], const int shmemsize)\
{\
	process_mmtgathering_CUDA_##BLOCKTYPE##PROCESSINGTYPE(vInfo, vBBPacks, nBlocks, vReadLUT, vWriteLUT , nSubspaces, p, pGPU, nThreads, shmemsize);\
}

#define _DECLARATION_CudaWrapper(BLOCKTYPE, PROCESSINGTYPE) \
extern "C" void process_mmt_CUDA_##BLOCKTYPE##PROCESSINGTYPE(BlockInfo * vInfo, const int nBlocks, PROCESSINGTYPE<BLOCKTYPE>& p, PROCESSINGTYPE<BLOCKTYPE>* pGPU, const int nThreads[3]); \
\
void process_mmt_CUDA(BlockInfo * vInfo, const int nBlocks, PROCESSINGTYPE<BLOCKTYPE>& p, PROCESSINGTYPE<BLOCKTYPE>* pGPU, const int nThreads[3])\
{\
	process_mmt_CUDA_##BLOCKTYPE##PROCESSINGTYPE(vInfo, nBlocks, p, pGPU, nThreads);\
}

//shitty macro for the implementation
#define _IMPLEMENTATION_CudaWrapper_Gathering(BLOCKTYPE, PROCESSINGTYPE) \
\
extern "C" void process_mmtgathering_CUDA_##BLOCKTYPE##PROCESSINGTYPE(BlockInfo * vInfo, BoundaryBlockPack** vBBPacks, const int nBlocks, UniformPartition * vReadLUT,  UniformPartition * vWriteLUT, const int nSubspaces, PROCESSINGTYPE<BLOCKTYPE> & p, PROCESSINGTYPE<BLOCKTYPE> * pGPU, const int nThreads[3], const int shmemsize)\
{\
	BlockProcessing_CUDACore<BLOCKTYPE>::process_mmt(vInfo, vBBPacks, nBlocks, vReadLUT, vWriteLUT, nSubspaces, p, pGPU, nThreads, shmemsize);\
}

#define _IMPLEMENTATION_CudaWrapper(BLOCKTYPE, PROCESSINGTYPE) \
\
extern "C" void process_mmt_CUDA_##BLOCKTYPE##PROCESSINGTYPE(BlockInfo * vInfo, const int nBlocks, PROCESSINGTYPE	<BLOCKTYPE> & p, PROCESSINGTYPE<BLOCKTYPE> * pGPU, const int nThreads[3])\
{\
	BlockProcessing_CUDACore<BLOCKTYPE>::process_mmt(vInfo, nBlocks, p, pGPU, nThreads);\
}


#ifdef _CUDA_SIDE
#define CudaWrapper _IMPLEMENTATION_CudaWrapper
#else
#define CudaWrapper _DECLARATION_CudaWrapper
#endif

#ifdef _CUDA_SIDE
#define CudaWrapper_Gathering _IMPLEMENTATION_CudaWrapper_Gathering
#else
#define CudaWrapper_Gathering _DECLARATION_CudaWrapper_Gathering
#endif


 
