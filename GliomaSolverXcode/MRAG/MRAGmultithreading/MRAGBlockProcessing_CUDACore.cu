/*
 *  MRAGBlockProcessing_CUDA.cu
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/10/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
 

#include <cutil.h>
#include <stdio.h>

#define _CUDA_SIDE
#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGBlock.h"
#include "MRAGcore/MRAGBoundaryBlockPack.h"
#undef _CUDA_SIDE

using namespace MRAG;
#include "MRAGBlockProcessing_CUDACore.h"

template<class B, class E, class Processing>
__global__  
void SimpleProcessing(Processing& processing, UniformPartition partition, BlockInfo * vInfo);


template< class B, class E, class Processing>
__global__
void GatheringProcessingInside(Processing& processing, UniformPartition partition, BlockInfo * vInfo);


template< class B, class E, class Processing>
__global__  void GatheringProcessingBoundary(Processing& processing, 
		UniformPartition* lutReadData, UniformPartition* lutWriteData,int nPatchesPerBlock,  
		BlockInfo * vInfo, BoundaryBlockPack** vBBPacks);
		

template<class BlockType_>
struct BlockProcessing_CUDACore						
{
	typedef BlockType_ B;
	typedef typename BlockType_::ElementType E;

	template< template <typename X> class Processing> 
	static void process_mmt(BlockInfo * vInfo,  const int nBlocks, Processing<B>& processing, Processing<B>* processingGPU, const int nThreads[3])
	{
		const int s[3] = {0,0,0};
		const int e[3] = {B::sizeX, B::sizeY, B::sizeZ};
		
		UniformPartition partition(s, e, nThreads);
		
		const int nOriginalAmount = sizeof(E)*(partition.work_per_thread[0]*partition.work_per_thread[1]*partition.work_per_thread[2])* nThreads[0]*nThreads[1]*nThreads[2];
		
		//limit SHMEM usage
		{
			const int KB = 1024;
			const int nMaxBytesSHMEM = 20*KB;
			
			
			if (nOriginalAmount>nMaxBytesSHMEM)
			{
				const bool not_one[3] = {
					partition.work_per_thread[0] > 1, 
					partition.work_per_thread[1] > 1,
					partition.work_per_thread[2] > 1
				};
				
				const int countNotOne = (int)not_one[0]+ (int)not_one[1]+ (int)not_one[2];

				const double proportionality_constant = pow(nMaxBytesSHMEM/(double)nOriginalAmount, 1./countNotOne);
				
				const int vMaximalAmountOfWork[3] = { //per thread, per pass
					not_one[0]? max(1,(int)floor(partition.work_per_thread[0]*proportionality_constant)) : 1,
					not_one[1]? max(1,(int)floor(partition.work_per_thread[1]*proportionality_constant)) : 1,
					not_one[2]? max(1,(int)floor(partition.work_per_thread[2]*proportionality_constant)) : 1
				};
				
				partition.setup(s, e, nThreads, vMaximalAmountOfWork);
			}
		}

			dim3 gridSize(nBlocks,1+0*partition.passes_CPU(),1);
		dim3 blockSize(nThreads[0]*nThreads[1]*nThreads[2], 1, 1);
		//printf("work per thread now: %d %d %d, passes %d\n", partition.work_per_thread[0], partition.work_per_thread[1], partition.work_per_thread[2], partition.passes_CPU());
		const int nSharedMemorySize = 
			(partition.work_per_thread[0]*partition.work_per_thread[1]*partition.work_per_thread[2])* 
			sizeof(E)*nThreads[0]*nThreads[1]*nThreads[2];
			
		Processing<BlockLab_CUDA<typename BlockType_::ElementType> >* labProcessing = (Processing<BlockLab_CUDA<typename BlockType_::ElementType> >*)processingGPU;
		
		//SimpleProcessing<BlockType_,  typename BlockType_::ElementType ><<< gridSize, blockSize, nSharedMemorySize >>>(*labProcessing, partition, vInfo);
		SimpleProcessing<BlockType_,  typename BlockType_::ElementType ><<< gridSize, blockSize >>>(*processingGPU, partition, vInfo);
	}
	
	template< template <typename X> class Processing> 
	static void process_mmt(BlockInfo * vInfo, BoundaryBlockPack** vBBPacks, const int nBlocks, 
		UniformPartition * ptrReadLutGPU, UniformPartition * ptrWriteLutGPU, const int nSubspaces, 
		Processing<BlockType_>& processing, Processing<BlockType_>*  processingGPU, const int nThreads[3], const int shmemsize)
	{
		//1. compute data-partitions for the inside, launch the kernels
		//2. compute data-partitions for the boundary, launch the kernels

		//1.
		{
			dim3 gridSize(nBlocks,1,1);
			dim3 blockSize(nThreads[0]*nThreads[1]*nThreads[2], 1, 1);
			
			const int inside_start[3] = {
				- processing.stencil_start[0],
				- processing.stencil_start[1],
				- processing.stencil_start[2]
			};
			
			const int inside_end[3] = {
				B::sizeX - processing.stencil_end[0]+ 1,
				B::sizeY - processing.stencil_end[1]+ 1,
				B::sizeZ - processing.stencil_end[2]+ 1
			};
			
			UniformPartition partition(inside_start, inside_end, nThreads);
			//printf("work per thread now: %d %d %d\n", partition.work_per_thread[0], partition.work_per_thread[1], partition.work_per_thread[2]);
			const int nOriginalAmountSHMEM = 
				sizeof(E)*
				(partition.work_per_thread[0]*partition.threads[0]  + processing.stencil_end[0] - 1 - processing.stencil_start[0])*
				(partition.work_per_thread[1]*partition.threads[1]  + processing.stencil_end[1] - 1 - processing.stencil_start[1])*
				(partition.work_per_thread[2]*partition.threads[2]  + processing.stencil_end[2] - 1 - processing.stencil_start[2]);
				
				//printf("AMOUNT SHMEM %f\n KB", nOriginalAmountSHMEM/1024.);

				
			//limit SHMEM usage
			{
				const int KB = 1024;
				const int nMaxBytesSHMEM = 8*KB;//16*KB;
				
				if (nOriginalAmountSHMEM>nMaxBytesSHMEM)
				{
					const bool not_one[3] = {
						partition.work_per_thread[0] > 1, 
						partition.work_per_thread[1] > 1,
						partition.work_per_thread[2] > 1
					};
					
					const int countNotOne = (int)not_one[0]+ (int)not_one[1]+ (int)not_one[2];
					
					const int nConstantAmount = 
						sizeof(E)* 
						(processing.stencil_end[0] - 1 - processing.stencil_start[0])*
						(processing.stencil_end[1] - 1 - processing.stencil_start[1])*
						(processing.stencil_end[2] - 1 - processing.stencil_start[2]);
					
					const int nVariableAmount = 
						sizeof(E)*
						(partition.work_per_thread[0]*partition.threads[0])*
						(partition.work_per_thread[1]*partition.threads[1])*
						(partition.work_per_thread[2]*partition.threads[2]);
						
					const double proportionality_constant = pow((nMaxBytesSHMEM - nConstantAmount)/(double)nVariableAmount, 1./countNotOne);
					
					const int vMaximalAmountOfWork[3] = { //per thread, per pass
						not_one[0]? max(1,(int)floor(partition.work_per_thread[0]*proportionality_constant)) : 1,
						not_one[1]? max(1,(int)floor(partition.work_per_thread[1]*proportionality_constant)) : 1,
						not_one[2]? max(1,(int)floor(partition.work_per_thread[2]*proportionality_constant)) : 1
					};
					
					partition.setup(inside_start, inside_end, nThreads, vMaximalAmountOfWork);
				}
			}
			
			const int nCorrectAmountSHMEM = 
			sizeof(E)*
			(partition.work_per_thread[0]*partition.threads[0]  + processing.stencil_end[0] - 1 - processing.stencil_start[0])*
			(partition.work_per_thread[1]*partition.threads[1]  + processing.stencil_end[1] - 1 - processing.stencil_start[1])*
			(partition.work_per_thread[2]*partition.threads[2]  + processing.stencil_end[2] - 1 - processing.stencil_start[2]);
			
			
			//printf("FIXED!\n");
			//printf("CORRECTED AMOUNT SHMEM %f\n KB", nCorrectAmountSHMEM/1024.);
			printf("work per thread now: %d %d %d, passes %d, SHMEM=%.2f KB\n", partition.work_per_thread[0], partition.work_per_thread[1], partition.work_per_thread[2], partition.passes_CPU(),nCorrectAmountSHMEM/1024.);
			//processing.dt = 0;
		
		
		//	exit(0);
			Processing<BlockLab_CUDA<typename BlockType_::ElementType> >* labProcessing = (Processing<BlockLab_CUDA<typename BlockType_::ElementType> >*)processingGPU;
			GatheringProcessingInside<BlockType_,  typename BlockType_::ElementType  ><<< gridSize, blockSize, nCorrectAmountSHMEM >>>(*labProcessing, partition, vInfo);
			//exit(0);
			
			
			//GatheringProcessingInside<BlockType_,  typename BlockType_::ElementType  ><<< gridSize, blockSize >>>(*processingGPU, partition, vInfo);
			//CUDA_SAFE_CALL(cudaMemcpy( &processing,processingGPU, sizeof(Processing<BlockType_>), cudaMemcpyDeviceToHost));
			//printf("TIMESTEP : %e\n", processing.dt);
		//	exit(0);
		}

		//2.
		{
			const int nPatches = nSubspaces;
			dim3 gridSize(nPatches*nBlocks,1,1);
			dim3 blockSize(nThreads[0]*nThreads[1]*nThreads[2], 1, 1);
			
			Processing<BlockLab_CUDA<typename BlockType_::ElementType> >* labProcessing = (Processing<BlockLab_CUDA<typename BlockType_::ElementType> >*)processingGPU;
			
			GatheringProcessingBoundary<BlockType_ , typename BlockType_::ElementType> <<< gridSize, blockSize, shmemsize >>>
				(*labProcessing, ptrReadLutGPU, ptrWriteLutGPU, nPatches, vInfo, vBBPacks);
			
		}
		int ii;
	//	scanf("%d\n", &ii);
		//cin>>ii;
	}
};

#define _CUDA_SIDE
#include "../DemoCompressibleFlow/CompressibleFlowTypes.h"
#undef _CUDA_SIDE
/*
template<class B, class E, class Processing>
__global__  
void SimpleProcessing(Processing& processing, UniformPartition partition, BlockInfo * vInfo)
{	
	Processing local_processing = processing;
	BlockInfo info = vInfo[blockIdx.x];
	B& block = *((B*)info.ptrBlock);
	
	int thread_index[3];
	partition.findThreadIndex(threadIdx.x, thread_index);
	
	const int pass = blockIdx.y;
	int s[3], e[3];
	if (partition.findWorkingInterval(thread_index, pass, s, e))
		local_processing(info, block, s, e);
}
*/
/*
template<class B, class E, class Processing>
__global__  
void SimpleProcessing(Processing& processing, UniformPartition partition, BlockInfo * vInfo)
{
	BlockInfo info = vInfo[blockIdx.x];
	B& block = *((B*)info.ptrBlock);
	
	Processing local_processing = processing;

	int thread_index[3];
	partition.findThreadIndex(threadIdx.x, thread_index);	
	const int nWorkPerThreadSize = 
		partition.work_per_thread[0]*
		partition.work_per_thread[1]*
		partition.work_per_thread[2];
	
	extern __shared__ float _cache[];
	E * cache = (E *)_cache;
	
	const char cache_start[3] = {0,0,0};
	const char cache_end[3] = {
		partition.work_per_thread[0],
		partition.work_per_thread[1],
		partition.work_per_thread[2]
	};
	BlockLab_CUDA<E> lab(cache_start, cache_end, cache + nWorkPerThreadSize*threadIdx.x);
	
	int s[3], e[3];
	for(int pass=0; pass<partition.passes(); pass++)
	{
		if(partition.findWorkingInterval(thread_index,pass, s, e))
		{
			int i[3];

			for(i[2]=s[2]; i[2]<e[2]; i[2]++)
			for(i[1]=s[1]; i[1]<e[1]; i[1]++)
			for(i[0]=s[0]; i[0]<e[0]; i[0]++)
				lab(i[0]-s[0], i[1]-s[1], i[2]-s[2]) = block(i[0], i[1], i[2]);
			
			 int local_start[3] = {0,0,0};
			 int local_end[3] = {
				partition.work_per_thread[0],
				partition.work_per_thread[1],
				partition.work_per_thread[2]
			};
			
			info.origin[0] += info.h[0]*s[0];
			info.origin[1] += info.h[1]*s[1];
			info.origin[2] += info.h[2]*s[2];
			
			local_processing(info, lab, local_start, local_end);
			
			info.origin[0] -= info.h[0]*s[0];
			info.origin[1] -= info.h[1]*s[1];
			info.origin[2] -= info.h[2]*s[2];
			
			for(i[2]=s[2]; i[2]<e[2]; i[2]++)
			for(i[1]=s[1]; i[1]<e[1]; i[1]++)
			for(i[0]=s[0]; i[0]<e[0]; i[0]++)
				block(i[0], i[1], i[2]) = lab(i[0]-s[0], i[1]-s[1], i[2]-s[2]);
		}
	}
}
*/

template<class B, class E, class Processing>
__global__  
void SimpleProcessing(Processing& processing, UniformPartition partition, BlockInfo * vInfo)
{	
	Processing local_processing = processing;
	BlockInfo info = vInfo[blockIdx.x];
	B& block = *((B*)info.ptrBlock);
	
	int thread_index[3];
	partition.findThreadIndex(threadIdx.x, thread_index);
	
	int s[3], e[3];
	for(int pass=0; pass<partition.passes(); pass++)
		if (partition.findWorkingInterval(thread_index, pass, s, e))
			local_processing(info, block, s, e);
}


template< class B, class E, class Processing>
__global__
void GatheringProcessingInside(Processing& processing, UniformPartition partition, BlockInfo * vInfo)
{
	Processing local_processing = processing;
	BlockInfo info = vInfo[blockIdx.x];
	B& block = *((B*)info.ptrBlock);

	int thread_index[3];
	partition.findThreadIndex(threadIdx.x, thread_index);	
			
	const int start_lab[3] = {
		processing.stencil_start[0],
		processing.stencil_start[1],
		processing.stencil_start[2] 
	};
	
	const int end_lab[3] = {
		partition.work_per_thread[0]*partition.threads[0]  + processing.stencil_end[0] - 1,
		partition.work_per_thread[1]*partition.threads[1]  + processing.stencil_end[1] - 1,
		partition.work_per_thread[2]*partition.threads[2]  + processing.stencil_end[2] - 1, 
	};
	
	extern __shared__ float _cache[];
	E * cache = (E *)_cache;
	BlockLab_CUDA<E> lab(start_lab, end_lab, (E *)cache);

	for(int pass=0; pass<partition.passes(); pass++)
	{
		int ref_s[3], ref_e[3];
		partition.findBlockWorkingArea(pass, ref_s, ref_e);

		int s[3], e[3], i[3];
		const bool bICanDoIt = partition.findWorkingInterval(thread_index,pass, s, e);

		const int is_first[3] = {
			(int)(thread_index[0] == 0),
			(int)(thread_index[1] == 0),
			(int)(thread_index[2] == 0),
		};

		const int is_last[3] = {
			(int)(s[0]<ref_e[0] && e[0]>=ref_e[0]),
			(int)(s[1]<ref_e[1] && e[1]>=ref_e[1]),
			(int)(s[2]<ref_e[2] && e[2]>=ref_e[2])
		};
		
		const int r_s[3] = {
			s[0] + is_first[0]*processing.stencil_start[0],
			s[1] + is_first[1]*processing.stencil_start[1],
			s[2] + is_first[2]*processing.stencil_start[2]
		};
		
		const int r_e[3] = {
			e[0] + is_last[0]*(processing.stencil_end[0]-1),
			e[1] + is_last[1]*(processing.stencil_end[1]-1),
			e[2] + is_last[2]*(processing.stencil_end[2]-1)
		};
			
		__syncthreads();

		if(bICanDoIt)
		{
			for(i[2]=r_s[2]; i[2]<r_e[2]; i[2]++)
			for(i[1]=r_s[1]; i[1]<r_e[1]; i[1]++)
			for(i[0]=r_s[0]; i[0]<r_e[0]; i[0]++)
				lab(i[0]-ref_s[0], i[1]-ref_s[1], i[2]-ref_s[2]) = block(i[0], i[1], i[2]);
		}
			
		const int start_processing[3] = {
			s[0]-ref_s[0], 
			s[1]-ref_s[1], 
			s[2]-ref_s[2]
		};
		
		const int end_processing[3] = {
			e[0]-ref_s[0], 
			e[1]-ref_s[1],
			e[2]-ref_s[2]
		};
			
		__syncthreads();
			
		if(bICanDoIt)
		{
			local_processing(info, lab, start_processing, end_processing);
		}
			
		if(bICanDoIt)
		{
			for(i[2]=s[2]; i[2]<e[2]; i[2]++)
			for(i[1]=s[1]; i[1]<e[1]; i[1]++)
			for(i[0]=s[0]; i[0]<e[0]; i[0]++)
				block(i[0], i[1], i[2]) = lab(i[0]-ref_s[0], i[1]-ref_s[1], i[2]-ref_s[2]) ;
		}
	}

}
/*


template< class BlockType_, class E, class Processing>
__global__
void GatheringProcessingInside(Processing& processing, UniformPartition partition, BlockInfo * vInfo)
{
	BlockInfo info = vInfo[blockIdx.x];
	BlockType_& block = *((BlockType_*)info.ptrBlock);

	int thread_index[3];
	partition.findThreadIndex(threadIdx.x, thread_index);	
	
	int s[3], e[3];
	for(int pass=0; pass<partition.passes(); pass++)
		if(partition.findWorkingInterval(thread_index,pass, s, e))
			processing(info, block, s, e);

}*/

/*
template< class BlockType_, class Processing>
__global__
void GatheringProcessingInside(Processing& processing, UniformPartition partition, BlockInfo * vInfo)
{
	BlockInfo& info = vInfo[blockIdx.x];
	BlockType_& block = *((BlockType_*)info.ptrBlock);

	int thread_index[3];
	partition.findThreadIndex(threadIdx.x, thread_index);	
	
	int s[3], e[3];
	for(int pass=0; pass<partition.passes(); pass++)
		if(partition.findWorkingInterval(thread_index,pass, s, e))
			processing(info, block, s, e);

}
*/

template< class B, class E, class Processing>
__global__  void GatheringProcessingBoundary(Processing& processing, 
		UniformPartition* lutReadData, UniformPartition* lutWriteData,int nPatchesPerBlock,  
		BlockInfo * vInfo, BoundaryBlockPack** vBBPacks)
{
	//1. setup
	//2. find thread idx, data to work on
	//3. construct the blocklab
	//4. work using the blocklab
	//5. write back the result
	
	//1.
	const int blockID = blockIdx.x / nPatchesPerBlock;
	const int patchID = blockIdx.x % nPatchesPerBlock;
	
	BlockInfo& info = vInfo[blockID];
	B& block = *((B*)info.ptrBlock);
	BoundaryBlockPack& bbpack = *vBBPacks[blockID];
	UniformPartition& read_partition = lutReadData[patchID];
	UniformPartition& write_partition = lutWriteData[patchID];
	
	int s[3], e[3];
	
	//2.
	int thread_index[3];
	read_partition.findThreadIndex(threadIdx.x, thread_index); //must be the same as write_partition.findThreadIndex
	
	//3.
	extern __shared__ float _cache[];
	E * cache = (E *)_cache;
	read_partition.getCoveredInterval(s,e);
	BlockLab_CUDA<E> lab(s, e, (E *)cache);
	
	for(int pass=0; pass<read_partition.passes(); pass++)
		if (read_partition.findWorkingInterval(thread_index, pass, s, e))
		{
			int i[3];

			for(i[2]=s[2]; i[2]<e[2]; i[2]++)
			for(i[1]=s[1]; i[1]<e[1]; i[1]++)
			for(i[0]=s[0]; i[0]<e[0]; i[0]++)
			{
				const bool bOutside =
					i[0]<0 || i[1]<0 || i[2]<0 ||
					i[0]>= B::sizeX || i[1]>= B::sizeY || i[2]>= B::sizeZ;  
				
				if (!bOutside)
					lab(i[0], i[1], i[2]) = block(i[0], i[1], i[2]);
				else
					constructGhosts<B,E>(bbpack, i, lab(i[0], i[1],i[2]));
			}
		}
	
	__syncthreads();

	//4.
	for(int pass=0; pass<write_partition.passes(); pass++)
		if (write_partition.findWorkingInterval(thread_index,pass, s, e))
			processing(info, lab, s, e);

	__syncthreads();
	
	//5.
	for(int pass=0; pass<write_partition.passes(); pass++)
		if (write_partition.findWorkingInterval(thread_index,pass, s, e))
		{
			int i[3];
			for(i[2]=s[2]; i[2]<e[2]; i[2]++)
			for(i[1]=s[1]; i[1]<e[1]; i[1]++)
			for(i[0]=s[0]; i[0]<e[0]; i[0]++)
				block(i[0], i[1], i[2]) = lab(i[0], i[1], i[2]);
		}
}
