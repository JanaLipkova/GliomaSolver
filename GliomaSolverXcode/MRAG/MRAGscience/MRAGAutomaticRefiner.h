/*
 *  MRAGAutomaticRefiner.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 12/1/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include "../MRAGcore/MRAGCommon.h"
#include "../MRAGcore/MRAGRefiner.h"
#include "../MRAGcore/MRAGProfiler.h"
#include "../MRAGcore/MRAGBlockFWT.h"

namespace MRAG
{
	
	template<typename Grid, typename BlockFWT >
	class AutomaticRefiner: public Refiner
	{
		static const bool bVerbose = false;
		
		double m_dMaxMemorySizeMB;
		int m_nMaxPasses;
		
		Profiler* m_refProfiler;
		Refiner * m_refRefiner;
		
		Grid& m_grid;
		BlockFWT& m_fwt;
		
		static float _computeAverageBlockSize(Grid& grid)
		{
			const int nBlocks = grid.getBlocksInfo().size();
			
			assert(nBlocks == grid.getBoundaryInfo().boundaryInfoOfBlock.size());
			
			float block_size = grid.getBlockCollection().getBlockSize()/1024./1024.;
			float bbinfo_size = grid.getBoundaryInfo().getMemorySize()/nBlocks;
			
			return block_size + bbinfo_size;
		}
		
	public:
		
		AutomaticRefiner(Grid& grid, BlockFWT& fwt, Refiner * refiner, double dMaxMemorySizeMB = 1000.0, Profiler* profiler = NULL):
		Refiner(refiner->getMaxLevelJump()), 
		m_grid(grid), m_fwt(fwt), m_refRefiner(refiner), m_refProfiler(profiler),
		m_dMaxMemorySizeMB(dMaxMemorySizeMB)
		{
			m_grid.setRefiner(this);
		}
		
		//inherited from refiner interface
		virtual RefinementPlan* createPlan(const HierarchyType& hierarchy, const NeighborhoodType& neighborhood, const bool vProcessingDirections[3], vector<NodeToRefine>& vRefinements)
		{
			RefinementPlan * plan = m_refRefiner->createPlan(hierarchy, neighborhood, vProcessingDirections, vRefinements);
			
			const int additionalBlocks = plan->nNewBlocks - plan->refinements.size();
			const double avgBlockSize = _computeAverageBlockSize(m_grid);
			const double currentGridSize = m_grid.getMemorySize();
			
			const double estimatedSize = currentGridSize + additionalBlocks*avgBlockSize;
			const double margin = 5.0*avgBlockSize;//MB
			const bool bExceeding = (estimatedSize + margin>= m_dMaxMemorySizeMB);
			
			if (bVerbose)
			{
				printf("down to the refinement: blocks to refine %d, new blocks %d\n",  plan->refinements.size(), additionalBlocks);
				printf("down to the refinement: avgBlockSize=%f, currentGridSize=%f, estimatedSize=%f\n", avgBlockSize,currentGridSize, estimatedSize);
			}
			
			if (bExceeding)
			{
				delete plan;
				return NULL;
			}
			
			return plan;
		}
		
		
		bool refine(double& tolerance, double& compression_tolerance, void (*fillGrid)(Grid& g)=NULL)
		{
			printf("START AutomaticRefiner::refine ========================================================\n");
			double maxDetailLeft;
			if (bVerbose)
				printf("FIRST ATTEMPT ========================================================\n");
			RefinementResult refinement_result = Science::AutomaticRefinementForLevelsets(m_grid, m_fwt, tolerance, -1, true, -1, m_refProfiler, &maxDetailLeft, fillGrid);
			
			if (!refinement_result.hasFailed())
			{
				if (bVerbose)
					printf("PASSED, maxDetailLeft=%e, MB = %f\n", maxDetailLeft, m_grid.getMemorySize());
				
				return false;
			}
			else
				printf("FAILED FIRST STEP  maxDetailLeft=%e ========================================================\n", maxDetailLeft);

			const double maxTolerance = maxDetailLeft;
			const double minTolerance = tolerance;
			
			assert(maxTolerance > minTolerance);
			
			double delta = maxTolerance -  minTolerance + 1e-3;
			double tol = minTolerance;
			int loop_counter = 0;
			const int max_loops = 10;
			
			do
			{
				const double candidate_tol = tol + delta*0.75;
				if (bVerbose)
					printf("LOOP ATTEMPT START candidate_tol=%e ========================================================\n", candidate_tol);
				RefinementResult refinement_result = Science::AutomaticRefinementForLevelsets(m_grid, m_fwt, candidate_tol, -1, true, -1, m_refProfiler, &maxDetailLeft, fillGrid);
				if (bVerbose)
					printf("LOOP ATTEMPT END candidate_tol=%e maxDetailLeft=%e ========================================================\n", candidate_tol, maxDetailLeft);
				if(maxDetailLeft != maxTolerance && maxDetailLeft<candidate_tol)
				{
					printf("qui sta andando tutto a puttane... %e %e\n", maxTolerance, maxDetailLeft);
					abort();
				}

				if(refinement_result.hasFailed())
				{
					tol = candidate_tol;
					delta *= 0.5;

					compression_tolerance = tol/8.;//min(5e-2, tol/20);
					if (bVerbose)
						printf("LOOP ATTEMPT  NOT PASSED, NOT COMPRESS with  compression_tolerance=%e ========================================================\n", compression_tolerance);
				//	Science::AutomaticCompressionForLevelsets(m_grid, m_fwt, compression_tolerance, true, -1, m_refProfiler); 
					
					if (bVerbose)
						printf("NEW TOL =%e\n", candidate_tol);
				}
				else
				{
					tolerance = candidate_tol;
					break;
				}
				
				if (loop_counter++ > max_loops)
				{
					printf("Probably something went wrong\n");
					abort();
				}
			}
			while(true);
			
			const float gridSize = m_grid.getMemorySize();
			if (gridSize > m_dMaxMemorySizeMB)
			{
				printf("AHI AHI memory constraint NOT SATISFIED! constraint: %f MB, grid = %f MB\n", 
						m_dMaxMemorySizeMB, m_grid.getMemorySize());
				abort();
			}

			printf("END AutomaticRefiner::refine ========================================================\n");
			return true;
		}
	};
}