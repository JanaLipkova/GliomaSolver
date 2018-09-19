/*
 *  MRAG_CompressionRefinementHelpers.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 9/4/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGProfiler.h"
#include "MRAGcore/MRAGBlockFWT.h"
#include "MRAGSimpleLevelsetBlock.h"

#include <vector>
#include <set>
using namespace std;
namespace MRAG
{
	namespace Science
	{
        /**
         * Runs an automatic refinement.
         * Uses Projector defined by MRAG::BlockFWT to compute "energy" for
         * each block. The Projector is evaluated with the integer template-
         * parameter set from iFirstChannel to iLastChannel.
         * The 
         * @param g                     THE MRAG::Grid.
         * @param fwt                   MRAG::BlockFWT instance for the wavelet transforms.
         * @param dAbsoluteTolerance    Tolerance in terms of "energy" for the refinement.
         * @param iMaxLevel             Maximal level of refinement.
         * @param iMaxLoops             Maximal number of refinement-iterations to be done.
         *                              (each block may be refined once per loop)
         * @param profiler              Optional profiler to monitor performance.
         * @param fillGrid(Grid&)       Optional function to fill the refined grid.
         *                              (otherwise wavelets used to interpolate values)
         * 
         * @see MRAG::BlockFWT, #make_projector
         */
        
        const bool bVerbose = false;

		template < int iFirstChannel, int iLastChannel, typename Grid, typename BlockFWT>
		int AutomaticRefinement(Grid& g, BlockFWT& fwt, const double dAbsoluteTolerance, 
								 const int iMaxLevel = -1, const int iMaxLoops=-1, 
								 MRAG::Profiler* profiler=NULL, void (*fillGrid)(Grid& g)=NULL)
		{
			int loopCounter = 0;
			
			if(bVerbose) printf("AutomaticRefinement\n");
			set<int> niceGuys;
			int nSkippedBlocks = 0;
			int nRefinedBlocks = 0;
			
			do
			{
				vector<BlockInfo> vInfo = g.getBlocksInfo(), vBlocksToFWT;
				for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++)
				{
					if (iMaxLevel>=0 && iMaxLevel <= it->level) continue;
					
					if (niceGuys.find(it->blockID) == niceGuys.end()) // maybe not really a nice guy
						vBlocksToFWT.push_back(*it);
					else
						nSkippedBlocks++;
				}
				if(bVerbose) printf("Blocks to FWT (levelset): %d\n", vBlocksToFWT.size());
				if (vBlocksToFWT.size() == 0) break;
				
				if (profiler !=NULL) profiler->getAgent("AutoRef::FWT").start();
				vector<FWTReport<iLastChannel - iFirstChannel+1> > vReports = 
					BlockFWT::template multichannel_fwt<iFirstChannel, iLastChannel>(vBlocksToFWT, g.getBlockCollection(), g.getBoundaryInfo());
				if (profiler !=NULL) profiler->getAgent("AutoRef::FWT").stop(vBlocksToFWT.size());
				
				set<int> shouldBeRefined;
				for(int i=0; i<vBlocksToFWT.size(); i++)
				{
					if(bVerbose) printf("getOverAll_DetailMaxMag()=%e >= dAbsoluteTolerance=%e -----> %s\n", vReports[i].getOverAll_DetailMaxMag(), dAbsoluteTolerance,vReports[i].getOverAll_DetailMaxMag() >= dAbsoluteTolerance ? "to refine!": "skept");
					
					if (vReports[i].getOverAll_DetailMaxMag() >= dAbsoluteTolerance)
						shouldBeRefined.insert(vBlocksToFWT[i].blockID); // i was right
					else
						niceGuys.insert(vBlocksToFWT[i].blockID); //i was wrong
				}
				if (shouldBeRefined.size() == 0) break;
				
				if (profiler !=NULL) profiler->getAgent("AutoRef::refine").start();
				
				RefinementResult result = g.refine(shouldBeRefined);
				nRefinedBlocks += result.nCollapsedParentBlocks;

				if (profiler !=NULL) profiler->getAgent("AutoRef::refine").stop(shouldBeRefined.size());
				
				if (profiler !=NULL) profiler->getAgent("boundaries").start();
				g.getBoundaryInfo();
				if (profiler !=NULL) profiler->getAgent("boundaries").stop();
				
				if (fillGrid != NULL) fillGrid(g);
				
				if(bVerbose) printf("AutomaticRefinement:: Refined %d, FWT skipped blocks:%d\n", shouldBeRefined.size(), nSkippedBlocks);
				
				loopCounter++;
				
				if (iMaxLoops>0 && loopCounter>=iMaxLoops) break;
				
			} while(true);

			return nRefinedBlocks;
		}
		
        /**
         * Runs an automatic compression.
         * Works just like MRAG::Science::AutomaticRefinement() but compressing (coarsening)
         * the grid instead of refining it. Also the arguments are essentially the same besides:
         * - dAbsoluteTolerance is now the tolerance in terms of "energy" for the compression.
         * - iMaxLevel is missing since the coarsest level is fix anyway.
         */
		template <int iFirstChannel, int iLastChannel, typename BlockFWT, typename Grid>
		int AutomaticCompression(Grid& g, BlockFWT& fwt, const double dAbsoluteTolerance, 
								 const int iMaxLoops=-1,
								 MRAG::Profiler* profiler=NULL, void (*fillGrid)(Grid& g)=NULL)
		{
			if(bVerbose) printf("AutomaticCompression\n");
			
			set<int> niceGuys;
			int nSkippedBlocks = 0, loopCounter = 0, nTotalCollapsed = 0; 
			
			do
			{ 
				vector<BlockInfo> vInfo = g.getBlocksInfo(), vBlocksToFWT;
				for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++)
					if (niceGuys.find(it->blockID) == niceGuys.end()) // maybe not really a nice guy
						vBlocksToFWT.push_back(*it);
					else
						nSkippedBlocks++;
				
				if (profiler !=NULL) profiler->getAgent("AutoCompr::FWT").start();
				vector<FWTReport<iLastChannel - iFirstChannel+1> > vReports = 
				BlockFWT::template multichannel_fwt<iFirstChannel, iLastChannel>(vBlocksToFWT, g.getBlockCollection(), g.getBoundaryInfo());
				if (profiler !=NULL) profiler->getAgent("AutoCompr::FWT").stop(vBlocksToFWT.size());
				
				set<int> shouldBeCompressed;
				
				for(int i=0; i<vBlocksToFWT.size(); i++)
					if (vReports[i].getOverAll_DetailMaxMag() < dAbsoluteTolerance) 
						shouldBeCompressed.insert(vBlocksToFWT[i].blockID);
					else
						niceGuys.insert(vBlocksToFWT[i].blockID);
							
	
				if (shouldBeCompressed.size() == 0) break;
				
				if (profiler !=NULL) profiler->getAgent("AutoCompr::compress").start();
				int nCollapsed;
				g.compress(shouldBeCompressed, nCollapsed);
				
				if (profiler !=NULL) profiler->getAgent("AutoCompr::compress").stop(nCollapsed);
				
				if (profiler !=NULL) profiler->getAgent("boundaries").start();
				g.getBoundaryInfo();
				if (profiler !=NULL) profiler->getAgent("boundaries").stop();
				
				if (fillGrid != NULL) fillGrid(g);
				
				if (nCollapsed == 0) break;
				
				if(bVerbose) printf("AutomaticCompression:: Compressed %d, FWT skipped blocks:%d\n", nCollapsed, nSkippedBlocks);
				
				nTotalCollapsed += nCollapsed;
				
				loopCounter++;
				
				if (iMaxLoops>0 && loopCounter>=iMaxLoops) break;
				
			} while(true);
			
			return nTotalCollapsed;
		}

		
      template <typename BlockFWT, typename Grid>
		int AutomaticCompressionForLevelsets(Grid& g, BlockFWT& fwt, const double dAbsoluteTolerance, const bool bKeepHUpdated = true,  int iMaxLoops=-1, MRAG::Profiler* profiler=NULL)
		{
			if(bVerbose) printf("AutomaticCompressionForLevelsets\n");
			
			set<int> niceGuys;
			int nSkippedBlocks = 0, loopCounter = 0, nTotalCollapsed = 0; 
			
			do
			{
				fwt.prepare(g.getBlockCollection(), g.getBoundaryInfo());

				vector<BlockInfo> vInfo = g.getBlocksInfo(), vBlocksToFWT;
				for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++)
					if (niceGuys.find(it->blockID) == niceGuys.end()) // maybe not really a nice guy
						vBlocksToFWT.push_back(*it);
					else
						nSkippedBlocks++;
				
				if (profiler !=NULL) profiler->getAgent("AutoCompr::FWT").start();

				set<int> shouldBeCompressed;
				for(vector<BlockInfo>::iterator it = vBlocksToFWT.begin(); it != vBlocksToFWT.end(); it++)
				{
					g.getBlockCollection().lock(it->blockID).setH(it->h[0]);
					
					fwt.template multichannel_fwt<0, 0>(*it);
					
					if (fwt.getReport().getOverAll_DetailMaxMag() < dAbsoluteTolerance) 
						shouldBeCompressed.insert(it->blockID);
					else
						niceGuys.insert(it->blockID);

					g.getBlockCollection().release(it->blockID);
				}
				if (profiler !=NULL) profiler->getAgent("AutoCompr::FWT").stop(vBlocksToFWT.size());
							
				if (shouldBeCompressed.size() == 0) break;
				
				if (profiler !=NULL) profiler->getAgent("AutoCompr::compress").start();
				int nCollapsed;
				g.compress(shouldBeCompressed, nCollapsed);
				if (profiler !=NULL) profiler->getAgent("AutoCompr::compress").stop(nCollapsed);
				
				if (profiler !=NULL) profiler->getAgent("boundaries").start();
				g.getBoundaryInfo();
				if (profiler !=NULL) profiler->getAgent("boundaries").stop();
				
				
				if (nCollapsed == 0) break;
				
				if(bVerbose) printf("AutomaticCompression:: Compressed %d, FWT skipped blocks:%d\n", nCollapsed, nSkippedBlocks);
				
				nTotalCollapsed += nCollapsed;
				
				loopCounter++;
				
				if (iMaxLoops>0 && loopCounter>=iMaxLoops) break;
				
			} while(true);
			
			
			if(bKeepHUpdated)
			{
				vector<BlockInfo> vInfo = g.getBlocksInfo();
				for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++)
				{
					typedef typename Grid::GridBlockType B;
					B& b = g.getBlockCollection().lock(it->blockID);
					b.setH(it->h[0]);
					b.bHIsSet = false;
					g.getBlockCollection().release(it->blockID);
				}
			}
			else {
				vector<BlockInfo> vInfo = g.getBlocksInfo();
				for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++)
				{
					g.getBlockCollection().lock(it->blockID).bHIsSet = false;
					g.getBlockCollection().release(it->blockID);
				}
				
			}
			
			return nTotalCollapsed;
		}

		template <typename BlockFWT, typename Grid>
		RefinementResult AutomaticRefinementForLevelsets(Grid& g, BlockFWT& fwt, const double dAbsoluteTolerance,  const int iMaxLevel = -1, const bool bKeepHUpdated = true, const int iMaxLoops=-1,
		MRAG::Profiler* profiler=NULL, double * dMaxDetailAlive= NULL, void (*fillGrid)(Grid& g)=NULL)
		{
		  RefinementResult refinement_result;
		  
		  if(bVerbose) printf("AutomaticRefinementForLevelsets\n");

		int loopCounter = 0;

		if(bVerbose) printf("AutomaticRefinement\n");
		set<int> niceGuys;
		int nSkippedBlocks = 0;
		//int nRefinedBlocks = 0;
			
			do
			{
				fwt.prepare(g.getBlockCollection(), g.getBoundaryInfo());

				vector<BlockInfo> vInfo = g.getBlocksInfo(), vBlocksToFWT;
				for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++)
				{
					if (iMaxLevel>=0 && iMaxLevel <= it->level) continue;
					
					if (niceGuys.find(it->blockID) == niceGuys.end()) // maybe not really a nice guy
						vBlocksToFWT.push_back(*it);
					else
						nSkippedBlocks++;
				}
				
				if (vBlocksToFWT.size() == 0) break;
				
				if (profiler !=NULL) profiler->getAgent("AutoRef::FWT").start();

				set<int> shouldBeRefined;
				//printf("Blocks to FWT (levelset): %d\n", vBlocksToFWT.size());
				for(vector<BlockInfo>::iterator it = vBlocksToFWT.begin(); it != vBlocksToFWT.end(); it++)
				{
					g.getBlockCollection().lock(it->blockID).setH(it->h[0]);
					
					fwt.template multichannel_fwt<0, 0>(*it);
					//printf("getOverAll_DetailMaxMag()=%e >= dAbsoluteTolerance=%e -----> %s\n", fwt.getReport().getOverAll_DetailMaxMag(), dAbsoluteTolerance,fwt.getReport().getOverAll_DetailMaxMag() >= dAbsoluteTolerance ? "to refine!": "skept");
					if (fwt.getReport().getOverAll_DetailMaxMag() >= dAbsoluteTolerance) 
						shouldBeRefined.insert(it->blockID);
					else
						niceGuys.insert(it->blockID);

					g.getBlockCollection().release(it->blockID);
				}
				if (profiler !=NULL) profiler->getAgent("AutoRef::FWT").stop(vBlocksToFWT.size());
				
				
				if (shouldBeRefined.size() == 0) break;
				
				if (profiler !=NULL) profiler->getAgent("AutoRef::refine").start();
				refinement_result  += g.refine(shouldBeRefined);
				if (profiler !=NULL) profiler->getAgent("AutoRef::refine").stop(shouldBeRefined.size());
				if (refinement_result.hasFailed()) break;
				
				if (profiler !=NULL) profiler->getAgent("boundaries").start();
				g.getBoundaryInfo();
				if (profiler !=NULL) profiler->getAgent("boundaries").stop();
				
				if (fillGrid != NULL) fillGrid(g);
				
				if(bVerbose) printf("AutomaticRefinement:: Refined %d, FWT skipped blocks:%d\n", shouldBeRefined.size(), nSkippedBlocks);
				
				loopCounter++;
				
				if (iMaxLoops>0 && loopCounter>=iMaxLoops) break;
				
			} while(true);

			if (dMaxDetailAlive !=NULL)
			{
				double maxVal = dAbsoluteTolerance;

				vector<BlockInfo> vInfo = g.getBlocksInfo();
				for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++)
					if (niceGuys.find(it->blockID) == niceGuys.end())
					{
						g.getBlockCollection().lock(it->blockID).setH(it->h[0]);
						
						fwt.template multichannel_fwt<0, 0>(*it);
						
						maxVal = max(maxVal, fwt.getReport().getOverAll_DetailMaxMag());
						
						g.getBlockCollection().release(it->blockID);
					}

				*dMaxDetailAlive = maxVal;
			}

			if(bKeepHUpdated)
			{
				vector<BlockInfo> vInfo = g.getBlocksInfo();
				for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++)
				{
					typedef typename Grid::GridBlockType B;
					B& b = g.getBlockCollection().lock(it->blockID);
					b.setH(it->h[0]);
					b.bHIsSet = false;
					g.getBlockCollection().release(it->blockID);
				}
			}
			else {
				vector<BlockInfo> vInfo = g.getBlocksInfo();
				for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++)
				{
					g.getBlockCollection().lock(it->blockID).bHIsSet = false;
					g.getBlockCollection().release(it->blockID);
				}
				
			}

			

			return refinement_result;
		}
	} /* namespace Science */
	

}
