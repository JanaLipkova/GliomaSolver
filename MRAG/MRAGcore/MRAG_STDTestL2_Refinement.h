/*
 *  MRAG_STDTestL2.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/25/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 *  $Log: MRAG_STDTestL2_Refinement.h,v $
 *  Revision 1.1.1.1  2011/03/08 14:49:55  bbayati
 *  Imported Sources
 *
 *  Revision 1.4  2010/01/06 10:27:39  hbabak
 *  *** empty log message ***
 *
 *  Revision 1.3  2008/07/28 10:18:20  michaebe
 *  Automated refinement, pasted reference results
 *
 *  Revision 1.2  2008/07/28 09:30:41  michaebe
 *  Does logging work?
 * 
 */

#pragma once

#include "MRAGCommon.h"
#include "MRAGrid.h"
#include "MRAGRefiner.h"
#include "MRAGCompressor.h"
#include "MRAGBlockLab.h"
#include "MRAGBlockFWT.h"

namespace MRAG
{
	template <typename Wavelets, typename Block>
	class MRAG_STDTestL2_Refinement //LEVEL 2: correct data representation, correct data-processing
		{
		protected:
			Grid<Wavelets, Block> grid;
			Refiner refiner;
			
			MRAG_STDTestL2_Refinement():
			grid(2,2), refiner()
			{
				grid.setRefiner(&refiner);
			}
			
			void _some_ic()
			{
				vector<BlockInfo> vInfo = grid.getBlocksInfo();
				
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					Block& block = grid.getBlockCollection()[info.blockID];
					 					
					for(int iy=0; iy<Block::sizeY; iy++)
						for(int ix=0; ix<Block::sizeX; ix++)
						{
							double x[3];
							info.pos(x,ix,iy);
							//printf("info.origin[0] = %e\n",info.origin[0]-0.5*info.h[0]);
							//exit(0);
							block(ix,iy) =  -1./(2.*M_PI)*(cos(2.*M_PI*(x[0]+0.5*info.h[0]))-cos(2.*M_PI*(x[0]-0.5*info.h[0])))/info.h[0];//sin((x[0])*2.0*M_PI);//*cos((x[1])*2.0*M_PI);
						}
				}
			}
			
			void _some_refinements() 
			{
				for(int iRef = 0; iRef<1; iRef++)
				{
					set<int> shouldBeRefined;
					vector<MRAG::BlockInfo> vInfo = grid.getBlocksInfo();
					
					for(int i=0; i<vInfo.size(); i++)
						shouldBeRefined.insert(vInfo[i].blockID); /* everything should be refined */
					
					grid.refine(shouldBeRefined);
				}
			}
			
			
			
			double _some_error_calculation()
			{
				vector<BlockInfo> vInfo = grid.getBlocksInfo();
				
				//lab.prepare(grid.getBlockCollection(), grid.getBoundaryInfo());
				double Linfty = -1.0e10;
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					Block& block = grid.getBlockCollection()[info.blockID];
					
			
					/* probably don't need the lab for this ? */
					//lab.load(info);
					for(int iy=0; iy<Block::sizeY; iy++)
						for(int ix=0; ix<Block::sizeX; ix++){
			
							double x[3];
							info.pos(x,ix,iy);
							double exact_value =  -1./(2.*M_PI)*(cos(2.*M_PI*(x[0]+0.5*info.h[0]))-cos(2.*M_PI*(x[0]-0.5*info.h[0])))/info.h[0];//sin((x[0])*2.0*M_PI);//*cos((x[1])*2.0*M_PI);
							double abs_error   = fabs(block(ix,iy) - exact_value);
							if(abs_error > Linfty) Linfty = abs_error;
							if(abs_error>0.1) {
								printf("Funny error at ix=%d,iy=%d,l=%d  e=%f, block=%f\n",ix,iy,info.level,exact_value,block(ix,iy));
							}
						}
					
					
					//lab.flush(info);
				}				
				std::cout << "Linfty = " << Linfty << std::endl;
				return Linfty;
			}
			
			
			virtual bool run()
			{
				/* Expected output for 4th order wavelets:
				 0	6.878051e-05
				 1	4.340549e-06
				 2	2.719393e-07
				 3	1.700645e-08
				 4	1.063063e-09
				 
				 Expected output for 2nd order wavelets:
				 0	9.515058e-03
				 1	2.401840e-03
				 2	6.019092e-04
				 3	1.505680e-04
				 4	3.764766e-05
				 
				
				UnExpected output for Average 5th order wavelets:
				
				 0	2.748648e-03
				 1	3.520020e-04
				 2	4.426953e-05
				 3	5.545773e-06
				 4	7.060442e-07
				 */
				const unsigned int nbase = 5;
				double linfty[nbase];
				for(unsigned int ibase = 0;ibase<nbase;++ibase){
					try
					{
						_some_ic();	/* Set some initial conditions */
						_some_refinements();
						linfty[ibase] = _some_error_calculation();
						printf("End Test\n");
					}
					catch(...)
					{
						return false;
					}
				}
				
				printf("\tbase\tlinfty\n");
				for(unsigned int ibase = 0;ibase<nbase;++ibase){
					printf("\t%d\t%e\n",ibase,linfty[ibase]);
				}
				return true;
			}
			
		public:
			static void runTests()
			{
				MRAG_STDTestL2_Refinement<Wavelets, Block> test;
				
				test.run();
			}
		};
	
	
}