/*
 *  MRAG_STDTestL2_Compression.h
 *  MRAG
 *
 *  Created by Babak Hejazi Alhosseini  on 7/28/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */


#pragma once

#include "MRAGCommon.h"
#include "MRAGrid.h"
#include "MRAGRefiner.h"
#include "MRAGCompressor.h"
#include "MRAGBlockLab.h"
#include "MRAGBlockFWT.h"
#include "MRAG_STDTestL2.h"
#include <iostream> 
#include <fstream>
#include <cmath>

#define N_MAX_LEVEL 2


namespace MRAG
{
	template <typename Wavelets, typename Block>
	class MRAG_STDTestL2_Compression:  public MRAG_STDTestL2<Wavelets, Block>//LEVEL 2: correct data representation, correct data-processing
	{
		protected:
			
			MRAG_STDTestL2_Compression():
				MRAG_STDTestL2<Wavelets, Block>()
			{
			}
		public:
			
			double compthresh;
			double comprule;
			double error;
			int whichlevel;
			int NpBeforeRef;
			vector<double> vMaxDetails;
			
			void _babak_ic(Grid<Wavelets, Block>& grid, bool bSetToOne = false)
			{
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			double deltah;	
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					Block& block = grid.getBlockCollection()[info.blockID];
					
				
					deltah = info.h[0];
					for(int iy=0; iy<Block::sizeY; iy++)
						for(int ix=0; ix<Block::sizeX; ix++)
						{
							double x[2];
							info.pos(x, ix, iy);
							//block(ix,iy) =  sin(x[0]+cOffset)*cos(x[1]+cOffset);
							//block(ix,iy) = bSetToOne? 1 : exp(-pow(2.0*sin(2.0*x[0]*M_PI),1.0));
							block(ix,iy) = bSetToOne? 1 : exp(-pow(2.0*sin(2.0*x[0]*M_PI)*sin(2.0*x[1]*M_PI),1.0));
							//block(ix,iy) = 1.0*pow(x[0],2.0);
							//block(ix,iy) = 1.0;
						}
					
				}
				printf("we r sarting at level: %d\n", vInfo[0].level);
				whichlevel = vInfo[0].level;
				printf("(delta_h)^2:%f\n",pow(deltah,2.0));

			}
			
			void _babak_fwt(Grid<Wavelets, Block>& grid)
			{
				BlockFWT<Wavelets, Block> blockfwt;
				
				blockfwt.prepare(grid.getBlockCollection(), grid.getBoundaryInfo());
				
				vector<MRAG::BlockInfo> vInfo = grid.getBlocksInfo();

				vMaxDetails = vector<double>(vInfo.size());
				vector<double> vMinDetails(vInfo.size());
				
				double maxDetail = 0;
				double globmax=0.0;
				double globmin=1000.0;

				for(int i=0; i<vInfo.size(); i++)
				{
					blockfwt.fwt(vInfo[i]);
					
					const FWTReport<>& report= blockfwt.getReport();
					maxDetail = std::max(maxDetail, report.getCoeffMaxMag(1,1));
				
					// get maximum and minimum coefficients not considering the (0,0)
					
					vMaxDetails[i] = max(report.getCoeffMaxMag(1,0),report.getCoeffMaxMag(0,1));
					vMaxDetails[i] = max(vMaxDetails[i],report.getCoeffMaxMag(1,1));
					vMinDetails[i] = min(report.getCoeffMinMag(1,0),report.getCoeffMinMag(0,1));
					vMinDetails[i] = min(vMinDetails[i],report.getCoeffMinMag(1,1));
					
					//printf("maxDetail: %f,%f \n",report.getCoeffMaxMag(1,0),report.getCoeffMinMag(1,0));
					//printf("coordinates:%d,%d,%d\n",vInfo[i].index[0],vInfo[i].index[1],vInfo[i].level);
					
					// get the global maximum and minimum
					if (vMaxDetails[i] > globmax) globmax = vMaxDetails[i];
					if (vMinDetails[i] < globmin) globmin = vMinDetails[i];
				}
						
				printf("global max and min are: %f, %f\n", globmax,globmin);
				
				//threshold is set
				compthresh = comprule*globmax + (1.0-comprule)*globmin;
			}
			
			vector<double> _fwt(Grid<Wavelets, Block>& grid)
			{
			BlockFWT<Wavelets, Block> blockfwt;
			
			blockfwt.prepare(grid.getBlockCollection(), grid.getBoundaryInfo());
			
			vector<MRAG::BlockInfo> vInfo = grid.getBlocksInfo();
			
			vector<double> vMaxDetail(vInfo.size());
			vector<double> vMinDetail(vInfo.size());
			
			double maxDetail = 0;
			double globmax=0.0;
			double globmin=1000.0;
			
			for(int i=0; i<vInfo.size(); i++)
			{
				blockfwt.template fwt<0>(vInfo[i]);
				
				const FWTReport<>& report= blockfwt.getReport();
				maxDetail = std::max(maxDetail, report.getCoeffMaxMag<0>(1,1));
				
				// get maximum and minimum coefficients
				vMaxDetail[i] = max(report.getCoeffMaxMag<0>(1,0),report.getCoeffMaxMag<0>(0,1));
				vMaxDetail[i] = max(vMaxDetail[i],report.getCoeffMaxMag<0>(1,1));
				vMinDetail[i] = min(report.getCoeffMinMag<0>(1,0),report.getCoeffMinMag<0>(0,1));
				vMinDetail[i] = min(vMinDetail[i],report.getCoeffMinMag<0>(1,1));
				//printf("maxDetail: %f,%f \n",report.getCoeffMaxMag(1,0),report.getCoeffMinMag(1,0));
				//printf("coordinates:%d,%d,%d\n",vInfo[i].index[0],vInfo[i].index[1],vInfo[i].level);
				
				// get the global maximum and minimum
				if (vMaxDetail[i] > globmax) globmax = vMaxDetail[i];
				if (vMinDetail[i] < globmin) globmin = vMinDetail[i];
			}
			

			//threshold
			compthresh = comprule*globmax + (1.0-comprule)*globmin;
			return vMaxDetail;
			}
			
			void _babak_computation(Grid<Wavelets, Block>& grid)  
			{
				FILE * file = fopen("error.dat","w");
				vector<MRAG::BlockInfo> vInfo = grid.getBlocksInfo();
				error = 0.0;
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					Block& block = grid.getBlockCollection()[info.blockID];
					
					const double h[2]= {pow(2.,-info.level), pow(2.,-info.level)};
					const double start[2] = {info.index[0]*h[0],info.index[1]*h[1]};
					const int n[2] = {Block::sizeX+1, Block::sizeY+1};
					const double d[2] = {h[0]/(n[0]-1),h[1]/(n[1]-1)};
					const double cOffset = Wavelets::CenteringOffset;
					
					for(int iy=0; iy<Block::sizeY; iy++)
						for(int ix=0; ix<Block::sizeX; ix++)
						{
							const double x[2] = {start[0]+(ix)*d[0], start[1]+(iy)*d[1]};
							const double f = error;
							double exactValue;
							exactValue = exp(-pow(2.0*sin(2.0*x[0]*M_PI)*sin(2.0*x[1]*M_PI),1.0));//sin(2.0*x[0]*M_PI);
							//exactValue = 1.0;
							error = max(error,fabs(exactValue-block(ix,iy)));
							fprintf(file,"%f %f %f %f\n",x[0],x[1],exactValue,block(ix,iy));
							if (f!=error) printf("%d %d %d %d -> error =%f\n", ix, iy, info.index[0], info.index[1], error);
							//error += pow(fabs(sin(1.0*x[0])-block(ix,iy)),2.0);
							//error = sqrt(error/(n[0]*n[1]));
							//error = error/(n[0]*n[1]);
						}
					
				}
				fclose(file);
				printf("L_inf error is: %f\n",error);
			}

			
			void _babak_compression(Grid<Wavelets, Block>& grid)  
			{
				for(int iCompressionStep = 0; iCompressionStep<1; iCompressionStep++)
				{
					vector<BlockInfo> vInfo = grid.getBlocksInfo();
					set<int> shouldBeCompressed;
					
					for(int i=0; i<vInfo.size(); i++)
						//thresholding criterion applied for compression
						if (vMaxDetails[i] < compthresh) shouldBeCompressed.insert(vInfo[i].blockID);
					
					int nRequested = shouldBeCompressed.size();
					int nCollapsed = 0;
					grid.compress(shouldBeCompressed, nCollapsed);
					printf("nCollapsed %d of %d requested\n", nCollapsed,nRequested);
				}
				
			}
			
			void _compress(Grid<Wavelets, Block>& grid)
			{
				unsigned int totalCollapsed = 0, totalRequested = 0;
				for(int iCompressionStep = 0; iCompressionStep<4; iCompressionStep++)
				{
					vector<BlockInfo> vInfo = grid.getBlocksInfo();
					set<int> shouldBeCompressed;
					vector<double> vMaxDetail = _fwt(grid);
					
					for(int i=0; i<vInfo.size(); i++)
						if (vMaxDetail[i] < compthresh) shouldBeCompressed.insert(vInfo[i].blockID);
					
					int nRequested = shouldBeCompressed.size();
					int nCollapsed = 0;
					grid.compress(shouldBeCompressed, nCollapsed);
					totalCollapsed+=nCollapsed;
					totalRequested+=nRequested;
				}
				
				printf("nCollapsed %d of %d requested\n", totalCollapsed,totalRequested);

		}
			
			
			void _babak_refinements(Grid<Wavelets, Block>& grid) 
			{
				for(int iRef = 0; iRef<1; iRef++)
				{
					set<int> shouldBeRefined;
					vector<MRAG::BlockInfo> vInfo = grid.getBlocksInfo();
					
					for(int i=0; i<vInfo.size(); i++){
								shouldBeRefined.insert(vInfo[i].blockID);
							
					}
					
					grid.refine(shouldBeRefined);
				}
			}
			
			void _refineUpto(Grid<Wavelets, Block>& grid, unsigned int level) 
			{
				while(true){
					set<int> shouldBeRefined;
					vector<MRAG::BlockInfo> vInfo = grid.getBlocksInfo();
					for(int i=0; i<vInfo.size(); i++){
						if(vInfo[i].level < level) 
						{
							shouldBeRefined.insert(vInfo[i].blockID);
							
						}
					}
					if(shouldBeRefined.size()==0) break;
					grid.refine(shouldBeRefined);
				}
						
			}

			void _babak_refineUpto(Grid<Wavelets, Block>& grid, unsigned int level) 
			{
				for(int iRef = 0; iRef<2; iRef++)
				{
					set<int> shouldBeRefined;
					vector<MRAG::BlockInfo> vInfo = grid.getBlocksInfo();
					
					for(int i=0; i<vInfo.size(); i++){
							if(vInfo[i].level < level) 
							{
								shouldBeRefined.insert(vInfo[i].blockID);
								
							}
					}
					
					grid.refine(shouldBeRefined);
				}
			}
			

			
		
			virtual bool run()
			{
				// OVERLOAD THIS
				try
					{
						
						
						FILE * file2 = fopen("L_inf_16_4th_new.txt","w");
						
						for (comprule=0.0; comprule<1.0; comprule+=0.1)
						{
							Grid<Wavelets, Block>  * myGrid = new Grid<Wavelets, Block>(8,8);
							myGrid->setCompressor(&this->compressor);
							myGrid->setRefiner(&this->refiner);
							

							{
								/* ic at level L
								fwt
								compress
								refine to level L
								compute error 
								*/
							}
														
							printf("Start Test This is Babak with comprule=%f\n", comprule);
							myGrid->getBlocksInfo();
							printf("Done BlockInfo\n");

							MRAG_STDTestL2_Compression<Wavelets, Block>::_babak_ic(*myGrid);
							printf("Done IC\n");
							
							//MRAG_STDTestL2_Compression<Wavelets, Block>::_babak_fwt(*myGrid);
							printf("Done FWT\n");
							printf("Threshold is: %f\n",compthresh);
							
							
							MRAG_STDTestL2_Compression<Wavelets, Block>::_compress(*myGrid);
							printf("Done Compression\n");
							// get the number of points beore refinement
							NpBeforeRef = myGrid->getBlocksInfo().size()*Block::sizeX*Block::sizeY;
							MRAG_STDTestL2_Compression<Wavelets, Block>::_refineUpto(*myGrid,whichlevel);
							
							printf("Done Refinement\n");
							MRAG_STDTestL2_Compression<Wavelets, Block>::_babak_computation(*myGrid);
							
							
							printf("End Test\n");
							
							if (comprule>0.0) fprintf(file2,"%d	%f\n",NpBeforeRef,error);
							
							delete myGrid;
						}
						
						//output the data using gnuplot
						{
						char output_name[1000]; 
						sprintf(output_name,"%s",file2);
						FILE * fgnuplot = popen("/usr/local/bin/gnuplot \n","w");
						assert(fgnuplot != NULL);
						fprintf(fgnuplot, "set xlabel \" number of points \" \n");
						fprintf(fgnuplot, "set ylabel \"Linf_error\" \n");
						fprintf(fgnuplot, "plot \"%s\" using 1:2 title 'L-inf postcomp' \n",file2);
						fprintf(fgnuplot, "set terminal postscript color \n");
						fprintf(fgnuplot, "set out \"%s.ps\"\n",output_name);
						fprintf(fgnuplot, "replot\n");
						fclose(fgnuplot);
						}
						
						
					}
				catch(...)
				{
					return false;
				}
				
				return true;
				
			}
			
		public:
			static void runTests()
			{
				MRAG_STDTestL2_Compression<Wavelets, Block> test;
				
				test.run();
			}
		};
	
	
}