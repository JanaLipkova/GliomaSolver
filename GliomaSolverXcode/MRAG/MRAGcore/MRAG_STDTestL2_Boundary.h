/*
 *  MRAG_STDTestL2_Boundary.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/25/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 *  //modified by mq: boundary representation test.
 *
 */

#pragma once

#include "MRAGCommon.h"
#include "MRAGrid.h"
#include "MRAGRefiner.h"
#include "MRAGCompressor.h"
#include "MRAGBlockLab.h"
#include "MRAGBlockFWT.h"
#include <string>
#include <fstream>

//const int sX=2;
//const int sY=2;
const int sSt[3]={-2,-2,0};
const int sNd[3]={3,3,1};

//#define mq_debug




namespace MRAG
{
	template <typename Wavelets, typename Block>
	class MRAG_STDTestL2_Boundary //LEVEL 2: correct data representation at boundaries
	{
		protected:
			Grid<Wavelets, Block> grid;
			BlockLab<Block> lab; //need this for the BCs
			BoundaryInfo* bInfo;
			Refiner refiner;

            //Return Type
			typedef float BaseType;
			
			//Constructor:
			MRAG_STDTestL2_Boundary(int NumOfBlocks):
			grid(NumOfBlocks,NumOfBlocks),refiner()
			{
			  grid.setRefiner(&refiner);
			  bInfo=0; //Set boundaryInfo to zero, initially.
			}
			
			//Destructor:
			~MRAG_STDTestL2_Boundary()
			{
				delete bInfo;
			}
			
		    
			/*put some initial condition on the domain*/
			/*also safe it to oldG, and initialG*/
			void _some_ic(const int* steStart, const int* steEnd)
			{
			   			
				vector<BlockInfo> vInfo = grid.getBlocksInfo();
			

				
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
							block(ix,iy).cur =  sin((x[0]+cOffset)*2.0*M_PI)*cos((x[1]+cOffset)*2.0*M_PI);
							block(ix,iy).setIC();
							block(ix,iy).update();

												
						}
					

				}
				
				//build the boundary_info: (?)
				bInfo=grid.createBoundaryInfo(sSt,sNd);
				//DELETE THIS afterwards. (...->ok doing this now in the destructor)
			}
		
			void _some_computation(int dt, int ux, int uy) //"integer-advection (obsolete)"
			{
				int steStart[3] ={-abs(ux*dt),-abs(uy*dt),0};
				int steEnd[3] = {abs(ux*dt)+1,abs(uy*dt)+1,1};
				
				
				#ifdef mq_debug
				std::cout<<"requiring stencils:"  <<std::endl;
				std::cout<<"start:" << steStart[0]<< " " <<  steStart[1]<< " " << steStart[2]<< " " <<std::endl;
				std::cout<<"end:" << steEnd[0] << " " << steEnd[1] <<" " << steEnd[2]  <<std::endl;
				std::cout <<"shifting data: lab(ix+"<< dt*ux << ",iy+"<<dt*uy<<")=oldlab(ix,iy)"<<std::endl;
			    #endif mq_debug
				
				vector<BlockInfo> vInfo = grid.getBlocksInfo();
				
				//lab.prepare(grid.getBlockCollection(), grid.getBoundaryInfo(),steStart,steEnd);
				lab.prepare(grid.getBlockCollection(), *bInfo,steStart,steEnd);
				
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
    				lab.load(info);
					
					for(int iy=0; iy<Block::sizeY; iy++)
						for(int ix=0; ix<Block::sizeX; ix++)
						{
	
					    #ifdef mq_debug
					    std::cout << "current,old,diff:" << lab(ix,iy).cur << " " << lab(ix,iy).old << " " << lab(ix,iy).cur-lab(ix,iy).old <<std::endl ;
					    #endif
						lab(ix,iy).cur = lab(ix-dt*ux,iy-dt*uy).old; 
						}
					lab.flush(info);
				 }
			 }
				 
				 
			void _some_advection_1stOrder(double dt, double ux, double uy) //upwinding advection
			{
					int steStart[3] ={-1,-1,0};
					int steEnd[3] = {2,2,1};
					
                    std::cout << "Doing in a 1st order upwinding step with: ux= " <<  ux << " uy= " <<uy<<std::endl;
					
					int selx=1;
					int sely=1;
					if (ux>0) selx=-1;
					if (uy>0) sely=-1;

													
					vector<BlockInfo> vInfo = grid.getBlocksInfo();

					
					//lab.prepare(grid.getBlockCollection(), grid.getBoundaryInfo(),steStart,steEnd);
					lab.prepare(grid.getBlockCollection(), *bInfo,steStart,steEnd);
					
					for(int i=0; i<vInfo.size(); i++)
					{
						BlockInfo& info = vInfo[i];
						lab.load(info);
					
						const double h[2]= {pow(2.,-info.level), pow(2.,-info.level)};
						
						for(int iy=0; iy<Block::sizeY; iy++)
							for(int ix=0; ix<Block::sizeX; ix++)
							{
		
							  lab(ix,iy).cur = lab(ix,iy).old-dt*(fabs(ux)/h[0]*(lab(ix+selx,iy).old - lab(ix+selx-1,iy).old)+fabs(uy)/h[1]*(lab(ix,iy+sely).old - lab(ix,iy-sely-1).old)); 
							}
						lab.flush(info);
					 }

				
				
				//update "oldG" (dont need a lab here, but direct access to BlockCollection)
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info2 = vInfo[i];
					Block& block = grid.getBlockCollection()[info2.blockID];
		
					for(int iy=0; iy<Block::sizeY; iy++)
						for(int ix=0; ix<Block::sizeX; ix++)
						{
							block(ix,iy).update();
						    #ifdef mq_debug
    					    std::cout << "[chkUpd:]current,old,diff:" << block(ix,iy).cur << " " << block(ix,iy).old << " " << block(ix,iy).cur-block(ix,iy).old <<std::endl ;
         				    #endif

						}
						
				 }

												
			}
			
			
			void _some_advection_2ndOrder(double dt, double ux, double uy) //upwinding advection
			{
			
			        double dudx,dudy;
					int steStart[3] ={-1,-1,0};
					int steEnd[3] = {2,2,1};
					
                    std::cout << "Doing in a 2ndorder upwinding step with: ux= " <<  ux << " uy= " <<uy<<std::endl;
					
					int selx=1;
					int sely=1;
					if (ux>double(0.0)) selx=-1;
					if (uy>double(0.0)) sely=-1;

													
					vector<BlockInfo> vInfo = grid.getBlocksInfo();

					
					//lab.prepare(grid.getBlockCollection(), grid.getBoundaryInfo(),steStart,steEnd);
					lab.prepare(grid.getBlockCollection(), *bInfo,steStart,steEnd);
					
					for(int i=0; i<vInfo.size(); i++)
					{
						BlockInfo& info = vInfo[i];
						lab.load(info);
					
						const double h[2]= {pow(2.,-info.level), pow(2.,-info.level)};
						
						for(int iy=0; iy<Block::sizeY; iy++)
							for(int ix=0; ix<Block::sizeX; ix++)
							{
		                      //2nd order, FD Upwinding
                              dudx=1.0/(2.0*h[0])*(3.0*lab(ix,iy).old-4.0*lab(ix+selx,iy).old+lab(ix-(selx*2),iy).old);
							  dudy=1.0/(2.0*h[1])*(3.0*lab(ix,iy).old-4.0*lab(ix,iy+sely).old+lab(ix,iy-(sely*2)).old);
							  lab(ix,iy).cur=lab(ix,iy).old-dt*(fabs(ux)*dudx+fabs(uy)*dudy);
							}
						lab.flush(info);
					 }

				
				
				//update "oldG" (dont need a lab here, but direct access to BlockCollection)
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info2 = vInfo[i];
					Block& block = grid.getBlockCollection()[info2.blockID];
		
					for(int iy=0; iy<Block::sizeY; iy++)
						for(int ix=0; ix<Block::sizeX; ix++)
						{
							block(ix,iy).update();
						    #ifdef mq_debug
    					    std::cout << "[chkUpd:]current,old,diff:" << block(ix,iy).cur << " " << block(ix,iy).old << " " << block(ix,iy).cur-block(ix,iy).old <<std::endl ;
         				    #endif

						}
						
				 }

												
			}
			
			
			BaseType _some_error_calculation()
			{
				vector<BlockInfo> vInfo = grid.getBlocksInfo();

				BaseType Linfty = BaseType(-1.0);
				BaseType abs_error;
				
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					Block& block = grid.getBlockCollection()[info.blockID];
					
					
					for(int iy=0; iy<Block::sizeY; iy++)
						{
							for(int ix=0; ix<Block::sizeX; ix++)
							{
								abs_error   = fabs(block(ix,iy).cur - block(ix,iy).ic); 
								if(abs_error > Linfty) Linfty = abs_error;
								#ifdef mq_debug
								if(abs_error>0.5 || abs_error < 0.0) 
									{
									
										printf("Funny error at ix=%d,iy=%d,l=%d  ini=%f, current=%f\n",ix,iy,info.level,block(ix,iy).ic,block(ix,iy).cur); 
									}
								#endif
							 }
						}
					
				}				
				std::cout << "Linfty = " << Linfty << std::endl;
				return Linfty;
			}
			
			
			void _dump2file(std::string fname) {

				std::ofstream outf(fname.c_str());
				
				vector<BlockInfo> vInfo = grid.getBlocksInfo();

				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					Block&block = grid.getBlockCollection()[info.blockID];
					
				
					const double h[2]= {pow(2.,-info.level), pow(2.,-info.level)};
					const double start[2] = {info.index[0]*h[0],info.index[1]*h[1]};
					const int n[2] = {Block::sizeX+1, Block::sizeY+1};
					const double d[2] = {h[0]/(n[0]-1),h[1]/(n[1]-1)};
				
					for(int iy=0; iy<Block::sizeY; iy++)
					   {
						for(int ix=0; ix<Block::sizeX; ix++)
							{
    							const double x[2] = {start[0]+(ix)*d[0], start[1]+(iy)*d[1]};
								outf << x[0] << " " << x[1]<< " " << block(ix,iy).cur << std::endl ;
							}
					    }
					
					
				 }
				
				
                outf.close();
			
			}
			
			
			void _some_refinements() 
			{
				for(int iRef = 0; iRef<2; iRef++)
				{
					set<int> shouldBeRefined;
					vector<MRAG::BlockInfo> vInfo = grid.getBlocksInfo();
					
					for(int i=0; i<vInfo.size(); i++)
						if (drand48() > 0.7) shouldBeRefined.insert(vInfo[i].blockID);
					
					grid.refine(shouldBeRefined);
				}
			}

		
					
			virtual bool run(int LevelOfRefinements, int NumOfBlocks, int OrderFD, int NumberOfSteps)
			{

				try
				{
                   double ux(1.0),uy(1.0);
				   int steStart[3] ={-2,-2,0};
				   int steEnd[3] = {3,3,1};
				   double CFL=.5;
				   double dt;
				   double err;
				   

					printf("Start Test\n");
					for (int i=0; i < LevelOfRefinements; ++i)
					{
					_some_refinements();
					}
					
				    //find finest level:
					int finest_level=-1;
	   				vector<BlockInfo> vInfo = grid.getBlocksInfo();
					for (int i=0;i<vInfo.size();++i)
						{
							finest_level=std::max(finest_level,(int)vInfo[i].level);
						}
					//Determine timestep, based on CFL
					const double h[2]= {pow(2.,-finest_level), pow(2.,-finest_level)};
					ux=h[0];
					uy=h[1];
					dt=CFL*h[0]/ux;
					dt=std::max(dt,CFL*h[1]/uy);
					std::cout <<" CFL= "<<CFL << std::endl;
					std::cout <<" dt = " <<dt <<std::endl;

				
					_some_ic(steStart,steEnd);
					_dump2file("iniC.dat");
					
					
					for (int i=0; i<NumberOfSteps;++i)
					{
						if (OrderFD==1)
						{
						   _some_advection_1stOrder(dt,ux,uy);
						}
						else
						{
							_some_advection_2ndOrder(dt,ux,uy);
						}
					}
					_dump2file("t00.dat");
					
					for (int i=0; i<NumberOfSteps;++i)
					{
						if (OrderFD==1)
						{
						   _some_advection_1stOrder(dt,-ux,-uy);
						}
						else
						{
							_some_advection_2ndOrder(dt,-ux,-uy);
						}
                    }
				
					_dump2file("t01.dat");
					
					err=_some_error_calculation();
					std::cout << "Error is: " << err << std::endl;
					std::cout << "CFL was: " << CFL << std::endl;
					std::cout <<" dt = " <<dt <<std::endl;
					std::cout << "hx,err: " << h[0] << " " << err <<std::endl;
					
					printf("End Test\n");
				}
				catch(...)
				{
					return false;
				}
				
				return true;
			}
			
		public:
			static void runTests(int LevelOfRefinements=2, int NumOfBlocks=2, int OrderFD=2, int NumberOfSteps=1)
			{
				MRAG_STDTestL2_Boundary<Wavelets, Block> test(NumOfBlocks);
				
				test.run(LevelOfRefinements,NumOfBlocks,OrderFD,NumberOfSteps);
                //usage example: 
			    //MRAG::MRAG_STDTestL2_Boundary<MRAG::Wavelets_Interp2Order, MRAG::Block<misc::memNum<double>,16,16,1> >::runTests(2,4,2,1);

			}
		};
	
	
}

namespace misc 
{



//use setIc to save IC based on current status
//use update to save current status on "old"
template <typename NType>
struct memNum
{
	NType cur;
	NType old;
	NType ic;

	void operator +=(memNum  v)
	{
		cur+=v.cur;
		old+=v.old;
		ic+=v.ic;
	}
	
	void operator=(float v)
	{
		cur= v;
		old= v;
		ic= v;
	}
	
	void setIC()
	{
	  ic=cur;
	}
	
	void update()
	{
	  old=cur;
	}

};
template <typename NType, typename A>
memNum<NType> operator*(const memNum<NType>& p, A v)
{
	memNum<NType> t = {p.cur*v, p.old*v,p.ic*v};
	return t;
}

template <typename NType>
memNum<NType> operator*(const memNum<NType>& p, float v)
{
	memNum<NType> t = {p.cur*v, p.old*v,p.ic*v};
	return t;
}




}