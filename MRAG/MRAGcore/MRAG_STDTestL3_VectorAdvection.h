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
#include "../MRAGScience/MRAGScienceCore.h"
#include "MRAG_SmartBoundaryInfo.h"
#include <string>
#include <fstream>
#include <sstream>

//const int sX=2;
//const int sY=2;
const int sSt[3]={-4,-4,0};
const int sNd[3]={5,5,1};
const int resJump=2;

//#define mq_debug



namespace misc 
{
	//use update to save current status on "old"
	//toDo: check loops for fast access. (which index is varying faster in c++ [][] arrays?)
	template <typename NType, int VecDim=3, int SpaceDim=2, int tStencil=2>
	struct tAdVector
	{
		
		typedef NType BType;
		
		BType data[VecDim][tStencil];
		BType velocity[SpaceDim];
		
		//Helpers, only for error checking purposes:
		BType ic[VecDim];
		
		void setIC()
		{
			for (int i=0;i<VecDim;++i)
			{
				ic[i]=data[i][0];
			}
		}
		
		
		// Operators, and Functions.
		
		
		typedef  tAdVector< BType,VecDim,SpaceDim,tStencil> thisType; 
		
		void operator +=(const thisType & v)
		{
			for (int i=0;i<VecDim;++i)
			{
				for (int j=0;j<tStencil; ++j)
				{
					data[i][j]+=v.data[i][j];
				}
			}
			
			for (int i=0;i<SpaceDim;++i)
			{
				velocity[i]+=v.velocity[i];	
			}
			
		}
		
		
		
		thisType& operator=(const thisType & v)
		{
			for (int i=0;i<VecDim;++i)
			{
				for (int j=0;j<tStencil; ++j)
				{
					data[i][j]=v.data[i][j];
				}
			}
			
			for (int i=0;i<SpaceDim;++i)
			{
				velocity[i]=v.velocity[i];	
			}
			return *this;
			
		}
		
		BType& operator()(int i, int j=0) //idea is to access with (vecComponent,t=0,t=-1...)
		{
			assert(i<VecDim);
			assert(-j<tStencil);
			return data[i][-j];
		}
		
		
		void update()
		{
			
			for (int j=tStencil-1;j>0;--j) //start at oldest value, and overwrite with newer values.
			{
				for (int i=0;i<VecDim;++i)
				{
					data[i][j]=data[i][j-1] ;  
				}
			}
			
		}
		
		
		
		
	};
	
	//Overloading Scalar Multiplication with someType (for convenience).
	template <typename NType, int VecDim, int SpaceDim, int tStencil, typename someType >
	tAdVector< NType,VecDim,SpaceDim,tStencil> operator*(const tAdVector< NType,VecDim,SpaceDim,tStencil> & p, someType v)
	{
		tAdVector< NType,VecDim,SpaceDim,tStencil> t;
		
		for (int i=0;i<VecDim;++i)
		{
			for (int j=0;j<tStencil; ++j)
		    {
				t.data[i][j]=v*p.data[i][j];
			}
		}
		
		for (int i=0;i<SpaceDim;++i)
		{
			t.velocity[i]=v*p.velocity[i];	
		}
		
		return t;
	}
	
	//Overloading Scalar Multiplication with float (required by MRAG).
	template <typename NType, int VecDim, int SpaceDim, int tStencil>
	tAdVector< NType,VecDim,SpaceDim,tStencil> operator*(const tAdVector< NType,VecDim,SpaceDim,tStencil>& p, float v)
	{
		tAdVector< NType,VecDim,SpaceDim,tStencil> t;
		
		for (int i=0;i<VecDim;++i)
		{
			for (int j=0;j<tStencil; ++j)
		    {
				t.data[i][j]=v*p.data[i][j];
			}
		}
		
		for (int i=0;i<SpaceDim;++i)
		{
			t.velocity[i]=v*p.velocity[i];	
		}
		
		return t;
	}
	
	
	
} //namespace misc




namespace MRAG
{

template <typename T, int i > inline Real burgers_projector_impl(const T&t)
{
	return (Real) t.data[0][0];
}
	
make_projector(burgers_projector, burgers_projector_impl)
	
	
template <typename Wavelets, typename Block>
class MRAG_STDTestL3_VectorAdvection //LEVEL 3: a "Real Life" Advection problem
	{
	protected:
		Grid<Wavelets, Block> grid;
		BlockLab<Block> lab; //need this for the BCs
		Refiner refiner;
		Compressor compressor;
		BlockFWT<Wavelets, Block, burgers_projector> blockfwt;
		//BoundaryInfo* bInfo;
		//bool bOwnerOfBoundaryInfo;
		SmartBoundaryInfo *myBoundaries;
		
		
		//Return Type (typename necessary to get templated type)
		typedef typename Block::ElementType::BType BaseType;
		
		//Vector Size (currently fixed, how about get int from type?):
		static const int vDim=3;
		static const int sDim=2;
		
		//Constructor:
		MRAG_STDTestL3_VectorAdvection(int NumOfBlocks):
		grid(NumOfBlocks,NumOfBlocks),refiner(resJump), compressor(resJump),
		lab(), blockfwt()//, myBoundaries()
		{
			grid.setRefiner(&refiner);
			grid.setCompressor(&compressor);
			// bInfo=0; //Set boundaryInfo to zero, initially.
			myBoundaries=new SmartBoundaryInfo();
		}
		
		//Destructor:
		~MRAG_STDTestL3_VectorAdvection()
		{
			delete myBoundaries;
		}
		
		
		
		static void _ic_func(BaseType res[vDim], BaseType x[sDim], int icType)
		{
		  bool hi;
		  int ds;
		  BaseType ofs[2]={0.25,0.45};
		  BaseType rDisc=0.15;
		  BaseType vAmp[3]={1.0,1.0,1.0};
		  
		  for (int ik=0; ik<vDim; ++ik)
		  {
		  switch (icType)
		  {
		  case 0: 
		     res[ik]=vAmp[ik]*sin(x[0]*2.0*M_PI)*cos(x[1]*2.0*M_PI);
		    break;
		  case 1: 
			res[ik]=vAmp[ik];
			break;
		  case 2:
		    const BaseType r = sqrt(pow(x[0]-ofs[0],2) + pow(x[1]-ofs[1], 2) );
		    r<BaseType(rDisc)? res[ik]=1: res[ik]=0; 
		    break;
			
		  case 3:
		     hi=true;
			 for (int d=0;d<sDim;++d)
			 {
			 if (x[d] < 0.5) hi=false;
			 }
			 hi? res[ik]=vAmp[ik]:res[ik]=0;
			 break;
		  case 4:
			 res[ik]= vAmp[ik]*sin(x[0]*2.0*M_PI);
			 break;
		  case 5:
		     hi=true;
			 for (int d=0;d<sDim;++d)
			 {
			 if ((x[d] < 1.0/3.0) || (x[d] > 2.0/3.0)) hi=false;
			 }
			 //if ((x[d] < 1.0/3.0) || (x[d] > 2.0/3.0)) res[ik]= -vAmp[ik];
			 //else res[ik]= vAmp[ik];
			 hi? res[ik]=vAmp[ik]:res[ik]=-vAmp[ik];
			 break;
		  case 6:
		     hi=true;
			 for (int d=0;d<sDim;++d)
			 {
			 if ((x[d] < 1.0/3.0) || (x[d] > 2.0/3.0)) hi=false;
			 }
			 hi? res[ik]=vAmp[ik]:res[ik]=0.0;
			 break;
			
		  case 7: //ND-Shock-and-Rarefaction
		     hi=true;
			 for (int d=0;d<sDim;++d)
			 {
			 if ((x[d] < 1.0/4.0) || (x[d] > 3.0/4.0)) hi=false;
			 }
			 hi? res[ik]=vAmp[ik]:res[ik]=0.0;
			 break;
		  case 8: //1D-Shock-and-Rarefaction
		     ds=0;
			 hi=true;
			 if ((x[ds] < 1.0/4.0) || (x[ds] > 3.0/4.0)) hi=false;
			 hi? res[ik]=vAmp[ik]:res[ik]=0.0;
			 break;
		  case 9: //1D-Shock-and-Rarefaction
 			 res[ik]=vAmp[ik]*_solShockNWave(x[0], 0.0);
			 break;
 		 //default: block(ix,iy)(ik,0) = 0;


		  }
		 } 
		  
		}
		
		//Solution to the Shock and Rarefaction-wave. mapped to space (0..1)
		static BaseType _solShockNWave(BaseType x, BaseType t)
		{
			BaseType res;
			BaseType shock_loc=BaseType(1.0)+BaseType(0.5)*t;
			BaseType xdash=BaseType(4.0)*x-BaseType(1.0); //mapping to -1..3
			if ((xdash<BaseType(0.0)) || xdash>shock_loc) res=BaseType(0.0);
			else if (t<BaseType(2.0) && t<=xdash && xdash <= shock_loc) res=1.0;
			else res=xdash/t;
			return res;

		}
		
		
		//put some initial condition on the domain
		//also safe it to oldG, and initialG
		void _some_ic(Grid<Wavelets, Block>& grid,int icType=0)
		{
			
			BaseType x[sDim];
			
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			
			
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				
				Block& block = grid.getBlockCollection()[info.blockID];
				
				
				
				for(int iy=0; iy<Block::sizeY; iy++)
					for(int ix=0; ix<Block::sizeX; ix++)
					{
						
	
						BaseType vals[vDim];
						info.pos(x,ix,iy);
						_ic_func(vals, x, icType);
						 for (int ik=0;ik<vDim;++ik)
						 {
						   block(ix,iy)(ik,0)=vals[ik];
						 }
						block(ix,iy).update();
						block(ix,iy).setIC();
						
					}
				}
			
		}

		
				
		//Update Velocities (currently just a constant velocity (yeah I know, its a memory killer))
		void _updateVelocities(float* uconst)
		{
			
			std::cout << "updating Velocity field: ... " ;
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				
				Block& block = grid.getBlockCollection()[info.blockID];
				
				
				
				for(int iy=0; iy<Block::sizeY; iy++)
					for(int ix=0; ix<Block::sizeX; ix++)
					{
						
						for (int sk=0; sk<sDim;++sk)
						{
							block(ix,iy).velocity[sk] =  uconst[sk];
						}
						
					}
				
				
			}
			
			std::cout << "...done." <<std::endl;
		}
		
		
		//Update Velocities
		//Idea is to use this before calling upwinding scheme. (case 1, 2) 
		BaseType _updateVelocities(vector<BlockInfo> & vInfo, BaseType forceCFL=BaseType(0.5))
		{
			float constvel[2]={0.25, 0.5};
			BaseType maxvel=0;
			int steStart[]={-1,-1,0};
			int steEnd[]={2,2,1};
			
			BaseType minDeltaT=BaseType(HUGE);
			BaseType h[sDim];
			short int finest_level=-1;
			lab.prepare(grid.getBlockCollection(), *myBoundaries,steStart,steEnd);
			
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				
				finest_level=std::max(finest_level,info.level);
				h[0]=pow(2.,-finest_level);
				h[1]=pow(2.,-finest_level);
				lab.load(info);					
				
				
				for(int iy=0; iy<Block::sizeY; iy++)
					for(int ix=0; ix<Block::sizeX; ix++)
					{
						
						//velocity set to: value of 1st vector component at old iteration step (-1). 
					    lab(ix,iy).velocity[0]=lab(ix,iy)(0,-1);
						lab(ix,iy).velocity[1]=lab(ix,iy)(0,-1); 
						minDeltaT=std::min(minDeltaT,forceCFL*h[0]/fabs(lab(ix,iy).velocity[0]));
						minDeltaT=std::min(minDeltaT,forceCFL*h[1]/fabs(lab(ix,iy).velocity[1]));
						
						maxvel=std::max(maxvel,fabs(lab(ix,iy).velocity[0]));
						maxvel=std::max(maxvel,fabs(lab(ix,iy).velocity[1]));
						
						
					}
				
				lab.flush(info);
				
			}
			
			//debug:
			//std::cout << "Calling updateVelocities, while forcing CFL="<< forceCFL << " mindT " << minDeltaT <<std::endl;
			std::cout << "Calling updateVelocities, maxVeL="<< maxvel <<std::endl;
			
			return minDeltaT;
			
		}
		
		
		
		
		
		BaseType _advectwFluxes(vector<BlockInfo> & vInfo, BaseType forceCFL=BaseType(0.5), int fType=0, int icType=9)
		{
			
			BaseType h[sDim];
			BaseType hpref;
			BaseType minDeltaT=BaseType(HUGE);
			BaseType tiny=numeric_limits<BaseType>::epsilon();
			BaseType x[2];
			
			//NOTE:THE MAPPING PARAMETER IS ONLY for ictyp=9
			if (icType==9) hpref=.25;
			else hpref=1.0;
            

			
			//determine mindeltaT
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info2 = vInfo[i];
				Block& block = grid.getBlockCollection()[info2.blockID];
				
				for (int sd=0;sd<sDim;++sd)
				{
					h[sd]=hpref*pow(2.,-info2.level);
				}
				
				for(int iy=0; iy<Block::sizeY; iy++)
					for(int ix=0; ix<Block::sizeX; ix++)
					{
						
						//set velocities to first component of vector (can change this later)
						block(ix,iy).velocity[0]=block(ix,iy)(0,-1);
						block(ix,iy).velocity[1]=block(ix,iy)(0,-1);
						minDeltaT=std::min(minDeltaT,forceCFL*h[0]/fabs(block(ix,iy).velocity[0]+tiny));
						minDeltaT=std::min(minDeltaT,forceCFL*h[1]/fabs(block(ix,iy).velocity[1]+tiny));
						
					}
			}
           
		    std::cout << "mindeltaTa=" << minDeltaT << std::endl;
			
			
			BaseType preFac[sDim];
			BaseType Fm[sDim];
			BaseType Fr[sDim];
			BaseType Fl[sDim];
			BaseType rsum,lsum;

			
			lab.prepare(grid.getBlockCollection(), *myBoundaries, myBoundaries->getStencilStart(), myBoundaries->getStencilEnd());
            //BaseType mxdif(0),mxvel(0),mnvel(HUGE);
			//the Real Loop: calculate new values.
			
			
			//diagonal matrix, used for advection in respective direction.
			int sel[sDim][sDim];
			for (int d=0;d<sDim;++d)
			{
				for (int e=0;e<sDim;++e)
				{
				  e==d? sel[d][e] = 1: sel[d][e]=0;
				}
			}
		    
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				lab.load(info);
				
			    for (int sd=0; sd<sDim; ++sd)
			    {
				   h[sd]=hpref*pow(2.0,-info.level);
				   preFac[sd]=minDeltaT/(2.0*h[sd]);
			    }
				
				
				for(int iy=0; iy<Block::sizeY; iy++)
					for(int ix=0; ix<Block::sizeX; ix++)
					{
	                        
							for (int ik=0;ik<vDim;++ik) //loop over vector components
							{
								rsum=0.0;
								lsum=0.0;
								for (int d=0; d<sDim; ++d)
								{
									bool r_gt_m = (lab(ix,iy).velocity[d] > lab(ix+sel[d][0],iy+sel[d][1]).velocity[d]) ;
									bool l_gt_m = (lab(ix-sel[d][0],iy-sel[d][1]).velocity[d] > lab(ix,iy).velocity[d]);
									Fl[d]=pow(lab(ix-sel[d][0],iy-sel[d][1]).velocity[d],2.0); //currently burgers-flux is hardcoded.
									Fm[d]=pow(lab(ix,iy).velocity[d],2.0); 
									Fr[d]=pow(lab(ix+sel[d][0],iy+sel[d][1]).velocity[d],2.0); 
									BaseType maxfl = l_gt_m? Fl[d]:Fm[d];
									BaseType minfl = Fl[d]+Fm[d]-maxfl;
							
									BaseType maxfr = r_gt_m? Fm[d]:Fr[d];
									BaseType minfr = Fm[d]+Fr[d]-maxfr;
									//a_gt_b? lab(ix,iy)(ik,0)=maxf:lab(ix,iy)(ik,0)=minf;
									r_gt_m? rsum+=maxfr:rsum+=minfr;
									l_gt_m? lsum+=maxfl:lsum+=minfl;
								 }
								 
								 
								 lab(ix,iy)(ik,0)=lab(ix,iy)(ik,-1)-preFac[0]*(rsum-lsum);

																	

							  }

															
					}
				lab.flush(info);
			}
			 
			//std::cout << "mxdif was: " << mxdif << "mxvel was: " << mxvel << " minvel was:" << mnvel << std::endl; 
			std::cout << "preFac was" << preFac[0] <<std::endl;

			//update	
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info3 = vInfo[i];
				Block& block = grid.getBlockCollection()[info3.blockID];
				for(int iy=0; iy<Block::sizeY; iy++)
					for(int ix=0; ix<Block::sizeX; ix++)
					{
						block(ix,iy).update();
					}
			}
			
		
		return minDeltaT;
			
		}
		
		
		void _some_advection_1stOrder(BaseType dt) //upwinding advection
		{
			int sel[sDim];
			BaseType h[sDim];
			BaseType dud[sDim];
			
			
			for (int s=0;s<sDim;++s)
			{
				sel[s]=0;
			}

			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			
			
			lab.prepare(grid.getBlockCollection(), *myBoundaries, myBoundaries->getStencilStart(), myBoundaries->getStencilEnd());
			//lab.prepare(grid.getBlockCollection(), grid.getBoundaryInfo(),steStart,steEnd);
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				lab.load(info);
				
				for (int sd=0;sd<sDim;++sd)
				{
					h[sd]=pow(2.,-info.level);
				}
				
				for(int iy=0; iy<Block::sizeY; iy++)
					for(int ix=0; ix<Block::sizeX; ix++)
					{
						
						//Stencil-Selector for upwinding:
						for (int sd=0;sd<sDim;++sd)
						{
							sel[sd]= lab(ix,iy).velocity[sd]>0? -1 : 1;
						}
						
						for (int ik=0;ik<vDim;++ik) //loop over vector components
						{
							
							dud[0]=1.0/h[0]*(lab(ix,iy)(ik,-1)-lab(ix+sel[0],iy)(ik,-1));
							dud[1]=1.0/h[1]*(lab(ix,iy)(ik,-1)-lab(ix,iy+sel[1])(ik,-1));
							//dud[2]=1.0/h[2]*(lab(ix,iy)(ik,-1)-lab(ix,iy+sel[2])(ik,-1));
							lab(ix,iy)(ik) = lab(ix,iy)(ik,-1)-dt*(fabs(lab(ix,iy).velocity[0])*dud[0]+fabs(lab(ix,iy).velocity[1])*dud[1]);
							
							
							// explanation:
							//lab(ix,iy)(ik) : Value of component ik, at current time
							//lab(ix,iy)(ik,-1) : Value of component ik, from previous timestep
							//lab(ix,iy).velocity[0] : X-compenent of local velocity field.
							//sel[0]: selector for upwinding scheme, based on local velocity field.			
							
							//does not work (would require some auxillary-arrays, so hardcoding is better) 
							
							//for (int sd=0; sd<sDim;++sd) //loop over space directions.
							// {
							//  lab(ix,iy)(ik)-=dt*
							//  fabs(lab(ix,iy).velocity[sd])/h[sd]*(lab(ix+sel[sd][0],iy+sel[sd][1])(ik,-1) - lab(ix+sel[sd][0]-1,iy+seld[sd][1]-1)(ik,-1));
							// }
							
							//2nd order, FD Upwinding
							//dudx=1.0/(2.0*h[0])*(3.0*lab(ix,iy).old-4.0*lab(ix+selx,iy).old+lab(ix-(selx*2),iy).old);
							//dudy=1.0/(2.0*h[1])*(3.0*lab(ix,iy).old-4.0*lab(ix,iy+sely).old+lab(ix,iy-(sely*2)).old);
							//lab(ix,iy).cur=lab(ix,iy).old-dt*(fabs(ux)*dudx+fabs(uy)*dudy);	 
							
							
						}
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
					}
			}
		}
		
		
		void _some_advection_2ndOrder(BaseType dt) //upwinding advection
		{
			
			int sel[sDim];
			BaseType h[sDim];
			BaseType dud[sDim];
			
			for (int s=0;s<sDim;++s)
			{
				sel[s]=0;
			}
			
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			std::cout << "[2ndOrderAdvection] DBG:" <<std::endl;
			std::cout <<  " s: " << int(myBoundaries->stencil_start[0]) << " " << int(myBoundaries->stencil_start[1]) <<  " " << int(myBoundaries->stencil_start[2])<< " " <<std::endl;
			std::cout <<  " e: "  <<  int(myBoundaries->stencil_end[0]) <<" " << int(myBoundaries->stencil_end[1]) <<  " " << int(myBoundaries->stencil_end[2])<< " " <<std::endl;
			lab.prepare(grid.getBlockCollection(), *myBoundaries,myBoundaries->getStencilStart(),myBoundaries->getStencilEnd());
			//lab.prepare(grid.getBlockCollection(), grid.getBoundaryInfo(),steStart,steEnd);
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				lab.load(info);
				
				for (int sd=0;sd<sDim;++sd)
				{
					h[sd]=pow(2.,-info.level);
				}
				
				for(int iy=0; iy<Block::sizeY; iy++)
					for(int ix=0; ix<Block::sizeX; ix++)
					{
						//Stencil-Selector for upwinding:
						for (int sd=0;sd<sDim;++sd)
						{
							sel[sd]= lab(ix,iy).velocity[sd]>0? -1 : 1;
						}
						
						for (int ik=0;ik<vDim;++ik) //loop over vector components
						{
							
							//2nd Order Upwinding Stencils:
							dud[0]=1.0/(2.0*h[0])*(3.0*lab(ix,iy)(ik,-1)-4.0*lab(ix+sel[0],iy)(ik,-1)+lab(ix-(sel[0]*2),iy)(ik,-1));
							dud[1]=1.0/(2.0*h[1])*(3.0*lab(ix,iy)(ik,-1)-4.0*lab(ix,iy+sel[1])(ik,-1)+lab(ix,iy-(sel[1]*2))(ik,-1));
							//dud[2]=1.0/(2.0*h[2])*(3.0*lab(ix,iy)(ik,-1)-4.0*lab(ix,iy+sel[2])(ik,-1)+lab(ix,iy-(sel[2]*2))(ik,-1));
							//calculation based on interim results. (the same for all (OrderN) stencils)
							lab(ix,iy)(ik) = lab(ix,iy)(ik,-1)-dt*(fabs(lab(ix,iy).velocity[0])*dud[0]+fabs(lab(ix,iy).velocity[1])*dud[1]);
							
							// explanation:
							//lab(ix,iy)(ik) : Value of component ik, at current time
							//lab(ix,iy)(ik,-1) : Value of component ik, from previous timestep
							//lab(ix,iy).velocity[0] : X-compenent of local velocity field.
							//sel[0]: selector for upwinding scheme, based on local velocity field.			
							
						}
					}
				lab.flush(info);
			} //end integrate-loop
			
			//update "oldG" (dont need a lab here, but direct access to BlockCollection)
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info2 = vInfo[i];
				Block& block = grid.getBlockCollection()[info2.blockID];
				for(int iy=0; iy<Block::sizeY; iy++)
					for(int ix=0; ix<Block::sizeX; ix++)
					{
						block(ix,iy).update();
					}
				
			} //end update-loop
			
		}	//end function	 _some_advection_2ndOrde	
		
	
		//error calculation based on all components of vector (i.e. RMS)
		
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
						
						
						abs_error=0.0;
						for (int ik=0;ik<vDim;++ik)
						{
							abs_error   += pow(fabs(block(ix,iy)(ik) - block(ix,iy).ic[ik]),2.0);
						} 
						abs_error =  sqrt(abs_error);
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
		
		
		BaseType getError(BaseType t, int solCase=0, double avgNumElems)
		{
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			BaseType vAmp[3]={1.0,1.0,1.0};
			BaseType Linfty = BaseType(-1.0);
			BaseType Ltwo(0);
			BaseType sumAbs(0);
			BaseType abs_error;
			BaseType exact[vDim];
			BaseType x[sDim];
			BaseType b1,b2,bs;
			b1=1.0/3.0-t;
			b2=1.0/3.0+t;
			bs=2.0/3.0;
			short int maxLevel=-1;
			
			
			std::ofstream exoutf;
			
			if (t==BaseType(0.0)) exoutf.open("exact.dat");
			if (t< BaseType(0.6) && t>BaseType(0.4)) exoutf.open("exact-t0.5.dat");
			
			
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				Block& block = grid.getBlockCollection()[info.blockID];
				maxLevel=std::max(maxLevel,info.level);
				
				
				for(int iy=0; iy<Block::sizeY; iy++)
				{
					for(int ix=0; ix<Block::sizeX; ix++)
					{
						
						info.pos(x,ix,iy);
						
						


						
						
						abs_error=0.0;
						for (int ik=0;ik<vDim;++ik)
						{
						
							for (int d=0; d < sDim; ++d)
							{
							switch(solCase)
								{
								  case 5:
								  {
								  if ((x[d] <= b1) || (x[d] >= bs)) exact[ik] = -1.0 ;
								  else if ( (x[d] > b1 ) && (x[d] <= b2)) exact[ik]=-1.0+2.0*(x[ik]-b1)/(b2-b1) ;
								  else if ( (x[d] > b1) && (x[d] < bs) ) exact[ik]=1.0;
								  else exact[ik]=5.0;
								  break;
								  }
								  case 7:
								  {
								  if ((x[d] <= b1) || (x[d] >= bs)) exact[ik] = -1.0 ;
								  else if ( (x[d] > b1 ) && (x[d] <= b2)) exact[ik]=-1.0+2.0*(x[ik]-b1)/(b2-b1) ;
								  else if ( (x[d] > b1) && (x[d] < bs) ) exact[ik]=1.0;
								  else exact[ik]=5.0;
								  break;
								  }
								  case 8:
								  {
								  int ds=0;
								  BaseType shock_loc=3.0/4.0+0.5*t;
								  if ((x[ds]<1.0/4.0) || x[ds]>shock_loc) exact[ik]=0.0;
								  else if (t<2.0 && t<(3.0/4.0)*x[ds] && x[ds] <= shock_loc) exact[ik]=1.0;
								  else exact[ik]=4.0*(x[ds]-(1.0/4.0))/t;
								  break;
								  }
								  case 9:
								  {
								  /*
								  BaseType shock_loc=BaseType(1.0)+BaseType(0.5)*t;
								  BaseType xdash=BaseType(4.0)*x[0]-BaseType(1.0); //mapping to -1..3
								  if ((xdash<BaseType(0.0)) || xdash>shock_loc) exact[ik]=BaseType(0.0);
								  else if (t<BaseType(2.0) && t<=xdash && xdash <= shock_loc) exact[ik]=1.0;
								  else exact[ik]=xdash/t;
								  */
								  exact[ik]=vAmp[ik]*_solShockNWave(x[0], t);
								  break;
								  }

								  default:
								  {
								  exact[ik]=block(ix,iy).ic[ik];
								  }
								 }
							}							
							
							abs_error   += pow((block(ix,iy)(ik) - exact[ik]),2.0);
					
						} 
						
						if (t==BaseType(0.0) || (t< BaseType(0.6) && t>BaseType(0.4)))
						{
						  exoutf << x[0] << " " << x[1] << " " << exact[0] << " " << block(ix,iy)(0,0) <<std::endl;
						}
						abs_error =  sqrt(abs_error);
						Ltwo+=1.0/pow(4.0,info.level)*abs_error; //(1/h^2)*Eloc
						sumAbs+=abs_error;
						if(abs_error > Linfty) Linfty = abs_error;
						if (abs_error > 0.0 && t==BaseType(0.0)) std::cout << "err. = " << abs_error << " at: " << x[0] << " " << x[1] <<std::endl;
						
					}
				}
				
			}
							
			std::cout << "hx, Linfty, Ltwo, abs/avgNumPoints = " << 1.0/(pow(2.0,maxLevel)) << " " << Linfty << " " << Ltwo << " " << sumAbs/avgNumElems << std::endl;
			exoutf.close();
			return Linfty;
		}
		
		//updated:
		void _dump2file(std::string fname) {
			
			std::ofstream outf(fname.c_str());
			std::cout << "requesting vInfo" <<std::endl;
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			std::cout << "starting loop over vInfo: " << vInfo.size() <<std::endl;
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				Block&block = grid.getBlockCollection()[info.blockID];
				
				
				float x[2];
				
				for(int iy=0; iy<Block::sizeY; iy++)
				{
					for(int ix=0; ix<Block::sizeX; ix++)
					{
						
						info.pos(x,ix,iy);
						outf << x[0] << " " << x[1]<< " " ;
						
						for (int ik=0; ik<vDim; ++ik)
						{			
							outf << block(ix,iy)(ik) << " " ;
						}
						
						for (int sd=0; sd<sDim; ++sd)
						{			
							outf << block(ix,iy).velocity[sd] << " " ;
						}
						outf <<  std::endl ;
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
				float uconst[2]={1.0,1.0};
				int steStart[3] ={-2,-2,0};
				int steEnd[3] = {3,3,1};
				double CFL=.5;
				double dt;
				float err;
				
				printf("Start Test\n");
				for (int i=0; i < LevelOfRefinements; ++i)
				{
					_some_refinements();
				}
				
				//find finest level:
				short int finest_level=-1;
				vector<BlockInfo> vInfo = grid.getBlocksInfo();
				for (int i=0;i<vInfo.size();++i)
				{
					finest_level=std::max(finest_level,vInfo[i].level);
				}
				
				for (int P=2; P>=1; --P)
				{
					std::cout << "Running test for order: "<< P <<std::endl; 
					
					//Determine timestep, based on CFL on finest grid
					const double h[2]= {pow(2.,-finest_level), pow(2.,-finest_level)};
					uconst[0]=h[0];
					uconst[1]=h[1];
					dt=CFL*h[0]/uconst[0];
					dt=std::max(dt,CFL*h[1]/uconst[1]);
					
					std::cout <<" CFL= "<<CFL << std::endl;
					std::cout <<" dt = " <<dt <<std::endl;
					
					//get correct stencil-width.
					switch (P)
					{
						case 1:
							OrderFD=1;
							steStart[0] =-1;
							steStart[1]=-1;
							steEnd[0] = 2;
							steEnd[1] =2 ;
							break;
						case 2:
							OrderFD=2;
							steStart[0] =-2;
							steStart[1]=-2;
							steEnd[0] = 3;
							steEnd[1] = 3 ;
							break;
						default:
							exit(1);
							break;
					}
					
					myBoundaries->init(grid,steStart,steEnd);
					std::cout << "orderFD is: " << OrderFD <<std::endl;
					
					_some_ic(grid);
					_dump2file("iniC.dat");
					_updateVelocities(uconst);
					
					for (int i=0; i<NumberOfSteps;++i)
					{
						switch (OrderFD)
						{
							case 1:  {_some_advection_1stOrder(dt); break;}
							case 2:  {_some_advection_2ndOrder(dt); break;}
							default: {_some_advection_1stOrder(dt); break;}
						}
					}
					_dump2file("t00.dat");
					
					uconst[0]=-uconst[0];
					uconst[1]=-uconst[1];
					_updateVelocities(uconst);
					
					for (int i=0; i<NumberOfSteps;++i)
					{
						switch (OrderFD)
						{
							case 1:  {_some_advection_1stOrder(dt); break;}
							case 2:  {_some_advection_2ndOrder(dt); break;}
							default: {_some_advection_1stOrder(dt); break;}
						}
					}
					
					_dump2file("t01.dat");
					
					err=_some_error_calculation();
					std::cout << "Error is: " << err << std::endl;
					std::cout << "CFL was: " << CFL << std::endl;
					std::cout <<" dt = " <<dt <<std::endl;
					std::cout << "hx,err: " << h[0] << " " << err <<std::endl;
					
					printf("End Test, with Order %i \n",P);
				}
			}
			catch(...)
			{
				return false;
			}
			
			return true;
		}
		
		
		
		//burgers2D-demo.
		//use ictype to choose IC, (recommended:9)
		//output: *.dat files in running directory. plot with gnuplot -> splot
		//use adapt to change Refinement-Frequency. 
		void burgers2D(int orderFD=4, bool WithSTR=false, float tmax=1.0, int maxLevel = 3)
		{
			std::stringstream outfname;
			vector<MRAG::BlockInfo> vInfo;
			int dumpfreq=5;
			int adapt=1e9; //never
			int iter(0),maxiter(20000);
			const BaseType toleranceIC=0.001;
			BaseType tcur(0.0),dt;
			BaseType forceCFL=0.5;
			int steStart[3] ={-2,-2,0};
			int steEnd[3] = {3,3,1};
			long int sumNumElements(0);
			double avgNumElements(0);
			double curAvg(0);
			int szX,szY,szZ;
			
			
	
			
			int ictype=9; //NOTE: if you use ictype 9: it is going to use a mapping from -1...3 to 0...1 to be compatible with an exact solution defined on that domain.
			
			
			std::cout << "putting initial Condition: " <<std::endl;
			_some_ic(grid,ictype);
			if (adapt<1e6)
			{
			Science::AutomaticRefinement< 0,0 >(grid, blockfwt, toleranceIC/2,maxLevel);
			std::cout << "refined initial Condition: " <<std::endl;
			_some_ic(grid,ictype);
			Science::AutomaticCompression< 0,0 >(grid, blockfwt, toleranceIC);
			_some_ic(grid,ictype);
			std::cout << "compressed initial Condition: " <<std::endl;
			}
			
			myBoundaries->init(grid,steStart,steEnd);
			myBoundaries->update(grid);
			
			vInfo=grid.getBlocksInfo();
			std::cout << "trying to dump file" <<std::endl;
			_dump2file("b2DiC.dat");
			
			std::cout << "myBoundaries->stest[0] " << int(myBoundaries->stencil_start[0]) << " end[0]: "<< int(myBoundaries->stencil_end[0]) <<std::endl;
			std::cout << "myBoundaries->_smarttest[0] " << myBoundaries->getStencilStart()[0] << " end[0] "<<myBoundaries->getStencilEnd()[0] << std::endl;
			
			
			szX=Block::sizeX;
			szY=Block::sizeY;
			szZ=Block::sizeZ;
			curAvg=vInfo.size()*szX*szY*szZ;
			std::cout << "Error at t= 0 (IC)"  << getError(tcur, ictype,curAvg) <<std::endl;
			
			
			iter=0;
			while (tcur<tmax && iter<maxiter)
			{ 
				vInfo = grid.getBlocksInfo();
				
				
				if (WithSTR) //SpaceTimeRefinement
				{ 
					//to be implemented 
				}
				
				else
				{
					
					if (iter>0 && iter%adapt==0) 
					{ 
						Science::AutomaticRefinement< 0,0 >(grid, blockfwt, toleranceIC/2,maxLevel);
						//bInfo=&grid.getBoundaryInfo(steStart,steEnd);
						myBoundaries->update(grid);
						vInfo.clear();
						vInfo=grid.getBlocksInfo();
						std::cout << "updated vInfo, size is: " << vInfo.size() <<std::endl;
					}
					
					//size after refinement
					curAvg=vInfo.size()*szX*szY*szZ;
															
					switch(orderFD)
					{
	
						
						case 1: 
						{dt=_updateVelocities(vInfo, forceCFL);
						 std::cout << "t= " << tcur << " iter= " << iter << " moving with timestep: " << dt << std::endl;
						 _some_advection_1stOrder(dt);
						 break;
						}
						case 2: 
						{ dt=_updateVelocities(vInfo, forceCFL);
						 std::cout << "t= " << tcur << " iter= " << iter << " moving with timestep: " << dt << std::endl;
						 _some_advection_2ndOrder(dt);
						break;
						
						case 4: //wFluxes
						dt=_advectwFluxes(vInfo,forceCFL,0,ictype);
						std::cout << "t= " << tcur << " iter= " << iter << " moving with timestep: " << dt << std::endl;
						break;
						}
					}
				}
				
				
				if (iter%dumpfreq==0) 
				{
					outfname.str("");
					outfname << "b2DI" << iter << ".dat";
					cout << "using filename: " << outfname.str();
					_dump2file(outfname.str().c_str());
				}
				
				if (iter>0 && iter%adapt==0) 
				{
					Science::AutomaticCompression< 0,0 >(grid, blockfwt, toleranceIC);
					myBoundaries->update(grid);
					vInfo.clear();
					vInfo=grid.getBlocksInfo();
				}
				
				curAvg+=vInfo.size()*szX*szY*szZ;
				curAvg=curAvg/2.0;
				sumNumElements+=curAvg;
				avgNumElements=double(sumNumElements)/double(iter+1);
				std::cout << "avgNumPoints: " << avgNumElements <<std::endl;

				tcur+=dt;
				std::cout << "Error at t= " << tcur << " " << getError(tcur, ictype,avgNumElements) <<std::endl;

				iter++;
				
			}
			
			//finishing up:
			
		} //end burgers2D
		

	public:
		static void runTests(int LevelOfRefinements=2, int NumOfBlocks=2, int OrderFD=2, int NumberOfSteps=1)
		{
			MRAG_STDTestL3_VectorAdvection<Wavelets, Block> test(NumOfBlocks);
			
			test.run(LevelOfRefinements,NumOfBlocks,OrderFD,NumberOfSteps);
			//usage example: 
			//	MRAG::MRAG_STDTestL3_VectorAdvection<MRAG::Wavelets_Interp2ndOrder, MRAG::Block<misc::tAdVector<double,3,2,2>,16,16,1> >::runTests(2,4,1,1);
			
			
		}
		
		
		//usage:
		//	MRAG::MRAG_STDTestL3_VectorAdvection<MRAG::Wavelets_Interp2ndOrder, MRAG::Block<misc::tAdVector<double,3,2,2>,16,16,1> >::runBurgers2D(8,2);
		static void runBurgers2D(int NumOfBlocks=2, int maxRefinementLevel)
		{
			

			MRAG_STDTestL3_VectorAdvection<Wavelets, Block> test(NumOfBlocks);
			//void burgers2D(int orderFD=4, bool WithSTR=false, float tmax=1.0, int maxLevel = 3)
			test.burgers2D(4,false,1.0,maxRefinementLevel);
		}
		
		
	};
	
	
}





