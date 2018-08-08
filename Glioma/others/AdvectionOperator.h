/*
 *  AdvectionOperator.h
 *  GliomaXcode
 *
 *  Created by Lipkova on 9/30/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#pragma once
#include "Glioma_Types.h"

using namespace MRAG;
using namespace Science;

template<int nDim = 3>
struct GetUMax
{
	map< int, Real>& local_max_velocities;
	
	GetUMax(map< int, Real>& local_max_velocities):local_max_velocities(local_max_velocities) 
	{ }
	
	GetUMax(const GetUMax& c): local_max_velocities(c.local_max_velocities)
	{ }
	
	template<typename BlockType>	
	inline void operator() (const BlockInfo& info, BlockType& b) const
	{
		Real maxVel[3] = {0,0,0};
		if(nDim == 2)
        {
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        maxVel[0] = max((Real)fabs( b(ix,iy).ux), maxVel[0]);
                        maxVel[1] = max((Real)fabs( b(ix,iy).uy), maxVel[1]);
                        
                    }
        }
        else
        {
		for(int iz=0; iz<BlockType::sizeZ; iz++)
			for(int iy=0; iy<BlockType::sizeY; iy++)
				for(int ix=0; ix<BlockType::sizeX; ix++)
				{
					maxVel[0] = max((Real)fabs( b(ix,iy,iz).ux), maxVel[0]);
					maxVel[1] = max((Real)fabs( b(ix,iy,iz).uy), maxVel[1]);
					maxVel[2] = max((Real)fabs( b(ix,iy,iz).uz), maxVel[2]);
					
				}
		
        }
        
		map< int, Real>::iterator it = local_max_velocities.find(info.blockID);
		assert(it != local_max_velocities.end());
		
		it->second = max(maxVel[0], max(maxVel[1], maxVel[2]));
        
	}
};


struct Upwind1rdOrder
{
	int stencil_start[3];
	int stencil_end[3];
	double dt;
	
	Upwind1rdOrder(double dt_)
	{
		stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
		stencil_end[0] = stencil_end[1] = stencil_end[2] = +2;
		dt		= dt_;
	}
	
	Upwind1rdOrder(const Upwind1rdOrder&, double dt_ )
	{
		stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
		stencil_end[0] = stencil_end[1] = stencil_end[2] = +2;
		dt		= dt_;
	}
	
	template<typename LabType, typename BlockType>
	inline void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
	{
		const Real factor = -1./(info.h[0]);
		
		
		for(int iz=0; iz<BlockType::sizeZ; iz++)
			for(int iy=0; iy<BlockType::sizeY; iy++)
				for(int ix=0; ix<BlockType::sizeX; ix++)
				{
					
					const Real dwdx[2] = {
						lab(ix,iy,iz).omega - lab(ix-1,iy,iz).omega,
					    lab(ix+1,iy,iz).omega - lab(ix,iy,iz).omega};
					
					const Real dwdy[2] = {
						lab(ix,iy,iz).omega - lab(ix,iy-1,iz).omega,
					    lab(ix,iy+1,iz).omega - lab(ix,iy,iz).omega};
					
					const Real dwdz[2] = {
						lab(ix,iy,iz).omega - lab(ix,iy,iz-1).omega,
					    lab(ix,iy,iz+1).omega - lab(ix,iy,iz).omega};
					
					const Real u[3] = {
						lab(ix, iy, iz).ux,
						lab(ix, iy, iz).uy,
						lab(ix, iy, iz).uz};
					
					o(ix, iy, iz).domegadt = factor*(
													 max(u[0], (Real)0.)*dwdx[0] + min(u[0], (Real)0.)*dwdx[1] + 
													 max(u[1], (Real)0.)*dwdy[0] + min(u[1], (Real)0.)*dwdy[1] + 
													 max(u[2], (Real)0.)*dwdz[0] + min(u[2], (Real)0.)*dwdz[1] ); 
					
					if (isnan(o(ix,iy,iz).domegadt))
					{
						printf("Nan in RHSUpwind1rdOrder operator, at block(%g, %g, %g ) \n", ix, iy, iz);
						printf("facotr=%g, dwdx[0]=%g, dwdx[1]=%g, dwdy[0]=%g, dwdy[1]=%g, dwdz[0]=%g, dwdz[1]=%g, \n", dwdx[0],dwdx[1],dwdy[0],dwdy[1],dwdz[0],dwdz[1]);
						abort();
					}
					
					
				}
	}	
};


struct Upwind3rdOrder
{
	int stencil_start[3];
	int stencil_end[3];
	double dt;
	
	Upwind3rdOrder(double dt_)
	{
		stencil_start[0] = stencil_start[1] = stencil_start[2] = -2;
		stencil_end[0] = stencil_end[1] = stencil_end[2] = +3;
		dt		= dt_;
	}
	
	Upwind3rdOrder(const Upwind3rdOrder&, double dt_ )
	{
		stencil_start[0] = stencil_start[1] = stencil_start[2] = -2;
		stencil_end[0] = stencil_end[1] = stencil_end[2] = +3;
		dt		= dt_;
	}
	
	template<typename LabType, typename BlockType>
	inline void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
	{
		const Real factor = -1./(6*info.h[0]);
				
		for(int iz=0; iz<BlockType::sizeZ; iz++)
			for(int iy=0; iy<BlockType::sizeY; iy++)
				for(int ix=0; ix<BlockType::sizeX; ix++)
				{

					const Real dwdx[2] = {
						lab(ix-2, iy, iz).omega - 6.*lab(ix-1, iy, iz).omega + 3.*lab(ix, iy, iz).omega + 2.*lab(ix+1,iy,iz).omega,
						-2.*lab(ix-1, iy, iz).omega - 3.*lab(ix, iy, iz).omega + 6.*lab(ix+1, iy, iz).omega - lab(ix+2, iy, iz).omega};
					

					const Real dwdy[2] = {
						lab(ix, iy-2, iz).omega - 6.*lab(ix, iy-1, iz).omega + 3.*lab(ix, iy, iz).omega + 2.*lab(ix,iy+1,iz).omega,
						-2.*lab(ix, iy-1, iz).omega - 3.*lab(ix, iy, iz).omega + 6.*lab(ix, iy+1, iz).omega - lab(ix, iy+2, iz).omega};
					
					const Real dwdz[2] = {
						lab(ix, iy, iz-2).omega - 6.*lab(ix, iy, iz-1).omega + 3.*lab(ix, iy, iz).omega + 2.*lab(ix,iy,iz+1).omega,
						-2.*lab(ix, iy, iz-1).omega - 3.*lab(ix, iy, iz).omega + 6.*lab(ix, iy, iz+1).omega - lab(ix, iy, iz+2).omega};
					

					const Real u[3] = {
						o(ix, iy, iz).ux,
						o(ix, iy, iz).uy,
						o(ix, iy, iz).uz};
					
					o(ix, iy, iz).domegadt = factor*(
													   max(u[0], (Real)0)*dwdx[0] + min(u[0], (Real)0)*dwdx[1] + 
													   max(u[1], (Real)0)*dwdy[0] + min(u[1], (Real)0)*dwdy[1] + 
													   max(u[2], (Real)0)*dwdz[0] + min(u[2], (Real)0)*dwdz[1] );
					
					if (isnan(o(ix,iy,iz).domegadt))
					{
						printf("Nan in RHSUpwind1rdOrder operator, at block(%g, %g, %g ) \n", ix, iy, iz);
						printf("facotr=%g, dwdx[0]=%g, dwdx[1]=%g, dwdy[0]=%g, dwdy[1]=%g, dwdz[0]=%g, dwdz[1]=%g, \n", dwdx[0],dwdx[1],dwdy[0],dwdy[1],dwdz[0],dwdz[1]);
						abort();
					}
				}
	}	
};

struct AdvectWeno5
{
	int stencil_start[3];
	int stencil_end[3];
	double dt;
	
	AdvectWeno5(double dt_)
	{
		stencil_start[0]=stencil_start[1]=stencil_start[2]= -3;  
		stencil_end[0]  =stencil_end[1]  =stencil_end[2]  = +4;
		dt		= dt_;
	}
	
	AdvectWeno5(const AdvectWeno5&, double dt_ )
	{
		
		stencil_start[0]=stencil_start[1]=stencil_start[2]= -3;
		stencil_end[0]  =stencil_end[1]  =stencil_end[2]  = +4;
		dt		= dt_;
	}
	
	template<typename LabType, typename BlockType>
	inline void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
	{
		Weno5<double> weno5;
		
		const double factor = 1./(info.h[0]);
		
		for(int iz=0; iz<BlockType::sizeZ; iz++)
			for(int iy=0; iy<BlockType::sizeY; iy++)
				for(int ix=0; ix<BlockType::sizeX; ix++)
				{
					
					// fluxes
					double omega_flux_x[2];   // omega_x- , omega_x+ 
					double omega_flux_y[2];
					double omega_flux_z[2];
					
					// in X
					weno5( factor * ( lab(ix-2,iy,iz).omega - lab(ix-3,iy,iz).omega ),
						  factor * ( lab(ix-1,iy,iz).omega - lab(ix-2,iy,iz).omega ),
						  factor * ( lab(ix  ,iy,iz).omega - lab(ix-1,iy,iz).omega ),
						  factor * ( lab(ix+1,iy,iz).omega - lab(ix  ,iy,iz).omega ),
						  factor * ( lab(ix+2,iy,iz).omega - lab(ix+1,iy,iz).omega ),
						  factor * ( lab(ix+3,iy,iz).omega - lab(ix+2,iy,iz).omega ),
						  omega_flux_x );
					
					// in Y
					weno5( factor * ( lab(ix,iy-2,iz).omega - lab(ix,iy-3,iz).omega ),
						  factor * ( lab(ix,iy-1,iz).omega - lab(ix,iy-2,iz).omega ),
						  factor * ( lab(ix,iy  ,iz).omega - lab(ix,iy-1,iz).omega ),
						  factor * ( lab(ix,iy+1,iz).omega - lab(ix,iy  ,iz).omega ),
						  factor * ( lab(ix,iy+2,iz).omega - lab(ix,iy+1,iz).omega ),
						  factor * ( lab(ix,iy+3,iz).omega - lab(ix,iy+2,iz).omega ),
						  omega_flux_y );
					
					
					// in Z
					weno5( factor * ( lab(ix,iy,iz-2).omega - lab(ix,iy,iz-3).omega ),
						  factor * ( lab(ix,iy,iz-1).omega - lab(ix,iy,iz-2).omega ),
						  factor * ( lab(ix,iy,iz  ).omega - lab(ix,iy,iz-1).omega ),
						  factor * ( lab(ix,iy,iz+1).omega - lab(ix,iy,iz  ).omega ),
						  factor * ( lab(ix,iy,iz+2).omega - lab(ix,iy,iz+1).omega ),
						  factor * ( lab(ix,iy,iz+3).omega - lab(ix,iy,iz+2).omega ),
						  omega_flux_z );
					
					// velocities
					Real velocity_plus[3] = {
						max(o(ix,iy,iz).ux, (Real) 0.),
						max(o(ix,iy,iz).uy, (Real) 0.),
						max(o(ix,iy,iz).uz, (Real) 0.)
					};
					
					Real velocity_minus[3] = {
						min(lab(ix,iy,iz).ux, (Real) 0.),
						min(lab(ix,iy,iz).uy, (Real) 0.),
						min(lab(ix,iy,iz).uz, (Real) 0.)
					};
					
					o(ix,iy,iz).domegadt = (-1.) * (
													velocity_plus[0] * omega_flux_x[0] + velocity_minus[0] * omega_flux_x[1] +
													velocity_plus[1] * omega_flux_y[0] + velocity_minus[1] * omega_flux_y[1] +
													velocity_plus[2] * omega_flux_z[0] + velocity_minus[2] * omega_flux_z[1]  );
					
                    
					if (isnan(o(ix,iy,iz).domegadt))
					{
						printf("Nan in Weno scheme \n");
						abort();
					}
					
					
				}
	}
};


