/*
 *  MRAG_STDTestL4_LevelsetAdvection.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 11/19/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGEnvironment.h"


#if _MRAG_OS == _MRAG_OS_APPLE
#ifdef _MRAG_GLUT_VIZ
#include "GLUT/glut.h"
#endif
#elif _MRAG_OS == _MRAG_OS_WINDOWS
#define _USE_MATH_DEFINES
#ifdef _MRAG_GLUT_VIZ
#include "GL/glew.h"
#include "GL/glut.h"
#endif
#endif

#undef min
#undef max

#include "MRAGcore/MRAGrid.h"
#include "MRAGcore/MRAGRefiner.h"
#include "MRAGcore/MRAGCompressor.h"
#include "MRAGcore/MRAGBlockLab.h"
#include "MRAGcore/MRAGBlockFWT.h"

#ifdef _MRAG_GLUT_VIZ
#include "MRAGvisual/GridViewer.h"
#endif

#include "MRAGscience/MRAGScienceCore.h"
#include "MRAGscience/MRAGAutomaticRefiner.h"
#include "MRAGscience/MRAGSimpleLevelsetBlock.h"
#include "MRAGscience/MRAGSpaceTimeSorter.h"
#include "MRAGscience/MRAGRefiner_SpaceExtension.h"

#include "MRAGmultithreading/MRAGBlockProcessing_SingleCPU.h"
#include "MRAGmultithreading/MRAGBlockProcessing_TBB.h"

#include "MRAGio/MRAG_IO_Native.h"
#include "MRAGio/MRAG_IO_VTK.h"

using namespace MRAG;

const int nDim = 2;

#pragma mark -
#pragma mark Element Type
struct PLS
{
	float rho,rho_0;
	float drho_dt;
	float tmp;
	PLS(): rho(0), drho_dt(0), tmp(0), rho_0(0) {}
	
	PLS(float rho_, float drho_dt_): rho(rho_), drho_dt(drho_dt_),tmp(0), rho_0(0) {}
	PLS(float rho_): rho(rho_), drho_dt(0),tmp(0), rho_0(0) {}
	
	Real levelset() const{
		return rho;
	}
	
	void operator += (PLS t) 
	{
		rho_0 +=t.rho_0;
		rho += t.rho;
		drho_dt += t.drho_dt;
	}
	
	operator Real() 
	
	{
		return (Real)rho;
	}
	
	void integrate(float dt)
	{
		rho += drho_dt*dt;
		drho_dt = 0;
		tmp = 0;
	}
	
	float evaluate_rho(double dt)
	{
		return rho + drho_dt*dt;
	}
};

PLS operator*(const PLS& p, Real v)
{
	PLS t(p.rho*v, p.drho_dt*v);
	return t;
}


template <typename T, int i> inline Real levelset_projector_impl(const T&t)
{
	return (Real)t.rho;
}


make_projector(levelset_projector, SimpleLevelsetBlock<PLS>::levelset_projector_impl)


template<typename RealType>
static void _velocity(const RealType x[3],  RealType  t, RealType  v[3])
{
	const double factor = 0.2;
	const float theta = atan2((x[2]-0.5),(x[1]-0.5));
	/*
	v[0] =  factor*2.0* pow(sin(M_PI*x[0]), 2) * sin(2*M_PI*x[1]) * sin(2*M_PI*x[2]);
	v[1] = -factor* pow(sin(M_PI*x[1]), 2) * sin(2*M_PI*x[0]) *  sin(2*M_PI*x[2]);
	v[2] = -factor* pow(sin(M_PI*x[2]), 2) * sin(2*M_PI*x[0]) * sin(2*M_PI*x[1]);
*/
	v[0] =      x[1]-0.5;//factor*sin(theta);
	v[1] =    -(x[0]-0.5);//-factor*cos(theta);
	if ((pow(v[0],2)+pow(v[1],2))>0.2025) {v[0]=0;v[1]=0;}
	v[2] = (nDim==3)? (-factor* pow(sin(M_PI*x[2]), 2) * sin(2*M_PI*x[0]) * sin(2*M_PI*x[1])):0.0;
}



const int maxStencil[2][3] = {
-3, -3, (nDim==3)? -3:0,
+4, +4, (nDim==3)? +4:1
};


namespace MRAG
{
	
#pragma mark -
#pragma mark Block Processing
	
	struct BP_UpdateRHS_TR
	{
		double dt;
		
		BP_UpdateRHS_TR(const double dt_): dt(dt_){}
		BP_UpdateRHS_TR(const BP_UpdateRHS_TR& c): dt(c.dt) {}
		
		template<typename B>
		inline void operator()(const BlockInfo& info, B& b) const
		{
			typedef typename B::ElementType E;
			
			const int n = B::sizeZ*B::sizeY*B::sizeX;
			
			E* ptrE = &(b[0]);
			
			for(int iE=0; iE<n; iE++, ptrE++)
			{
				ptrE->drho_dt = ptrE->tmp;
				ptrE->rho -= dt*ptrE->tmp;
			}		
		}
	};
	
	struct BP_Integrate_TR
	{
		double dt;
		
		BP_Integrate_TR(const double dt_): dt(dt_){}
		BP_Integrate_TR(const BP_Integrate_TR& c): dt(c.dt) {}
		
		template<typename B>
		inline void operator()(const BlockInfo& info, B& b) const
		{
			typedef typename B::ElementType E;
			
			const int n = B::sizeZ*B::sizeY*B::sizeX;
			
			E* ptrE = &(b[0]);
			
			for(int iE=0; iE<n; iE++, ptrE++)
			{
				ptrE->rho += dt*ptrE->drho_dt;
				ptrE->drho_dt = 0;
			}		
		}
	};

	struct BP_Smoothing
	{
		int stencil_start[3], stencil_end[3];
		
		BP_Smoothing(const BP_Smoothing& c)
		{
			memcpy(this, &c, sizeof(BP_Smoothing) );
		}
		
		BP_Smoothing()
		{
			stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
			stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
		}
		
		template<typename LabType, typename BlockType>
		inline void operator()(LabType& i, const BlockInfo& info, BlockType& o) const
		{	
			typedef BlockType B;
			typedef typename BlockType::ElementType E;
			
			Real w[3]= {0.25,0.5,0.25};

			for(int iz=0; iz<B::sizeZ; iz++)
				for(int iy=0; iy<B::sizeY; iy++)
					for(int ix=0; ix<B::sizeX; ix++)
					{		
						Real val = 0;
						for(int cz =-1; cz<2; cz++)
						for(int cy =-1; cy<2; cy++)
						for(int cx =-1; cx<2; cx++)
							val += i(ix+cx, iy+cy, iz+cz).rho*w[cx+1]*w[cy+1]*w[cz+1];

						 o(ix,iy,iz).tmp = val;
					}
		}
	};

	struct BP_SmoothingUpdate
	{
		template<typename B>
		inline void operator()(const BlockInfo& info, B& b) const
		{
			typedef typename B::ElementType E;
			
			const int n = B::sizeZ*B::sizeY*B::sizeX;
			
			E* ptrE = &(b[0]);
			
			for(int iE=0; iE<n; iE++, ptrE++)
				ptrE->rho = ptrE->tmp;	
		}
	};

	struct BP_ComputeRHS_TR_WENO5
	{
		int stencil_start[3], stencil_end[3];
		Real dt, sign_val;
		
		BP_ComputeRHS_TR_WENO5(const BP_ComputeRHS_TR_WENO5& c): dt(0.),sign_val(1.) 
		{
			memcpy(this, &c, sizeof(BP_ComputeRHS_TR_WENO5) );
		}
		
		BP_ComputeRHS_TR_WENO5(Real dt_, Real sign_): dt(dt_), sign_val(sign_)
		{
			stencil_start[0] = stencil_start[1] = -3;
			stencil_start[2] = (nDim==3)? -3:0;
			stencil_end[0] = stencil_end[1] = 4;
			stencil_end[2] = (nDim==3)? 4:1;
		}
		
		template<typename LabType, typename BlockType>
		inline void operator()(LabType& i, const BlockInfo& info, BlockType& o) const
		{	
			typedef BlockType B;
			typedef typename BlockType::ElementType E;
			
			const Real h = info.h[0];
			
			Real x[3], v[3];
			for(int iz=0; iz<B::sizeZ; iz++)
				for(int iy=0; iy<B::sizeY; iy++)
					for(int ix=0; ix<B::sizeX; ix++)
					{		
						info.pos(x,ix,iy,iz);
						_velocity(x, 0.0f, v);
						
						v[0] *= sign_val;
						v[1] *= sign_val;
						v[2] *= ((nDim==3)? sign_val:0);
						
						o(ix,iy,iz).tmp = -mainRHS(i, h , ix, iy, iz, v, dt);
					}
		}
		
		inline const Real Upwind_Flux(const Real u_plus, const Real u_minus, const Real v_plus, const Real v_minus, const Real w_plus, const Real w_minus, const Real v[3]) const
		{
			const Real term[3] = {
				(v[0]<0)? u_plus : ((v[0]>0)? u_minus : 0.5*(u_minus+u_plus)),
				(v[1]<0)? v_plus : ((v[1]>0)? v_minus : 0.5*(v_minus+v_plus)),
				(v[2]<0)? w_plus : ((v[2]>0)? w_minus : 0.5*(w_minus+w_plus)),
			};
			
			return v[0]*term[0] + v[1]*term[1] + ( (nDim==3)? v[2]*term[2]:0.0 );
		}
		
		inline const Real phi_WENO(const Real a, const Real b, const Real c, const Real d) const
		{
			const Real eps = 1e-6;
			
			const Real IS[3] = {
				13*(a-b)*(a-b) + 3*(a-3*b)*(a-3*b),
				13*(b-c)*(b-c) + 3*(b+c)*(b+c),
				13*(c-d)*(c-d) + 3*(3*c-d)*(3*c-d)
			};
			
			const Real alpha[3] = {
				1./((eps + IS[0])*(eps + IS[0])),
				6./((eps + IS[1])*(eps + IS[1])),
				3./((eps + IS[2])*(eps + IS[2])),
			};
			
			const Real sum_alpha = alpha[0] + alpha[1] + alpha[2];
			
			const Real w[2] = {
				alpha[0]/sum_alpha,
				alpha[2]/sum_alpha
			};
			
			return 1./3*w[0]*(a-2*b+c) + 1./6*(w[1]-0.5)*(b-2*c+d);
		}
		
		template <int index, typename Field> inline const Real Dplus(Field& f, const int ix, const int iy, const int iz) const
		{ 
			//return (f(ix + (int)(index==0), iy + (int)(index==1), iz + (int)(index==2)).rho - f(ix, iy, iz).rho);
			return (f(ix + (int)(index==0), iy + (int)(index==1), iz + (int)(index==2)).evaluate_rho(dt) - f(ix, iy, iz).evaluate_rho(dt));
		}
		
		template <int index, typename Field> inline const Real Dminus(Field& f, const int ix, const int iy, const int iz) const
		{ 
			//return (f(ix, iy, iz).rho - f(ix - (int)(index==0), iy - (int)(index==1), iz - (int)(index==2)).rho );
			return (f(ix, iy, iz).evaluate_rho(dt) - f(ix - (int)(index==0), iy - (int)(index==1), iz - (int)(index==2)).evaluate_rho(dt) );
		}
		
		template <int index, typename Field> inline const Real DplusDminus(Field& f, const int ix, const int iy, const int iz) const
		{ 
			return (f(ix + (int)(index==0), iy + (int)(index==1), iz + (int)(index==2)).evaluate_rho(dt) +
					f(ix - (int)(index==0), iy - (int)(index==1), iz - (int)(index==2)).evaluate_rho(dt) +
					-2*f(ix, iy,iz).evaluate_rho(dt));
		}
		
		template <typename Field>
		inline const Real mainRHS(Field& f, const Real spacing, const int ix, const int iy, const int iz, const Real v[3], const Real dt) const
		{
			const Real u_common_term = 1./(12*spacing)*
			(-Dplus<0>(f, ix-2, iy, iz) + 7*Dplus<0>(f,ix-1, iy, iz) + 7*Dplus<0>(f,ix, iy, iz) -Dplus<0>(f, ix+1, iy, iz) );
			
			const Real v_common_term = 1./(12*spacing)*
			(-Dplus<1>(f, ix, iy-2, iz) + 7*Dplus<1>(f,ix, iy-1, iz) + 7*Dplus<1>(f,ix, iy, iz) -Dplus<1>(f, ix, iy+1, iz) );
			
			const Real w_common_term = (nDim==3)? 1./(12*spacing)*
			(-Dplus<2>(f, ix, iy, iz-2) + 7*Dplus<2>(f,ix, iy, iz-1) + 7*Dplus<2>(f,ix, iy, iz) -Dplus<2>(f, ix, iy, iz+1) ) : 0.0;
			
			const Real phi_weno_x[2] = {
				phi_WENO(DplusDminus<0>(f, ix-2, iy, iz)/spacing,	DplusDminus<0>(f, ix-1, iy, iz)/spacing, 
						 DplusDminus<0>(f, ix, iy, iz)/spacing,		DplusDminus<0>(f, ix+1, iy, iz)/spacing),
				phi_WENO(DplusDminus<0>(f, ix+2, iy, iz)/spacing,	DplusDminus<0>(f, ix+1, iy, iz)/spacing, 
						 DplusDminus<0>(f, ix, iy, iz)/spacing,		DplusDminus<0>(f, ix-1, iy, iz)/spacing)
				
			};
			
			const Real phi_weno_y[2] = {
				phi_WENO(DplusDminus<1>(f, ix, iy-2, iz)/spacing,	DplusDminus<1>(f, ix, iy-1, iz)/spacing, 
						 DplusDminus<1>(f, ix, iy, iz)/spacing,		DplusDminus<1>(f, ix, iy+1, iz)/spacing),
				phi_WENO(DplusDminus<1>(f, ix, iy+2, iz)/spacing,	DplusDminus<1>(f, ix, iy+1, iz)/spacing, 
						 DplusDminus<1>(f, ix, iy, iz)/spacing,		DplusDminus<1>(f, ix, iy-1, iz)/spacing)
			};
			
			
			const Real phi_weno_z[2] = { (nDim==3)?
				phi_WENO(DplusDminus<2>(f, ix, iy, iz-2)/spacing,	DplusDminus<2>(f, ix, iy, iz-1)/spacing, 
						 DplusDminus<2>(f, ix, iy, iz)/spacing,		DplusDminus<2>(f, ix, iy, iz+1)/spacing) : 0.0, (nDim==3)?
				phi_WENO(DplusDminus<2>(f, ix, iy, iz+2)/spacing,	DplusDminus<2>(f, ix, iy, iz+1)/spacing, 
						 DplusDminus<2>(f, ix, iy, iz)/spacing,		DplusDminus<2>(f, ix, iy, iz-1)/spacing) : 0.0
			
			};
			const Real u_minus = u_common_term - phi_weno_x[0];
			const Real u_plus  = u_common_term + phi_weno_x[1];
			const Real v_minus = v_common_term - phi_weno_y[0];
			const Real v_plus  = v_common_term + phi_weno_y[1];
			const Real w_minus = w_common_term - phi_weno_z[0];
			const Real w_plus  = w_common_term + phi_weno_z[1];
			
			return Upwind_Flux( u_plus, u_minus, v_plus, v_minus, w_plus, w_minus, v);
		}
	};
	
	struct BP_Integrate_REINIT
	{
		public:
			double dt;
			BP_Integrate_REINIT(double dt_): dt(dt_) {}
			BP_Integrate_REINIT(const BP_Integrate_REINIT& i): dt(i.dt){}
			
			template <typename BlockType>
			inline void operator() (const BlockInfo& info, BlockType& b) const
			{
				typedef BlockType B;
				
				const int n = B::sizeZ*B::sizeY*B::sizeX;
				
				//LevelsetPoint* ptrE = &(b(0));
				PLS* ptrE = &(b(0));
				for(int iE=0; iE<n; iE++, ptrE++)
				{
					//ptrE->phi += dt*ptrE->dphidt;
					ptrE->rho += dt*ptrE->drho_dt;
					ptrE->drho_dt = 0;
				}		
			}
		};
	
	struct BP_ComputeRHS_WENO5_REINIT
	{
		inline const Real phi_WENO(const Real a, const Real b, const Real c, const Real d) const
		{
			const Real eps = 1e-6;
			
			const Real IS[3] = {
				13*(a-b)*(a-b) + 3*(a-3*b)*(a-3*b),
				13*(b-c)*(b-c) + 3*(b+c)*(b+c),
				13*(c-d)*(c-d) + 3*(3*c-d)*(3*c-d)
			};
			
			const Real alpha[3] = {
				1./((eps + IS[0])*(eps + IS[0])),
				6./((eps + IS[1])*(eps + IS[1])),
				3./((eps + IS[2])*(eps + IS[2]))
			};
			
			const Real sum_alpha = alpha[0] + alpha[1] + alpha[2];
			
			const Real w[2] = {
				alpha[0]/sum_alpha,
				alpha[2]/sum_alpha
			};
			
			return 1./3*w[0]*(a-2*b+c) + 1./6*(w[1]-0.5)*(b-2*c+d);
		}
		
		inline const Real H_GODUNOV(const Real spacing, const Real phi, const Real u_plus, const Real u_minus, const Real v_plus, const Real v_minus, const Real w_plus, const Real w_minus) const
		{
			const Real s = phi/sqrt(phi*phi + spacing);
			
			if (phi>=0)
			{
				const Real term1 = max( -min(u_plus, (Real)0.0), max(u_minus, (Real)0.0) );
				const Real term2 = max( -min(v_plus, (Real)0.0), max(v_minus, (Real)0.0) );
				const Real term3 = max( -min(w_plus, (Real)0.0), max(w_minus, (Real)0.0) );
				
				return s*(sqrt(term1*term1 + term2*term2 + ( (nDim==3)? term3*term3 : 0.0 ) )-1);
			}
			else
			{
				const Real term1 = max( -min(u_minus, (Real)0.0), max(u_plus, (Real)0.0) );
				const Real term2 = max( -min(v_minus, (Real)0.0), max(v_plus, (Real)0.0) );
				const Real term3 = max( -min(w_minus, (Real)0.0), max(w_plus, (Real)0.0) );
				
				return s*(sqrt(term1*term1 + term2*term2 + ( (nDim==3)? term3*term3 : 0.0 ) )-1);
			}
		}
		
		template <int index, typename Field> inline const Real Dplus(Field& f, const int ix, const int iy, const int iz) const
		{ 
			return (f(ix + (int)(index==0), iy + (int)(index==1), iz + (int)(index==2)).rho - f(ix, iy, iz).rho);
		}
		
		template <int index, typename Field> inline const Real Dminus(Field& f, const int ix, const int iy, const int iz) const
		{ 
			return (f(ix, iy, iz).rho - f(ix - (int)(index==0), iy - (int)(index==1), iz - (int)(index==2)).rho );
		}
		
		template <int index, typename Field> inline const Real DplusDminus(Field& f, const int ix, const int iy, const int iz) const
		{ 
			return (f(ix + (int)(index==0), iy + (int)(index==1), iz + (int)(index==2)).rho +
					f(ix - (int)(index==0), iy - (int)(index==1), iz - (int)(index==2)).rho +
					-2*f(ix, iy,iz).rho );
		}
		
		template <typename Field>
		const Real mainRHS(Field& f, const Real spacing, const int ix, const int iy, const int iz) const
		{
			const Real u_common_term = 1./(12*spacing)*
			(-Dplus<0>(f, ix-2, iy, iz) + 7*Dplus<0>(f,ix-1, iy, iz) + 7*Dplus<0>(f,ix, iy, iz) -Dplus<0>(f, ix+1, iy, iz) );
			
			const Real v_common_term = 1./(12*spacing)*
			(-Dplus<1>(f, ix, iy-2, iz) + 7*Dplus<1>(f,ix, iy-1, iz) + 7*Dplus<1>(f,ix, iy, iz) -Dplus<1>(f, ix, iy+1, iz) );
			
			const Real w_common_term = (nDim==3)? 1./(12*spacing)*
			(-Dplus<2>(f, ix, iy, iz-2) + 7*Dplus<2>(f,ix, iy, iz-1) + 7*Dplus<2>(f,ix, iy, iz) -Dplus<2>(f, ix, iy, iz+1) ) : 0.0;
			
			const Real phi_weno_x[2] = {
				phi_WENO(DplusDminus<0>(f, ix-2, iy, iz)/spacing,	DplusDminus<0>(f, ix-1, iy, iz)/spacing, 
						 DplusDminus<0>(f, ix, iy, iz)/spacing,		DplusDminus<0>(f, ix+1, iy, iz)/spacing),
				phi_WENO(DplusDminus<0>(f, ix+2, iy, iz)/spacing,	DplusDminus<0>(f, ix+1, iy, iz)/spacing, 
						 DplusDminus<0>(f, ix, iy, iz)/spacing,		DplusDminus<0>(f, ix-1, iy, iz)/spacing)
				
			};
			
			const Real phi_weno_y[2] = {
				phi_WENO(DplusDminus<1>(f, ix, iy-2, iz)/spacing,	DplusDminus<1>(f, ix, iy-1, iz)/spacing, 
						 DplusDminus<1>(f, ix, iy, iz)/spacing,		DplusDminus<1>(f, ix, iy+1, iz)/spacing),
				phi_WENO(DplusDminus<1>(f, ix, iy+2, iz)/spacing,	DplusDminus<1>(f, ix, iy+1, iz)/spacing, 
						 DplusDminus<1>(f, ix, iy, iz)/spacing,		DplusDminus<1>(f, ix, iy-1, iz)/spacing)
			};
			
			const Real phi_weno_z[2] = { (nDim==3)?
				phi_WENO(DplusDminus<2>(f, ix, iy, iz-2)/spacing,	DplusDminus<2>(f, ix, iy, iz-1)/spacing, 
						 DplusDminus<2>(f, ix, iy, iz)/spacing,		DplusDminus<2>(f, ix, iy, iz+1)/spacing) : 0.0, (nDim==3)?
				phi_WENO(DplusDminus<2>(f, ix, iy, iz+2)/spacing,	DplusDminus<2>(f, ix, iy, iz+1)/spacing, 
						 DplusDminus<2>(f, ix, iy, iz)/spacing,		DplusDminus<2>(f, ix, iy, iz-1)/spacing) : 0.0
			};
			
			const Real u_minus = u_common_term - phi_weno_x[0];
			const Real u_plus = u_common_term + phi_weno_x[1];
			const Real v_minus = v_common_term - phi_weno_y[0];
			const Real v_plus = v_common_term + phi_weno_y[1];
			const Real w_minus = w_common_term - phi_weno_z[0];
			const Real w_plus = w_common_term + phi_weno_z[1];
			
			return -H_GODUNOV(spacing, f(ix,iy, iz).rho_0, u_plus, u_minus, v_plus, v_minus, w_plus, w_minus);
		}
		
	public:
		int stencil_start[3], stencil_end[3];
		
		BP_ComputeRHS_WENO5_REINIT()
		{
			stencil_start[0] = -3;
			stencil_start[1] = -3;
			stencil_start[2] = (nDim==3)? -3:0;
			
			stencil_end[0] = 4;
			stencil_end[1] = 4;
			stencil_end[2] = (nDim==3)? 4:1;
		}
		
		BP_ComputeRHS_WENO5_REINIT(const BP_ComputeRHS_WENO5_REINIT& c)
		{
			stencil_start[0] = -3;
			stencil_start[1] = -3;
			stencil_start[2] = (nDim==3)? -3:0;
			
			stencil_end[0] = 4;
			stencil_end[1] = 4;
			stencil_end[2] = (nDim==3)? 4:1;
		}
		
		template<typename LabType, typename BlockType>
		inline void operator()(LabType& i, const BlockInfo& info, BlockType& o) const
		{
			typedef BlockType B;
			typedef typename BlockType::ElementType E;
			
			const Real h = info.h[0];
			
			for(int iz=0; iz<B::sizeZ; iz++)
				for(int iy=0; iy<B::sizeY; iy++)
					for(int ix=0; ix<B::sizeX; ix++)
						o(ix,iy,iz).drho_dt = mainRHS(i, h , ix, iy, iz);
		}
		
	};
	
	struct ReadFromFile
		{
			/*int m_vSize[3];
			char * m_data;
			string m_sFileName;*/
			
			int m_vSize[3];
			int m_iCurrent;
			char * m_data;
			bool m_bCollocated;
			bool m_bWritten;
			
	
			ReadFromFile():
			m_data(NULL)
			{
				
			}
			
			~ReadFromFile()
			{
				if (m_data != NULL) delete [] m_data;
			}
			void check_file(string m_sFileName)
			{
				FILE * f = fopen(m_sFileName.data(), "rb");
				assert(f!=NULL);
				
				string sTrueString = m_sFileName;
				fread(this, sizeof(ReadFromFile), 1, f);
				printf("Reading... size= %dx%dx%d, collocated:%s\n", m_vSize[0], m_vSize[1], m_vSize[2], m_bCollocated?"yes":"no");
				
				m_data = new char[m_vSize[0]*m_vSize[1]*m_vSize[2]];
				assert(m_data!=NULL);
				fread(m_data, sizeof(char), m_vSize[0]*m_vSize[1]*m_vSize[2], f);
				fclose(f);
			}
			
			const float operator()(const int ix, const int iy, const int iz)
			{
				assert(ix>=0 && ix<m_vSize[0]);
				assert(iy>=0 && iy<m_vSize[1]);
				assert(iz>=0 && iz<m_vSize[2]);
				
				const int index = ix + iy*m_vSize[0] + iz*m_vSize[0]*m_vSize[1];
				
				return (float)m_data[index];
			}
		};
	
	template <typename Wavelets>
	class STDTestL4_LevelsetAdvection //LEVEL 4: correct 3D simulations, correct data-processing
	{
	public:
		enum TestMode{TestMode_Advection, TestMode_IC, TestMode_PostProcessing};
		static const TestMode test_mode = TestMode_Advection;
		//static const TestMode test_mode = TestMode_IC;
		static const bool bRestartFromFile = false;
		static const bool bPostProcessData = false;
		static const bool bUseSpaceTimeSorter = true;
		static const int nStepsPerSerialization = 1;
		static const int blocksPerDimension = 1;
		static const int blockSize = 32;
		static const int maxLevel = 2;
		static const int resJump = 2;
		static const int narrowBandWidth = 6;
		
		typedef SimpleLevelsetBlock< PLS, narrowBandWidth, blockSize, blockSize, (nDim==3)? blockSize:1 > B;
		//typedef Wavelets_Interp2ndOrder W;
		typedef Wavelets_Interp4thOrder W;
		Grid<W, B> grid;
		
#ifndef _MRAG_TBB
		typedef Multithreading::BlockProcessing_SingleCPU<B> BlockProcessing;
#else
		//typedef Multithreading::BlockProcessing_TBB<B> BlockProcessing;
		typedef Multithreading::BlockProcessing_Pipeline_TBB<B, BlockLab, _MRAG_TBB_NTHREADS_HINT> BlockProcessing;
#endif
		
	protected:
		BlockProcessing blockProcessing;
		Refiner_SpaceExtension refiner;
		Compressor compressor;
		BlockFWT<W, B, levelset_projector> blockfwt;
		AutomaticRefiner<Grid<W, B>, BlockFWT<W, B, levelset_projector> > automatic_refiner;
#ifdef _MRAG_GLUT_VIZ
		GridViewer viewer;
#endif
		SpaceTimeSorter stSorter;
		Profiler profiler;
		string m_sFormat, m_sToRenderFormat;
		int m_iCurrentDumpTime, serialize_counter;
		float t;
		double refinement_tolerance, compression_tolerance;
#pragma mark -
#pragma mark Class Management
		
	public:
		
		STDTestL4_LevelsetAdvection():
#ifdef _MRAG_GLUT_VIZ
			viewer(true,true),
#endif
		grid(blocksPerDimension,blocksPerDimension, (nDim==3)? blocksPerDimension:1, maxStencil), 
			blockProcessing(),
			refiner(resJump), compressor(resJump), profiler(), 
			blockfwt(), t(0), m_iCurrentDumpTime(0), serialize_counter(0),
			stSorter(),
			m_sToRenderFormat("Grid_At_Time%03d.%s"), 
			m_sFormat("./data/SERIALIZED_%03d"),
			refinement_tolerance(1e-4), compression_tolerance(1e-5),
			automatic_refiner(grid, blockfwt, &refiner, 200.0, &profiler)
		{
			grid.setCompressor(&compressor);
			grid.setRefiner(&automatic_refiner);
			
			stSorter.connect(grid);
			
			if (test_mode == TestMode_Advection)
			{
				if (bRestartFromFile)
					_RestartFromFile("E:\\Smoke\\VAI_CAZZO!\\SERIALIZED_099",99, t, grid, m_iCurrentDumpTime);
				else
				{
					_ic(grid);
					 //automatic_refiner.refine(refinement_tolerance, compression_tolerance, _ic);
					Science::AutomaticRefinementForLevelsets(grid, blockfwt, refinement_tolerance);
					 _ic(grid);
					 Science::AutomaticCompressionForLevelsets(grid, blockfwt, compression_tolerance);         
				 
					 _ic(grid);
				}
				
					//_RestartFromFile("./data/ICDATA",0, t, grid, m_iCurrentDumpTime);
			}
			else if (test_mode == TestMode_IC)
			{
		
				_icFromBinFile(grid);
				/*automatic_refiner.refine(refinement_tolerance, compression_tolerance, _icFromBinFile);

				_icFromBinFile(grid);
				Science::AutomaticCompressionForLevelsets(grid, blockfwt, compression_tolerance);         

				_icFromBinFile(grid);*/
			}
		}
		
	protected:
		
		void _RestartFromFile(string sFilePath, int stepid, float& output_time, Grid<W,B>& output_grid, int& out_currentDumpTime) 
		{
			//1. compute the actual time t
			//2. fill the grid from the restart file
			//3. fill the dump time, serialize counter
			
			//1.
			const double dt = computeDeltaTime();
			output_time = dt*stepid*nStepsPerSerialization;
			
			//2.
			IO_Native<W, B, levelset_projector> serializer;
			serializer.Read(output_grid, sFilePath);
			output_grid.getBoundaryInfo();
			printf("READ, size =%f MB (boundary size=%f MB)\n",  output_grid.getMemorySize(), output_grid.getBoundaryInfo().getMemorySize());
			printf("#BLOCKS = %d\n", output_grid.getBlocksInfo().size());
			
			//3.
			out_currentDumpTime = stepid;
			serialize_counter = stepid*nStepsPerSerialization + 1;
		}
		
	public:
		
		static void idle(void)
		{
#ifdef	_MRAG_GLUT_VIZ
			glutPostRedisplay();
#endif
		}

		static STDTestL4_LevelsetAdvection<Wavelets> * singletone; 
		
		static void runTests(int argc, char ** argv, bool bVisual=true)
		{
#ifdef _MRAG_GLUT_VIZ
			if (bVisual)
			{
				glutInit(&argc, argv);
				glutInitWindowSize(800,800);
				glutInitWindowPosition(0, 0);
				glutInitDisplayMode(GLUT_DEPTH| GLUT_STENCIL |GLUT_RGBA | GLUT_DOUBLE );
				
				glutCreateWindow("MRAG Refinement Test");
				
				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();
				
				glOrtho(-0.2, 1.2, -0.2, 1.2, -1, 1);
				glMatrixMode(GL_MODELVIEW);
				
				glEnableClientState(GL_VERTEX_ARRAY);
				glEnableClientState(GL_TEXTURE_COORD_ARRAY);
				glEnable(GL_TEXTURE_2D);
				
				glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
				glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
				
				glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
				
				glutDisplayFunc(display);
				glutIdleFunc(idle);
			}
#endif			
			singletone = new STDTestL4_LevelsetAdvection<Wavelets>();

#ifdef _MRAG_GLUT_VIZ			
			if (bVisual) 
				glutMainLoop();
			else
#endif
				for(int i=0; i<15; i++)
					singletone->run();
		}

	
#pragma mark -
#pragma mark Computing
		
	protected:
		
		static float _ic_func4(float x[3]) 
		{
			const float r = sqrt(pow(x[0]-0.5,2) + pow(x[1]-0.75, 2) + ( (nDim==3)? pow(x[2]-0.35,2):0.0 ) );
			return 0.15-r;
		}
		
		static inline double _Bspline4(double x) 
		{
			const double t = fabs(x);
			
			if (t>2) return 0;
			
			if (t>1) return pow(2-t,3)/6;
			
			return (1 + 3*(1-t)*(1 + (1-t)*(1 - (1-t))))/6;
			
			return 0.;
		}

		template<int dir>
		static inline const double HeavySide(double x0, double x)
		{
			//if (dir == 1) return (int)(x>x0);
			//return (int)(x<x0);
			const double eps = 0.1/5;

			const double alpha = M_PI*min(1., max(0., (x-x0+0.5*eps)/eps));
			if (dir==1) return 0.5+0.5*sin(alpha - M_PI/2);
			else return 0.5+0.5*cos(alpha);
		}
		
		static inline double _sample(ReadFromFile& rtf, double x, double y, double z, const int n)
		{
			const int anchor_point[3] = {
				(int)floor(x),
				(int)floor(y),
				(int)floor(z)
			};
			
			double val = 0;
			double wsum = 0;
			
			//return sqrt(pow(x-n/2,2)+pow(y-n/2,2)+pow(z-n/2,2)) - 0.25*n;

			for(int iz=-2; iz<3; iz++)
			for(int iy=-2; iy<3; iy++)
			for(int ix=-2; ix<3; ix++)
			{
				const int src[3] = {
					(anchor_point[0]+ ix + n) % n, 
					(anchor_point[1]+ iy + n) % n,
					(anchor_point[2]+ iz + n) % n
				};
				
				const double w = 
						_Bspline4(anchor_point[0]+ix+0.5- x)*
						_Bspline4(anchor_point[1]+iy+0.5- y)*
						_Bspline4(anchor_point[2]+iz+0.5- z);
				
				val += rtf(src[0], src[1],src[2])*w;
				
				wsum+= w;
			}
			
			assert(fabs(wsum-1.0)<1e-4);
			
			return val;
		}
		
		static void _icFromBinFile(Grid<W, B>& grid)
		{
			vector<BlockInfo> vInfo = grid.getBlocksInfo();

			ReadFromFile rtf;
			
			rtf.check_file( "./data/horse.bin");
			 
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				
				block.setH(info.h[0]);
				
				const int n_per_dim = rtf.m_vSize[0];
				const double h = 1./rtf.m_vSize[0];
				
				for(int iz=0; iz<B::sizeZ; iz++)
					for(int iy=0; iy<B::sizeY; iy++)
						for(int ix=0; ix<B::sizeX; ix++)
						{
							double x[3];
							info.pos(x, ix, iy,iz);
							/*const double val = _sample(rtf,
										x[0]*n_per_dim, 
										x[1]*n_per_dim,
										x[2]*n_per_dim, n_per_dim*/
							const double hwidth = 0.1;
							
							const double val =  HeavySide<1>(0.5-hwidth,x[0]) * HeavySide<-1>(0.5+hwidth,x[0]) *
												HeavySide<1>(0.5-hwidth,x[1]) * HeavySide<-1>(0.5+hwidth,x[1]) *
												HeavySide<1>(0.5-hwidth,x[2]) * HeavySide<-1>(0.5+hwidth,x[2]);
							//);
							
							block(ix,iy,iz).rho_0 = block(ix,iy,iz).rho = val * (-2.) + 1.;
							//rtf(ix+info.index[0]*B::sizeX, iy+info.index[1]*B::sizeY ,iz+info.index[2]*B::sizeZ);
							//block(ix,iy,iz).rho;
						}
				grid.getBlockCollection().release(info.blockID);
			}
		}
		
		static void _ic(Grid<W, B>& grid) 
		{
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection().lock(info.blockID);
				
				block.setH(info.h[0]);
				
				for(int iz=0; iz<B::sizeZ; iz++)
					for(int iy=0; iy<B::sizeY; iy++)
						for(int ix=0; ix<B::sizeX; ix++)
						{
							float x[3];
							
							info.pos(x, ix, iy, iz);
							
							block(ix,iy,iz).rho = _ic_func4(x);
						}
				
				grid.getBlockCollection().release(info.blockID);
			}
		}

		double computeDeltaTime() const
		{
			const double dx = 1.0/blockSize;
			const double CFL = 0.25;
			const double dt = dx*CFL/0.2;//8e-1;
			
			return dt;
		}

		void _smoothing_step(const int nSteps, BoundaryInfo* boundaryInfo, const int nParallelGranularity)
		{
			BP_Smoothing smoothing;
			BP_SmoothingUpdate update;
			
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			profiler.getAgent("computing").start();
			
			for(int i=0; i<nSteps; i++)
			{
				//BlockProcessing::process<BlockLab>(vInfo, grid.getBlockCollection(), grid.getBoundaryInfo(), computeRhs, -1);
				blockProcessing.pipeline_process(vInfo, grid.getBlockCollection(), *boundaryInfo, smoothing);
				
				BlockProcessing::process(vInfo, grid.getBlockCollection(),  update, nParallelGranularity);
			}
			profiler.getAgent("computing").stop();
			
			cout<<"Done smoothing"<<endl;
		}
		void _advection_step(float t, float dt, BoundaryInfo* boundaryInfo, const int nParallelGranularity)
		{
			stSorter.startSession(dt, 2, 0);
			
			if(t>9.0) exit(0);
		
			while(true)
			{
				double currTime, currDeltaT;
				int level;
				SpaceTimeSorter::ETimeInterval type;
				vector<BlockInfo> vInfo;
				
				const bool bContinue = stSorter.getBlocks(level, currDeltaT, currTime, vInfo, type);
				
				if (type == SpaceTimeSorter::ETimeInterval_Start)
				{
					BP_ComputeRHS_TR_WENO5 task1(currTime, /*(t+currTime<4.5)? 1 : -1*/ 1);
					BP_UpdateRHS_TR task2(currTime);
					
					profiler.getAgent("computing").start();
					//BlockProcessing::process< BlockLab >(vInfo, grid.getBlockCollection(), *boundaryInfo, task1, nParallelGranularity);
					blockProcessing.pipeline_process(vInfo, grid.getBlockCollection(), *boundaryInfo, task1);
					BlockProcessing::process(vInfo, grid.getBlockCollection(), task2, nParallelGranularity);
					profiler.getAgent("computing").stop();
				}
				else
				{
					BP_Integrate_TR task(currTime);
					profiler.getAgent("computing").start();
					BlockProcessing::process(vInfo, grid.getBlockCollection(), task, nParallelGranularity);
					profiler.getAgent("computing").stop();
				}
				
				if (!bContinue) break;
			}
			stSorter.endSession();		
		}
		
		void _reinitialization_step(const int nSteps, BoundaryInfo* boundaryInfo, const int nParallelGranularity)
		{
			BP_ComputeRHS_WENO5_REINIT computeRhs;
			BP_Integrate_REINIT integrate(1e-4);
			
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			profiler.getAgent("computing").start();
			
			for(int i=0; i<nSteps; i++)
			{
				//BlockProcessing::process<BlockLab>(vInfo, grid.getBlockCollection(), grid.getBoundaryInfo(), computeRhs, -1);
				blockProcessing.pipeline_process(vInfo, grid.getBlockCollection(), *boundaryInfo, computeRhs);
				
				BlockProcessing::process(vInfo, grid.getBlockCollection(),  integrate, nParallelGranularity);
			}
			profiler.getAgent("computing").stop();
			
			cout<<"Done re-init"<<endl;
		}
		
		void _step(float t, float dt)
		{
			const int nParallelGranularity = (grid.getBlocksInfo().size()<=8 ? 1 : 4);
			BoundaryInfo* boundaryInfo = &grid.getBoundaryInfo();
			
			
			if (test_mode == TestMode_Advection)
				_advection_step(t, dt, boundaryInfo, nParallelGranularity);
			
			if(test_mode == TestMode_Advection){
				_reinitialization_step(40, boundaryInfo, nParallelGranularity);
			}
			else if (test_mode == TestMode_IC)
			{
				//_smoothing_step(2, boundaryInfo, nParallelGranularity);

				_reinitialization_step(4*200, boundaryInfo, nParallelGranularity);
			}
		}
		
		bool run()
		{
			if (test_mode == TestMode_Advection)
			{
				if (serialize_counter++ % nStepsPerSerialization == 0)
				{
					profiler.printSummary();

					{
						IO_Native<W, B, levelset_projector> serializer;
						char buf[300];
						sprintf(buf, m_sFormat.data(), m_iCurrentDumpTime);
						serializer.Write(grid, buf);
					}

					/*{
						char buf[300];
						IO_VTK<W, B, levelset_projector, 0> vtkserializer;
						sprintf(buf, (m_sFormat+"VTK").data(), m_iCurrentDumpTime);
						vtkserializer.Write(grid, buf);
					}*/
					m_iCurrentDumpTime++;
				}

				const int nStep = 1;
				const float dt = computeDeltaTime();
				
				//soft refinement for memory reasons
				const double old_refinement_tolerance = refinement_tolerance;
				const double old_compression_tolerance = compression_tolerance;
				
				//const bool bTolChanged = automatic_refiner.refine(refinement_tolerance, compression_tolerance);
				
				Science::AutomaticRefinementForLevelsets(grid, blockfwt, refinement_tolerance);
				
				/*if (bTolChanged)
				{
					printf("refinement_tolerance changed! from %e to %e\n", old_refinement_tolerance, refinement_tolerance);
					printf("compression_tolerance changed! from %e to %e\n", old_compression_tolerance, compression_tolerance);
				}*/
				cout<<"I am at time= "<<t<<" and my refinement_tolerance is "<< refinement_tolerance <<endl;
			
				printf("/////////// WHILE COMPUTING MEM MB: %f  (%d blocks, BI size=%f MB) ////////////////// \n", grid.getMemorySize(), grid.getBlocksInfo().size(), grid.getBoundaryInfo().getMemorySize());
				
				for(int i=0; i<nStep; i++, t+=dt)
					_step(t, dt);
				
				Science::AutomaticCompressionForLevelsets(grid, blockfwt, compression_tolerance, true, -1, &profiler);  
				
				printf("/////////// AFTER COMPUTING MEM MB: %f  (%d blocks, BI size=%f MB) ////////////////// \n", grid.getMemorySize(), grid.getBlocksInfo().size(), grid.getBoundaryInfo().getMemorySize());
				
				//serialize the data
				//static int serialize_counter = 0;
				
				
				return false;
			}
			else if (test_mode == TestMode_IC)
			{
				IO_VTK<W, B, levelset_projector, 0> vtkserializer;
			//	vtkserializer.Write(grid, "./data/ICcube_noreinit");
				_step(t, computeDeltaTime());
				vtkserializer.Write(grid, "./data/ICcube_reinit_800");

			//	Science::AutomaticCompressionForLevelsets(grid, blockfwt, compression_tolerance, true, -1, &profiler); 
				
				IO_Native<W, B, levelset_projector> serializer;
				serializer.Write(grid, "./data/ICDATA");	
				//IO_VTK<W, B, levelset_projector, 0> vtkserializer;
				//vtkserializer.Write(grid, "./data/ICDATAVTK");

				exit(0);
			}
			else if (test_mode == TestMode_PostProcessing)
			{
				IO_Native<W, B, levelset_projector> serializer;
				
				const int nGrids = 2;
				for(int iGrid=0; iGrid<=nGrids; iGrid++)
				{
					printf("POSTPROCESS PHASE \n");
					char buf[300];
					sprintf(buf, m_sFormat.data(), iGrid);
					serializer.Read(grid, buf);
					
					DumpData();
				}
				exit(0);
				return true;
			}
		}
		
#pragma mark -
#pragma mark Postprocessing
		
#ifdef _MRAG_GLUT_VIZ
		void _drawGridPoint3f(const int block_index[3],  int level, const float point_index[3], const float *vColorBoundary = NULL, double factor=1.)
		{
			const double h[2]= {pow(2.,-level), pow(2.,-level)};
			const bool bCollocated = true;
			const double start[2] = {block_index[0]*h[0],block_index[1]*h[1]};
			//const double end[2] = {(block_index[0]+1)*h[0],(block_index[1]+1)*h[1]};
			const int n[2] = {B::sizeX+1, B::sizeY+1};
			const double d[2] = {h[0]/(n[0]-1),h[1]/(n[1]-1)};
			const double point_center[2] = {start[0] + (point_index[0]+ (bCollocated?0: 0.5))*d[0], start[1] + (point_index[1]+(bCollocated?0: 0.5))*d[1]};
			
			glPointSize(600*d[0]/2*factor);
			glBegin(GL_POINTS);
			if (vColorBoundary != NULL)
				glColor3fv(vColorBoundary);
			else
				glColor3f(1,0,0);
			
			
			glVertex2f(point_center[0], point_center[1]);
			glEnd();
		}		
		
		void _drawLevelsetIntersections()
		{
			const float vContent[3] = {1.0,0.0,0.573};
			
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				const int level = info.level;
				const int block_index[] = {
					info.index[0],info.index[1],info.index[2]
				};
				for(int iz=0; iz<B::sizeZ; iz++)
					for(int iy=0; iy<B::sizeY-1; iy++)
						for(int ix=0; ix<B::sizeX-1; ix++)
						{
							if(block(ix,iy,iz).rho*block(ix+1,iy,iz).rho<0.0){
								const float x[] = {
									(float)ix + block(ix,iy,iz)/(block(ix,iy,iz)-block(ix+1,iy,iz)),
									(float)iy,
									0.0
								};
								_drawGridPoint3f(block_index, level, x,vContent,1.0);
							} 
							if(block(ix,iy,iz).rho*block(ix,iy+1,iz).rho<0.0){
								const float x[] = {
									(float)ix,
									(float)iy + block(ix,iy,iz)/(block(ix,iy,iz)-block(ix,iy+1,iz)),
									0.0
								};
								_drawGridPoint3f(block_index, level, x,vContent,1.0);
							} 
						}
			}
		}
		
		static void display(void)
		{
			glClear(GL_COLOR_BUFFER_BIT);
			
			const bool bStop = singletone->run();
			
			singletone->viewer.drawSketch(singletone->grid,false);
			singletone->_drawLevelsetIntersections();
			
			if (bStop) exit(0);
			
			glutSwapBuffers();
		}
#endif
		
		void DumpData()
		{
			//1. create suitable filenames
			//2. dump the information of the blocks
			//3. dump the information about the data
			
			//1.
			profiler.getAgent("dumping").start();
			
			char buf[300], buf2[300];
			sprintf(buf, m_sToRenderFormat.data(), m_iCurrentDumpTime, "txt");
			sprintf(buf2, m_sToRenderFormat.data(), m_iCurrentDumpTime, "grid");
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			
			
			//2.
			{
				FILE * file = fopen(buf, "w");
				
				assert(file!=NULL);
				fprintf (file, "Wavelets: %s\n", typeid(W).name()); 
				fprintf(file, "Cell-centered? %s\n", W::bIsCellCentered ? "yes" : "no");
				fprintf(file, "Blocks: %d\n", vInfo.size());
				fprintf(file, "Block size: %d %d %d\n", blockSize, blockSize, blockSize);
				fprintf(file, "Ghosts at 0: [%d, %d[\n", W::HsSupport[0], W::HsSupport[1]);
				
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					fprintf(file, "Block %d: Tree Index: %d %d %d, %d\n", i,
							info.index[0], info.index[1], info.index[2], info.level); 
				}
				
				fclose(file);
			}
			
			//3.
			{
				FILE * file = fopen(buf2, "wb");
				assert(file!=NULL);
				
				const int steStart[3] ={ W::HsSupport[0], W::HsSupport[0], W::HsSupport[0]};
				const int steEnd[3] ={  W::HsSupport[1], W::HsSupport[1], W::HsSupport[1]};
				
				BlockLab<B> lab;
				lab.prepare(grid.getBlockCollection(), grid.getBoundaryInfo(),steStart,steEnd);
				
				const int sX = W::HsSupport[0];
				const int sY = W::HsSupport[0];
				const int sZ = W::HsSupport[0];
				
				const int eX = blockSize + W::HsSupport[1] - 1;
				const int eY = blockSize + W::HsSupport[1] - 1;
				const int eZ = blockSize + W::HsSupport[1] - 1;
				
				Matrix3D<float> * matData = new Matrix3D<float>(eX - sX, eY - sY, eZ - sZ);
				
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					
					lab.load(info);
					
					for(int iz=sZ; iz<eZ; iz++)
						for(int iy=sY; iy<eY; iy++)
							for(int ix=sX; ix<eX; ix++)
								matData->Access(ix-sX, iy-sY, iz-sZ) = (float)lab(ix, iy, iz).rho;
					
					matData->Serialize(file);
				}
				
				fclose(file);
			}
			
			m_iCurrentDumpTime ++;
			profiler.getAgent("dumping").stop();
		}		
		
	
	};
			
	template <typename Wavelets>
	STDTestL4_LevelsetAdvection<Wavelets> * STDTestL4_LevelsetAdvection<Wavelets>::singletone; 
}
