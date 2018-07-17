/*
 *  MRAG_STDTestL4_LevelsetAdvection.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 11/19/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include "../MRAGcore/MRAGCommon.h"
#include "../MRAGcore/MRAGEnvironment.h"


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

#include "../MRAGcore/MRAGrid.h"
#include "../MRAGcore/MRAGRefiner.h"
#include "../MRAGcore/MRAGCompressor.h"
#include "../MRAGcore/MRAGBlockLab.h"
#include "../MRAGcore/MRAGBlockFWT.h"

#ifdef _MRAG_GLUT_VIZ
#include "../MRAGvisual/GridViewer.h"
#endif

#include "../MRAGscience/MRAGScienceCore.h"
#include "../MRAGscience/MRAGSimpleLevelsetBlock.h"
#include "../MRAGscience/MRAGSpaceTimeSorter.h"
#include "../MRAGscience/MRAGRefiner_SpaceExtension.h"

#include "../MRAGmultithreading/MRAGBlockProcessing_SingleCPU.h"
#include "../MRAGmultithreading/MRAGBlockProcessing_TBB.h"

#include "../MRAGio/MRAG_IO_Native.h"

using namespace MRAG;

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
	
	v[0] =  factor*2.0* pow(sin(M_PI*x[0]), 2) * sin(2*M_PI*x[1]) * sin(2*M_PI*x[2]);
	v[1] = -factor* pow(sin(M_PI*x[1]), 2) * sin(2*M_PI*x[0]) *  sin(2*M_PI*x[2]);
	v[2] = -factor* pow(sin(M_PI*x[2]), 2) * sin(2*M_PI*x[0]) * sin(2*M_PI*x[1]);

}


template<typename RealType>
static void _velocity_RBRot(const RealType x[3],  RealType  t, RealType  v[3])
{
	const double factor = 0.2;
	
	v[0] =  factor* M_PI * (x[1] + x[2] - 1.);
	v[1] =  factor* M_PI * (0.5 - x[0]);
	v[2] =  factor* M_PI * (0.5 - x[0]);
	
}

template<typename RealType>
static void _velocity_RBTrans( RealType x[3],  RealType  t, RealType  v[3])
{
	const double factor = 0.2;
	
	v[0] =  factor;
	v[1] =  factor;
	v[2] =  -factor;
	
}


//float toleranceIC = 0.001;

const int maxStencil[2][3] = {
	-3, -3, -3,
	+4, +4, +4
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
			stencil_start[0] = stencil_start[1] = stencil_start[2] = -3;
			stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
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
						//_velocity(x, 0.0f, v);
						_velocity_RBRot(x, 0.0f, v);
						//_velocity_RBTrans(x, 0.0f, v);
						
						v[0] *= sign_val;
						v[1] *= sign_val;
						v[2] *= sign_val;
						
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
			
			return v[0]*term[0] + v[1]*term[1] + v[2]*term[2];
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
			
			const Real w_common_term = 1./(12*spacing)*
			(-Dplus<2>(f, ix, iy, iz-2) + 7*Dplus<2>(f,ix, iy, iz-1) + 7*Dplus<2>(f,ix, iy, iz) -Dplus<2>(f, ix, iy, iz+1) );
			
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
			
			const Real phi_weno_z[2] = {
				phi_WENO(DplusDminus<2>(f, ix, iy, iz-2)/spacing,	DplusDminus<2>(f, ix, iy, iz-1)/spacing, 
						 DplusDminus<2>(f, ix, iy, iz)/spacing,		DplusDminus<2>(f, ix, iy, iz+1)/spacing),
				phi_WENO(DplusDminus<2>(f, ix, iy, iz+2)/spacing,	DplusDminus<2>(f, ix, iy, iz+1)/spacing, 
						 DplusDminus<2>(f, ix, iy, iz)/spacing,		DplusDminus<2>(f, ix, iy, iz-1)/spacing)
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
				3./((eps + IS[2])*(eps + IS[2])),
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
				
				return s*(sqrt(term1*term1 + term2*term2 + term3*term3)-1);
			}
			else
			{
				const Real term1 = max( -min(u_minus, (Real)0.0), max(u_plus, (Real)0.0) );
				const Real term2 = max( -min(v_minus, (Real)0.0), max(v_plus, (Real)0.0) );
				const Real term3 = max( -min(w_minus, (Real)0.0), max(w_plus, (Real)0.0) );
				
				return s*(sqrt(term1*term1 + term2*term2 + term3*term3)-1);
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
			
			const Real w_common_term = 1./(12*spacing)*
			(-Dplus<2>(f, ix, iy, iz-2) + 7*Dplus<2>(f,ix, iy, iz-1) + 7*Dplus<2>(f,ix, iy, iz) -Dplus<2>(f, ix, iy, iz+1) );
			
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
			
			const Real phi_weno_z[2] = {
				phi_WENO(DplusDminus<2>(f, ix, iy, iz-2)/spacing,	DplusDminus<2>(f, ix, iy, iz-1)/spacing, 
						 DplusDminus<2>(f, ix, iy, iz)/spacing,		DplusDminus<2>(f, ix, iy, iz+1)/spacing),
				phi_WENO(DplusDminus<2>(f, ix, iy, iz+2)/spacing,	DplusDminus<2>(f, ix, iy, iz+1)/spacing, 
						 DplusDminus<2>(f, ix, iy, iz)/spacing,		DplusDminus<2>(f, ix, iy, iz-1)/spacing)
			};
			
			const Real u_minus = u_common_term - phi_weno_x[0];
			const Real u_plus = u_common_term + phi_weno_x[1];
			const Real v_minus = v_common_term - phi_weno_y[0];
			const Real v_plus = v_common_term + phi_weno_y[1];
			const Real w_minus = w_common_term - phi_weno_z[0];
			const Real w_plus = w_common_term + phi_weno_z[1];
			
			return -H_GODUNOV(spacing, f(ix,iy).rho_0, u_plus, u_minus, v_plus, v_minus, w_plus, w_minus);
		}
		
	public:
		int stencil_start[3], stencil_end[3];
		
		BP_ComputeRHS_WENO5_REINIT()
		{
			stencil_start[0] = -3;
			stencil_start[1] = -3;
			stencil_start[2] = -3;
			
			stencil_end[0] = 4;
			stencil_end[1] = 4;
			stencil_end[2] = 4;
		}
		
		BP_ComputeRHS_WENO5_REINIT(const BP_ComputeRHS_WENO5_REINIT& c)
		{
			stencil_start[0] = -3;
			stencil_start[1] = -3;
			stencil_start[2] = -3;
			
			stencil_end[0] = 4;
			stencil_end[1] = 4;
			stencil_end[2] = 4;
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
	
	template <typename Wavelets>
	class STDTestL4_LevelsetAdvection //LEVEL 4: correct 3D simulations, correct data-processing
	{
	public:
		static const bool bRestartFromFile = false;
		static const bool bPostProcessData = false;
		static const bool bUseSpaceTimeSorter = true;
		static const int nStepsPerSerialization = 4;
		static const int blocksPerDimension = 1;
		static const int blockSize = 26;
		static const int maxLevel = 6;
		static const int resJump = 2;
		static const int narrowBandWidth = 6;
		
		float toleranceIC;
		
		typedef SimpleLevelsetBlock< PLS, narrowBandWidth, blockSize, blockSize, blockSize> B;
		//typedef Wavelets_Interp2ndOrder W;
		typedef Wavelets_Interp4thOrder W;
		Grid<W, B> grid;
		
#ifndef _MRAG_TBB
		typedef Multithreading::BlockProcessing_SingleCPU<B> BlockProcessing;
#else
		typedef Multithreading::BlockProcessing_TBB<B> BlockProcessing;
//		typedef Multithreading::BlockProcessing_SingleCPU<B> BlockProcessing;
		//typedef Multithreading::BlockProcessing_Pipeline_TBB<B, BlockLab, 6> BlockProcessing;
#endif
		
	protected:
		Refiner_SpaceExtension refiner;
		BlockProcessing block_processing;
		Compressor compressor;
		BlockFWT<W, B, levelset_projector> blockfwt;
#ifdef _MRAG_GLUT_VIZ
		GridViewer viewer;
#endif
		SpaceTimeSorter stSorter;
		Profiler profiler;
		string m_sFormat, m_sToRenderFormat;
		int m_iCurrentDumpTime;
		float t;

#pragma mark -
#pragma mark Class Management
		
	public:
		
		STDTestL4_LevelsetAdvection(float tol= 1e-3):
#ifdef _MRAG_GLUT_VIZ
			viewer(true,true),
#endif
			grid(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil), 
			refiner(resJump), compressor(resJump), profiler(),toleranceIC(tol), block_processing(),
			blockfwt(), t(0), m_iCurrentDumpTime(0),
			stSorter(),
			m_sToRenderFormat("Grid_At_Time%03d.%s"), 
			m_sFormat("SERIALIZED_%03d")
		{
			grid.setCompressor(&compressor);
			grid.setRefiner(&refiner);
			
			stSorter.connect(grid);
			
			if (bRestartFromFile)
			{
				_RestartFromFile("SERIALIZED_116",116, t, grid, m_iCurrentDumpTime);
			}
			else if (!bPostProcessData)
			{
				_ic(grid);
			//	s(Grid& g, BlockFWT& fwt, const double dAbsoluteTolerance,  const int iMaxLevel = -1, const bool bKeepHUpdated = true, const int iMaxLoops=-1,
			//	  MRAG::Profiler* profiler=NULL, double * dMaxDetailAlive= NULL, void (*fillGrid)(Grid& g)=NULL)
				printf("*/*/*/*/*/*/ Initializing!!! %e\n", tol);
				Science::AutomaticRefinementForLevelsets(grid, blockfwt, tol, maxLevel, true, -1,  NULL, NULL, _ic);
				_ic(grid);
				Science::AutomaticCompressionForLevelsets(grid, blockfwt, tol);
				_ic(grid);
				
			}
		}
		
	protected:
		
		void _RestartFromFile(string sFilePath, int stepid, float& output_time, Grid<W,B>& output_grid, int& out_currentDumpTime) const
		{
			//1. compute the actual time t
			//2. fill the grid from the restart file
			//3. fill the dump time
			
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
		}
		
	public:
		
		static void idle(void)
		{
#ifdef	_MRAG_GLUT_VIZ
			glutPostRedisplay();
#endif
		}

		static STDTestL4_LevelsetAdvection<Wavelets> * singletone; 
		
		static void runTests(int argc, char ** argv, bool bVisual)
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
			

#ifdef _MRAG_GLUT_VIZ			
			if (bVisual) 
			{
				
				singletone = new STDTestL4_LevelsetAdvection<Wavelets>();
				glutMainLoop();
			}
			else
#endif				
			{
				int nRuns = 10;
				float tols[10] = {1e-1, 5e-2, 1e-2,  5e-3, 1e-3,  5e-4, 1e-4, 5e-5, 1e-5, 5e-6};
							
				for(int iRun=0; iRun<nRuns; iRun++)
				{
					
					singletone = new STDTestL4_LevelsetAdvection<Wavelets>(tols[iRun]);
						
					
					
					bool bStop = false;
					do
					{
						bStop = singletone->run();
					}
					while(!bStop);
					
					//singletone->RefineOnceEverywhere(resJump);
					//double error = singletone->_computeError();
					
					delete singletone;
					singletone = NULL;
					

				}
			}
		}

	
#pragma mark -
#pragma mark Computing
		
	protected:
		
		float _ic_func4(float x[3]) const
		{
			const float r = sqrt(pow(x[0]-0.35,2) + pow(x[1]-0.35, 2) + pow(x[2]-0.35,2));
			return 0.15-r;
		}
		
		static float _ic_func3(float x[3])
		{
			const float r = sqrt(pow(x[0]-0.5,2) + pow(x[1]-0.5, 2) + pow(x[2]-0.5,2));
			return 0.15-r;
		}
		
		template<typename Grid>
		static void _ic(Grid& grid) 
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
							
							//block(ix,iy,iz).rho = _ic_func4(x);
							block(ix,iy,iz).rho = _ic_func3(x);
						}
				
				grid.getBlockCollection().release(info.blockID);
			}
		}
		
		
		static inline double _eta(double r_, double eps) 
		{
			const double r = r_/eps;
			//return 1/(2*M_PI*eps)*(4-r*r)*exp(-r*r);
			return 1/(eps*sqrt(M_PI))*exp(-r*r); 
		}
		
		double getAverageResolutionSize()
		{
			return 1.0/sqrt(1.0*grid.getBlocksInfo().size()*B::sizeX*B::sizeY);
		}
		
		//static double _computeError(Grid<W, B>& grid, double t)
		double _computeError()
		{
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			double dTotalError = 0.0;
			double avgResolution = getAverageResolutionSize();	
			
			const double smallest_h = pow(2.0,-grid.getCurrentMaxLevel());
			const double smallest_scale = smallest_h*sqrt(2.);
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				
				double dLocalError = 0.0;
				
				for(int iz=0; iz<B::sizeZ; iz++)
					for(int iy=0; iy<B::sizeY; iy++)
						for(int ix=0; ix<B::sizeX; ix++)
						{
							float x[3];

							info.pos(x, ix, iy, iz);
							
							const double dError = _eta(block(ix,iy,iz).rho, smallest_scale)*fabs(block(ix,iy,iz).rho - _ic_func3(x));
							
							dLocalError += dError;
						}
				
				dTotalError += dLocalError*info.h[0]*info.h[1]*info.h[2];
			}
			
			printf("dTotalError = %f\n", dTotalError);
			
			static int counter = 0;
			FILE * f = fopen("Error.dat", counter++==0?"w":"a");
			fprintf(f,"%e\t%e\t%e\n", toleranceIC, avgResolution, dTotalError);
			fclose(f);
			

			
			if (counter>=2 && counter %1 == 0)
			{
				
				FILE * fgnuplot = popen("/usr/local/bin/gnuplot -persist\n", "w");
				assert(fgnuplot != NULL);
				
				fprintf(fgnuplot, "set title \"Error vs time (%d points)\"\n", vInfo.size()*B::sizeX*B::sizeY*B::sizeZ);
				fprintf(fgnuplot, "set xlabel \" time \" \n");
				fprintf(fgnuplot, "set ylabel \"error\" \n");
				fprintf(fgnuplot, "set grid\n");
				fprintf(fgnuplot, "set logscale\n");
				fprintf(fgnuplot, "plot \"Error.dat\" using 1:2  title 'error of the computed solution' with linespoint\n");
				fprintf(fgnuplot, "set terminal postscript color\n");
				fprintf(fgnuplot, "set out \"Error.ps\"\n");
				fprintf(fgnuplot, "replot\n");
				
				fclose(fgnuplot);
			}
			
			return dTotalError;
		}
		

		double computeDeltaTime() const
		{
			const double dx = 1.0/blockSize;
			const double CFL = 4.0;
			const double dt = dx/CFL;//8e-1;
			
			return dt;
		}

		bool _step(float t, float dt)
		{
			const int nParallelGranularity = (grid.getBlocksInfo().size()<8 ? 1 : 4);
			BoundaryInfo* boundaryInfo = &grid.getBoundaryInfo();
			
			stSorter.startSession(dt, 2, 0);
			
			if(t>4.) 
			{				
				RefineOnceEverywhere(0);
				double error = _computeError();
				return true;
			}			
			
						
			//advection here
			while(true)
			{
				double currTime, currDeltaT;
				int level;
				SpaceTimeSorter::ETimeInterval type;
				vector< BlockInfo> vInfo;
				
				const bool bContinue = stSorter.getBlocks(level, currDeltaT, currTime, vInfo, type);
				
				if (type == SpaceTimeSorter::ETimeInterval_Start)
				{
					BP_ComputeRHS_TR_WENO5 task1(currTime, (t+currTime<=4.)? 1 : -1);
					BP_UpdateRHS_TR task2(currTime);
					
					BlockProcessing::process< BlockLab >(vInfo, grid.getBlockCollection(), *boundaryInfo, task1, nParallelGranularity);
					//block_processing.pipeline_process(vInfo, grid.getBlockCollection(), *boundaryInfo, task1);
					BlockProcessing::process(vInfo, grid.getBlockCollection(), task2, nParallelGranularity);
				}
				else
				{
					BP_Integrate_TR task(currTime);
					BlockProcessing::process(vInfo, grid.getBlockCollection(), task, nParallelGranularity);
				}
				
				if (!bContinue) break;
			}
			stSorter.endSession();
			
			
			// reinitialize here
			{
				//double time=0.0;
	
				BP_ComputeRHS_WENO5_REINIT computeRhs;
				BP_Integrate_REINIT integrate(5e-4);
				
				vector<BlockInfo> vInfo = grid.getBlocksInfo();
				
				for(int i=0; i<20; i++)
				{
					BlockProcessing::process<BlockLab>(vInfo, grid.getBlockCollection(), grid.getBoundaryInfo(), computeRhs, nParallelGranularity);
					//blockProcessing.pipeline_process(vInfo, grid.getBlockCollection(), *boundaryInfo, computeRhs);
					BlockProcessing::process(vInfo, grid.getBlockCollection(),  integrate, nParallelGranularity);
				}
				
				cout<<"Done re-init"<<endl;
			}
			
			
			return false;
		}
		
		bool run()
		{
			if (!bPostProcessData)
			{
				const int nStep = 1;
				const float dt = computeDeltaTime();
				bool bStop;
				
				//soft refinement for memory reasons
				{					
					float memSize = grid.getMemorySize();
					while(memSize>300.0)
					{
						toleranceIC *= 2.0;

						Science::AutomaticCompressionForLevelsets(grid, blockfwt, toleranceIC, true, -1, &profiler); 

						memSize = grid.getMemorySize();
					}
				
					Science::AutomaticRefinementForLevelsets(grid, blockfwt, toleranceIC, maxLevel, true, -1, &profiler);
				}
				
				cout<<"I am at time= "<<t<<" and my tolerance is "<< toleranceIC <<endl;
				printf("/////////// WHILE COMPUTING MEM MB: %f  (%d blocks, BI size=%f MB) ////////////////// \n", grid.getMemorySize(), grid.getBlocksInfo().size(), grid.getBoundaryInfo().getMemorySize());
				
				profiler.getAgent("computing").start();
				
				for(int i=0; i<nStep; i++, t+=dt)
					bStop = _step(t, dt);
				
				profiler.getAgent("computing").stop();
				
				Science::AutomaticCompressionForLevelsets(grid, blockfwt, toleranceIC, true, -1, &profiler);  
				
				printf("/////////// AFTER COMPUTING MEM MB: %f  (%d blocks, BI size=%f MB) ////////////////// \n", grid.getMemorySize(), grid.getBlocksInfo().size(), grid.getBoundaryInfo().getMemorySize());
				
				//serialize the data
				static int serialize_counter = 0;
				if (serialize_counter++ % nStepsPerSerialization == 0)
				{
					profiler.printSummary();

					IO_Native<W, B, levelset_projector> serializer;
					char buf[300];
					sprintf(buf, m_sFormat.data(), m_iCurrentDumpTime);
					//serializer.Write(grid, buf);	
					m_iCurrentDumpTime++;
				}
				
				return bStop;
			}
			else
			{
				IO_Native<W, B, levelset_projector> serializer;
				
				const int nGrids = 234;
				for(int iGrid=0; iGrid<=nGrids; iGrid++)
				{
					printf("POSTPROCESS PHASE \n");
					char buf[300];
					sprintf(buf, m_sFormat.data(), iGrid);
					serializer.Read(grid, buf);
					
					DumpData();
				}
				
				return true;
			}
		}
		
#pragma mark -
#pragma mark Postprocessing
		
		void RefineOnceEverywhere(int maxResJump)
		{
			for(int i=0; i<=maxResJump; i++)
			{
				blockfwt.prepare(grid.getBlockCollection(), grid.getBoundaryInfo());
				vector<BlockInfo> vInfo = grid.getBlocksInfo();
				
				set<int> shouldBeRefined;
				for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++)
				{
					blockfwt.fwt<0>(*it);
					
					if (blockfwt.getReport().getOverAll_DetailMaxMag()>0)
						shouldBeRefined.insert(it->blockID);
				}
				
				grid.refine(shouldBeRefined);
			}
			
			
			/*	{
			 vector<BlockInfo> vInfo = grid.getBlocksInfo();
			 blockfwt.prepare(grid.getBlockCollection(), grid.getBoundaryInfo());
			 set<int> shouldBeRefined;
			 for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++)
			 {
			 blockfwt.fwt<0>(*it);
			 assert(blockfwt.getReport().getOverAll_DetailMaxMag() < 1e-15);
			 }
			 }*/
			
		}
		
		
		
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

