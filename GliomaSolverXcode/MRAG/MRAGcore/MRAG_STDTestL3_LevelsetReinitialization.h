/*
 *  MRAG_STDTestL3_LevelsetReinitialization.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/2/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */



#ifndef _MRAG_GLUT_VIZ
#error Please do not include this file if you are not using glut, or if you do not define _MRAG_GLUT_VIZ. Diego.
#else
#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGEnvironment.h"

#if _MRAG_OS == _MRAG_OS_APPLE
#include "GLUT/glut.h"
#elif _MRAG_OS == _MRAG_OS_WINDOWS
#define _USE_MATH_DEFINES
#include "GL/glew.h"
#include "GL/glut.h"
#endif

#undef min
#undef max

#include "MRAGcore/MRAGrid.h"
#include "MRAGcore/MRAGRefiner.h"
#include "MRAGcore/MRAGCompressor.h"
#include "MRAGcore/MRAGBlockLab.h"
#include "MRAGcore/MRAGBlockFWT.h"

#include "MRAGvisual/GridViewer.h"

#pragma once

#include "MRAGscience/MRAGScienceCore.h"
#include "MRAGscience/MRAGSimpleLevelsetBlock.h"
#include "MRAGmultithreading/MRAGBlockProcessing_SingleCPU.h"
#include "MRAGmultithreading/MRAGBlockProcessing_TBB.h"


struct LevelsetPoint
{
	Real phi_0, phi, dphidt;
	Real tmp;
	
	void operator += (LevelsetPoint t)
	{
		phi_0 += t.phi_0;
		phi += t.phi;
		dphidt += t.dphidt;
	}

	operator Real()  { return phi_0; }
	
	Real levelset() const { return phi_0; }
};

LevelsetPoint operator*(const LevelsetPoint& p, Real v)
{
	LevelsetPoint t = {p.phi_0*v, p.phi*v, p.dphidt*v};
	return t;
}

template <typename T, int i> inline Real lvlstinit_projector_impl(const LevelsetPoint&t)
{
	return (Real)1e-3/max(fabs(t.phi), 1e-3);
}

make_projector(levelset_projector, lvlstinit_projector_impl)


RGBA convertToRGBA( LevelsetPoint& LevelsetPoint)
{
	RGBA c;

	double mag = fabs(LevelsetPoint.phi);
	
	c.g = mag*max((double)LevelsetPoint.phi*1000.0, 0.0);
	c.r = -mag*min((double)LevelsetPoint.phi*1000.0, 0.0);
	
	c.b =  0;//LevelsetPoint.phi>0?0.0:0.5;
	return c;
}


class ComputeRHS
{
	inline const double phi_WENO(const double a, const double b, const double c, const double d) const
	{
		const double eps = 1e-6;
		
		const double IS[3] = {
			13*(a-b)*(a-b) + 3*(a-3*b)*(a-3*b),
			13*(b-c)*(b-c) + 3*(b+c)*(b+c),
			13*(c-d)*(c-d) + 3*(3*c-d)*(3*c-d)
		};
		
		const double alpha[3] = {
			1./((eps + IS[0])*(eps + IS[0])),
			6./((eps + IS[1])*(eps + IS[1])),
			3./((eps + IS[2])*(eps + IS[2])),
		};
		
		const double sum_alpha = alpha[0] + alpha[1] + alpha[2];
		
		const double w[2] = {
			alpha[0]/sum_alpha,
			alpha[2]/sum_alpha
		};
		
		return 1./3*w[0]*(a-2*b+c) + 1./6*(w[1]-0.5)*(b-2*c+d);
	}
	
	inline const double compute_sign1(const double phi, const double h, const double grad_phi_mag) const
	{
		return phi/sqrt(phi*phi + h*h);
	}
	
	inline const double compute_sign2(const double phi, const double h, const double grad_phi_mag) const
	{
		return phi/sqrt(phi*phi + h*h*grad_phi_mag*grad_phi_mag);
	}
	
	inline const double H_h(const double phi, const double h) const
	{
		if (phi<-h) return 0;
		if (phi>h) return 1;
		
		return 0.5*(1+phi/h+1/M_PI*sin(M_PI*phi/h));
	}
	
	inline const double compute_sign3(const double phi, const double h, const double grad_phi_mag) const
	{
		return 2*H_h(phi, h) - 1;
	}
	
	inline const double H_GODUNOV(const double spacing, const double phi, const double u_plus, const double u_minus, const double v_plus, const double v_minus) const
	{
		if (phi>=0)
		{
			const double term1 = max( -min(u_plus, 0.0), max(u_minus, 0.0) );
			const double term2 = max( -min(v_plus, 0.0), max(v_minus, 0.0) );
			
			const double grad_phi_mag = sqrt(term1*term1 + term2*term2);
			const double s = compute_sign2(phi, spacing, grad_phi_mag);
			return s*(grad_phi_mag-1);
		}
		else
		{
			const double term1 = max( -min(u_minus, 0.0), max(u_plus, 0.0) );
			const double term2 = max( -min(v_minus, 0.0), max(v_plus, 0.0) );
			const double grad_phi_mag = sqrt(term1*term1 + term2*term2);
			const double s = compute_sign2(phi, spacing, grad_phi_mag);
			return s*(grad_phi_mag-1);

		}
	}
	
	template <int index, typename Field> inline const double Dplus(Field& f, const int ix, const int iy) const
	{ 
		return (f(ix + (int)(index==0), iy + (int)(index==1),0).phi - f(ix, iy,0).phi);
	}
	
	template <int index, typename Field> inline const double Dminus(Field& f, const int ix, const int iy) const
	{ 
		return (f(ix, iy,0).phi - f(ix - (int)(index==0), iy - (int)(index==1),0).phi );
	}
	
	template <int index, typename Field> inline const double DplusDminus(Field& f, const int ix, const int iy) const
	{ 
		return (f(ix + (int)(index==0), iy + (int)(index==1),0).phi +
				f(ix - (int)(index==0), iy - (int)(index==1),0).phi  
				-2*f(ix, iy,0).phi );
	}
	
	template <typename Field>
	const double mainRHS(Field& f, const double spacing, const int ix, const int iy) const
	{
		const double u_common_term = 1./(12*spacing)*
			(-Dplus<0>(f, ix-2, iy) + 7*Dplus<0>(f,ix-1, iy) + 7*Dplus<0>(f,ix, iy) -Dplus<0>(f, ix+1, iy) );
		
		const double v_common_term = 1./(12*spacing)*
			(-Dplus<1>(f, ix, iy-2) + 7*Dplus<1>(f,ix, iy-1) + 7*Dplus<1>(f,ix, iy) -Dplus<1>(f, ix, iy+1) );
		
		const double phi_weno_x[2] = {
			phi_WENO(DplusDminus<0>(f, ix-2, iy)/spacing,	DplusDminus<0>(f, ix-1, iy)/spacing, 
					 DplusDminus<0>(f, ix, iy)/spacing,		DplusDminus<0>(f, ix+1, iy)/spacing),
			phi_WENO(DplusDminus<0>(f, ix+2, iy)/spacing,	DplusDminus<0>(f, ix+1, iy)/spacing, 
					 DplusDminus<0>(f, ix, iy)/spacing,		DplusDminus<0>(f, ix-1, iy)/spacing)
		};
		
		const double phi_weno_y[2] = {
			phi_WENO(DplusDminus<1>(f, ix, iy-2)/spacing,	DplusDminus<1>(f, ix,iy-1)/spacing, 
					 DplusDminus<1>(f, ix, iy)/spacing,		DplusDminus<1>(f, ix, iy+1)/spacing),
			phi_WENO(DplusDminus<1>(f, ix, iy+2)/spacing,	DplusDminus<1>(f, ix, iy+1)/spacing, 
					 DplusDminus<1>(f, ix, iy)/spacing,		DplusDminus<1>(f, ix, iy-1)/spacing)
		};
		
		const double u_minus = u_common_term - phi_weno_x[0];
		const double u_plus = u_common_term + phi_weno_x[1];
		const double v_minus = v_common_term - phi_weno_y[0];
		const double v_plus = v_common_term + phi_weno_y[1];
	
		return -H_GODUNOV(spacing, f(ix,iy).phi_0, u_plus, u_minus, v_plus, v_minus);
	}
	
public:
	int stencil_start[3], stencil_end[3];
	
	ComputeRHS()
	{
		stencil_start[0] = -3;
		stencil_start[1] = -3;
		stencil_start[2] = 0;
		
		stencil_end[0] = 4;
		stencil_end[1] = 4;
		stencil_end[2] = 1;
	}
	
	ComputeRHS(const ComputeRHS& c)
	{
		stencil_start[0] = -3;
		stencil_start[1] = -3;
		stencil_start[2] = 0;
		
		stencil_end[0] = 4;
		stencil_end[1] = 4;
		stencil_end[2] = 1;
	}
	
	template<typename LabType, typename BlockType>
	inline void operator()(LabType& i, const BlockInfo& info, BlockType& o) const
	{
		typedef BlockType B;
		typedef typename BlockType::ElementType E;
		
		const double h = info.h[0];
		
		for(int iz=0; iz<B::sizeZ; iz++)
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
					o(ix,iy,iz).dphidt = mainRHS(i, h , ix, iy);
	}
	
};

class Integrate
{
public:
	double dt;
	Integrate(double dt_): dt(dt_) {}
	Integrate(const Integrate& i): dt(i.dt){}
	
	template <typename BlockType>
	inline void operator() (const BlockInfo& info, BlockType& b) const
	{
		typedef BlockType B;
				
		const int n = B::sizeZ*B::sizeY*B::sizeX;
		
		LevelsetPoint* ptrE = &(b(0));
		for(int iE=0; iE<n; iE++, ptrE++)
		{
			ptrE->phi += dt*ptrE->dphidt;
			ptrE->dphidt = 0;
		}		
	}
};

namespace MRAG
{
	template <typename Wavelets>
	class STDTestL3_LevelsetReinitialization //LEVEL 3: correct 2D simulations, correct data-processing
	{
	protected:
		typedef SimpleLevelsetBlock<LevelsetPoint,4,32,32,1> BlockType;
	
		double time;
		Grid<Wavelets, BlockType> grid;
		Refiner refiner;
		Compressor compressor;
		BlockLab<BlockType> lab;
		BlockFWT<Wavelets, BlockType, levelset_projector> blockfwt;
		GridViewer viewer;
#ifdef _MRAG_TBB
		typedef Multithreading::BlockProcessing_TBB<BlockType> BlockProcessing;
#else
		typedef Multithreading::BlockProcessing_SingleCPU<BlockType> BlockProcessing;
#endif
		
		STDTestL3_LevelsetReinitialization():
		grid(8,8), refiner(), compressor(), lab(), 
		blockfwt(), viewer(true, true), time(0)
		{
			grid.setCompressor(&compressor);
			grid.setRefiner(&refiner);
			
			_ic(grid);
			Science::AutomaticRefinement<0,0>(grid, blockfwt, 0.01, 5, -1, NULL, _ic);
		}
		
		template<int dir>
		static inline const double HeavySide(double x0, double x)
		{
			const double eps = 0.1/2;
			
			/*const double alpha = M_PI*min(1., max(0., (x-x0+0.5*eps)/eps));
			if (dir==1) return 0.5+0.5*sin(alpha - M_PI/2);
			else return 0.5+0.5*cos(alpha);*/
			const double alpha = min(1., max(0., (x-x0+0.5*eps)/eps));
			if (dir==1) return alpha;
			else return 1-alpha;
		}
		
		template <bool bPerturbation>
		static float _ic_func(float x[3])
		{
			/*const double r[2] = {x[0]-  0.5, x[1]-0.5};
			
			const double radius = 0.25;
			
			const double d = sqrt(r[0]*r[0] + r[1]*r[1]) - radius;
			
			const double epsilon = 0.2;
			
			if (!bPerturbation || fabs(d)>epsilon) return d;
			
			const double teta = atan2(x[1], x[0]);
			
			const double perturbation = 4*epsilon/(16*M_PI)*sin(4*M_PI*d*sin(5*teta)/epsilon);
			
			return d+perturbation;*/
			const double hw = 0.1;
			if (bPerturbation)
			{
			
			const double val =  HeavySide<1>(0.5-hw,x[0]) * HeavySide<-1>(0.5+hw,x[0]) *
			HeavySide<1>(0.5-hw,x[1]) * HeavySide<-1>(0.5+hw,x[1]); 
			
			return val * (-2.) + 1.;
			}
			else
			{
				double p[2] = {
					fabs(x[0]-0.5), fabs(x[1]-0.5)
				};
				
				int c[2] = {
					(int)(p[0]>=hw), (int)(p[1]>=hw)
				};
				
				const bool bInside = (c[0]+c[1] == 0);
				const double dInside = -min(hw-p[0], hw-p[1]);
				const double dOutside = sqrt(c[0]*pow(hw-p[0],2) + c[1]*pow(hw-p[1],2));
				return bInside? dInside:dOutside;
			}
		}
		
		/*
		template<int dir>
		static inline const double HeavySide(double x0, double x)
		{
			//if (dir == 1) return (int)(x>x0);
			//return (int)(x<x0);
			const double eps = 0.1/10.;
			
			const double alpha = M_PI*min(1., max(0., (x-x0+0.5*eps)/eps));
			if (dir==1) return 0.5+0.5*sin(alpha - M_PI/2);
			else return 0.5+0.5*cos(alpha);
		}*/
		
		
		static void _ic(Grid<Wavelets, BlockType>& grid)
		{
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				BlockType& block = grid.getBlockCollection()[info.blockID];
				
				for(int iz=0; iz<BlockType::sizeZ; iz++)
					for(int iy=0; iy<BlockType::sizeY; iy++)
						for(int ix=0; ix<BlockType::sizeX; ix++)
						{
							float x[3];
							
							info.pos(x, ix, iy, iz);
							
							//block(ix,iy,iz).phi = block(ix,iy,iz).phi_0 = _ic_func<true>(x);
							
							
							const double hwidth = 0.1;
							
							const double val =  HeavySide<1>(0.5-hwidth,x[0]) * HeavySide<-1>(0.5+hwidth,x[0]) *
												HeavySide<1>(0.5-hwidth,x[1]) * HeavySide<-1>(0.5+hwidth,x[1]);												
							
							block(ix,iy,iz).phi = block(ix,iy,iz).phi_0 = val * (-2.) + 1.;
							
						}
			}
		}
		
		static inline double _eta(double r_, double eps) 
		{
			const double r = r_/eps;
			return 1/(2*M_PI*eps)*(4-r*r)*exp(-r*r); 
		}
		
		static double _computeError(Grid<Wavelets, BlockType>& grid, double time)
		{
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			double dTotalError = 0.0;
			
			const double smallest_h = pow(2.0,-grid.getCurrentMaxLevel());
			
			const double smallest_scale = smallest_h*sqrt(2.);
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				BlockType& block = grid.getBlockCollection()[info.blockID];
				
				double dLocalError = 0.0;
				for(int iz=0; iz<BlockType::sizeZ; iz++)
					for(int iy=0; iy<BlockType::sizeY; iy++)
						for(int ix=0; ix<BlockType::sizeX; ix++)
						{
							float x[3];
							
							info.pos(x, ix, iy, iz);

							const double dError = _eta(block(ix,iy,iz).phi, smallest_scale)*fabs(block(ix,iy,iz).phi - _ic_func<false>(x));
							
							dLocalError += dError;
						}
				
				dTotalError += dLocalError*info.h[0]*info.h[1]*info.h[2];
			}
			
			printf("dTotalError = %f\n", dTotalError);
				
			static int counter = 0;
			FILE * f = fopen("Error.dat", counter++==0?"w":"a");
			fprintf(f,"%e\t%e\n", time, dTotalError);
			fclose(f);
			
			if (counter>=2 && counter %1 == 0)
			{
				/*FILE * fgnuplot = popen("/usr/local/bin/gnuplot -persist\n", "w");
				assert(fgnuplot != NULL);
				
				fprintf(fgnuplot, "set title \"Error vs time (%d points)\"\n", vInfo.size()*BlockType::sizeX*BlockType::sizeY*BlockType::sizeZ);
				fprintf(fgnuplot, "set xlabel \" time \" \n");
				fprintf(fgnuplot, "set ylabel \"error\" \n");
				fprintf(fgnuplot, "set grid\n");
				fprintf(fgnuplot, "set logscale\n");
				fprintf(fgnuplot, "plot \"Error.dat\" using 1:2  title 'error of the computed solution' with linespoint\n");
				fprintf(fgnuplot, "set terminal postscript color\n");
				fprintf(fgnuplot, "set out \"Error.ps\"\n");
				fprintf(fgnuplot, "replot\n");
				
				fclose(fgnuplot);*/
			}
			
			return dTotalError;
		}
		
		static void idle(void)
		{
			glutPostRedisplay();
		}
		
		bool run()
		{
			if (time==0)
				_computeError(grid, time);
			
			ComputeRHS computeRhs;
			Integrate integrate(1e-4);
			
			BoundaryInfo* boundaryInfo = grid.createBoundaryInfo(computeRhs.stencil_start, computeRhs.stencil_end);
			
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			
			for(int i=0; i<100; i++)
			{
				BlockProcessing::process<BlockLab>(vInfo, grid.getBlockCollection(), *boundaryInfo, computeRhs);
				BlockProcessing::process(vInfo, grid.getBlockCollection(),  integrate);
				time += integrate.dt;
			}
				
			const double error = _computeError(grid, time);
			
			delete boundaryInfo;
			
			//if (error<0.00002) exit(0);
			
			return true;
		}
		
		static void display(void)
		{
			glClear(GL_COLOR_BUFFER_BIT);
			
			singletone->run();
			singletone->viewer.drawContent(singletone->grid, singletone->grid.getBlockCollection());
			singletone->viewer.drawSketch(singletone->grid,false);
			
			glutSwapBuffers();
		}
		
	public:
		
		static STDTestL3_LevelsetReinitialization<Wavelets> * singletone; 
		
		static void runTests(int argc, char ** argv, bool bVisual=true)
		{
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
			
			singletone = new STDTestL3_LevelsetReinitialization<Wavelets>();
			
			if (bVisual) glutMainLoop();
			else
			{
				for(int i=0; i<300; i++)
					singletone->run();
			}
		}
	};
	
	template <typename Wavelets>
	STDTestL3_LevelsetReinitialization<Wavelets> * STDTestL3_LevelsetReinitialization<Wavelets>::singletone; 
}

#endif