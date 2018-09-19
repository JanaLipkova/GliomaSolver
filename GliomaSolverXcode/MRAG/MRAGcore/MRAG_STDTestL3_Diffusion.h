/*
 *  MRAG_STDTestL3_Diffusion.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/29/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include "../MRAGcore/MRAGCommon.h"
#include "../MRAGcore/MRAGEnvironment.h"

#if _MRAG_OS == _MRAG_OS_APPLE
#include "GLUT/glut.h"
#elif _MRAG_OS == _MRAG_OS_WINDOWS
#define _USE_MATH_DEFINES
#include "GL/glew.h"
#include "GL/glut.h"
#endif

#undef min
#undef max

#include "../MRAGcore/MRAGWavelets_Interp2ndOrder.h"
#include "../MRAGcore/MRAGWavelets_Interp4thOrder.h"
#include "../MRAGcore/MRAGBlock.h"
#include "../MRAGcore/MRAGrid.h"
#include "../MRAGcore/MRAGBlockFWT.h"
#include "../MRAGscience/MRAGScienceCore.h"
#include "../MRAGscience/MRAGSpaceTimeSorter.h"
#include "../MRAGvisual/GridViewer.h"
#include "../MRAGvisual/MRAGVisualTypes.h"

using namespace MRAG;

const Real cazzo_viscosity = 0.05; 

struct P
{
	Real rho;
	Real drho_dt;
	Real tmp;
	
	P(): rho(0), drho_dt(0), tmp(0){}
	
	P(Real rho_, Real drho_dt_): rho(rho_), drho_dt(drho_dt_) {}
	
	void operator += (P t)
	{
		rho += t.rho;
		drho_dt += t.drho_dt;
	}
	
	operator Real()
	{
		return (Real)rho;
	}
	
	void partA(Real t, Real dt)
	{
		drho_dt = tmp;
		rho -= tmp*t;
	}
	
	void partC(Real t, Real dt)
	{
		rho = evaluate_rho(t-dt*0.5) + dt*0.5*tmp;
		drho_dt = 0;
	}
	
	Real evaluate_rho(Real t)
	{
		return rho + drho_dt*t;
	}
};

P operator*(const P& p, Real v)
{
	P t;
	
	t.rho = v*p.rho;
	t.drho_dt = v*p.drho_dt;

	return t;
}

RGBA convertToRGBA(P& p)
{
	RGBA c(max(0.0, 2.0*(p.rho+2.223)/4.445 - 0.5), 1.0 - 2.0*fabs((p.rho+2.223)/4.445 - 0.5),max(0.0, 1.0-2.0*(p.rho+2.223)/4.445),0);
	//RGBA c(max(0.0, 2.0*(p.rho+2.223)/4.445 - 0.5),  1.0 - 2.0*fabs((p.rho+2.223)/4.445 - 0.5),min(1.0, 0.0+2.0*(p.rho+2.223)/4.445),0);
	return c;
}

template <typename T, int i> inline Real diffusion_projector_impl(const T&t)
{
	return (Real)0.001/(0.1+t.rho);
}

make_projector(diffusion_projector, diffusion_projector_impl)

const Real beta = 2.0;

namespace MRAG
{
	typedef Wavelets_Interp4thOrder W;
	
	template <typename WaveletsAAA>
	class STDTestL3_Diffusion //LEVEL 3: correct 2D simulations, correct data-processing
	{
		static const bool bUseSpaceTimeSorter = true;
		static const bool b3D = false;
		static const int blocksPerDimension = 4;
		static const int blockSize = 20;
		static const int maxLevel = 3;
		static const int resJump = 1;	
		
		Real toleranceIC, mu, ti, beta2;
		
		typedef Block< P, blockSize, blockSize, b3D? blockSize: 1> B;
		
		Grid<W, B> grid;
		Refiner refiner;
		Compressor compressor;
		BlockLab<B> lab;
		BlockFWT<W, B, diffusion_projector> blockfwt;
		GridViewer viewer;
		SpaceTimeSorter stSorter;
		Real t;
		
#pragma mark IC
		static float _ic_circle(float x[3])
		{
			abort();
			const float r = sqrt(pow(x[0]-0.5,2) + pow(x[1]-0.5, 2));
			return r<0.2? 2: -2;
		}
		
		float _ic_square(float x[3])
		{
			abort();
			const float r = max(fabs(x[0]-0.5),fabs(x[1]-0.7));
			return r<0.2? 1: 0;
		}
		
		float _ic_gauss(float x[3])
		{
			abort();
			const float nu = 0.5;
			const float ro = 0.1;
			const float r2 = pow(x[0]-nu,2) + pow(x[1]-nu, 2);
			return 1.0/(ro*sqrt(2.0*M_PI))*exp(-r2/(3.0*pow(ro,2)))*10.0;
		}
		
		static Real _ic_test1(Real x[3],Real time)
		{
			//return sin(2.0*M_PI*x[0]*100);
			return (exp(-4.0*M_PI*M_PI*time)*sin(2.0*M_PI*x[0])*sin(2.0*M_PI*x[1])
					+
					exp(-4.0*M_PI*M_PI*time*beta*beta)*sin(2.0*M_PI*beta*x[0])*sin(2.0*beta*M_PI*x[1]));
		}
		
		static Real _ic_test2(Real x[3],Real time)
		{
			return (exp(-8.0*M_PI*M_PI*time)*(sin(2.0*M_PI*x[0])+sin(2.0*M_PI*x[1])));
		}
		
		static Real _ic_test3(Real x[3],Real time)
		{
			return (sin(2.0*M_PI*x[0])*sin(2.0*M_PI*x[1])*exp(-8*M_PI*M_PI*time));//exp(-time);//
		}
		
		template<typename Grid>
		static void _ic(Grid& grid)
		{
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				
				for(int iz=0; iz<B::sizeZ; iz++)
					for(int iy=0; iy<B::sizeY; iy++)
						for(int ix=0; ix<B::sizeX; ix++)
						{
							Real x[3];
							
							info.pos(x, ix, iy, iz);
							
							block(ix,iy,iz) = P();
							block(ix,iy,iz).rho = _ic_test3(x,0);
							//block(ix,iy,iz).rho = _ic_circle(x);							;
						}
			}			
			
		}

		
#pragma mark Simulation
		void _step_computeRHS(vector<BlockInfo>& vInfo)
		{
			const int steStart[3] ={ -1,-1,0};
			const int steEnd[3] ={ +2,+2,+1};
			
			lab.prepare(grid.getBlockCollection(), grid.getBoundaryInfo(),steStart,steEnd);
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				
				lab.load(info);
				
				const Real h = info.h[0];
				const Real factor = 1.0/(h*h);
				
				for(int iy=0; iy<B::sizeY; iy++)
					for(int ix=0; ix<B::sizeX; ix++) //-block(ix,iy).rho;
						block(ix,iy).tmp = -8*M_PI*M_PI*lab(ix,iy).rho*(1-cazzo_viscosity) + cazzo_viscosity*factor*(
												   lab(ix-1 ,iy).rho + 
												   lab(ix+1,iy).rho +
												   lab(ix,iy-1).rho + 
												   lab(ix,iy+1).rho - 4.0*lab(ix,iy).rho );
											
			}
		}

		
		void _step_computeRHS_TR(vector<BlockInfo>& vInfo,  double t, double dt)
		{
			const int steStart[3] ={ -2,-2,0};
			const int steEnd[3] ={ +3,+3,+1};
			
			lab.prepare(grid.getBlockCollection(), grid.getBoundaryInfo(),steStart,steEnd);
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				
				lab.load(info);
				const Real A = -1./12;
				const Real B = 16./12;
				const Real C = -30./12;
				
				const Real h = info.h[0];
				const Real factor = 1.0/(h*h);
				assert(mu==1);
				for(int iy=0; iy<B::sizeY; iy++)
					for(int ix=0; ix<B::sizeX; ix++)
						block(ix,iy).tmp = -8*M_PI*M_PI*lab(ix,iy).evaluate_rho(t)*(1-cazzo_viscosity) + 
						cazzo_viscosity*factor*(
							A*lab(ix-2 ,iy).evaluate_rho(t) +
							B*lab(ix-1 ,iy).evaluate_rho(t) + 
							B*lab(ix+1,iy).evaluate_rho(t) +
							A*lab(ix+2 ,iy).evaluate_rho(t) +
							A*lab(ix,iy-2).evaluate_rho(t) +
							B*lab(ix,iy-1).evaluate_rho(t) + 
							B*lab(ix,iy+1).evaluate_rho(t) +
							A*lab(ix,iy+2).evaluate_rho(t) +
							2*C*lab(ix,iy).evaluate_rho(t) );
			}
		}
		
		template <int mode>
		void _step_part(vector<BlockInfo>& vInfo, Real t, Real dt)
		{
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				
				switch(mode)
				{
					case 0:
						for(int iy=0; iy<B::sizeY; iy++)
							for(int ix=0; ix<B::sizeX; ix++)
								block(ix,iy).partA(t, dt);
						break;
					case 1:
						abort();
						break;
					case 2:
						for(int iy=0; iy<B::sizeY; iy++)
							for(int ix=0; ix<B::sizeX; ix++)
								block(ix,iy).partC(t, dt);
						break;
					default:
						abort();
				}
			}		
		}
		
		
		/*void _prepare()
		{
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				
				for(int iy=0; iy<B::sizeY; iy++)
					for(int ix=0; ix<B::sizeX; ix++)
						block(ix,iy).prepare();
						
			}		
		}*/
		
		void _step(Real t, Real dt)
		{
		/*	vector<BlockInfo> vInfo = grid.getBlocksInfo();
			
			_step_computeRHS(vInfo);
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				
				for(int iy=0; iy<B::sizeY; iy++)
					for(int ix=0; ix<B::sizeX; ix++)
						block(ix,iy).rho += dt*block(ix,iy).tmp;
			}	*/
			
			//
		//	_prepare();
			
			vector<BlockInfo> vInfo;
			
			double currTime, currDeltaT;
			int level;
			stSorter.startSession(dt, 4, 0);
			
			while(true)
			{
				SpaceTimeSorter::ETimeInterval type;
				const bool bContinue = stSorter.getBlocks(level, currDeltaT, currTime, vInfo, type);
			
				if (type == SpaceTimeSorter::ETimeInterval_Start)
				{
					_step_computeRHS_TR(vInfo, currTime, currDeltaT);
					_step_part<0>(vInfo, currTime, currDeltaT);
				}
				else if (type == SpaceTimeSorter::ETimeInterval_End)
				{
					_step_computeRHS_TR(vInfo, currTime, currDeltaT);
					_step_part<2>(vInfo, currTime, currDeltaT);
				}
				else
					abort();
				
				if (!bContinue) break;
			}
			
			stSorter.endSession();
		}
		
	public:
		
		STDTestL3_Diffusion(double tol= 1e-3):
			toleranceIC(tol), mu(1.0), ti(0.0), beta2(3.0),
			grid(blocksPerDimension,blocksPerDimension, b3D?blocksPerDimension:1), 
			refiner(resJump), compressor(resJump), 
			lab(), blockfwt(), viewer(false,true), t(0),
			stSorter()
		{
			grid.setCompressor(&compressor);
			grid.setRefiner(&refiner);
			
			stSorter.connect(grid);
			
			_ic(grid);

			Science::AutomaticRefinement< 0,0 >(grid, blockfwt, toleranceIC/2, maxLevel,-1, NULL, _ic);
		
			_ic(grid);
			
			Science::AutomaticCompression< 0,0 >(grid, blockfwt, toleranceIC, -1, NULL, _ic);
		}
							
		static void idle(void)
		{
			glutPostRedisplay();
		}
		
		static void display(void)
		{
			glClear(GL_COLOR_BUFFER_BIT);
			
			singletone->Step();
			singletone->Render();
			
			glutSwapBuffers();
		}
			
		public:
		
		double getAverageElementSize()
		{
			double sum = 0;
			
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++)
				sum += it->h[0];
		
			return sum/vInfo.size();
		}
		
		double getAverageResolutionSize()
		{
			return 1.0/sqrt(1.0*grid.getBlocksInfo().size()*B::sizeX*B::sizeY);
		}
		
		
		void Step(double factorDT=1.0)
		{
			const int nStep = 1;
			const double dx = 1.0/blockSize;
			const double F = 1./8;
			const double dt = factorDT*F*dx*dx/cazzo_viscosity;//(4.0*mu);//8e-1;
			
			printf("Time: %f\n",getTime());
			
		//	Science::AutomaticRefinement< 0,0 >(grid, blockfwt, toleranceIC/2,maxLevel);
			
			for(int i=0; i<nStep; i++)
			{
				printf("------------- Step %d over (%d). t is %e \n", i+1, nStep,t);
				_step(t, dt);
				t+= dt;
				printf("------------- END Step %d over (%d) t is %e \n", i+1, nStep,t );
			}
			
		//	Science::AutomaticCompression< 0,0 >(grid, blockfwt, toleranceIC);
		}
		
				
		//}
		
		double getTime() { return t;}
		
		void calculateError(double outputError[2])
		{
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			double err[2] = {0,0};
			double ext;
			double diff;
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				for(int iz=0; iz<B::sizeZ; iz++)
					for(int iy=0; iy<B::sizeY; iy++)
						for(int ix=0; ix<B::sizeX; ix++)
						{
							assert(iz==0);
							Real x[3];
							
							info.pos(x, ix, iy, iz);
							
							ext = _ic_test3(x,getTime());
							//assert(block(0,0,0).rho == block(ix,iy,iz).rho);
							
							diff = (ext-(double)block(ix,iy,iz).rho);
							//printf("error is (%d %d) = %e\n", ix, iy,diff);
							err[0]=err[0]+(diff*diff);
							err[1]=std::max(err[1],fabs(diff));
						}
			}
			err[0] = err[0]/(vInfo.size()*B::sizeY*B::sizeX);
			
			outputError[0] = err[0];
			outputError[1] = err[1];
		}
		
		
		void Render()
		{
			viewer.drawContent(grid, grid.getBlockCollection());
			viewer.drawSketch(grid,false);
		}
		
		static STDTestL3_Diffusion<W> * singletone; 
		
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
			

			if (bVisual)
			{
				singletone = new STDTestL3_Diffusion<W>();
				glutMainLoop();
			}
			else
			{
				/* do here the convergence study */
				const bool bDoTimestepAnalysis = true;
				const int nSteps = bDoTimestepAnalysis? 4 : 100;
				const int nRuns = bDoTimestepAnalysis? 1: 11;
				double tols[11] = {(bDoTimestepAnalysis? 1e-4: 1e-1) , 5e-2, 1e-2,  5e-3, 1e-3, 5e-4,1e-4, 5e-5, 1e-5, 5e-6, 1e-6};//, 1e-3, 5e-4};//, 1e-4};

				for(int iTimeStep=0; iTimeStep<(bDoTimestepAnalysis?6:1); iTimeStep++)
				for(int iRun=0; iRun<nRuns; iRun++)
				{
					const double factorDT = pow(1.333, -iTimeStep);
					singletone = new STDTestL3_Diffusion<W>(tols[iRun]);
					
					double avgLengthScale = singletone->getAverageElementSize();
					double avgResolution = singletone->getAverageResolutionSize();
					
					//for(int iStep=0; iStep<nSteps*pow(2, iTimeStep); iStep++)
					int iStep = 0;
					while(singletone->getTime()<0.01)
					{
						singletone->Step(factorDT);
						iStep++;
					}
					
					printf("****** Simulation Finished (t=%e)*********\n", singletone->getTime());
					//singletone->RefineOnceEverywhere(resJump);
					
					double error[2] = {0,0};
					singletone->calculateError(error);
					printf("=========================================\n");
					printf("Error(%d): %e %e\n", iRun, error[0], error[1]); 
					
					{
						static int counter = 0;
						
						FILE * f = fopen("Error.dat", counter++==0?"w":"a");
						
						if (bDoTimestepAnalysis)
							fprintf(f,"%e\t%e\t%e\t%e\t%e\t%e\n", factorDT, error[1],  factorDT,  avgLengthScale, avgResolution, error[1]);
						else
							fprintf(f,"%e\t%e\t%e\t%e\t%e\t%e\n", tols[iRun],  avgLengthScale, avgResolution, error[1], log(avgResolution), log(error[1]));
						
						fclose(f);
					}
					
					delete singletone;
					singletone = NULL;
				}
				
				{				
					FILE * fgnuplot = popen("/sw/bin/gnuplot -persist\n", "w");
					assert(fgnuplot != NULL);
					
					fprintf(fgnuplot, "set title \"Error (after %d sim steps)\n", nSteps);
					fprintf(fgnuplot, "set xlabel \" threshold \" \n");
					fprintf(fgnuplot, "set ylabel \"error\" \n");
					fprintf(fgnuplot, "set grid\n");
					fprintf(fgnuplot, "set logscale\n");
					
					if (bDoTimestepAnalysis)
						fprintf(fgnuplot, "plot \"Error.dat\" using 1:2  title 'cazz0' with linespoint\n");
					else
						fprintf(fgnuplot, "plot \"Error.dat\" using 1:4  title 'error vs compression threshold' with linespoint, \"Error.dat\" using 2:4 title 'error vs average element length' with linespoint , \"Error.dat\" using 3:4 title 'error vs 1./sqrt(N)' with linespoint\n");
					

					fprintf(fgnuplot, "set terminal postscript color\n");
					fprintf(fgnuplot, "set out \"Error.ps\"\n");
					fprintf(fgnuplot, "replot\n");
					
					 usleep(50000);
					 fclose(fgnuplot);
					 usleep(50000);
				}
			}
		}
	};
	
	template <typename Wavelets>
	STDTestL3_Diffusion< W> * STDTestL3_Diffusion< Wavelets>::singletone; 
}

