/*
 *  MRAG_STDTestL2.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/25/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#include "MRAGEnvironment.h"

#ifndef _MRAG_GLUT_VIZ
#error Please do not include this file if you are not using glut, or if you do not define _MRAG_GLUT_VIZ. Diego.
#else

#include "GLUT/glut.h"
#include "../MRAGvisual/GridViewer.h"

#pragma once

#include "MRAGScienceCore.h"
#include "MRAGCommon.h"
#include "MRAGrid.h"
#include "MRAGRefiner.h"
#include "MRAGCompressor.h"
#include "MRAGBlockLab.h"
#include "MRAGBlockFWT.h"

namespace MRAG
{
	template <typename Wavelets, typename Block>
	class MRAG_STDTestL2 //LEVEL 2: correct data representation, correct data-processing
	{
		protected:
			Grid<Wavelets, Block> grid;
			Refiner refiner;
			Compressor compressor;
			BlockLab<Block> lab;
			BlockFWT<Wavelets, Block> blockfwt;
			GridViewer viewer;
		
			MRAG_STDTestL2():
				grid(4,4), refiner(), compressor(), lab(), 
				blockfwt(), viewer()
			{
				grid.setCompressor(&compressor);
				grid.setRefiner(&refiner);
				
				// OVERLOAD THIS
				try
				{
					printf("Start Test\n");
					_some_ic();
					_some_fwt();
					_some_refinements();
					_some_compression();
					_some_computation();
					
					printf("End Test\n");
				}
				catch(...)
				{
					
				}
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
							double x[2];
							info.pos(x, ix, iy);
							
							block(ix,iy) = sin(x[0]*2*3.1415)*cos(x[1]*2*3.1415);
						}
				}
			}
		
			void _some_computation()
			{
				const int steStart[3] ={ -1,-1,0};
				const int steEnd[3] ={ +2,+2,+1};
				
				vector<BlockInfo> vInfo = grid.getBlocksInfo();
				
				lab.prepare(grid.getBlockCollection(), grid.getBoundaryInfo(),steStart,steEnd);
				
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					
					lab.load(info);
					
					for(int iy=0; iy<Block::sizeY; iy++)
						for(int ix=0; ix<Block::sizeX; ix++)
							lab(ix,iy) += 0.5*lab(-1+ix,iy+1);
					
					lab.flush(info);
				}				
			}
		
			void _some_refinements() 
			{
				Science::AutomaticRefinement< 0,0 >(grid, blockfwt, 0.002);
			}
		
			
			void _some_compression()  
			{
				Science::AutomaticCompression< 0,0 >(grid, blockfwt, 0.002);
			}
		
			void _some_fwt()
			{

				BlockFWT<Wavelets, Block>::template multichannel_fwt<0,0>(grid.getBlocksInfo(),grid.getBlockCollection(), grid.getBoundaryInfo());
				
				double maxDetail = blockfwt.getReport().getOverAll_DetailMaxMag();
			}
		
		static void idle(void)
		{
			glutPostRedisplay();
		}
		
		bool run()
		{
			viewer.drawContent(grid, grid.getBlockCollection());
			viewer.drawSketch(grid,false);
			
			return true;
		}
		
		static void display(void)
		{
			glClear(GL_COLOR_BUFFER_BIT);
			
			singletone->run();
			
			glutSwapBuffers();
		}
			
		public:
		
		static MRAG_STDTestL2<Wavelets, Block> * singletone; 
		
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
			
			
			singletone = new MRAG_STDTestL2<Wavelets, Block>();
			
			if (bVisual) glutMainLoop();
			else
			{
				for(int i=0; i<300; i++)
				{
					delete singletone;
					singletone = new MRAG_STDTestL2<Wavelets, Block>(); 
				}
				printf("Binfo MB: %f\n" , singletone->grid.getBoundaryInfo().getMemorySize());
			}
			
		}
		
		};
	
		template <typename Wavelets, typename Block>
		MRAG_STDTestL2<Wavelets, Block> * MRAG_STDTestL2<Wavelets, Block>::singletone; 
}

#endif