/*
 *  MRAG_STDTestL1_BoundaryInfo.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/23/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */



#include "MRAG_STDTestL1.h"

#ifndef _MRAG_GLUT_VIZ
#error Please do not include this file if you are not using glut, or if you do not define _MRAG_GLUT_VIZ. Diego.
#else

#include "GLUT/glut.h"
#include "../MRAGvisual/GridViewer.h"

#pragma once
namespace MRAG
{
	template <typename Wavelets, typename Block>
	class MRAG_STDTestL1_BoundaryInfo: public MRAG_STDTestL1<Wavelets, Block>
	{
	private:
		GridViewer viewer;
		Compressor compressor;
		Refiner refiner;
		MRAG::Grid<Wavelets, Block> grid;
		
		static MRAG_STDTestL1_BoundaryInfo<Wavelets, Block> * singletone;
		
		static void display(void)
		{
			glClear(GL_COLOR_BUFFER_BIT);
			
			singletone->run();
			
			glutSwapBuffers();
		}
		
		static void idle(void)
		{
			//int i;
			//cin>>i;
			glutPostRedisplay();
		}
		
		bool run()
		{
			viewer.drawSketch(grid,true);
			viewer.drawGhosts(grid, grid.getBlockCollection(),grid.getBoundaryInfo());
			
			return true;
		}
		
		void checkJumpResolution()
		{
			for(HierarchyType::const_iterator it=grid.m_hierarchy.begin(); it!=grid.m_hierarchy.end(); it++)
			{
				GridNode * node = it->first;
				vector<GridNode *>& neighbors = grid.m_neighborhood[node];
				
				for(int n=0; n<neighbors.size(); n++)
					assert(abs(neighbors[n]->level - node->level)<=2);
			}
		}
		
		template<typename Grid>
		int _getID(int ix, int iy, int level, Grid& g)
		{
			vector<BlockInfo> vInfo = g.getBlocksInfo();
			
			for(vector<BlockInfo>::iterator it = vInfo.begin(); it!=vInfo.end(); it++)
				if (ix == it->index[0]  && iy == it->index[1] && it->level == level)
					return it->blockID;
			
			abort();
			
			return -1;
		}
		
		MRAG_STDTestL1_BoundaryInfo(int nBlocks):
		viewer(true), compressor(), refiner(2), grid(Block::sizeX>1?2:1,
										 Block::sizeY>1?2:1,
										 Block::sizeZ>1?2:1)
		{
			const bool bStandardTest = true;
			
			if (bStandardTest)
			{
				for(int iCompressionStep = 0; iCompressionStep<3; iCompressionStep++)
				{
					grid.setCompressor(&compressor);
					
					set<int> shouldBeCompressed;
					vector<MRAG::BlockInfo> vInfo = grid.getBlocksInfo();
					for(int i=0; i<vInfo.size(); i++)
						if (drand48() < 0.85)
							shouldBeCompressed.insert(vInfo[i].blockID);
					
					int nCollapsed = 0;
					grid.compress(shouldBeCompressed, nCollapsed);
					printf("Collapsed = %d\n", nCollapsed);
				}
				
				for(int iRef = 0; iRef<3; iRef++)
				{
					grid.setRefiner(&refiner);
					
					set<int> shouldBeRefined;
					vector<MRAG::BlockInfo> vInfo = grid.getBlocksInfo();
					for(int i=0; i<vInfo.size(); i++)
						if (drand48() > 0.7)
							shouldBeRefined.insert(vInfo[i].blockID);
					
					grid.refine(shouldBeRefined);
				}
			}
			else {
				grid.setRefiner(&refiner);
				
				set<int> ref1;
				ref1.insert(_getID(1,0,1, grid));
				ref1.insert(_getID(0,1,1, grid));
				ref1.insert(_getID(1,1,1, grid));
				grid.refine(ref1);

				set<int> ref2;
				ref2.insert(_getID(0,2,2, grid));
				ref2.insert(_getID(1,2,2, grid));
				ref2.insert(_getID(2,2,2, grid));
				grid.refine(ref2);
			}

		//	grid.getBoundaryInfo();
			//exit(0);
			//checkJumpResolution();
		}
		
	public:
		
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
			
			
			singletone = new MRAG_STDTestL1_BoundaryInfo<Wavelets, Block>(4);
			
			if (bVisual) glutMainLoop();
			else
			{
				for(int i=0; i<300; i++)
				{
					delete singletone;
					singletone = new MRAG_STDTestL1_BoundaryInfo<Wavelets, Block>(8); 
				}
				printf("Binfo MB: %f\n" , singletone->grid.getBoundaryInfo().getMemorySize());
			}
			
		}
	};
	
	template <typename Wavelets, typename Block> 
	MRAG_STDTestL1_BoundaryInfo<Wavelets, Block> * MRAG_STDTestL1_BoundaryInfo<Wavelets, Block>::singletone = NULL;
}
#endif
