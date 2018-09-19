/*
 *  MRAG_STDTestL1_Compression.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/22/08.
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
	class MRAG_STDTestL1_Compression: public MRAG_STDTestL1<Wavelets, Block>
	{
	private:
		GridViewer viewer;
		Compressor compressor;
		MRAG::Grid<Wavelets, Block> grid;
		
		static MRAG_STDTestL1_Compression<Wavelets, Block> * singletone;
		
		static void display(void)
		{
			glClear(GL_COLOR_BUFFER_BIT);
			
			singletone->run();
			
			glutSwapBuffers();
		}
		
		static void idle(void)
		{
			glutPostRedisplay();
		}
		
		bool run()
		{
			viewer.drawSketch(grid,true);
			
			return true;
		}
		
		void checkJumpResolution()
		{
			for(HierarchyType::const_iterator it=grid.m_hierarchy.begin(); it!=grid.m_hierarchy.end(); it++)
			{
				GridNode * node = it->first;
				vector<GridNode *>& neighbors = grid.m_neighborhood[node];
				
				for(int n=0; n<neighbors.size(); n++)
					assert(abs(neighbors[n]->level - node->level)<=1);
			}
		}
		
		
		MRAG_STDTestL1_Compression(int nBlocks):
		viewer(), compressor(), grid(nBlocks,nBlocks)
		{
			for(int iCompressionStep = 0; iCompressionStep<8; iCompressionStep++)
			{
				grid.setCompressor(&compressor);
				
				set<int> shouldBeCompressed;
				vector<MRAG::BlockInfo> vInfo = grid.getBlocksInfo();
				for(int i=0; i<vInfo.size(); i++)
					if (drand48() < 0.75)
						shouldBeCompressed.insert(vInfo[i].blockID);
				
				int nCollapsed = 0;
				grid.compress(shouldBeCompressed, nCollapsed);
				printf("Collapsed = %d\n", nCollapsed);
			}
			
			checkJumpResolution();
		}
		
	public:
		
		static void runTests(int argc, char ** argv)
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
			
			singletone = new MRAG_STDTestL1_Compression<Wavelets, Block>(16);
			glutMainLoop();
		}
	};
	
	template <typename Wavelets, typename Block> 
	MRAG_STDTestL1_Compression<Wavelets, Block> * MRAG_STDTestL1_Compression<Wavelets, Block>::singletone = NULL;
	
}
#endif
