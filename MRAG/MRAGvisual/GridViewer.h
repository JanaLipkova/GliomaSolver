/*
 *  GridViewer.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 4/24/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once
#undef max
#undef min
#include <stdio.h>
#undef max
#undef min


#include "MRAGcore/MRAGBlock.h"
#include "MRAGcore/MRAGrid.h"
#include "Texture2D.h"
#include "MRAGVisualTypes.h"
using namespace MRAG;
using namespace MRAG::Visual;

/**
 * Used to visualize 2D grids.
 */
class GridViewer
{
private:
	Texture2D * m_texBlock;
	bool m_bCollocated;
	bool m_bRGBA;
	
	template <typename Block>
	void _Setup()
	{
		if (m_texBlock != NULL) return;
		
		m_texBlock = new Texture2D(Block::sizeX, Block::sizeY, m_bRGBA? GL_RGBA32F_ARB : GL_LUMINANCE32F_ARB, GL_NEAREST);
	}
	
	template<typename Block>
	void _drawBlockContent(Block& block, const short int index[3],  short int level)
	{
		const double h[2]= {pow(2.,-level), pow(2.,-level)};
		
		const double start[2] = {index[0]*h[0],index[1]*h[1]};
		const double end[2] = {(index[0]+1)*h[0],(index[1]+1)*h[1]};

		
		typedef typename Block::ElementType ElementType;
/*#if  ElementType != float
#error wrong Block::ElementType for this method
#endif*/
		const int n = Block::sizeX*Block::sizeY;
		
		glPushAttrib(GL_ENABLE_BIT);
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, m_texBlock->handle);
		if (!m_bRGBA)
		{
			/*abort();
			float * buf = new float[n];
			for(int i=0; i<n; i++) buf[i] = (float)block[i].s.r;
			glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, Block::sizeX,Block::sizeY, GL_RED, GL_FLOAT, buf);
			delete [] buf;*/
		}
		else
		{
			RGBA * buf = new RGBA[n];
			for(int i=0; i<n; i++) buf[i] = convertToRGBA(block[i]);
			glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, Block::sizeX,Block::sizeY, GL_RGBA, GL_FLOAT, buf);
			delete [] buf;
		}
		
		
		//m_texBlock->ShowContent();
		if (!m_bCollocated)
		{
			glBegin(GL_QUADS);
				glTexCoord2f(0, 0);
				glVertex2f(start[0],start[1]);
				glTexCoord2f(Block::sizeX/(double)m_texBlock->size[0], 0);
				glVertex2f(end[0],start[1]);
				glTexCoord2f(Block::sizeX/(double)m_texBlock->size[0], Block::sizeY/(double)m_texBlock->size[1]);
				glVertex2f(end[0],end[1]);
				glTexCoord2f(0, Block::sizeY/(double)m_texBlock->size[1]);
				glVertex2f(start[0],end[1]);
			glEnd();
		}
		else
		{
			//const int n[2] = {Block::sizeX+1, Block::sizeY+1};
			//const double d[2] = {h[0]/(n[0]-1),h[1]/(n[1]-1)};
			
			//const double offsetX = 0.5/Block::sizeX;
			//const double offsetY = 0.5/Block::sizeY;
			
			glBegin(GL_QUADS);
			glTexCoord2f(0, 0);
			glVertex2f(start[0],start[1]);
			glTexCoord2f(Block::sizeX/(double)m_texBlock->size[0], 0);
			glVertex2f(end[0],start[1]);
			glTexCoord2f(Block::sizeX/(double)m_texBlock->size[0], Block::sizeY/(double)m_texBlock->size[1]);
			glVertex2f(end[0],end[1]);
			glTexCoord2f(0, Block::sizeY/(double)m_texBlock->size[1]);
			glVertex2f(start[0],end[1]);
			
			glEnd();
		}
		
		glBindTexture(GL_TEXTURE_2D, 0);
		glPopAttrib();
	}
	
	template<typename Block>
	void _drawBlockSkeleton(const short int index[3],  short int level, const bool bDrawBlockElements, const float *vColorBoundary = NULL ,
					  const float *vColorContent = NULL)
	{
		const double h[2]= {pow(2.,-level), pow(2.,-level)};
		
		const double start[2] = {index[0]*h[0],index[1]*h[1]};
		const double end[2] = {(index[0]+1)*h[0],(index[1]+1)*h[1]};
		
		if (bDrawBlockElements)
		{
			const int n[2] = {Block::sizeX + 1, Block::sizeY+1};
			const double d[2] = {h[0]/(n[0]-1),h[1]/(n[1]-1)};
			
			glBegin(GL_LINES);
			
			if (vColorBoundary != NULL)
				glColor3fv(vColorContent);
			else
				glColor3f(0.2,0.2,0.2);
				
			for(int ix = 0; ix<n[0]; ix++)
			{
				glVertex2f(start[0] + ix*d[0], start[1]);
				glVertex2f(start[0] + ix*d[0], end[1]);
			}
			
			for(int iy = 0; iy<n[1]; iy++)
			{
				glVertex2f(start[0], start[1] + iy*d[1]);
				glVertex2f(end[0], start[1] + iy*d[1]);
			}
			glEnd();
			
			if (m_bCollocated)
			{
				//abort();
				const double start[2] = {index[0]*h[0],index[1]*h[1]};
			
				const double ptsize = 8.0;
				glPointSize(ptsize*pow(2,-(level-1)*0.5));
				glBegin(GL_POINTS);
				glColor3f(0.2,0.2,0.2);

				//glColor3f(drand48()+0.2,drand48()+0.2, drand48()+0.2);
				for(int iy = 0; iy<n[1]-1; iy++)
					for(int ix = 0; ix<n[0]-1; ix++)
					{
						glVertex2f(start[0] + ix*d[0], start[1] + iy*d[1]);
					}
				
				glEnd();
			}
		}
		
		glBegin(GL_LINES);
					
		if (vColorBoundary != NULL)
			glColor3fv(vColorBoundary);
		else
			glColor3f(0.8,0.8,0.8);	
					
		glVertex2f(start[0],start[1]);
		glVertex2f(start[0],end[1]);
		
		glVertex2f(end[0],start[1]);
		glVertex2f(end[0],end[1]);
		
		glVertex2f(start[0],start[1]);
		glVertex2f(end[0],start[1]);
		
		glVertex2f(start[0],end[1]);
		glVertex2f(end[0],end[1]);
		glEnd();
	}
	
	template<typename Block, typename A,typename B_, typename C>
	void _drawGridPoint(A block_index[3],  B_ level, C point_index[3], const float *vColorBoundary = NULL, double factor=1.)
	{
		const double h[2]= {pow(2.,-level), pow(2.,-level)};
		
		const double start[2] = {block_index[0]*h[0],block_index[1]*h[1]};
		//const double end[2] = {(block_index[0]+1)*h[0],(block_index[1]+1)*h[1]};
		const int n[2] = {Block::sizeX+1, Block::sizeY+1};
		const double d[2] = {h[0]/(n[0]-1),h[1]/(n[1]-1)};
		const double point_center[2] = {start[0] + (point_index[0]+ (m_bCollocated?0: 0.5))*d[0], start[1] + (point_index[1]+(m_bCollocated?0: 0.5))*d[1]};

		glPointSize(600*d[0]/2*factor);
		glBegin(GL_POINTS);
		if (vColorBoundary != NULL)
			glColor3fv(vColorBoundary);
		else
			glColor3f(1,0,0);
		
		
			glVertex2f(point_center[0], point_center[1]);
		glEnd();
	}
	
	template<typename Block, typename A,typename B_, typename C>
	void _drawConnection(A block_index[3],  B_ level, C point_index[3], 
						 A block_index2[3],  B_ level2, C point_index2[3], const float *vColorBoundary = NULL)
	{
		
		
				
		
		glBegin(GL_LINES);
		if (vColorBoundary != NULL)
			glColor3fv(vColorBoundary);
		else
			glColor3f(1,0,0);
		
			{
				const double h[2]= {pow(2.,-level), pow(2.,-level)};
				const double start[2] = {block_index[0]*h[0],block_index[1]*h[1]};
			///	const double end[2] = {(block_index[0]+1)*h[0],(block_index[1]+1)*h[1]};
				const int n[2] = {Block::sizeX+1, Block::sizeY+1};
				const double d[2] = {h[0]/(n[0]-1),h[1]/(n[1]-1)};
				const double point_center[2] = {start[0] + (point_index[0]+(m_bCollocated?0: 0.5))*d[0], start[1] + (point_index[1]+(m_bCollocated?0: 0.5))*d[1]};
				
				glVertex2f(point_center[0], point_center[1]);
			}
			{
				const double h[2]= {pow(2.,-level2), pow(2.,-level2)};
				const double start[2] = {block_index2[0]*h[0],block_index2[1]*h[1]};
			//	const double end[2] = {(block_index2[0]+1)*h[0],(block_index2[1]+1)*h[1]};
				const int n[2] = {Block::sizeX+1, Block::sizeY+1};
				const double d[2] = {h[0]/(n[0]-1),h[1]/(n[1]-1)};
				const double point_center[2] = {start[0] + (point_index2[0]+(m_bCollocated?0: 0.5))*d[0], start[1] + (point_index2[1]+(m_bCollocated?0: 0.5))*d[1]};
				
				glColor3f(1,1,1);
				glVertex2f(point_center[0], point_center[1]);
			}
		
		glEnd();
	}
	
	
public:

    /**
     * Default constructor for a GridViewer.
     * @param bCollocated   Optional (def = false): true if we are working on a collocated grid
     * @param bRGBA         Optional (def = false): true if convertToRGBA function used to get color for each element in grid.
     * @see GridViewer::drawContent()
     */
	GridViewer(bool bCollocated = false, bool bRGBA = false):m_texBlock(NULL), m_bCollocated(bCollocated), m_bRGBA(bRGBA) {}
	
    /**
     * Draws content of all blocks in the grid.
     * Requires:
     * - function MRAG::Visual::RGBA convertToRGBA(Block::ElementType & p) returning color of a grid-element.
     *   (only if bRGBA set to true in constructor)
     * - operator Real() returning floating point number representing a grid-element
     */
	template <typename Wavelets, typename Block>
	void drawContent(Grid<Wavelets, Block>& grid, const BlockCollection<Block>& collection)
	{
		_Setup<Block>();
		
		vector<BlockInfo> vInfo = grid.getBlocksInfo();
		
		for(int i=0; i<vInfo.size(); i++)
		{
			BlockInfo& info = vInfo[i];
			collection.lock(info.blockID);
			Block& block = collection[info.blockID];
			//Block& block = collection.lock(info.blockID);
			_drawBlockContent<Block>(block, info.index, info.level);
			collection.release(info.blockID);
		}
		
	}
	
    /**
     * Draws an outline of the grid-blocks.
     */
	template <typename Wavelets, typename Block>
	void drawSketch(Grid<Wavelets, Block>& grid, bool bDrawBlockElements=false)
	{
		glBegin(GL_LINES);
		glColor3f(0,1,0);
		glVertex2f(0,0);
		glVertex2f(0,1);
		
		glVertex2f(1,0);
		glVertex2f(1,1);
		
		glVertex2f(0,0);
		glVertex2f(1,0);
		
		glVertex2f(0,1);
		glVertex2f(1,1);
		glEnd();

		
		vector<BlockInfo> vInfo = grid.getBlocksInfo();
		
		for(int i=0; i<vInfo.size(); i++)
		{
			 BlockInfo& info = vInfo[i];
			_drawBlockSkeleton<Block>(info.index, info.level,  bDrawBlockElements);
		}
	}
	
	
	template <typename Wavelets, typename Block>
	void drawNeighborsOfABlock(MRAG::Grid<Wavelets, Block>& grid)
	{
		static int counter = 0;
		NeighborhoodType::const_iterator it = grid.m_neighborhood.begin();
		for(int i=0; i<counter; i++) it++;
		counter++;
		counter = counter%grid.m_neighborhood.size();
		GridNode& node = *(it->first);
	
		printf("NNNEgihhvbors %d\n", it->second.size());
		for(int i=0; i<it->second.size(); i++)
		{
			const float vContent[3] = {1,0,0};
			const float vBoundary[3] = {1,0,0};
			const GridNode& node = *(it->second[i]);
			_drawBlockSkeleton<Block>(node.index, node.level, true, vContent, vBoundary);
		}
		
		const float vContent[3] = {0,0,1};
		const float vBoundary[3] = {0,0,1};
		
		_drawBlockSkeleton<Block>(node.index, node.level, true, vContent, vBoundary);
	}
	
	template <typename Wavelets, typename Block>
	void drawGhosts(MRAG::Grid<Wavelets, Block>& grid, const BlockCollection<Block>& collection, BoundaryInfo & boundaryInfo)
	{
		static int counter = 10;
		counter++;
		
		vector<BlockInfo> vInfo = grid.getBlocksInfo();
		map<int, BlockInfo*> mapBlockIDToBlockInfo;
		for(int i=0; i<vInfo.size(); i++)
			mapBlockIDToBlockInfo[vInfo[i].blockID] = &vInfo[i];
		
		counter = counter % vInfo.size();
		
		BlockInfo& blockinfo =vInfo[counter];//*mapBlockIDToBlockInfo[10];// ////;////
		printf("BLOCKID = ======= %d\n", blockinfo.blockID);
		const float vBlock[3] = {1.,1.,1.};
		const float vContent[3] = {0.3,0.7,0.7};
		//const float vBoundary[3] = {1,0,0};
		const float vConnecions[3] = {1,1,0};
	//	const float vGhost[3] = {1,0,0};
		
		_drawBlockSkeleton<Block>(blockinfo.index, blockinfo.level, true, vBlock, vContent);
		
		const int nX = Block::sizeX;
		const int nY = Block::sizeY;
		const int nZ = Block::sizeZ;
	
		for(int iz=0; iz<nZ; iz++)
			for(int iy=0; iy<nY; iy++)
				for(int ix=0; ix<nX; ix++)
				{
					const int idx[3] = {ix,iy,iz};
					_drawGridPoint<Block>(blockinfo.index, blockinfo.level, idx, vContent);
				}
		
		BoundaryInfoBlock& bbinfo = *boundaryInfo.boundaryInfoOfBlock[blockinfo.blockID];
		bbinfo.lock();
		const vector< vector<IndexWP> >& ghosts = bbinfo.getGhosts();
		int ghost_counter = 0;
		for(int icode=0; icode<27; icode++)
		{
			if (icode == 1*1 + 3*1 + 9*1) continue;
			
			const int code[3] = { icode%3-1, (icode/3)%3-1, (icode/9)%3-1};
			
			const int s[3] = { 
				code[0]<1? (code[0]<0 ? boundaryInfo.stencil_start[0]:0 ) : nX, 
				code[1]<1? (code[1]<0 ? boundaryInfo.stencil_start[1]:0 ) : nY, 
				code[2]<1? (code[2]<0 ? boundaryInfo.stencil_start[2]:0 ) : nZ };
			
			const int e[3] = {
				code[0]<1? (code[0]<0 ? 0:nX ) : nX+boundaryInfo.stencil_end[0]-1, 
				code[1]<1? (code[1]<0 ? 0:nY ) : nY+boundaryInfo.stencil_end[1]-1, 
				code[2]<1? (code[2]<0 ? 0:nZ ) : nZ+boundaryInfo.stencil_end[2]-1};
			
			
			assert(bbinfo.boundary[icode].nGhosts == (e[2]-s[2])*(e[1]-s[1])*(e[0]-s[0]));
			
			int currentghost = bbinfo.boundary[icode].start;
			
			for(int iz=s[2]; iz<e[2]; iz++)
				for(int iy=s[1]; iy<e[1]; iy++)
					for(int ix=s[0]; ix<e[0]; ix++, currentghost++)
					{
						const vector<IndexWP>& vWP= ghosts[currentghost];
						const int n = vWP.size();
						
						double sum = 0;
						for(int i=0; i<n; i++)
							sum +=  bbinfo.weightsPool[vWP[i].weights_index[0]]* bbinfo.weightsPool[vWP[i].weights_index[1]]* bbinfo.weightsPool[vWP[i].weights_index[2]];
						printf("Ghost number: %d (start: %d nGhosts: %d, code=%d)\n",ghost_counter, currentghost, bbinfo.boundary[icode].nGhosts, icode);
						ghost_counter++;
						const bool bOk =  fabs(sum-1.0)<0.000001;
						
						if (!bOk)
						{
							printf("Failed: block %d %d %d l=%d, ghost: %d %d %d weight: %f (Contributions: %d)\n", 
								   blockinfo.index[0],  blockinfo.index[1],  blockinfo.index[2],  blockinfo.level, 
								  ix, iy,  iz, sum, n);
						}
						else
							{
						/*		printf("PASSED!! block %d %d %d l=%d, ghost: %d %d %d weight: %f (Contributions: %d)\n", 
									   blockinfo.index[0],  blockinfo.index[1],  blockinfo.index[2],  blockinfo.level, 
									   ix, iy,  iz, sum, n);*/
							}
							assert(bOk);	   
						for(int i=0; i<n; i++)
						{
							const int blockID = bbinfo.getIndexPool()[vWP[i].point_index].blockID;
							const int pointID = bbinfo.getIndexPool()[vWP[i].point_index].index;
							const int point_index[3]= {
							pointID % Block::sizeX, 
							(pointID/Block::sizeX )% Block::sizeY, 
							pointID/(Block::sizeX*Block::sizeY)
							};
							
							const int pt[3]= {ix, iy,iz};
							//const int * const block_index = mapBlockIDToBlockInfo[blockID]->index;
							_drawGridPoint<Block>(mapBlockIDToBlockInfo[blockID]->index, mapBlockIDToBlockInfo[blockID]->level, point_index, vConnecions,1.0);
							_drawConnection<Block>(mapBlockIDToBlockInfo[blockID]->index, mapBlockIDToBlockInfo[blockID]->level, point_index, 
											blockinfo.index, blockinfo.level, pt, vConnecions);
							
							if (!bOk)
							{
								printf("(Failed) Contribution %d: block %d %d %d l=%d, point %d %d %d, w=%f %f %f\n", 
									   i,
									   mapBlockIDToBlockInfo[blockID]->index[0],  mapBlockIDToBlockInfo[blockID]->index[1],  mapBlockIDToBlockInfo[blockID]->index[2],  mapBlockIDToBlockInfo[blockID]->level, 
									  point_index[0], point_index[1], point_index[2], 
									  bbinfo.weightsPool[vWP[i].weights_index[0]], bbinfo.weightsPool[vWP[i].weights_index[1]], bbinfo.weightsPool[vWP[i].weights_index[2]]);
							}
							
							
							//const unsigned char* const iw =vWP[i].weights_index;
							//const double w = (weight_pool[iw[0]]*weight_pool[iw[1]]*weight_pool[iw[2]]);
							//ghost += (value_pool[vWP[i].point_index]*w);//((weight_pool[iw[0]]*weight_pool[iw[1]]*weight_pool[iw[2]])*value_pool[vWP[i].point_index]);
						}
						
						
						const int pt[3]= {ix, iy,iz};
						const float vColor[3] = {bOk?1:0.5,bOk?0:0.5,bOk?0:0.5};
						_drawGridPoint<Block>(blockinfo.index, blockinfo.level, pt, vColor,1.0);
					//	const int point_index[3]= {ix, iy, iz};
					//	_drawGridPoint<Block>(blockinfo.index, blockinfo.level, point_index, vBoundary);
					 
					}
			
		}
		bbinfo.release();
		
//counter++;
		
	}
};
















