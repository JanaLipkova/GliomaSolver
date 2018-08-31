/*
 *  DumpScalarToVP.h
 *  XcodeGliomaProject
 *
 *  Created by Lipkova on 4/25/14, based on Wim's code from IF3D MRAG, adapted for glioma code
 *  Copyright 2014 Grid. All rights reserved.
 *
 */

#pragma once
#include <vector>
#include <stack>
#include <string>
#include <fstream>

using namespace std;
using namespace MRAG;


template <typename BlockLabType, typename W,  typename B>
class DumpScalarToVP
{
    Grid<W,B> * grid_ptr;
    
public:
    DumpScalarToVP(Grid<W,B> * grid_ptr): grid_ptr(grid_ptr)
    {
    }
    
    void Write(std::string filename)
    {
        //1.
        vector<BlockInfo> vInfo = grid_ptr->getBlocksInfo();
        const BlockCollection<B>& coll = grid_ptr->getBlockCollection();
        const BoundaryInfo& binfo=grid_ptr->getBoundaryInfo();
        
        /* VP needs at least 2 ghosts on each side, since it used cubic interploation, therefore here one need to modify the stencilStart/End */
        const int stencilStart[3] ={ -2, -2, -2};
        const int stencilEnd[3] ={  +3, +3, +3};
        
        
        const int ghostsize[3] = {stencilEnd[0]-stencilStart[0]-1,stencilEnd[1]-stencilStart[1]-1,stencilEnd[2]-stencilStart[2]-1};
        //2. dump header
        {
            FILE * file = fopen((filename+".txt").c_str(), "w");
            
            assert(file!=NULL);
            fprintf(file, "Blocks: %d\n", vInfo.size());
            fprintf(file, "Block size: %d %d %d\n", blockSize, blockSize, blockSize);
            fprintf(file, "Ghostsize: %d %d %d\n", ghostsize[0], ghostsize[1], ghostsize[2]);
            fprintf(file, "sizeof(Real): %lu\n", sizeof(Real));
            
            for(int i=0; i<vInfo.size(); i++)
            {
                const BlockInfo& info = vInfo[i];
                
                fprintf(file, "Block %d: Tree Index: %d %d %d, %d\n", i,
                        info.index[0], info.index[1], info.index[2], info.level);
            }
            fclose(file);
        }
        
        // 3. dump data
        {
            FILE * file = fopen((filename+".dat").c_str(), "wb");
            assert(file!=NULL);
            
            const size_t pointsPerBlock = (blockSize+ghostsize[0])*(blockSize+ghostsize[1])*(blockSize+ghostsize[2]);
            const size_t bytesPerBlock = pointsPerBlock*sizeof(Real);
            
            struct Dummy{Real t;Dummy():t(0){}};
            
            vector<Real> tmp(pointsPerBlock,0.0);
            
#pragma omp parallel
            {
                BlockLabType m_blockLab;
                m_blockLab.prepare(coll, binfo, stencilStart, stencilEnd);
                Dummy dummy;
                //m_blockLab.inspect(dummy);
                
#pragma omp for
                for(int i=0;i<vInfo.size();i++)
                {
                    const size_t offset = bytesPerBlock*i;
                    m_blockLab.load(vInfo[i]);
                    
                    
                    int Nx = B::sizeX + stencilEnd[0]-1;
                    int Ny = B::sizeY + stencilEnd[1]-1;
                    int Nz = B::sizeZ + stencilEnd[2]-1;
                    
                    int Xsize = B::sizeX + ghostsize[0];
                    int Ysize = B::sizeY + ghostsize[1];
                    int Zsize = B::sizeZ + ghostsize[2];
                    
                    Real tmp_phi, tmp_w, tmp_g;
                    
                    
                    for (int iz = stencilStart[2]; iz < Nz; iz++)
                        for (int iy = stencilStart[1]; iy < Ny; iy++)
                            for (int ix = stencilStart[0]; ix < Nx; ix++)
                                tmp[ix - stencilStart[0] + (iy-stencilStart[1]) * Xsize + (iz - stencilStart[2]) * Xsize * Ysize] =m_blockLab(ix,iy,iz).vp;
                    
                    
#pragma omp critical
                    {
                        // put position indicator
                        fseek(file, offset, SEEK_SET);
                        const size_t nwrites = fwrite((void*) &tmp[0], sizeof(Real), pointsPerBlock, file);
                        
                        assert(nwrites == pointsPerBlock);
                    }
                }
            }
            fclose(file);
        }
    }    
};

