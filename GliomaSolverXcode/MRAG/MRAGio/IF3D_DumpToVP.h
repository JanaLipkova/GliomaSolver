//
//  IF3D_DumpToVP.h
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 5/15/13.
//
//

#ifndef IncompressibleFluids3D_IF3D_DumpToVP_h
#define IncompressibleFluids3D_IF3D_DumpToVP_h


template <typename BlockLabType>
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
        
        const int stencilStart[3] ={ -2, -2, -2};
        const int stencilEnd[3] ={  +3, +3, +3};
        const int ghostsize[3] = {stencilEnd[0]-stencilStart[0]-1,stencilEnd[1]-stencilStart[1]-1,stencilEnd[2]-stencilStart[2]-1};
        //2. dump header
        {
            FILE * file = fopen((filename+".txt").c_str(), "w");
            
            assert(file!=NULL);
            fprintf(file, "Blocks: %d\n", vInfo.size());
            fprintf(file, "Block size: %d %d %d\n", _BLOCKSIZE_, _BLOCKSIZE_, _BLOCKSIZE_);
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
            
            const size_t pointsPerBlock = (_BLOCKSIZE_+ghostsize[0])*(_BLOCKSIZE_+ghostsize[1])*(_BLOCKSIZE_+ghostsize[2]);
            const size_t bytesPerBlock = pointsPerBlock*sizeof(Real);
            
            struct Dummy{Real t;Dummy():t(0){}};
            
            
#pragma omp parallel
            {
                BlockLabType m_blockLab;
                
                m_blockLab.prepare(coll, binfo, stencilStart, stencilEnd);
                Dummy dummy;
                m_blockLab.inspect(dummy);
                
#pragma omp for
                for(int i=0;i<vInfo.size();i++)
                {
                    const size_t offset = bytesPerBlock*i;
                    m_blockLab.load(vInfo[i]);
#pragma omp critical
                    {
                        // put position indicator
                        fseek(file, offset, SEEK_SET);
                        const size_t nwrites = fwrite((void*) &m_blockLab.read(stencilStart[0],stencilStart[1],stencilStart[2]), sizeof(Real), pointsPerBlock, file);
                        assert(nwrites == pointsPerBlock);
                    }
                }
            }
            fclose(file);
        }
    }    
};

#endif


