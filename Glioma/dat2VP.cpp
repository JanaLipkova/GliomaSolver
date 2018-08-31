//
//  dat2VP.cpp
//  GliomaSolverXcode
//
//  Created by Lipkova on 31/08/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

#include "dat2VP.h"

// need biger stencil for the VP
static int maxStencil[2][3] = {
    -2, -2, -2,
    +3, +3, +3
};

dat2VP::dat2VP(int argc, const char ** argv): parser(argc, argv)
{
    bVerbose = parser("-verbose").asBool();
    
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("//////////////////                 DUMP 2 VP                    ////////////////\n");
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    
    if(bVerbose) printf("suggested commands:\n");
    if(bVerbose) printf("mv test test_t%d_b%d_w%s\n", nThreads, blockSize, "w");
    if(bVerbose) printf("RD INIT! nThreads=%d, blockSize=%d Wavelets=w%s (blocksPerDimension=%d, maxLevel=%d)\n", nThreads, blockSize, "w", blocksPerDimension, maxLevel);
    
    refiner		= new Refiner_SpaceExtension(resJump,maxLevel);
    compressor	= new Compressor(resJump);
    Environment::setup();
    
    grid = new Grid<W,B>(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil);
    grid->setCompressor(compressor);
    grid->setRefiner(refiner);
    stSorter.connect(*grid);
    
    inputFileName = parser("-inFileName").asString();
    
    _ic(*grid,inputFileName);
    isDone              = false;
}

dat2VP::~dat2VP()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions
void dat2VP::_ic(Grid<W,B>& grid, std::string inputFileName)
{
    
    printf("Reading anatomy from: %s \n", inputFileName.c_str());
    MatrixD3D DATA(inputFileName.c_str());

    
    /* Read in binary matrices */
    
    // Tumour
    int brainSizeX = (int) DATA.getSizeX();
    int brainSizeY = (int) DATA.getSizeY();
    int brainSizeZ = (int) DATA.getSizeZ();
    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    
    std::cout<<"DATA: "<<brainSizeX<<","<<brainSizeY<<","<<brainSizeZ<<std::endl;
    
    double brainHx = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    int mappedBrainX = (int)floor( x[0] / brainHx  );
                    int mappedBrainY = (int)floor( x[1] / brainHy  );
                    int mappedBrainZ = (int)floor( x[2] / brainHz  );
                    
                    // aspect ratio correction
                    mappedBrainX -= (int) ( (brainSizeMax - brainSizeX) * 0.5);
                    mappedBrainY -= (int) ( (brainSizeMax - brainSizeY) * 0.5);
                    mappedBrainZ -= (int) ( (brainSizeMax - brainSizeZ) * 0.5);
                    
                    
                    if ( (mappedBrainX > 0 && mappedBrainX < brainSizeX-1) && (mappedBrainY > 0 && mappedBrainY-1 < brainSizeY-1) && (mappedBrainZ > 0 && mappedBrainZ < brainSizeZ-1) )
                        block(ix,iy,iz).vp = DATA(mappedBrainX,mappedBrainY,mappedBrainZ);
                }
        
        grid.getBlockCollection().release(info.blockID);

    }
}


#pragma mark DumpingOutput
void dat2VP:: _dumpVTK(int counter)
{
    
    if (parser("-vtk").asBool())
    {
        printf("dumping VTK data \n");
        
        char filename[256];
        sprintf(filename,"Visualisation%04i", counter);
        
        IO_VTKNative3D<W,B, 1,0 > vtkdumper;
        vtkdumper.Write(*grid, grid->getBoundaryInfo(), filename);
    }
}

void dat2VP:: _dumpVP(int counter)
{
    if (parser("-vp").asBool())
    {
        printf("dumping VP data \n");
        char filename[256];
        sprintf(filename,"Visualisation%04i", counter);
        
        DumpScalarToVP < BlockLab<B>, W, B>	vpdumper(grid);
        vpdumper.Write(filename);
    }
}

void dat2VP::run()
{
    int count = 0;
    _dumpVTK(count);
    _dumpVP(count);
    
    
    if(bVerbose) printf("**** Dumping done\n");
    if(bVerbose) printf("\n\n Run Finished \n\n");
}
