//
//  Glioma_ComputePFF_CahnHilliard.cpp
//  GliomaSolverXcode
//
//  Created by Lipkova on 21/09/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

#include "Glioma_ComputePFF_CahnHilliard.h"

static int maxStencil[2][3] = {
    -3, -3, -3,
    +4, +4, +4
};


Glioma_ComputePFF_CahnHilliard::Glioma_ComputePFF_CahnHilliard(int argc, const char ** argv): parser(argc, argv)
{
    bVerbose    = parser("-verbose").asBool(1);
    bVTK        = parser("-vtk").asBool();
    bProfiler   = parser("-profiler").asBool(1);
    bAdaptivity = parser("-adaptive").asBool();
    
    printf("////////////////////////////////////////////////////////////////////////////////\n");
    printf("///////   Computing Phase Filed Function using Cahn-Hilliard Operator    ///////\n");
    printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("Set up: blockSize=%d Wavelets=w%s (blocksPerDimension=%d, maxLevel=%d)\n", blockSize, "w", blocksPerDimension, maxLevel);
    
    refiner		= new Refiner_SpaceExtension();
    compressor	= new Compressor();
    grid        = new Grid<W,B>(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil);
    grid->setCompressor(compressor);
    grid->setRefiner(refiner);
    stSorter.connect(*grid);
    
    PatientFileName = parser("-PatFileName").asString();
    L = 1;
    
    _ic(*grid, PatientFileName, L);
    
    if(bAdaptivity)
    {
        Science::AutomaticRefinement	<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);
        Science::AutomaticCompression	<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler);
    }
    
    if(parser("-bDumpIC").asBool(0)) _dump(0);
    
    isDone              = false;
    whenToWriteOffset	= parser("-dumpfreq").asDouble();
    whenToWrite			= whenToWriteOffset;
    numberOfIterations	= 0;
}

Glioma_ComputePFF_CahnHilliard::~Glioma_ComputePFF_CahnHilliard()
{
    std::cout << "------Adios muchachos------" << std::endl;
}

#pragma mark InitialCondition
void Glioma_ComputePFF_CahnHilliard::_ic(Grid<W,B>& grid, string PatientFileName, Real&L)
{
    printf("Reading data from file: %s \n", PatientFileName.c_str());
    
    char anatomy      [200];
    sprintf(anatomy, "%s_mask.dat", PatientFileName.c_str());
    MatrixD3D MASK(anatomy);
    
    int brainSizeX = (int) MASK.getSizeX();
    int brainSizeY = (int) MASK.getSizeY();
    int brainSizeZ = (int) MASK.getSizeZ();
    
    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    L    = brainSizeMax * 0.1;   // voxel spacing 1mm, convert from mm to cm  // L = 25.6 cm
    
    printf("brainSizeX=%i, brainSizeY=%i, brainSizeZ= %i \n", brainSizeX, brainSizeY, brainSizeZ);
    
    double brainHx = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
#pragma omp paraller for
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
                    
                    /* Anatomy */
                    int mappedBrainX = (int)floor( x[0] / brainHx  );
                    int mappedBrainY = (int)floor( x[1] / brainHy  );
                    int mappedBrainZ = (int)floor( x[2] / brainHz  );
                    
                    // aspect ratio correction
                    mappedBrainX -= (int) ( (brainSizeMax - brainSizeX) * 0.5);
                    mappedBrainY -= (int) ( (brainSizeMax - brainSizeY) * 0.5);
                    mappedBrainZ -= (int) ( (brainSizeMax - brainSizeZ) * 0.5);
                    
                    if ( (mappedBrainX >= 0 && mappedBrainX < brainSizeX) & (mappedBrainY >= 0 && mappedBrainY < brainSizeY) && (mappedBrainZ >= 0 && mappedBrainZ < brainSizeZ) )
                        block(ix,iy,iz).phi  =  MASK( mappedBrainX,mappedBrainY,mappedBrainZ);
                    
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}



void Glioma_ComputePFF_CahnHilliard::_CahnHilliardStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, double dt, const int w)
{
    CahnHilliardPotential<_DIM>	potential(w);
    CahnHilliardEquation <_DIM> rhs;
    CahnHilliardUpdate   <_DIM> update(dt);
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, potential);
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
    BlockProcessing::process(vInfo, collecton, update, nParallelGranularity);
}

#pragma mark DumpingOutput
void Glioma_ComputePFF_CahnHilliard:: _dump(int counter)
{
    if (parser("-vtk").asBool())
    {
        char filename[256];
        sprintf(filename,"PFF_%04d", counter);
        
        IO_VTKNative3D<W,B, 1,0 > vtkdumper2;
        vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
    }
}


// dump the uniform resolution data into binary matrix
void Glioma_ComputePFF_CahnHilliard::_dumpBinary(int counter)
{
    int gpd    = blocksPerDimension * blockSize;
    double hf  = 1./gpd;
    double eps = hf*0.5;
    
    MatrixD3D PFF(gpd,gpd,gpd);
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
#pragma omp parallel for
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        double h = info.h[0];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    //mapped coordinates
                    int mx = (int)round( (x[0]) / hf  );
                    int my = (int)round( (x[1]) / hf  );
                    int mz = (int)round( (x[2]) / hf  );
                    
                    if( h < hf + eps)
                        PFF(mx,my,mz) = block(ix,iy,iz).phi;
                    else if(h < 2.*hf + eps)
                    {
                        for(int cz=0; cz<2; cz++)
                            for(int cy=0; cy<2; cy++)
                                for(int cx=0; cx<2; cx++)
                                    PFF(mx+cx,my+cy,mz+cz) = block(ix,iy,iz).phi;
                    }
                    else if (h < 3.*hf + eps)
                    {
                        for(int cz=0; cz<3; cz++)
                            for(int cy=0; cy<3; cy++)
                                for(int cx=0; cx<3; cx++)
                                    PFF(mx+cx,my+cy,mz+cz) = block(ix,iy,iz).phi;
                    }
                    else
                    {
                        for(int cz=0; cz<4; cz++)
                            for(int cy=0; cy<4; cy++)
                                for(int cx=0; cx<4; cx++)
                                    PFF(mx+cx,my+cy,mz+cz) = block(ix,iy,iz).phi;
                    }
                }
        
    }
    
    char filename[256];
    sprintf(filename,"PFF_%04d.dat",counter);
    PFF.dump(filename);
    
}

void Glioma_ComputePFF_CahnHilliard::run()
{
    const int nParallelGranularity	= (grid->getBlocksInfo().size()<=8 ? 1 : 4);
    BoundaryInfo* boundaryInfo		= &grid->getBoundaryInfo();
    
    double h            = 1./(blockSize*blocksPerDimension);
    double dt           = h*h*h / (2*_DIM);        // 4th order method should be h^4
    int Niter           = parser("-Niter").asInt(100);
    int iCounter        = 0;
    int w               = parser("-width").asInt(9); // smoothen over w grid points
    
    while (iCounter <= Niter )
    {
        _CahnHilliardStep(boundaryInfo, nParallelGranularity, dt, w);
        iCounter++;
        
        if ( iCounter >= ((double)(whenToWrite)) )
        {
            if(bAdaptivity)
            {
                Science::AutomaticRefinement	<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);
                Science::AutomaticCompression	<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler);
            }
            
            _dump(iCounter);
            _dumpBinary(iCounter);
            whenToWrite = whenToWrite + whenToWriteOffset;
        }
    }
        
    printf("**** Dumping done\n");
    printf("\n\n Run Finished \n\n");
}