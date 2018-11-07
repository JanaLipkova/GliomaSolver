//
//  Glioma_ReactionDiffusion.cpp
//  GliomaSolver
//
//  Created by Lipkova on 10/06/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#include "Glioma_ReactionDiffusion.h"

static int maxStencil[2][3] = {
    -1, -1, -1,
    +2, +2, +2
};

Glioma_ReactionDiffusion::Glioma_ReactionDiffusion(int argc, const char ** argv): parser(argc, argv)
{
    bVerbose    = parser("-verbose").asBool(1);
    bVTK        = parser("-vtk").asBool(1);
    bUQ         = parser("-UQ").asBool(0);
    bAdaptivity = parser("-adaptive").asBool(1);
    PatientFileName = parser("-PatFileName").asString();
    
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("//////////////////          Glioma Reaction Diffusion           ////////////////\n");
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("Set up: blockSize=%d Wavelets=w%s (blocksPerDimension=%d, maxLevel=%d)\n", blockSize, "w", blocksPerDimension, maxLevel);
    
    refiner		= new Refiner_SpaceExtension(resJump,maxLevel);
    compressor	= new Compressor(resJump);
    grid        = new Grid<W,B>(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil);
    grid->setCompressor(compressor);
    grid->setRefiner(refiner);
    stSorter.connect(*grid);
    
    L = 1;
    
    if(bUQ){
        ifstream mydata("TumorIC.txt");
        
        if (mydata.is_open()){
            mydata >> tumor_ic[0];
            mydata >> tumor_ic[1];
            mydata >> tumor_ic[2];
            mydata.close();
        }
        else{
            printf("Aborting: missing input file TumorIC.txt \n");
            abort();
        }
    }
    else{
        tumor_ic[0] = 0.315;
        tumor_ic[1] = 0.71;
        tumor_ic[2] = 0.4;
    }
    
    _ic(*grid, PatientFileName, L, tumor_ic);
    
    if(!bUQ) _dump(0);
    _dump(0);
    
    isDone              = false;
    whenToWriteOffset	= parser("-dumpfreq").asDouble();
    whenToWrite			= whenToWriteOffset;
    numberOfIterations	= 0;
}

Glioma_ReactionDiffusion::~Glioma_ReactionDiffusion()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialCondition
// 1) Read in anatomies - rescaled to [0,1]^3
// 2) Initialize tumor
// 3) Set the characteristic length L as the length of the data
void Glioma_ReactionDiffusion::_ic(Grid<W,B>& grid, string PatientFileName, Real& L, Real tumor_ic[3])
{
    printf("Reading data from file: %s \n", PatientFileName.c_str());
   
    char anatomy      [200];
    sprintf(anatomy, "%sGM.dat", PatientFileName.c_str());
    MatrixD3D GM(anatomy);
    sprintf(anatomy, "%sWM.dat", PatientFileName.c_str());
    MatrixD3D WM(anatomy);
    sprintf(anatomy, "%sCSF.dat", PatientFileName.c_str());
    MatrixD3D CSF(anatomy);
    
    int brainSizeX = (int) GM.getSizeX();
    int brainSizeY = (int) GM.getSizeY();
    int brainSizeZ = (int) GM.getSizeZ();
    printf("brainSizeX=%i, brainSizeY=%i, brainSizeZ= %i \n", brainSizeX, brainSizeY, brainSizeZ);
    
    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    L   = brainSizeMax * 0.1;   // voxel spacing 1mm, convert from mm to cm  // L = 25.6 cm
    printf("Characteristic Lenght L=%f \n", L);
    
    double brainHx = 1.0 / ((double)(brainSizeMax)); //  w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax)); //  w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax)); //  w.r.t. longest dimension for correct aspect ratio
    
    
    // Tumor set up
    const Real tumorRadius = 0.005;
    const Real smooth_sup  = 2.;		// suppor of smoothening
    const Real h           = 1./128;    // use fixed h, for same IC smoothening at all resolutions
    const Real iw          = 1./(smooth_sup * h); // widht of smoothening
    
    Real pGM, pWM, pCSF;
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
#pragma omp parallel for schedule(static)
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
                    
                    int mappedBrainX = (int)floor( x[0] / brainHx  );
                    int mappedBrainY = (int)floor( x[1] / brainHy  );
                    int mappedBrainZ = (int)floor( x[2] / brainHz  );
                    
                    // aspect ratio correction
//                    mappedBrainX -= (int) ( (brainSizeMax - brainSizeX) * 0.5);
//                    mappedBrainY -= (int) ( (brainSizeMax - brainSizeY) * 0.5);
//                    mappedBrainZ -= (int) ( (brainSizeMax - brainSizeZ) * 0.5);
                    
                    if ( (mappedBrainX >= 0 && mappedBrainX < brainSizeX) & (mappedBrainY >= 0 && mappedBrainY < brainSizeY) && (mappedBrainZ >= 0 && mappedBrainZ < brainSizeZ) )
                    {
                        // anatomy
                        pGM     =  GM( mappedBrainX,mappedBrainY,mappedBrainZ);
                        pWM     =  WM( mappedBrainX,mappedBrainY,mappedBrainZ);
                        pCSF    =  CSF(mappedBrainX,mappedBrainY,mappedBrainZ);
                        
                        // separat tissue and fluid + normalise
                        pCSF = (pCSF > 0.1)  ? 1. : pCSF;
                        pWM  = (pCSF > 0.1)  ? 0. : pWM;
                        pGM  = (pCSF > 0.1)  ? 0. : pGM;
                        
                        block(ix,iy,iz).p_csf = pCSF;
                        
                        double tissue = pWM + pGM;
                        block(ix,iy,iz).p_w = (tissue > 0.) ? (pWM / tissue) : 0.;
                        block(ix,iy,iz).p_g = (tissue > 0.) ? (pGM / tissue) : 0.;
                        
                        // tumor
                        const Real p[3] = {x[0] - tumor_ic[0], x[1] - tumor_ic[1], x[2] - tumor_ic[2]};
                        const Real dist = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);    // distance of curent voxel from tumor center
                        const Real psi = (dist - tumorRadius)*iw;
                        
                        if ((psi < -1) && (pGM + pWM >0.001))		// we are in tumor
                            block(ix,iy,iz).phi = 1.0;
                        else if(( (-1 <= psi) && (psi <= 1) )&& (pGM + pWM >0) )
                            block(ix,iy,iz).phi = 1.0 * 0.5 * (1 - psi - sin(M_PI*psi)/(M_PI));
                        else
                            block(ix,iy,iz).phi = 0.0;
                        
                    }
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}


#pragma mark ReactionDiffusion
void Glioma_ReactionDiffusion::_reactionDiffusionStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, const Real Dw, const Real Dg, const Real rho, double dt)
{
    
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    ReactionDiffusionOperator<_DIM>  rhs(Dw,Dg,rho);
    UpdateTumor              <_DIM>  updateTumor(dt);
    
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
    BlockProcessing::process(vInfo, collecton, updateTumor, nParallelGranularity);
}


#pragma mark DumpingOutput
void Glioma_ReactionDiffusion:: _dump(int counter)
{
    if(bVTK)
    {
        char filename[256];
        sprintf(filename,"Data_%04d", counter);
        
        IO_VTKNative3D<W,B, 5,0 > vtkdumper2;
        vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
    }
}


/* Dump output for UQ likelihood. Requirements:
 - dump at the uniform finest resolution
 - use 3D Matrix structure to dump data in binary format
 - assume 3D simulation */
void Glioma_ReactionDiffusion::_dumpUQoutput()
{
    int gpd    = blocksPerDimension * blockSize;
    double hf  = 1./gpd;
    double eps = hf*0.5;
    
    MatrixD3D tumor(gpd,gpd,gpd);
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
#pragma omp parallel for schedule(static)
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
                    int mx = (int)floor( (x[0]) / hf  );
                    int my = (int)floor( (x[1]) / hf  );
                    int mz = (int)floor( (x[2]) / hf  );
                    
                    if(h < hf + eps)
                        tumor(mx,my,mz) = block(ix,iy,iz).phi;
                    else if(h < 2.*hf + eps)
                    {
                        for(int cz=0; cz<2; cz++)
                            for(int cy=0; cy<2; cy++)
                                for(int cx=0; cx<2; cx++)
                                    tumor(mx+cx,my+cy,mz+cz) = block(ix,iy,iz).phi;
                    }
                    else if (h < 3.*hf + eps)
                    {
                        for(int cz=0; cz<3; cz++)
                            for(int cy=0; cy<3; cy++)
                                for(int cx=0; cx<3; cx++)
                                    tumor(mx+cx,my+cy,mz+cz) = block(ix,iy,iz).phi;
                    }
                    else
                    {
                        for(int cz=0; cz<4; cz++)
                            for(int cy=0; cy<4; cy++)
                                for(int cx=0; cx<4; cx++)
                                    tumor(mx+cx,my+cy,mz+cz) = block(ix,iy,iz).phi;
                    }
                }
        
    }
    
    char filename2[256];
    sprintf(filename2,"HGG_data.dat");
    tumor.dump(filename2);
    
}



#pragma mark RUN
void Glioma_ReactionDiffusion::run()
{
    const int nParallelGranularity	= (grid->getBlocksInfo().size()<=8 ? 1 : 4);
    BoundaryInfo* boundaryInfo		= &grid->getBoundaryInfo();
    
    /* Tumor growth parameters*/
    Real Dw, Dg, rho, tend;
    
    if(bUQ){
        ifstream mydata("InputParameters.txt");
        if (mydata.is_open())
        {
            mydata >> Dw;
            mydata >> rho;
            mydata >> tend;
            mydata.close();
        }
        else{
            printf("Aborting: missing input file InputParameters.txt \n");
            abort();
        }
    }
    else
    {
        Dw   = (Real) parser("-Dw").asDouble(0.0013);
        rho  = (Real) parser("-rho").asDouble(0.025);
        tend = (Real) parser("-Tend").asDouble(300);
    }
    
    /*rescale for correct space dimension*/
    Dw = Dw/(L*L);
    Dg = 0.1*Dw;
    
    Real t			= 0.0;
    Real h          = 1./(blockSize*blocksPerDimension);
    Real dt         = 0.99 * h*h / ( 2.* _DIM * max(Dw, Dg) );
    int iCounter    = 1;
    if(bVerbose) printf("Dg=%e, Dw=%e, dt= %f, rho=%f , h=%f\n", Dg, Dw, dt, rho,h);
    
    while (t <= tend)
    {
        _reactionDiffusionStep(boundaryInfo, nParallelGranularity, Dw, Dg, rho, dt);
        t                   += dt   ;
        numberOfIterations  ++      ;
        
        if ( t >= ((double)(whenToWrite)) )
        {
            if(bAdaptivity){
                Science::AutomaticRefinement	<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);
                Science::AutomaticCompression	<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler);
            }
            
            _dump(iCounter++);
            whenToWrite = whenToWrite + whenToWriteOffset;
            if( (bVerbose) && (bVTK)) printf("Dumping data at time t=%f\n", t);
        }
    }
    
    
    // Refine + save the last one
    if(bAdaptivity)
        Science::AutomaticRefinement<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);
    
    _dump(iCounter);
    if(bUQ) _dumpUQoutput();
    
    if(bVerbose) printf("**** Dumping done\n");
    if(bVerbose) printf("\n\n Run Finished \n\n");
}
