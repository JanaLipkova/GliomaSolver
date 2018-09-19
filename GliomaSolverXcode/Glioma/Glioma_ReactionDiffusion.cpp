//
//  Glioma_ReactionDiffusion.cpp
//  GliomaXcode
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
    bVerbose  = parser("-verbose").asBool(1);
    bProfiler = parser("-dumpfreq").asBool(1);
    
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("//////////////////          Glioma Reaction Diffusion           ////////////////\n");
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("RD INIT! nThreads=%d, blockSize=%d Wavelets=w%s (blocksPerDimension=%d, maxLevel=%d)\n", nThreads, blockSize, "w", blocksPerDimension, maxLevel);
    
    refiner		= new Refiner_SpaceExtension(resJump,maxLevel);
    compressor	= new Compressor(resJump);
    Environment::setup();
    
    grid = new Grid<W,B>(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil);
    grid->setCompressor(compressor);
    grid->setRefiner(refiner);
    stSorter.connect(*grid);
    
    bVTK        = parser("-vtk").asBool();
    bAdaptivity = parser("-adaptive").asBool();
    pID         = parser("-pID").asInt();
    L = 1;
    
    _ic(*grid, pID, L);
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
// Patient 01 data
// 1) read in anatomies - rescaled to [0,1]^3
// 2) read in tumor center of mass + initialize tumor around
// 3) set the characteristic length L as the length of the data
void Glioma_ReactionDiffusion::_ic(Grid<W,B>& grid, int pID, Real& L)
{
    char dataFolder   [200];
    char patientFolder[200];
    char anatomy      [200];
    
#ifdef KRAKEN
    sprintf(dataFolder,"/home/jana/Work/GliomaAdvance/source/Anatomy/");
#elif defined(LRZ_CLUSTER)
    sprintf(dataFolder,"/home/hpc/txh01/di49zin/GliomaAdvance/GliomaSolver/Anatomy/");
#elif defined(JANA)
    sprintf(dataFolder,"/Users/lipkova 1/WORK/GliomaSolver/Anatomy/");
#else
    sprintf(dataFolder,"../../Anatmoy/");
#endif
    
    sprintf(patientFolder, "%sPatient%02d/P%02d",dataFolder,pID,pID);
    printf("Reading anatomy from: %s \n", patientFolder);
    
    sprintf(anatomy, "%s_GM.dat", patientFolder);
    MatrixD3D GM(anatomy);
    sprintf(anatomy, "%s_WM.dat", patientFolder);
    MatrixD3D WM(anatomy);
    sprintf(anatomy, "%s_CSF.dat", patientFolder);
    MatrixD3D CSF(anatomy);
    
    int brainSizeX = (int) GM.getSizeX();
    int brainSizeY = (int) GM.getSizeY();
    int brainSizeZ = (int) GM.getSizeZ();
    
    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    L    = brainSizeMax * 0.1;   // voxel spacing 1mm, convert from mm to cm  // L = 25.6 cm
    
    printf("brainSizeX=%i, brainSizeY=%i, brainSizeZ= %i \n", brainSizeX, brainSizeY, brainSizeZ);
    std::cout<<"brainSizeX="<<brainSizeX<<" brainSizeY="<<brainSizeY<<" brainSizeZ="<<brainSizeZ<<std::endl;
    
    double brainHx = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
    /* Tumor Set UP */
    vector<Real> tumor_ic(_DIM);
    tumor_ic[0] = 0.315;
    tumor_ic[1] = 0.67;
    tumor_ic[2] = 0.5;
    
    const Real tumorRadius = 0.005;
    const Real smooth_sup  = 2.;		// suppor of smoothening, over how many gp to smooth
    
    Real pGM, pWM, pCSF;
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
#pragma omp paraller for
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        //        const Real h = vInfo[0].h[0];
        const Real h = 1./128;
        const Real iw = 1./(smooth_sup * h);   // width of smoothening => now it is over two grid points
        
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
                    {
                        pGM     =  GM( mappedBrainX,mappedBrainY,mappedBrainZ);
                        pWM     =  WM( mappedBrainX,mappedBrainY,mappedBrainZ);
                        pCSF    =  CSF(mappedBrainX,mappedBrainY,mappedBrainZ);
                        
                        // Anatomy
                        double all = pWM + pGM + pCSF;
                        
                        if(all > 0)
                        {
                            // normalize
                            pGM    = pGM  / all;
                            pWM    = pWM  / all;
                            pCSF   = pCSF / all;
                            
                            pCSF = ( pCSF > 0.1 ) ? 1. : pCSF;  // enhance csf for hemisphere separation
                            block(ix,iy,iz).p_csf = pCSF;
                            
                            if(pCSF  < 1.)
                            {
                                block(ix,iy,iz).p_csf = pCSF / (pCSF + pWM + pGM); // low level mixing for nicer visualisation
                                block(ix,iy,iz).p_w   = pWM  / (pCSF + pWM + pGM);
                                block(ix,iy,iz).p_g   = pGM  / (pCSF + pWM + pGM);
                            }
                            
                        }
                        
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
    UpdateTumor                     <_DIM>  updateTumor(dt);
    
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
    BlockProcessing::process(vInfo, collecton, updateTumor, nParallelGranularity);
}


#pragma mark DumpingOutput
void Glioma_ReactionDiffusion:: _dump(int counter)
{
    if(bVTK)
    {
        char filename[256];
        sprintf(filename,"P%02d_data_%04d",pID, counter);
        
        IO_VTKNative3D<W,B, 5,0 > vtkdumper2;
        vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
    }
}


void Glioma_ReactionDiffusion::run()
{
    const int nParallelGranularity	= (grid->getBlocksInfo().size()<=8 ? 1 : 4);
    BoundaryInfo* boundaryInfo		= &grid->getBoundaryInfo();
    
    /* Tumor growth parameters*/
    Real Dw     = (Real) parser("-Dw").asDouble(0.0013);
    Real rho    = (Real) parser("-rho").asDouble(0.025);
    double tend = parser("-Tend").asDouble(300);
    
    /*rescale*/
    Dw = Dw/(L*L);
    Real Dg = 0.1*Dw;
    
    double t			= 0.0;
    int iCounter        = 1;
    double h            = 1./(blockSize*blocksPerDimension);
    double dt           = 0.99 * h*h / ( 2.* _DIM * max(Dw, Dg) );
    if(bVerbose)  printf("Dg=%e, Dw=%e, dt= %f, rho=%f , h=%f\n", Dg, Dw, dt, rho,h);
    
    while (t <= tend)
    {
        if(bProfiler) profiler.getAgent("RD_Step").start();
        _reactionDiffusionStep(boundaryInfo, nParallelGranularity, Dw, Dg, rho, dt);
        if(bProfiler) profiler.getAgent("RD_Step").stop();
        
        t                   += dt   ;
        numberOfIterations  ++      ;
        
        if ( t >= ((double)(whenToWrite)) )
        {
            if(bAdaptivity)
            {
                Science::AutomaticRefinement	<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);
                Science::AutomaticCompression	<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler);
            }
            
            _dump(iCounter);
            iCounter++;
            whenToWrite = whenToWrite + whenToWriteOffset;
            if( (bVerbose) && (bVTK)) printf("Dumping data at time t=%f\n", t);
        }
    }
    
    
    // Refine + save the last one
    if(bAdaptivity)
        Science::AutomaticRefinement	<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);
    
    _dump(iCounter);
    
    if(bVerbose) profiler.printSummary();
    if(bVerbose) printf("**** Dumping done\n");
    if(bVerbose) printf("\n\n Run Finished \n\n");
}
