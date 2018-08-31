//
//  HelmholtzTest.cpp
//  GliomaSolver
//
//  Created by Lipkova on 31/03/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#include "HelmholtzTest.h"
#include <mpi.h>

static int maxStencil[2][3] = {
    -1, -1, -1,
    +2, +2, +2
};


HelmholtzTest::HelmholtzTest(int argc, const char ** argv): parser(argc, argv), helmholtz_solver3D(argc,argv), helmholtz_solver2D(argc,argv)
{
    bVerbose  = parser("-verbose").asBool(1);
    bProfiler = parser("-profiler").asBool(1);
    bVTK      = parser("-vtk").asBool(1);
    
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("//////////////////        HELMHOLTZ TEST: USING HYPRE           ////////////////\n");
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("RD INIT! nThreads=%d, blockSize=%d Wavelets=w%s (blocksPerDimension=%d, maxLevel=%d)\n", nThreads, blockSize, "w", blocksPerDimension, maxLevel);
    
    refiner		= new Refiner_SpaceExtension();
    compressor	= new Compressor();
    Environment::setup();
    
    grid = new Grid<W,B>(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil);
    grid->setCompressor(compressor);
    grid->setRefiner(refiner);
    stSorter.connect(*grid);
    
    _ic_Square(*grid);
    _dump(0);
    
    MPI_Init(&argc, (char ***)&argv);
    
    isDone              = false;
    numberOfIterations	= 0;
}

HelmholtzTest::~HelmholtzTest()
{
    MPI_Finalize();
    std::cout << "------Adios muchachos------" << std::endl;
}

#pragma mark InitialConditions

/* TEST CASE 0 - SQUARE
 -----------------------
 -grad(psi grad p) + psi*p = f*psi  in Ω1 = [0.25,0.75]^2 in Ω2=[0,1]^2
 p = 0 on ∂Ω2
 -----------------------
 p = cos(c*PI*x) * cos(c*PI*y)
 f =  (dim. * c*c * PI*PI + 1.) * cos(c*PI*x) * cos(c*PI*y); */
void HelmholtzTest::_ic_Square(Grid<W,B>& grid)
{
    std::cout <<" Test case on a SQUARE"<< std::endl;
    
    const float radius      = 0.25;
    const double tau        = 1.e-10;    // cut of phase field function on LHS
    const float  c          = 4;         // mod of solution
    const float center[3]   = {0.5, 0.5, 0.5};
    Real r;
    
    if(_DIM == 2)
    {
        vector<BlockInfo> vInfo = grid.getBlocksInfo();;
        
        for(int i=0; i<vInfo.size(); i++)
        {
            BlockInfo& info = vInfo[i];
            B& block = grid.getBlockCollection()[info.blockID];
            
            double h           = info.h[0];
            double eps        = 1.1 * h;

            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double x[3];
                    info.pos(x, ix, iy);
                    
                    const double rx = x[0] - center[0];
                    const double ry = x[1] - center[0];
                    
                    if(fabs(rx)>=fabs(ry))
                        r = fabs(rx) - radius;
                    else
                        r = fabs(ry) - radius;
                    
                    const double rhs   = (2. * c*c * M_PI*M_PI + 1.) * cos(c * M_PI * x[0]) * cos(c * M_PI * x[1]) ;
                    const double pff   = 0.5 * (1. - tanh(3.*r / eps)) ;
                    
                    block(ix,iy).pff    = max(pff,tau);                                                     // phase field fun.
                    block(ix,iy).chi    = (r<= 0) ? 1. : 0.;                                                //domain characteristic function
                    block(ix,iy).exact  = cos(c * M_PI * x[0]) * cos(c * M_PI * x[1]) * block(ix,iy).chi ;  // exact. sol
                    block(ix,iy).p      = 0;                                                                // init. guess for pressure
                    block(ix,iy).f      = rhs * pff;
                    
                }
            
            grid.getBlockCollection().release(info.blockID);
        }
    }
    else
    {
        vector<BlockInfo> vInfo = grid.getBlocksInfo();;
        
        for(int i=0; i<vInfo.size(); i++)
        {
            BlockInfo& info = vInfo[i];
            B& block = grid.getBlockCollection()[info.blockID];
            
            double h          = info.h[0];
            double eps        = 2.*h;//1.1 * h;

            for(int iz=0; iz<B::sizeZ; iz++)
                for(int iy=0; iy<B::sizeY; iy++)
                    for(int ix=0; ix<B::sizeX; ix++)
                    {
                        double x[3];
                        info.pos(x, ix, iy,iz);
                        
                        const double rx = x[0] - center[0];
                        const double ry = x[1] - center[1];
                        const double rz = x[2] - center[2];
                        
                        if( (fabs(rx)>=fabs(ry)) && (fabs(rx)>=fabs(rz) ))
                            r = fabs(rx) - radius;
                        else if((fabs(ry)>=fabs(rx)) && (fabs(ry)>=fabs(rz)))
                            r = fabs(ry) - radius;
                        else
                            r = fabs(rz) - radius;
                        
                        const double rhs = (3. * c*c * M_PI*M_PI + 1.) * cos(c * M_PI * x[0]) * cos(c * M_PI * x[1]) * cos(c * M_PI * x[2]) ;
                        const double pff = 0.5 * (1. - tanh(3.*r / eps)) ;
                        
                        block(ix,iy,iz).pff     = max(pff,tau);                 // phase field fun
                        block(ix,iy,iz).chi     = (r<= 0) ? 1. : 0.;            //domain characteristic function
                        block(ix,iy,iz).exact   = cos(c * M_PI * x[0]) * cos(c * M_PI * x[1]) * cos(c * M_PI * x[2]) * block(ix,iy,iz).chi ;  // exact sol.
                        block(ix,iy,iz).p       = 0;                            // init. guess for pressure
                        block(ix,iy,iz).f       = rhs * pff;
                        
                    }
            
            grid.getBlockCollection().release(info.blockID);
        }
    }
}


#pragma mark DumpingOutput
void HelmholtzTest:: _dump(int counter)
{
    if (bVTK){
        char filename[256];
        sprintf(filename,"Data_%04d",counter);
        
        if( _DIM == 2){
            IO_VTKNative<W,B, 5,0 > vtkdumper2;
            vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
        }
        else{
            IO_VTKNative3D<W,B, 5,0 > vtkdumper2;
            vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
        }
    }
}


#pragma mark ErrorComputatation
void HelmholtzTest::_computeError()
{
    
    printf("Computing error \n");
    
    Real L1 = 0.0;
    Real L2 = 0.0;
    Real LI = 0.0;
    
    
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    if (_DIM ==2)
    {
        for(int i=0; i<vInfo.size(); i++)
        {
            BlockInfo info = vInfo[i];
            B& block = grid->getBlockCollection()[info.blockID];
            
            const Real dv = vInfo[i].h[0]*vInfo[i].h[1];	 // 2D case
            
                for(int iy=0; iy<B::sizeY; iy++)
                    for(int ix=0; ix<B::sizeX; ix++)
                    {
                        const Real err = fabs( block(ix,iy).p*block(ix,iy).chi - block(ix,iy).exact );
                        
                        L1 += err * dv;
                        L2 += err * err * dv;
                        LI = std::max(LI,err);
                    }
        }
    }
    else
    {
        for(int i=0; i<vInfo.size(); i++)
        {
            BlockInfo info = vInfo[i];
            B& block = grid->getBlockCollection()[info.blockID];
            
            const Real dv = vInfo[i].h[0]*vInfo[i].h[1];	 // 2D case
            
            for(int iz=0; iz<B::sizeZ; iz++)
                for(int iy=0; iy<B::sizeY; iy++)
                    for(int ix=0; ix<B::sizeX; ix++)
                    {
                        const Real err = fabs( block(ix,iy,iz).p*block(ix,iy,iz).chi - block(ix,iy,iz).exact );
                        
                        L1 += err * dv;
                        L2 += err * err * dv;
                        LI = std::max(LI,err);
                    }
        }
    }
    
    L2 =std::sqrt(L2);
    
    const int res = blockSize * blocksPerDimension;
    
    printf("========= PRESSURE ERRORS %d ========\n",res);
    printf("L1, L2, LI: %e %e %e\n",L1,L2,LI);       	
    printf("========= END PRESSURE ===========\n");
    
}


void HelmholtzTest::run()
{
    bool bCG=1;

    if (_DIM == 2)
        helmholtz_solver2D(*grid, bVerbose, bCG);
    else
        helmholtz_solver3D(*grid, bVerbose, bCG);

    
    numberOfIterations++;
    _dump(numberOfIterations);
    _computeError();

    isDone = 1;
    
    printf("**** Dumping done\n");
    printf("\n\n Run Finished \n\n");
}
