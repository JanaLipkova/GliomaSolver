//
//  Glioma_BrainDeformation.cpp
//  DeformationsXcode
//
//  Created by Lipkova on 17/07/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

/*  Dimension-less analysis
 ---------------------------
 Characteristic Units:
 T*   = 1 day
 L*   = brain lenght [cm]
 M*   = 1.04 * L*^3 [g]  , where concentration of brain tissue = 1.04 g/cm^3
 
 Model parametres:
 phi   = tumour volume faction         [ ]
 wm    = tissue concentration          [ M* / L*^3 ]
 p_w   = percantage of wm              [ ]
 rho   = proliferation rate            [ 1/T* ]
 p     = pressure induced by tumour    [Pa] = [kg / (m s^2) ]  ( 1mmHg = 133.32 Pa)
 kappa = relaxation                    [ 1/(Pa.s)] = [kg / (m.s)]
 M     = mobility/hydraulic conduct.   [ m^2 / (Pa.s) ] = [ m^3 s / kg ]
 */

#include "Glioma_BrainDeformation.h"
#include <mpi.h>


static int maxStencil[2][3] = {
    -3, -3, -3,
    +4, +4, +4
};


Glioma_BrainDeformation::Glioma_BrainDeformation(int argc, const char ** argv): parser(argc, argv), helmholtz_solver3D(argc,argv)
{
    bVerbose  = parser("-verbose").asBool(1);
    bProfiler = parser("-profiler").asBool(1);
    
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("//////////////////       Glioma Brain Deformations Model        ////////////////\n");
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("RD INIT! nThreads=%d, blockSize=%d Wavelets=w%s (blocksPerDimension=%d, maxLevel=%d)\n", nThreads, blockSize, "w", blocksPerDimension, maxLevel);
    
    refiner		= new Refiner_SpaceExtension();
    compressor	= new Compressor();
    Environment::setup();
    
    grid = new Grid<W,B>(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil);
    grid->setCompressor(compressor);
    grid->setRefiner(refiner);
    stSorter.connect(*grid);
    
    pID =  parser("-pID").asInt();
    L = 1;
    
    _ic(*grid, pID, L);
    
    if(parser("-bDumpIC").asBool(0))
        _dump(0);
    
    MPI_Init(&argc, (char ***)&argv);
    
    isDone              = false;
    whenToWriteOffset	= parser("-dumpfreq").asDouble();
    whenToWrite			= whenToWriteOffset;
    numberOfIterations	= 0;
}

Glioma_BrainDeformation::~Glioma_BrainDeformation()
{
    MPI_Finalize();
    std::cout << "------Adios muchachos------" << std::endl;
}



#pragma mark InitialConditions
void Glioma_BrainDeformation:: _ic(Grid<W,B>& grid, int pID, Real& L)
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
    sprintf(anatomy, "%s_PFF.dat", patientFolder);
    MatrixD3D PFF(anatomy);
    
    
    int brainSizeX = (int) GM.getSizeX();
    int brainSizeY = (int) GM.getSizeY();
    int brainSizeZ = (int) GM.getSizeZ();
    printf("brainSizeX=%d, brainSizeY=%d, brainSizeZ= %d \n", brainSizeX, brainSizeY, brainSizeZ);

    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    L                = brainSizeMax * 0.1;   // voxel spacing 1mm, convert from mm to cm  // L = 22.9 cm

    double brainHx = 1.0 / ((double)(brainSizeMax-1)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax-1)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax-1)); // should be w.r.t. longest dimension for correct aspect ratio
    
    // tumour parameters
    const Real tumorRadius = 0.01;
    const Real smooth_sup  = 3.;		// suppor of smoothening, over how many gp to smooth
    const Real c[3] = { 0.6, 0.7, 0.5};
    
    double pGM, pWM, pCSF, pPFF;
    const double tau = 1.e-10;

#pragma mark Initialization
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        const float h = 1./128; //vInfo[0].h[0]; i.e.ersusres same initialisaition for grids of different resolution
        const float iw = 1./(smooth_sup * h);   // width of smoothening => now it is over two grid points
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    Real x[3];
                    info.pos(x, ix, iy, iz);
                    
                    /* Anatomy */
                    int mappedBrainX = (int)floor( x[0] / brainHx  );
                    int mappedBrainY = (int)floor( x[1] / brainHy  );
                    int mappedBrainZ = (int)floor( x[2] / brainHz  );
                    
                    // aspect ratio correction
                    mappedBrainX -= (int) ( (brainSizeMax - brainSizeX) * 0.5);
                    mappedBrainY -= (int) ( (brainSizeMax - brainSizeY) * 0.5);
                    mappedBrainZ -= (int) ( (brainSizeMax - brainSizeZ) * 0.5);
                    
                    if ( (mappedBrainX >= 0 && mappedBrainX < brainSizeX) && (mappedBrainY >= 0 && mappedBrainY < brainSizeY) && (mappedBrainZ >= 0 && mappedBrainZ < brainSizeZ) )
                    {
                        pGM     = GM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        pWM     = WM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        pCSF    = CSF(mappedBrainX,mappedBrainY,mappedBrainZ);
                        pPFF    = PFF(mappedBrainX,mappedBrainY,mappedBrainZ);
                        
                        // remove background signal
                        pGM =  (pGM  < 1e-05) ? 0. : pGM;
                        pWM  = (pWM  < 1e-05) ? 0. : pWM;
                        pCSF = (pCSF < 1e-05) ? 0. : pCSF;
                        pPFF = (pPFF < 1e-05) ? 0. : pPFF;
                        
                        double all = pGM + pWM + pCSF;
                        if(all > 0)
                        {
                            // normalize
                            pGM    = pGM  / all;
                            pWM    = pWM  / all;
                            pCSF   = pCSF / all;
                            
                            //pCSF = ( pCSF > 0.1 ) ? 1. : pCSF;  // threashold to ensure hemispehre separation
                            //block(ix,iy,iz).p_csf = pCSF;
                            
                            //if(pCSF  < 1.)
                            {
                                block(ix,iy,iz).p_g   = pGM  / (pCSF + pWM + pGM);
                                block(ix,iy,iz).p_w   = pWM  / (pCSF + pWM + pGM);
                                block(ix,iy,iz).p_csf = pCSF / (pCSF + pWM + pGM);
                            }
                        }
                        
                        //tissue concetration
                        block(ix,iy,iz).wm  = block(ix,iy,iz).p_w;
                        block(ix,iy,iz).gm  = block(ix,iy,iz).p_g;
                        block(ix,iy,iz).csf = block(ix,iy,iz).p_csf;
                        
                        // tumor
                        const Real p[3] = {x[0] - c[0], x[1] - c[1], x[2] - c[2]};
                        const Real dist = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);    // distance of curent voxel from tumor center
                        const Real psi  = (dist - tumorRadius)*iw;
                        
                        bool bTissue = ( block(ix,iy,iz).p_w + block(ix,iy,iz).p_g > 0. ) ? 1 : 0 ;
                        
                        if ((psi < -1) && bTissue )		// we are in tumor
                            block(ix,iy,iz).phi = 1.0;
                        else if(( (-1 <= psi) && (psi <= 1) )&& (bTissue) )
                            block(ix,iy,iz).phi = 1.0 * 0.5 * (1 - psi - sin(M_PI*psi)/(M_PI));
                        else
                            block(ix,iy,iz).phi = 0.0;
                        
                        // auxiliary functions
                        block(ix,iy,iz).pff = max(tau, pPFF);            //phase field func.
                        block(ix,iy,iz).chi = (pPFF >= 0.5) ? 1. : 0.;   // domain char. func
                    }
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}



#pragma mark TimeStepEstimation
double Glioma_BrainDeformation::_estimate_dt(double Diff_dt, double CFL)
{
    const Real maxvel = _compute_maxvel();
    const Real max_dx = (1./B::sizeX)*pow(0.5,grid->getCurrentMinLevel());
    const Real min_dx = (1./B::sizeX)*pow(0.5,grid->getCurrentMaxLevel());
    
    const double largest_dt = CFL * max_dx/(maxvel * _DIM);
    const double smallest_dt = CFL * min_dx/(maxvel * _DIM);
    
    assert(largest_dt >= smallest_dt);
    
    if (largest_dt <= Diff_dt)
        printf("Advection dominated time step Adts = %f, ADLdt=%f, Ddt=%f, returning=%f \n", smallest_dt, largest_dt, Diff_dt, (largest_dt < 1.e-10) ? Diff_dt : min(largest_dt, Diff_dt) );
    
    return (smallest_dt < 1.e-10) ? Diff_dt : min(largest_dt, Diff_dt);
}

Real Glioma_BrainDeformation::_compute_maxvel()
{
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    map<int, Real> velocities;
    for(vector<BlockInfo>::iterator it=vInfo.begin(); it!=vInfo.end(); it++)
        velocities[it->blockID] = 0;
    
    
    GetUMax<_DIM> get_velocities(velocities);
    BlockProcessing::process(vInfo, collecton, get_velocities);
    
    Real maxvel = 0;
    for(map< int, Real>::iterator it=velocities.begin(); it!=velocities.end(); it++)
        maxvel = max(maxvel, it->second);
    
    return maxvel;
}

#pragma mark MainOperators
void Glioma_BrainDeformation::_computePressureSource(const int nParallelGranularity,const Real rho)
{
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    PressureSourceOperator<_DIM> pressureSource(rho);
    BlockProcessing::process(vInfo, collecton, pressureSource, nParallelGranularity);
}

void Glioma_BrainDeformation::_computeVelocities(BoundaryInfo* boundaryInfo, const bool bMobility, std::vector<Real>* mobility=NULL)
{
    if(bMobility)
        assert (mobility!=NULL);

    Real  Mcsf   = (bMobility) ? (*mobility)[0] : 1.;
    Real  Mwm    = (bMobility) ? (*mobility)[1] : 1.;
    Real  Mgm    = (bMobility) ? (*mobility)[2] : 1.;

    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();

    PressureGradient<_DIM> grad(bMobility, Mwm, Mgm, Mcsf );
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, grad);
}

void Glioma_BrainDeformation::_reactionDiffusionStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, const Real Dw, const Real Dg, const Real rho)
{
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();

    ReactionDiffusionOperator<_DIM>  rhs(Dw,Dg,rho);
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
}

void Glioma_BrainDeformation::_advectionConvectionStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, double dt)
{
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();

    TissueTumorAdvectionWeno5 <_DIM> rhsA;      // v grad(E) & v grad(tumor)
    TissueConvection          <_DIM> rhsB;      // E div(v)
    TimeUpdate                <_DIM> update(dt);

    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhsA);
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhsB);
    BlockProcessing::process(vInfo, collecton, update, nParallelGranularity);
}


#pragma mark ingOutput
void Glioma_BrainDeformation:: _dump(int counter)
{
    if (parser("-vtk").asBool())
    {
        char filename[256];
        sprintf(filename,"%P%02d_data_%04d",pID,counter);
        
        IO_VTKNative3D<W,B, 16,0 > vtkdumper2;
        vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
    }
}


void Glioma_BrainDeformation::run()
{
    const int nParallelGranularity	= (grid->getBlocksInfo().size()<=8 ? 1 : 4);
    BoundaryInfo* boundaryInfo		= &grid->getBoundaryInfo();
    
    // model parameters
    const double CFL    = parser("-CFL").asDouble(0.8);
    const double tend   = parser("-tend").asDouble();  // [day]
    const double rho    = parser("-rho").asDouble();   // [1/day]
    double Dw           = parser("-Dw").asDouble();
    Dw = Dw/(L*L);
    double Dg = 0.1*Dw;
    
    double  h            = 1./(blockSize*blocksPerDimension);
    double  Diff_dt     = (h*h)/( 2.*_DIM * max(Dw, Dg) );
    double  dt          = Diff_dt;
    int     iCounter    = 1;
    double  t           = 0.;
    
    printf("Dg=%f, Dw=%f, dt= %f, rho=%f \n", Dg, Dw, dt, rho);
    
    // relaxation parameters
    const bool bCG          = 1;
    const bool bRelaxation  = 1;
    
    const Real kCSF = parser("-kCSF").asDouble(); //  [m.s / kg] * Kstar = [L* T* / M* ]
    const Real kWM  = parser("-kWM").asDouble();  //  [m.s / kg] * Kstar = [L* T* / M* ]
    const Real kGM  = parser("-kGM").asDouble();  //  [m.s / kg] * Kstar = [L* T* / M* ]
    
    vector<Real> kappa;
    kappa.push_back(kCSF);
    kappa.push_back(kWM);
    kappa.push_back(kGM);
    
    printf("kCSF=%f, kWM=%f, kGM=%f, bVerbose=%i \n", kCSF, kWM, kGM, bVerbose);

    // mobility parameters
    const bool bMobility = parser("-bMobility").asBool(0);
    vector<Real> mobility;   // CSF, WM, GM
//
//    if(bMobility)
//    {
//        const Real mCSF = parser("-mTissue").asDouble(); // [m^3 s/ kg] * Mstar = [L*^3 T* / M*]
//        const Real mWM = mCSF;
//        const Real mGM = mCSF;
//        
//        mobility.push_back(mCSF);
//        mobility.push_back(mWM);
    //        mobility.push_back(mGM);
    //
    //        printf("mCSF=%f, mWM=%f, mGM=%f\n", mCSF, mWM, mGM);
    //    }
    
    
    while (t <= tend)
    {
        dt = _estimate_dt(Diff_dt, CFL);
        _computePressureSource(nParallelGranularity,rho);
        helmholtz_solver3D(*grid, bVerbose, bCG, bRelaxation, &kappa, bMobility, &mobility);
        _computeVelocities(boundaryInfo, bMobility, &mobility);
        _reactionDiffusionStep(boundaryInfo, nParallelGranularity, Dw, Dg, rho);
        _advectionConvectionStep(boundaryInfo, nParallelGranularity, dt);
        
        t                   += dt   ;
        numberOfIterations  ++      ;
        
        if ( t >= ((double)(whenToWrite)) )
        {
            _dump(iCounter);
            whenToWrite = whenToWrite + whenToWriteOffset;
            iCounter++;
        }
    }
    
    _dump(iCounter);
    
    if(bVerbose) profiler.printSummary();
    if(bVerbose) printf("**** Dumping done\n");
    if(bVerbose) printf("\n\n Run Finished \n\n");
}

