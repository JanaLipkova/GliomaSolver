//
//  Glioma_BrainDeformationTimeRelaxation.cpp
//  GliomaSolverXcode
//
//  Created by Lipkova on 08/08/18.
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
 beta  = compresibitliy                [ 1/(Pa.s)] = [kg / (m.s)]
 M     = mobility/hydraulic conduct.   [ m^2 / (Pa.s) ] = [ m^3 s / kg ]
 */

#include "Glioma_BrainDeformationTimeRelaxation.h"

// max stencil size from Weno5
static int maxStencil[2][3] = {
    -3, -3, -3,
    +4, +4, +4
};



Glioma_BrainDeformationTimeRelaxation::Glioma_BrainDeformationTimeRelaxation(int argc, const char ** argv): parser(argc, argv)
{
    bVerbose  = parser("-verbose").asBool(1);
    bProfiler = parser("-profiler").asBool(1);
    
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("/////   Glioma Brain Deformation Model with Time Pressure Relaxation   /////////\n");
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("RD INIT! nThreads=%d, blockSize=%d Wavelets=w%s (blocksPerDimension=%d, maxLevel=%d)\n", nThreads, blockSize, "w", blocksPerDimension, maxLevel);
    
    refiner		= new Refiner_SpaceExtension();
    compressor	= new Compressor();
    Environment::setup();
    
    grid = new Grid<W,B>(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil);
    grid->setCompressor(compressor);
    grid->setRefiner(refiner);
    stSorter.connect(*grid);
    
    bAdaptivity = parser("-adaptive").asBool();
    bVTK        = parser("-vtk").asBool();
    pID         = parser("-pID").asInt();
    L = 1;
    
    if(pID == 100)
        _icSphere3Parts(*grid, L);
    else
        _ic(*grid, pID, L);
    
    if(parser("-bDumpIC").asBool(0))
        _dump(0);
    
    isDone              = false;
    whenToWriteOffset	= parser("-dumpfreq").asDouble();
    whenToWrite			= whenToWriteOffset;
    numberOfIterations	= 0;
}

Glioma_BrainDeformationTimeRelaxation::~Glioma_BrainDeformationTimeRelaxation()
{
    std::cout << "------Adios muchachos------" << std::endl;
}



#pragma mark InitialConditions
/* 3 componten structure: WM, GM, CSF + Tumor*/
void Glioma_BrainDeformationTimeRelaxation::_icSphere3Parts(Grid<W,B>& grid, Real& L)
{
    std::cout <<" Test case: Sphere with 3 components"<< std::endl;
    
    const double AnatmoyRadius	= 0.4;
    const double tumorRad		= 0.1;    // tumor radius
    const double smooth_sup		= 3.;     // support, over how many grid points to smooth
    
    const double tau        = 1.e-10;     // cut of phase field function on LHS
    const Real center[3]   = {0.5, 0.5, 0.5};
    L = 20;

    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        double h    = vInfo[0].h[0];
        double eps  = 1.1 * h;             // phase field fun. smoothening
        double iw = 1./(smooth_sup * h);   // width of tumor smoothening
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    Real x[3];
                    info.pos(x, ix, iy,iz);
                    
                    const Real p[3]  = { x[0] - center[0], x[1] - center[1], x[2] - center[2]};
                    const Real dist  = sqrt( p[0]*p[0] + p[1]*p[1] +p[2]*p[2] );
                    const Real r     = dist - AnatmoyRadius;    // sign distance function
                    const Real tmp   = (dist - tumorRad) * iw;
                    
                    // tumor
                    if (tmp < -1)
                        block(ix,iy,iz).phi = 1.0;
                    else if( (-1 <= tmp) && (tmp <= 1) )
                        block(ix,iy,iz).phi = 0.5 * (1 - tmp - sin(M_PI * tmp) / (M_PI) );
                    
                    // compontents of pressure eq.
                    const double pff = 0.5 * (1. - tanh(3.*r / eps)) ;
                    block(ix,iy,iz).pff = max(pff,tau);
                    block(ix,iy,iz).chi = (r<= 0) ? 1. : 0.;
                    block(ix,iy,iz).p   = 0.;
                    
                    // anatomy
                    const Real theta = atan2((x[1]-0.5), (x[0]-0.5) ) * 180. / (M_PI);
                    
                    if (block(ix,iy,iz).pff > tau)
                    {
                        if ((0<=theta)&&(theta < 120))
                            block(ix,iy,iz).p_w = 1.;
                        else if ((-120 <= theta)&&(theta < 0))
                            block(ix,iy,iz).p_g = 1.;
                        else
                            block(ix,iy,iz).p_csf = 1.;
                        
                        block(ix,iy,iz).wm  = block(ix,iy,iz).p_w;
                        block(ix,iy,iz).gm  = block(ix,iy,iz).p_g;
                        block(ix,iy,iz).csf = block(ix,iy,iz).p_csf;
                        
                    }
                }
        
        grid.getBlockCollection().release(info.blockID);
    }
    
}



void Glioma_BrainDeformationTimeRelaxation:: _ic(Grid<W,B>& grid, int pID, Real& L)
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
    sprintf(anatomy, "%s_Mask.dat", patientFolder);
    MatrixD3D Mask(anatomy); // brain mask
    
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
    
    double pGM, pWM, pCSF, pMask;
    
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
                        pMask   = Mask(mappedBrainX,mappedBrainY,mappedBrainZ);
                        
                        // remove background signal
                        pGM  =  (pGM  < 1e-05)  ? 0. : pGM;
                        pWM   = (pWM  < 1e-05)  ? 0. : pWM;
                        pCSF  = (pCSF < 1e-05)  ? 0. : pCSF;
                        pMask = (pMask < 1e-05) ? 0. : pMask;
                        
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
                        block(ix,iy,iz).chi = pMask;   // domain char. func
                    }
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}


#pragma mark TimeStepEstimation
double Glioma_BrainDeformationTimeRelaxation::_estimate_dt(double dt, double CFL)
{
    const Real maxvel = _compute_maxvel();
    const Real max_dx = (1./B::sizeX)*pow(0.5,grid->getCurrentMinLevel());
    const Real min_dx = (1./B::sizeX)*pow(0.5,grid->getCurrentMaxLevel());
    
    const double largest_dt  = CFL * max_dx/(maxvel * _DIM);
    const double smallest_dt = CFL * min_dx/(maxvel * _DIM);
    
    assert(largest_dt >= smallest_dt);
    
    
    if (largest_dt <= dt)
        printf("Advection dominated time step Adts = %f, ADLdt=%f, Ddt=%f \n", smallest_dt, largest_dt, dt);
    
    return  (maxvel < 1.0e-10 ) ? dt : min(largest_dt, dt);
}

Real Glioma_BrainDeformationTimeRelaxation::_compute_maxvel()
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

#pragma mark Operators
void Glioma_BrainDeformationTimeRelaxation::_computeVelocities(BoundaryInfo* boundaryInfo, const bool bMobility, std::vector<Real>* mobility=NULL)
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

#pragma mark Operators
void Glioma_BrainDeformationTimeRelaxation::_pressureStep(BoundaryInfo* boundaryInfo,std::vector<Real> mobolity, std::vector<Real> beta, const Real rho)
{
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    PressureTimeRelaxationOperator<_DIM>  rhs(mobolity[0],mobolity[1],mobolity[2],beta[0], beta[1], beta[2], rho);
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
}


void Glioma_BrainDeformationTimeRelaxation::_reactionDiffusionStep(BoundaryInfo* boundaryInfo, const Real Dw, const Real Dg, const Real rho)
{
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    ReactionDiffusionOperator<_DIM>  rhs(Dw,Dg,rho);
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
}

void Glioma_BrainDeformationTimeRelaxation::_advectionConvectionStep(BoundaryInfo* boundaryInfo)
{
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    TissueTumorAdvectionWeno5 <_DIM> rhsA;      // v grad(E) & v grad(tumor)
    TissueConvection          <_DIM> rhsB;      // E div(v)
    
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhsA);
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhsB);
}

void Glioma_BrainDeformationTimeRelaxation::_timeUpdate(const int nParallelGranularity, const Real dt )
{
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    TimeUpdate        <_DIM> update(dt);
    PressureTimeUpdate<_DIM> pressureUpdate(dt);
    
    BlockProcessing::process(vInfo, collecton, update, nParallelGranularity);
    BlockProcessing::process(vInfo, collecton, pressureUpdate, nParallelGranularity);
}

#pragma mark ingOutput
void Glioma_BrainDeformationTimeRelaxation:: _dump(int counter)
{
    if (bVTK)
    {
        char filename[256];
        sprintf(filename,"%P%02d_data_%04d",pID,counter);
        
        IO_VTKNative3D<W,B, 14,0 > vtkdumper2;
        vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
    }
}


void Glioma_BrainDeformationTimeRelaxation::run()
{
    const int nParallelGranularity	= (grid->getBlocksInfo().size()<=8 ? 1 : 4);
    BoundaryInfo* boundaryInfo		= &grid->getBoundaryInfo();
    
    // model parameters
    const double CFL    = parser("-CFL").asDouble(0.8);
    const double tend   = parser("-Tend").asDouble(300);  // [day]
    const double rho    = parser("-rho").asDouble(0.025);   // [1/day]
    double Dw           = parser("-Dw").asDouble(0.0013);
    Dw = Dw/(L*L);
    double Dg = 0.1*Dw;
    
    // compresibility parameters
    vector<Real> beta;
    beta.push_back(parser("-betaCSF").asDouble(1.));
    beta.push_back(parser("-betaWM").asDouble(1.));
    beta.push_back(parser("-betaGM").asDouble(1.));
    
    bool bMobility = 1;
    vector<Real> mobility;
    mobility.push_back(parser("-mobCSF").asDouble(1.) / (L*L*L) );
    mobility.push_back(parser("-mobWM").asDouble(1.)  / (L*L*L)) ;
    mobility.push_back(parser("-mobGM").asDouble(1.)  / (L*L*L)) ;
    
    
    
    
    double  h         = 1./(blockSize*blocksPerDimension);
    double  dt_dif      = (h*h)/( 2.*_DIM * max(Dw, Dg) );
    double  dt_mob      = (h*h)/( 2.*_DIM * max( max(mobility[0],mobility[1]) ,mobility[2]) );
    double  dt          = min(dt_dif,dt_mob);
    int     iCounter    = 1;
    double  t           = 0.;
    
    if(bVerbose) printf("Dg=%e, Dw=%e, rho=%e, dt_D= %e\n", Dg, Dw, rho, dt_dif);
    if(bVerbose) printf("mCSF=%e, mWM=%e, mGM=%e, dt_p=%e \n",  mobility[0], mobility[1], mobility[2], dt_mob);
    if(bVerbose) printf("bCSF=%e, bWM=%e, bGM=%e \n", beta[0], beta[1], beta[2]);
    
    
    while (t <= tend)
    {
        if(bProfiler) profiler.getAgent("velocity").start();
        _computeVelocities(boundaryInfo, bMobility, &mobility);
        if(bProfiler) profiler.getAgent("velocity").stop();
        
        if(bProfiler) profiler.getAgent("pressure").start();
        _pressureStep(boundaryInfo, mobility, beta, rho);
        if(bProfiler) profiler.getAgent("pressure").stop();
        
        if(bProfiler) profiler.getAgent("RD").start();
        _reactionDiffusionStep(boundaryInfo, Dw, Dg, rho);
        if(bProfiler) profiler.getAgent("RD").stop();
        
        if(bProfiler) profiler.getAgent("AdvectConv").start();
        _advectionConvectionStep(boundaryInfo);
        if(bProfiler) profiler.getAgent("AdvectConv").stop();
        
        if(bProfiler) profiler.getAgent("dt").start();
        _estimate_dt( min(dt_dif, dt_mob), CFL);
        if(bProfiler) profiler.getAgent("dt").stop();
        
        if(bProfiler) profiler.getAgent("TimeIntegr").start();
        _timeUpdate(nParallelGranularity, dt );
        if(bProfiler) profiler.getAgent("TimeIntegr").stop();
        
        t                   += dt   ;
        numberOfIterations  ++      ;
        
        if ( t >= ((double)(whenToWrite)) )
        {
            if(bAdaptivity)
            {
                Science::AutomaticRefinement	<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);
                Science::AutomaticCompression	<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler);
            }
            
            if(bProfiler) profiler.getAgent("Dump").start();
            _dump(iCounter);
            if(bProfiler) profiler.getAgent("Dump").stop();
            
            if(bVerbose) printf("Dumping data at time t=%f\n", t);
            whenToWrite = whenToWrite + whenToWriteOffset;
            iCounter++;
        }
    }
    
    _dump(iCounter);
    
    if(bVerbose) profiler.printSummary();
    if(bVerbose) printf("**** Dumping done\n");
    if(bVerbose) printf("\n\n Run Finished \n\n");
}

