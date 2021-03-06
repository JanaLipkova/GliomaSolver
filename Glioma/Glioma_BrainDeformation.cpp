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


Glioma_BrainDeformation::Glioma_BrainDeformation(int argc, const char ** argv): parser(argc, argv)
{
    bVerbose  = parser("-verbose").asBool(1);
    bProfiler = parser("-profiler").asBool(1);
    bVTK      = parser("-vtk").asBool();
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    if(rank==0){
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("//////////////////       Glioma Brain Deformations Model        ////////////////\n");
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("Set up: blockSize=%d Wavelets=w%s (blocksPerDimension=%d, maxLevel=%d)\n", blockSize, "w", blocksPerDimension, maxLevel);
    }
    
    refiner		= new Refiner_SpaceExtension();
    compressor	= new Compressor();
    grid        = new Grid<W,B>(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil);
    grid->setCompressor(compressor);
    grid->setRefiner(refiner);
    stSorter.connect(*grid);
    
    
    PatientFileName = parser("-PatFileName").asString();
    L               = 1;
        
    if(parser("-ICtype").asInt(1)==0)
        _icSphere3Parts(*grid,rank,L);
    else
    {
        tumor_ic[0] = parser("-icx").asDouble(0.6);
        tumor_ic[1] = parser("-icy").asDouble(0.7);
        tumor_ic[2] = parser("-icz").asDouble(0.5);
        _ic(*grid, rank, PatientFileName, L, tumor_ic);
    }
    
    if((parser("-bDumpIC").asBool(0))&&(rank==0)) _dump(0);
    
    isDone              = false;
    whenToWriteOffset	= parser("-dumpfreq").asDouble();
    whenToWrite			= whenToWriteOffset;
    numberOfIterations	= 0;
}

Glioma_BrainDeformation::~Glioma_BrainDeformation()
{
    if(rank==0) std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions
/* 3 componten structure: WM, GM, CSF + Tumor*/
void Glioma_BrainDeformation::_icSphere3Parts(Grid<W,B>& grid, int rank, Real& L)
{
    if(rank==0) std::cout <<" Test case: Sphere with 3 components \n"<< std::endl;
    
    const double AnatmoyRadius	= 0.4;
    const double TumorRadius	= 0.1;    // tumor radius
    const double smooth_sup		= 3.;     // support, over how many grid points to smooth
    
    const double tau        = 1.e-10;     // cut of phase field function on LHS
    const Real AnatCenter[3]   = {0.5, 0.5, 0.5};
    const Real TumCenter[3]    = {0.62, 0.5, 0.5};

    L = 20;
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
#pragma omp parallel for
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        double h    = vInfo[0].h[0];
        double eps  = 3.*h;//1.1 * h;             // phase field fun. smoothening
        double iw = 1./(smooth_sup * h);   // width of tumor smoothening
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    Real x[3];
                    info.pos(x, ix, iy,iz);
                    
                    const Real pAnat[3]  = { x[0] - AnatCenter[0], x[1] - AnatCenter[1], x[2] - AnatCenter[2]};
                    const Real distAnat  = sqrt( pAnat[0]*pAnat[0] + pAnat[1]*pAnat[1] +    pAnat[2]*pAnat[2]);
                    const Real rAnat     = distAnat - AnatmoyRadius;    // sign distance function
                    
                    const Real pTum[3]  = { x[0] - TumCenter[0], x[1] - TumCenter[1], x[2] - TumCenter[2]};
                    const Real distTum  = sqrt( pTum[0]*pTum[0] + pTum[1]*pTum[1] + pTum[2]*pTum[2] );
                    const Real tmp      = (distTum - TumorRadius) * iw;
                    
                    // tumor
                    if (tmp < -1)
                        block(ix,iy,iz).phi = 1.0;
                    else if( (-1 <= tmp) && (tmp <= 1) )
                        block(ix,iy,iz).phi = 0.5 * (1 - tmp - sin(M_PI * tmp) / (M_PI) );
                    
                    // compontents of pressure eq.
                    const double pff = 0.5 * (1. - tanh(3.*rAnat / eps)) ;
                    block(ix,iy,iz).pff = max(pff,tau);
                    block(ix,iy,iz).chi = (rAnat<= 0) ? 1. : 0.;
                    
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
                    }
                    
                    //tissue concetration
                    block(ix,iy,iz).wm  = block(ix,iy,iz).p_w;
                    block(ix,iy,iz).gm  = block(ix,iy,iz).p_g;
                    block(ix,iy,iz).csf = block(ix,iy,iz).p_csf;
                }
        
        grid.getBlockCollection().release(info.blockID);
    }
    
}




void Glioma_BrainDeformation:: _ic(Grid<W,B>& grid, int rank, string PatientFileName, Real& L, Real tumor_ic[3] )
{
    if(rank==0) printf("Reading data from file: %s \n", PatientFileName.c_str());
    
    char anatomy      [200];
    sprintf(anatomy, "%sGM.dat", PatientFileName.c_str() );
    MatrixD3D GM(anatomy);
    sprintf(anatomy, "%sWM.dat", PatientFileName.c_str());
    MatrixD3D WM(anatomy);
    sprintf(anatomy, "%sCSF.dat", PatientFileName.c_str());
    MatrixD3D CSF(anatomy);
    sprintf(anatomy, "%sPFF.dat", PatientFileName.c_str());
    MatrixD3D PFF(anatomy);
    
    int brainSizeX = (int) GM.getSizeX();
    int brainSizeY = (int) GM.getSizeY();
    int brainSizeZ = (int) GM.getSizeZ();
    if(rank==0) printf("brainSizeX=%d, brainSizeY=%d, brainSizeZ= %d \n", brainSizeX, brainSizeY, brainSizeZ);
    if(rank==0) printf("PFFSizeX=%d, PFFSizeY=%d, PFFSizeZ= %d \n", (int) PFF.getSizeX(), (int) PFF.getSizeY(), (int)  PFF.getSizeZ());
    assert(GM.getSizeX() == PFF.getSizeX());
    assert(GM.getSizeY() == PFF.getSizeY());
    assert(GM.getSizeZ() == PFF.getSizeZ());

    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    L                = brainSizeMax * 0.1;   // voxel spacing 1mm, convert from mm to cm  // L = 22.9 cm
    
    double brainHx = 1.0 / ((double)(brainSizeMax-1)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax-1)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax-1)); // should be w.r.t. longest dimension for correct aspect ratio
    
    // tumour parameters
    const Real tumorRadius = 3./(blockSize * blocksPerDimension);//0.01;//0.005;//0.01;
    const Real smooth_sup  = 4.;		// suppor of smoothening, over how many gp to smooth
    const Real h           = 1./(blockSize*blocksPerDimension);    // use fixed h, for same IC smoothening at all resolutions
    const Real iw          = 1./(smooth_sup * h); // widht of smoothening
    
    double pGM, pWM, pCSF, pPFF;
    const double tau = 1.e-10;
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
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
                    }
                    else
                        pGM = pWM = pCSF = pPFF = 0.;
                    
                    double all = pGM + pWM + pCSF;
                    if(all > 0.1)
                    {
                        // enhance fluid -> ensures hemisphere separaiton
                        pCSF = ( pCSF > 0.1 ) ? 1. : pCSF;
                        pWM  = ( pCSF > 0.1 ) ? 0. : pWM;
                        pGM  = ( pCSF > 0.1 ) ? 0. : pGM;
                        
                        // alows low level tissue mixing -> smooth interface
                        block(ix,iy,iz).p_g   = pGM  / (pCSF + pWM + pGM);
                        block(ix,iy,iz).p_w   = pWM  / (pCSF + pWM + pGM);
                        block(ix,iy,iz).p_csf = pCSF / (pCSF + pWM + pGM);
                        
                    }
                    
                    //tissue concetration
                    block(ix,iy,iz).wm  = block(ix,iy,iz).p_w;
                    block(ix,iy,iz).gm  = block(ix,iy,iz).p_g;
                    block(ix,iy,iz).csf = block(ix,iy,iz).p_csf;
                    
                    // tumor
                    const Real p[3] = {x[0] - tumor_ic[0], x[1] - tumor_ic[1], x[2] - tumor_ic[2]};
                    const Real dist = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);    // distance of curent voxel from tumor center
                    const Real psi  = (dist - tumorRadius)*iw;
                    
                    bool bTissue = ( block(ix,iy,iz).p_w + block(ix,iy,iz).p_g > 0. ) ? 1 : 0 ;
                    
                    if ((psi < -1) && bTissue )		// we are in tumor
                        block(ix,iy,iz).phi = 1.0;
                    else if(( (-1 <= psi) && (psi <= 1) )&& (bTissue) )
                        block(ix,iy,iz).phi = 1.0 * 0.5 * (1 - psi - sin(M_PI*psi)/(M_PI));
                    else
                        block(ix,iy,iz).phi = 0.0;
                    
                    // avoid rounding errors
                    block(ix,iy,iz).phi = (block(ix,iy,iz).phi < 0.) ? 0. : block(ix,iy,iz).phi;
                    
                    // auxiliary functions
                    block(ix,iy,iz).pff = max(tau, pPFF);
                    block(ix,iy,iz).chi = (pPFF >= 0.5) ? 1. : 0.;   // domain char. func
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
    
    if ((largest_dt <= Diff_dt))
        printf("Advection dominated time ADdt=%f, Ddt=%f \n", largest_dt, Diff_dt );
    
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

void Glioma_BrainDeformation::_reactionDiffusionStep(BoundaryInfo* boundaryInfo, const Real Dw, const Real Dg, const Real rho)
{
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    ReactionDiffusionOperator<_DIM>  rhs(Dw,Dg,rho);
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
}

void Glioma_BrainDeformation::_advectionConvectionStep(BoundaryInfo* boundaryInfo)
{
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    TissueTumorAdvectionWeno5 <_DIM> rhsA;      // v grad(E) & v grad(tumor)
    TissueConvection          <_DIM> rhsB;      // E div(v)
    
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhsA);
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhsB);
}


void Glioma_BrainDeformation::_timeUpdate(const int nParallelGranularity, double dt)
{
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();

    TimeUpdate                <_DIM> update(dt);
    BlockProcessing::process(vInfo, collecton, update, nParallelGranularity);
}

#pragma mark ingOutput
void Glioma_BrainDeformation:: _dump(int counter)
{
    if (bVTK)
    {
        if(bVerbose) printf("dumping data %04d \n", counter);
        char filename[256];
        sprintf(filename,"Data_%04d",counter);
        
        IO_VTKNative3D<W,B, 15,0 > vtkdumper2;
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
    const double rho    = parser("-rho").asDouble(1);   // [1/day]
    double Dw           = parser("-Dw").asDouble();
    Dw = Dw/(L*L);
    double Dg = 0.1*Dw;
    
    double  h            = 1./(blockSize*blocksPerDimension);
    double  Diff_dt     = (h*h)/( 2.*_DIM * max(Dw, Dg) );
    double  dt          = Diff_dt;
    int     iCounter    = 1;
    double  t           = 0.;
    
    if(rank==0) printf("Dg=%e, Dw=%e, dt= %e, rho=%e, L=%e \n", Dg, Dw, dt, rho, L);
    
    // relaxation and mobility parameters
    const bool bRelaxation  = 1;
    const bool bMobility    = parser("-bMobility").asBool(0);

    vector<Real> kappa(3,1);
    vector<Real> mobility(3,1);   // CSF, WM, GM

    if(bRelaxation){
        kappa[0] = (Real) parser("-kCSF").asDouble(0.2) * L * L;
        kappa[1] = (Real) parser("-kWM").asDouble(0.02)  * L * L;
        kappa[2] = (Real) parser("-kGM").asDouble(0.002)  * L * L;
       if(rank==0) printf("kCSF=%f, kWM=%f, kGM=%f \n", kappa[0], kappa[1], kappa[2]);
    }
    
    if(bMobility){
        mobility[0] = (Real) parser("-mCSF").asDouble(1);  //[m^3 s/ kg] * Mstar = [L*^3 T* / M*]
        mobility[1] = (Real) parser("-mWM").asDouble(1);
        mobility[2] = (Real) parser("-mGM").asDouble(1);
        if(rank==0) printf("mCSF=%f, mWM=%f, mGM=%f\n", mobility[0], mobility[1], mobility[2]);
    }
    
    /* Set up Hypre solver */
    bool bCG = 1;
    HelmholtzSolver3D_Hypre_MPI helmholtz_solver;
    helmholtz_solver.setup_hypre(*grid, rank, nprocs, bVerbose, bCG, bRelaxation, &kappa, bMobility, &mobility);
    
    while (t <= tend)
    {
        if(!bProfiler){
            dt = _estimate_dt(Diff_dt, CFL);
            _computePressureSource(nParallelGranularity,rho);
            helmholtz_solver.solve();
            _computeVelocities(boundaryInfo, bMobility, &mobility);
            _reactionDiffusionStep(boundaryInfo, Dw, Dg, rho);
            _advectionConvectionStep(boundaryInfo);
            _timeUpdate(nParallelGranularity, dt);
            
            t                   += dt   ;
            numberOfIterations  ++      ;
            
            if ( t >= ((double)(whenToWrite)) ){
                if(rank==0) _dump(iCounter++);
                if(rank==0) printf(" dumping at time t=%f, as Data_%04d \n", t, iCounter-1);
                whenToWrite = whenToWrite + whenToWriteOffset;
            }
        }
        else{
            profiler.getAgent("TimeStep").start();
            dt = _estimate_dt(Diff_dt, CFL);
            profiler.getAgent("TimeStep").stop();
            
            profiler.getAgent("PressureSource").start();
            _computePressureSource(nParallelGranularity,rho);
            profiler.getAgent("PressureSource").stop();
            
            profiler.getAgent("Helmholtz").start();
            helmholtz_solver.solve();
            profiler.getAgent("Helmholtz").stop();
            
            profiler.getAgent("Velocity").start();
            _computeVelocities(boundaryInfo, bMobility, &mobility);
            profiler.getAgent("Velocity").stop();
            
            profiler.getAgent("RD").start();
            _reactionDiffusionStep(boundaryInfo, Dw, Dg, rho);
            profiler.getAgent("RD").stop();
            
            profiler.getAgent("Advect-Conv").start();
            _advectionConvectionStep(boundaryInfo);
            profiler.getAgent("Advect-Conv").stop();
            
            profiler.getAgent("TimeUpdate").start();
            _timeUpdate(nParallelGranularity, dt);
            profiler.getAgent("TimeUpdate").stop();
            
            t                   += dt   ;
            numberOfIterations  ++      ;
            
            if ( t >= ((double)(whenToWrite))) {
                profiler.getAgent("I/O").start();
                if(rank==0) _dump(iCounter++);
                profiler.getAgent("I/O").stop();
                whenToWrite = whenToWrite + whenToWriteOffset;
            }
            
        }
    }
    
   //if(rank == 0)      _dump(iCounter);
    
    helmholtz_solver.clean();
    
    if(rank==0){
    if(bVerbose) profiler.printSummary();
    if(bVerbose) printf("**** Dumping done\n");
    if(bVerbose) printf("\n\n Run Finished \n\n");
    isDone = 1;
    }
}

