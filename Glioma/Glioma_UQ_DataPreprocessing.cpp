//
//  Glioma_UQ_DataPreprocessing.cpp
//  GliomaSolverXcode
//
//  Created by Lipkova on 01/11/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

#include "Glioma_UQ_DataPreprocessing.h"


static int maxStencil[2][3] = {
    -1, -1, -1,
    +2, +2, +2
};

Glioma_UQ_DataPreprocessing::Glioma_UQ_DataPreprocessing(int argc, const char ** argv): parser(argc, argv)
{
    bVerbose    = parser("-verbose").asBool(1);
    bVTK        = parser("-vtk").asBool();
    
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("//////////////////          Glioma UQ Data Preprocessing        ////////////////\n");
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("Set up: blockSize=%d Wavelets=w%s (blocksPerDimension=%d, maxLevel=%d)\n", blockSize, "w", blocksPerDimension, maxLevel);
    
    refiner		= new Refiner_SpaceExtension(resJump,maxLevel);
    compressor	= new Compressor(resJump);
    grid        = new Grid<W,B>(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil);
    grid->setCompressor(compressor);
    grid->setRefiner(refiner);
    stSorter.connect(*grid);
    
    L = 1;
    PatientFileName = parser("-PatFileName").asString();
    int bSynt = parser("-bSynt").asInt(0);
    
    switch (bSynt) {
            
        case 0:
            _ic(*grid, PatientFileName, L);
            break;

        case 1:
            _ic_Synthetic(*grid, PatientFileName, L);
            break;
            
        default:
            break;
    }
    
    isDone              = false;
}

Glioma_UQ_DataPreprocessing::~Glioma_UQ_DataPreprocessing()
{
    if(bVerbose) std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialCondition
// Patient 01 data
// 1) read in anatomies - rescaled to [0,1]^3
// 2) read in tumor center of mass + initialize tumor around
// 3) set the characteristic length L as the length of the data
void Glioma_UQ_DataPreprocessing::_ic(Grid<W,B>& grid, string PatientFileName, Real &L)
{
    printf("Reading data from file: %s \n", PatientFileName.c_str());
    
    char anatomy      [200];
    sprintf(anatomy, "%sGM.dat", PatientFileName.c_str());
    MatrixD3D GM(anatomy);
    sprintf(anatomy, "%sWM.dat", PatientFileName.c_str());
    MatrixD3D WM(anatomy);
    sprintf(anatomy, "%sCSF.dat", PatientFileName.c_str());
    MatrixD3D CSF(anatomy);
    
    sprintf(anatomy, "%sTum_T1c.dat", PatientFileName.c_str());
    MatrixD3D TumT1c(anatomy);
    sprintf(anatomy, "%sTum_FLAIR.dat", PatientFileName.c_str());
    MatrixD3D TumFLAIR(anatomy);
    sprintf(anatomy, "%sTum_FET.dat", PatientFileName.c_str());
    MatrixD3D TumFET(anatomy);
    
    sprintf(anatomy, "%sFLAIR.dat", PatientFileName.c_str());
    MatrixD3D FLAIR(anatomy);
    sprintf(anatomy, "%sT1c.dat", PatientFileName.c_str());
    MatrixD3D T1c(anatomy);
    sprintf(anatomy, "%sFET.dat", PatientFileName.c_str());
    MatrixD3D FET(anatomy);
    
    int brainSizeX = (int) GM.getSizeX();
    int brainSizeY = (int) GM.getSizeY();
    int brainSizeZ = (int) GM.getSizeZ();
    printf("dataSizeX=%i, dataSizeY=%i, dataSizeZ= %i \n", brainSizeX, brainSizeY, brainSizeZ);
    
    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    L   = brainSizeMax * 0.1;   // voxel spacing 1mm, convert from mm to cm  // L = 25.6 cm
    
    double brainHx = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
    
    Real pGM, pWM, pCSF, pTumT1, pTumFLAIR;
    
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
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    int mappedBrainX = (int)floor( x[0] / brainHx  );
                    int mappedBrainY = (int)floor( x[1] / brainHy  );
                    int mappedBrainZ = (int)floor( x[2] / brainHz  );
                    
                    if ( (mappedBrainX >= 0 && mappedBrainX < brainSizeX) & (mappedBrainY >= 0 && mappedBrainY < brainSizeY) && (mappedBrainZ >= 0 && mappedBrainZ < brainSizeZ) )
                    {
                        // Anatomy
                        pGM        =  GM( mappedBrainX,mappedBrainY,mappedBrainZ);
                        pWM        =  WM( mappedBrainX,mappedBrainY,mappedBrainZ);
                        pCSF       =  CSF(mappedBrainX,mappedBrainY,mappedBrainZ);
                        
                        // separate tissue from fluid
                        pWM  = (pCSF > 0.2) ? 0. : pWM;
                        pGM  = (pCSF > 0.2) ? 0. : pGM;
                        pCSF = (pCSF > 0.2) ? 1. : pCSF;
                        double all = pWM + pGM + pCSF;
                        
                        // normalise
                        block(ix,iy,iz).p_w   = (all > 0.) ? (pWM  / all) : 0.;
                        block(ix,iy,iz).p_g   = (all > 0.) ? (pGM  / all) : 0.;
                        block(ix,iy,iz).p_csf = (all > 0.) ? (pCSF / all) : 0.;
                        
                        // Tumor Observations + normalise
                        pTumT1     =  TumT1c(mappedBrainX,mappedBrainY,mappedBrainZ);
                        pTumFLAIR  =  TumFLAIR(mappedBrainX,mappedBrainY,mappedBrainZ);
                        
                        block(ix,iy,iz).wm  = (pTumT1    > 0.1) ? 1. : 0.;
                        block(ix,iy,iz).gm  = (pTumFLAIR > 0.1) ? 1. : 0.;
                        block(ix,iy,iz).phi = TumFET(mappedBrainX,mappedBrainY,mappedBrainZ);
                        
                        // MRI + PET scnas
                        block(ix,iy,iz).csf = T1c(mappedBrainX,mappedBrainY,mappedBrainZ);
                        block(ix,iy,iz).ux  = FLAIR(mappedBrainX,mappedBrainY,mappedBrainZ);
                        block(ix,iy,iz).uy  = FET(mappedBrainX,mappedBrainY,mappedBrainZ);
                        
                    }
                    
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}



#pragma mark InitialCondition
// Patient 01 data
// 1) read in anatomies - rescaled to [0,1]^3
// 2) read in tumor center of mass + initialize tumor around
// 3) set the characteristic length L as the length of the data
void Glioma_UQ_DataPreprocessing::_ic_Synthetic(Grid<W,B>& grid, string PatientFileName, Real &L)
{
    printf("Reading data from file: %s \n", PatientFileName.c_str());
    
    char anatomy      [200];
    sprintf(anatomy, "%sHGG_data_GT.dat", PatientFileName.c_str());
    MatrixD3D GT(anatomy);
    
    int brainSizeX = (int) GT.getSizeX();
    int brainSizeY = (int) GT.getSizeY();
    int brainSizeZ = (int) GT.getSizeZ();
    printf("dataSizeX=%i, dataSizeY=%i, dataSizeZ= %i \n", brainSizeX, brainSizeY, brainSizeZ);
    
    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    L   = brainSizeMax * 0.1;   // voxel spacing 1mm, convert from mm to cm  // L = 25.6 cm
    
    double brainHx = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
    Real data;
    Real ucFLAIR = 0.25;
    Real ucT1    = 0.7;
    
    double xi_1, xi_2;
    double alpha = 0.02;  // 5% noise
    int n = 0;
    
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
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    int mappedBrainX = (int)floor( x[0] / brainHx  );
                    int mappedBrainY = (int)floor( x[1] / brainHy  );
                    int mappedBrainZ = (int)floor( x[2] / brainHz  );
                    
//                    // aspect ratio correction
//                    mappedBrainX -= (int) ( (brainSizeMax - brainSizeX) * 0.5);
//                    mappedBrainY -= (int) ( (brainSizeMax - brainSizeY) * 0.5);
//                    mappedBrainZ -= (int) ( (brainSizeMax - brainSizeZ) * 0.5);
                    
                    if ( (mappedBrainX >= 0 && mappedBrainX < brainSizeX) & (mappedBrainY >= 0 && mappedBrainY < brainSizeY) && (mappedBrainZ >= 0 && mappedBrainZ < brainSizeZ) )
                    {
                        data =  GT( mappedBrainX,mappedBrainY,mappedBrainZ);
                        
                        // transfer GT to synthetic tumor observations
                        block(ix,iy,iz).wm  = (data > ucT1)    ? 1. : 0.;
                        block(ix,iy,iz).gm  = (data > ucFLAIR) ? 1. : 0.;
                        block(ix,iy,iz).phi = ( block(ix,iy,iz).wm + block(ix,iy,iz).gm >= 0.9 ) ? data : 0.;
                        
                        // add noise to pet signal
                        // Box-Muller
                        if(  (n % 2) == 0 ){
                            double u1 = drand48();
                            double u2 = drand48();
                            
                            xi_1 = sqrt( - 2.0*log( u1 ) ) * cos( 2.0*M_PI*u2 );
                            xi_2 = sqrt( - 2.0*log( u1 ) ) * sin( 2.0*M_PI*u2 );
                        }
                        
                        if( block(ix,iy,iz).phi > 0.  ){
                            block(ix,iy,iz).phi = max( 0., block(ix,iy,iz).phi + alpha * xi_1);
                            xi_1 = xi_2;
                            n++;
                        }
                    
                    }
                    
                }
        
        
        grid.getBlockCollection().release(info.blockID);
    }
}



#pragma mark Preprocessing
void Glioma_UQ_DataPreprocessing::_normalisePET()
{
    Real maxPET = 0.;
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++){
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                    maxPET = max(maxPET, block(ix,iy,iz).phi);
    }
    
    
    Real scale = 1./ maxPET;
    
#pragma omp parallel for schedule(static)
    for(int i=0; i<vInfo.size(); i++){
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                    block(ix,iy,iz).phi = block(ix,iy,iz).phi * scale;
        
    }
}

/* Case specific prior:
   1) Load generic prior range,
   2) Compute patient specific parts
       - IC prior = Box centered at CM of T1c tumor, side = diameter of sphere with same volume as T1c lesion
       - T        = [Tmin, Tmax + 300], where Tmin = TumRad/v_max, Tmax = TumRad + 2cm / v_min, where
                     v = 2*sqrt(D*rho), get v_max, v_min from generic prior 
                     TumRad is a radius of sphere with volume same as FLAIR-enhancing lesion
   3) Write prior and LogPrior  */
void Glioma_UQ_DataPreprocessing::_computePriorRange()
{
    // 1) Load generic prior range:
    // Input order: D, rho, Tend, icx, icy, icz, PETsigma, ucT1, ucT2, Tisigma2
    vector<float> prior;
    float tmp;
    
    ifstream myfile("GenericPrior.txt");
    if (myfile.is_open()){
         while( myfile.good()){
             myfile >> tmp;
             if(myfile.eof()) break;  // since EOF is after the last element, without this the last element would be written to prior 2x
             prior.push_back(tmp);
          }
        myfile.close();
    }
    else{
        printf("Aborting: missing input file GenericPrior.txt \n");
        abort();
    }
    
    
    // 2) Compute IC and T prior
    Real massT1 = 0.;
    Real massT2 = 0.;
    Real cmx = 0.;
    Real cmy = 0.;
    Real cmz = 0.;
    Real h;
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        h = info.h[0];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    massT1 += block(ix,iy,iz).wm;
                    massT2 += block(ix,iy,iz).gm;
                    
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    cmx += x[0] * block(ix,iy,iz).wm;
                    cmy += x[1] * block(ix,iy,iz).wm;
                    cmz += x[2] * block(ix,iy,iz).wm;
                }
    }
    
    cm[0] = cmx / massT1;
    cm[1] = cmy / massT1;
    cm[2] = cmz / massT1;
    
    assert(!isnan(cm[0]) || !isnan(cm[1]) || !isnan(cm[2]));

    Real h3 = h*h*h;
    Real volT1 = massT1 * h3;
    volFLAIR   = massT2 * h3;
    
    // Radius of sphere with volume same as T1c lesion
    Real T1radius3 = 3. * volT1 / (4. * M_PI);
    Real T1radius  = pow(T1radius3, 1./3.);
    T1radius = T1radius * 1.25;  // add extra 25% just in case
    
    Real ICprior[6] = {
        cm[0] - T1radius, cm[0] + T1radius,
        cm[1] - T1radius, cm[1] + T1radius,
        cm[2] - T1radius, cm[2] + T1radius };
    
    Real FLAIRradius3 = 3.* volFLAIR / (4.* M_PI);
    Real FLAIRradius  = pow(FLAIRradius3, 1./3.);
    FLAIRradius       = 1.25 * FLAIRradius * L; // add 25% to be sure[cm]
   
    Real vmin = 2*pow(prior[0]*prior[2], 1./2.);
    Real vmax = 2*pow(prior[1]*prior[3], 1./2.);
    
    Real Tmin = FLAIRradius/vmax;
    Real Tmax = (FLAIRradius + 2.) / vmin + 300.; // add 2cm infiltration + 300 days
    
    if(bVerbose) printf("Tmin =%f, Tmax=%f \n",  Tmin, Tmax);

    // Restrict prior range to min 30 days and max 5 years
    Tmin = max( (Real) 30., Tmin);
    Tmax = min( (Real) 1800, Tmax);
    if(bVerbose) printf("Restricted range: Tmin =%f, Tmax=%f \n",  Tmin, Tmax);

    Real cmin = Tmin * Tmin * prior[3];
    Real cmax = Tmax * Tmax * prior[2];
    
    // Time
    prior[4]  = cmin;
    prior[5]  = cmax;
    prior[6]  = ICprior[0];
    prior[7]  = ICprior[1];
    prior[8]  = ICprior[2];
    prior[9]  = ICprior[3];
    prior[10] = ICprior[4];
    prior[11] = ICprior[5];

    if(bVerbose) printf("FLAIRradius = %f, vmin =%f, vmax=%f \n", FLAIRradius, vmin, vmax);
    
    FILE * pFile;
    pFile = fopen ("PriorRange.txt",    "w");
    for(int i = 0; i < prior.size(); i=i+2)
        fprintf(pFile, "%.4f     %.4f\n", prior[i], prior[i+1]);
    fclose (pFile);

    pFile = fopen ("LogPriorRange.txt",    "w");
    for(int i = 0; i < prior.size(); i=i+2)
        fprintf(pFile, "%.4f     %.4f\n", log(prior[i]), log(prior[i+1]));
    fclose (pFile);
    
}


#pragma mark OUTPUT
void Glioma_UQ_DataPreprocessing:: _dump(int counter)
{
    if(bVTK)
    {
        char filename[256];
        sprintf(filename,"Data_%04d", counter);
        
        IO_VTKNative3D<W,B, 10,0 > vtkdumper2;
        vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
    }
}


void Glioma_UQ_DataPreprocessing::_dumpUQdata()
{
    int gpd = blocksPerDimension * blockSize;
    MatrixD3D PET(gpd,gpd,gpd);
    MatrixD3D T1(gpd,gpd,gpd);
    MatrixD3D T2(gpd,gpd,gpd);
    
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
#pragma omp parallel for schedule(static)
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    PET(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1], iz + B::sizeZ * info.index[2] ) = block(ix,iy,iz).phi;
                    T1(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1], iz + B::sizeZ * info.index[2] ) = block(ix,iy,iz).wm;
                    T2(ix + B::sizeX * info.index[0], iy + B::sizeY * info.index[1], iz + B::sizeZ * info.index[2] ) = block(ix,iy,iz).gm;
                }
        
    }
    
    char filename[256];
    sprintf(filename,"tumPET.dat");
    PET.dump(filename);
    
    sprintf(filename,"tumT1c.dat");
    T1.dump(filename);
    
    sprintf(filename,"tumFLAIR.dat");
    T2.dump(filename);
}


/* ROI = regions used for likelihood computation
 (exclude points outside brain and far from tumor since they do not inform likelihood)
 ROI = sphere cenered at tumor CM,
 with radius given as radius of sphere with same volume as FLAIR lesion,
 -> the radius is increased by 100% to make sure whole tumor + infiltration around is included */
void Glioma_UQ_DataPreprocessing::_dumpInferenceROI()
{
    // radius of sphere with same volume as FLAIR lesion
    Real radius3 = 3. * volFLAIR / (4. * M_PI);
    Real radius  = pow(radius3, 1./3.);
    radius       = radius * 2.;
    
    if(bVerbose) printf("Center of mass =[ %f, %f, %f]\n", cm[0],cm[1],cm[2]);
    if(bVerbose) printf("ROI radius %f\n", radius);
    
    int points         = 0;
    vector<double>     ROI;
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++ )
            for(int iy=0; iy<B::sizeY; iy++ )
                for(int ix=0; ix<B::sizeX; ix++ )
                {
                    Real x[3];
                    info.pos(x, ix, iy, iz);
                    
                    const Real p[3] = {x[0] - cm[0], x[1] - cm[1], x[2] - cm[2]};
                    const Real dist = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
                    
                    Real tissue = block(ix,iy,iz).p_w + block(ix,iy,iz).p_g;
                    
                    if(( tissue > 0 )&&(dist <= radius))
                    {
                        const int gix = ix + info.index[0] * B::sizeX;
                        const int giy = iy + info.index[1] * B::sizeY;
                        const int giz = iz + info.index[2] * B::sizeZ;
                        
                        ROI.push_back(gix);
                        ROI.push_back(giy);
                        ROI.push_back(giz);
                        points++;
                    }
                }
    }
    
    MatrixD2D out(points,3);
    for (int i = 0; i<points; i++)
    {
        out(i,0) = ROI[i*3    ];
        out(i,1) = ROI[i*3 + 1];
        out(i,2) = ROI[i*3 + 2];
    }
    
    char filename[256];
    sprintf(filename,"ROI.dat");
    out.dump(filename);
}



void Glioma_UQ_DataPreprocessing::run()
{
    _normalisePET();
    _computePriorRange();
    _dumpUQdata();
    _dumpInferenceROI();
    _dump(0);
    
    if(bVerbose)printf("\n\n Run Finished \n\n");
}
