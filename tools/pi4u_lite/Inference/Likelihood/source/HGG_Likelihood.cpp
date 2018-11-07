//
//  HGG_Likelihood.cpp
//  GliomaXcode
//
//  Created by Lipkova on 15/06/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#include "HGG_Likelihood.h"

HGG_Likelihood::HGG_Likelihood(const int argc, const char ** argv):parser(argc, argv)
{
    ifstream mydata("LikelihoodInput.txt");
    
    if (mydata.is_open())
    {
        mydata >> PETsigma2;
        mydata >> PETscale;
        mydata >> T1uc;
        mydata >> T2uc;
        mydata >> slope;
        mydata.close();
    }
    else
    {
        printf("Aborting: missing input file LikelihoodInput.txt \n");
        abort();
    }
    
    bSelectivePoints = parser("-bROI").asBool(0);
    stepPET = parser("-stepPET").asInt(1);
    stepTi  = parser("-stepMRI").asInt(1);
    
    printf("PET: PETsigma2=%f, PETscale=%f \n",PETsigma2,PETscale );
    printf("MRI: T1uc=%f, T2uc =%f, slope=%f \n", T1uc, T2uc, slope);
    printf("stepPET=%d, stepTi =%d \n", stepPET, stepTi);

}

long double HGG_Likelihood::_computePETLogLikelihood(MatrixD3D model)
{
    /* 1) Read in PET data */
    char filename[256];
    sprintf(filename,"tumPET.dat");
    MatrixD3D PETdata(filename);   // pet signal, u_real = PETscale*pet_signal
    
    int dataX = PETdata.getSizeX();
    int dataY = PETdata.getSizeY();
    int dataZ = PETdata.getSizeZ();
    
    int N = 0;
    long double sum = 0.;
    
    
    for (int iz = 0; iz < dataZ; iz++)
        for (int iy = 0; iy < dataY; iy++)
            for (int ix = 0; ix < dataX; ix++)
            {
                if ( (ix % stepPET == 0) &&  (iy % stepPET == 0) && (iz % stepPET == 0) )
                {
                    if(PETdata(ix,iy,iz) > 0.)
                    {
                        sum += ( model(ix,iy,iz) - PETscale * PETdata(ix,iy,iz) )*( model(ix,iy,iz) - PETscale * PETdata(ix,iy,iz) );
                        N++;
                    }
                }
            }
    
    
    long double p1 = -0.5 * N * log(2.* M_PI * PETsigma2);
    long double p2 = -0.5 * (1./PETsigma2) * sum;
    
    printf("Points=%i, stepPET=%i \n", N, stepPET);
    return (p1 + p2);
}

long double HGG_Likelihood::_computeTiLogLikelihood(MatrixD3D model, int Ti)
{
    /* 1) Read in T1 data */
    char filename[256];
    
    if(Ti ==1 )
        sprintf(filename,"tumT1c.dat");
    else
        sprintf(filename,"tumFLAIR.dat");
    
    MatrixD3D data(filename);
    int dataX = data.getSizeX();
    int dataY = data.getSizeY();
    int dataZ = data.getSizeZ();
    
    long int N = dataX*dataY*dataZ;
    assert(N == model.getSizeX() * model.getSizeY() * model.getSizeZ() );
    
    long double sum = 0.;
    
    if(!bSelectivePoints)
    {
        //#pragma omp parallel for reduction(+:sum)
        for (int iz = 0; iz < dataZ; iz++)
            for (int iy = 0; iy < dataY; iy++)
                for (int ix = 0; ix < dataX; ix++)
                    sum += _computeLogBernoulli(model(ix,iy,iz), data(ix,iy,iz), Ti);
    }
    else
    {
        sprintf(filename,"ROI.dat");
        MatrixD2D Points(filename);
        int Npoints = Points.getSizeX();
        
        for(int i=0; i<Npoints; i++)
        {
            int ix = Points(i,0);
            int iy = Points(i,1);
            int iz = Points(i,2);            

            if ( (ix % stepTi == 0) &&  (iy % stepTi == 0) && (iz % stepTi == 0) )
                sum += _computeLogBernoulli(model(ix,iy,iz), data(ix,iy,iz), Ti);
            
        }
        
        printf("Ti Npoints = %i \n",Npoints);

    }
    
    printf("LogLike of T%i = %Lf \n",Ti, sum);
    return sum;
}


/* Likelihood based on Bernoulli distr. of double logsitic sigmoid
 1) compute alpha
 - should be in (0,1)
 - rounding errros can make alpha = 0 or =1
 - if that is the case, correct it since it will be used in log() */
long double HGG_Likelihood::_computeLogBernoulli(double u, double y, int Ti)
{
    double uc, is2;
    
    if(Ti == 1){
        uc = T1uc;
        is2 = 1./slope;
    }
    else{
        uc = T2uc;
        is2 = 1./slope;
    }
    
    double diff = u - uc;
    
    long double omega2 = (diff > 0.) ? 1. : diff*diff;
    long double alpha = 0.5 + 0.5 * sgn(diff) * (1. - exp( -omega2 * is2));
    
    // OLD version: long double alpha = 0.5 + 0.5 * sgn(diff) * (1. - exp( -diff * diff * is2));
    return  (y == 1 ) ? log(alpha) : log(1.-alpha);
}

int HGG_Likelihood::sgn(double d)
{
    double eps = 0.0;
    if (d < -eps) { return -1; }
    else { return d > eps; }
}

void HGG_Likelihood::_writeToFile(long double output)
{
    long double MinusLogLikelihood = - output;
    
    FILE* myfile = fopen("Likelihood.txt", "w");
    if (myfile!=NULL)
        fprintf(myfile, "%Lf \n", MinusLogLikelihood);
    
    fclose(myfile);
}

void HGG_Likelihood::run()
{
    char filename[256];
    sprintf(filename,"HGG_data.dat");
    MatrixD3D model(filename);
    
    long double Lpet = _computePETLogLikelihood(model);
    long double Lt1  = _computeTiLogLikelihood(model, 1);
    long double Lt2  = _computeTiLogLikelihood(model, 2);
    
    long double costFunction = Lpet + Lt1 + Lt2;
   
    printf("L_Pet=%Lf, L_T1=%Lf, L_T2 = %Lf \n", Lpet, Lt1, Lt2);
    printf("LogLike = %Lf \n", costFunction);
    _writeToFile(costFunction);

}
