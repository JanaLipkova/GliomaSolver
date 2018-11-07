//
//  HGG_Likelihood.h
//  GliomaXcode
//
//  Created by Lipkova on 15/06/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//


/* CODE DESCRIPTION:
 This function evaluate likelihood for Glioma paramter estimation
 reads in data and results of the simulation, and write to the file the final likelihood so it can be then read in CMA/TMCMC
 */

#ifndef __Glioma_Likelihood__
#define __Glioma_Likelihood__


#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <map>
#include <vector>
#include <math.h>
#include <ctime>
#include <limits>
#include <sys/time.h>
using namespace std;

#include "ArgumentParser.h"
#include "Matrix.h"

typedef Matrix::D3D<double> MatrixD3D;
typedef Matrix::D2D<double> MatrixD2D;


class HGG_Likelihood
{
private:
    ArgumentParser		 parser;
    long double          _computePETLogLikelihood(MatrixD3D model);
    long double          _computeTiLogLikelihood(MatrixD3D model, int Ti);
    long double          _computeLogBernoulli(double y, double u, int Ti);
    void                 _writeToFile(long double output);
    int                  sgn(double d);

    double               PETsigma2, PETscale;
    double               slope;  //same for T1 & T2(sigma2 double sigmoid, k singel sigmoid)
    double               T1uc,    T2uc;
    bool                 bSelectivePoints;
    int                  stepPET, stepTi;
    
    
public:
    HGG_Likelihood(const int argc, const char ** argv);
    ~HGG_Likelihood(){};
    void run();

};

#endif /* defined(__Glioma_Likelihood__) */
