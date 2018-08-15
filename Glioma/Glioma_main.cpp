/*
 *  mainGlioma.cpp
 *  Untitled
 *
 *  Created by Lipkova on 9/14/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */


#include <iostream>
#include <xmmintrin.h>

#include "../Tests/Test.h"
#include "../Tests/HelmholtzTest.h"


#include "Glioma_ReactionDiffusion.h"
#include "Glioma_BrainDeformationTimeRelaxation.h"

#ifdef HYPRE
#include "Glioma_BrainDeformation.h"
#endif

using namespace std;
using namespace MRAG;

int main(int argc,const char ** argv)
{
    std::cout << std::endl << "MRAG Launched" << std::endl << std::endl;
    
    ArgumentParser parser(argc, argv);
    Environment::setup(max(1, parser("-nthreads").asInt()));
    
    Glioma * s = NULL;

    
    printf("INPUT IS %s\n", parser("-model").asString().data());
    
    if(parser("-model").asString() == "RD")
        s = new Glioma_ReactionDiffusion(argc, (const char **)argv);
    else if(parser("-model").asString() == "deformTime")
        s = new Glioma_BrainDeformationTimeRelaxation(argc, (const char **)argv);
#ifdef HYPRE
    else if(parser("-model").asString() == "deform")
        s = new Glioma_BrainDeformation(argc, (const char **)argv);
    else if(parser("-model").asString() == "helmholtzTest")
        s = new HelmholtzTest(argc, (const char **)argv);
#endif
    else
        s = new Test(argc, (const char **)argv);
  
    tbb::tick_count t1,t0;
    {
        t0=tbb::tick_count::now();
        s->run();
        t1=tbb::tick_count::now();
    }
    
    printf("we spent: %2.2f \n",(t1-t0).seconds());	
    delete s;
    
    std::cout << std::endl << "MRAG Terminated" << std::endl;
}



