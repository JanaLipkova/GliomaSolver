//
//  Likelihood_main.cpp
//  GliomaXcode
//
//  Created by Lipkova on 15/06/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#include "HGG_Likelihood.h"

double my_gettime()
{
   struct timeval t;
   gettimeofday(&t, NULL);
   return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
}

int main(int argc,const char ** argv)
{
    printf("==========================================\n");
    printf("          Computing likelihood            \n");
    printf("==========================================\n");
    
    ArgumentParser parser(argc, argv);
       double t_begin = my_gettime();

    HGG_Likelihood* l = new HGG_Likelihood(argc, (const char **)argv);
    l -> run();
    delete l;
    
      double t_end = my_gettime();
   printf("Likelihood: elapsed time=%f (sec) \n",t_end - t_begin);

    
    return 0;
};
