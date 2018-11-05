/*
 *  fitfun.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <math.h>

#define _USE_ROSENBROCK_

double fitfun(double /*const*/ *x, int N, void *output, int *info)
{
	double f;
	int i;

	f = 0.0;
	for (i=0; i<N-1; i++)	/* rosenbrock */
		f = f + 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);
//	f = -f;
	f = -log(f);	// peh xxx: logval = 1

	return f;
}

