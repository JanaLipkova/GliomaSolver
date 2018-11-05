/* --------------------------------------------------------- */
/* --------------- A Very Short Example -------------------- */
/* --------------------------------------------------------- */

#define _XOPEN_SOURCE 500
#define _BSD_SOURCE

#include <stdio.h>
#include <stdlib.h> /* free() */
#include "cmaes_interface.h"
#include <mpi.h>
#include <torc.h>
#include <unistd.h>
#include <math.h>

#define _STEALING_
//#define _RESTART

#include "fitfun.c" 

/* the objective (fitness) function to be minimized */
void taskfun(double *x, int *pn, double *res, int *info)
{
	int n = *pn;
//	int gen, chain, step, task;
//	gen = info[0]; chain = info[1]; step = info[2]; task = info[3];
//	printf("executing task (%d,%d,%d,%d)\n", gen, chain, step, task);
	
	double f = -fitfun(x, n, (void *)NULL, info);	/* CMA-ES needs this minus sign */

	*res = f;
	return;
}

double lower_bound[] = {-6.0, -6.0};
double upper_bound[] = {+6.0, +6.0};

int is_feasible(double *pop, int dim)
{
	int i, good;
	// printf("is_feasible %d\n", dim);
	for (i = 0; i < dim; i++) {
		good = (lower_bound[i] <= pop[i]) && (pop[i] <= upper_bound[i]);
		if (!good) {
			//printf("%d not good\n", i);
			//usleep(1000);
			return 0;
		}
	}
	return 1;
}


/* the optimization loop */
int main(int argn, char **args)
{
	cmaes_t evo; /* an CMA-ES type struct or "object" */
	double *arFunvals, *const*pop, *xfinal;
	int i; 
	/* peh - start */
	int lambda, dim;
	double gt0, gt1, gt2, gt3;
	double tt0, tt1, stt = 0.0; 
	int step = 0;
	int info[4];	/* gen, chain, step, task */
	/* peh - end */

	torc_register_task(taskfun);

	/* Initialize everything into the struct evo, 0 means default */
	torc_init(argn, args, MODE_MS);

	gt0 = torc_gettime();
	arFunvals = cmaes_init(&evo, 0, NULL, NULL, 0, 0, "initials.par"); 
	printf("%s\n", cmaes_SayHello(&evo));
	cmaes_ReadSignals(&evo, "signals.par");  /* write header and initial values */

	/* Iterate until stop criterion holds */
	gt1 = torc_gettime();
	while(!cmaes_TestForTermination(&evo))
	{ 
		/* generate lambda new search points, sample population */
		pop = cmaes_SamplePopulation(&evo); /* do not change content of pop */

		/* Here you may resample each solution point pop[i] until it
		   becomes feasible, e.g. for box constraints (variable
		   boundaries). function is_feasible(...) needs to be
		   user-defined.  
		   Assumptions: the feasible domain is convex, the optimum is
		   not on (or very close to) the domain boundary, initialX is
		   feasible and initialStandardDeviations are sufficiently small
		   to prevent quasi-infinite looping.
		*/

		lambda = cmaes_Get(&evo, "lambda"); 
		dim = cmaes_Get(&evo, "dim"); 

		for (i = 0; i < cmaes_Get(&evo, "popsize"); ++i) 
			while (!is_feasible(pop[i], dim)) 
				cmaes_ReSampleSingle(&evo, i); 

		/* evaluate the new search points using fitfun from above */ 
		tt0 = torc_gettime();
		for (i = 0; i < lambda; ++i) {
#if 0
			/*arFunvals[i] =*/ fitfun(pop[i], (int) dim, &arFunvals[i]);
#endif
			info[0] = 0; info[1] = 0; info[2] = step; info[3] = i; 	/* gen, chain, step, task */
			torc_create(-1, taskfun, 4,
				dim, MPI_DOUBLE, CALL_BY_VAL,
				1, MPI_INT, CALL_BY_COP,
				1, MPI_DOUBLE, CALL_BY_RES,
				4, MPI_INT, CALL_BY_COP,
				pop[i], &dim, &arFunvals[i], info);
		}
#if defined(_STEALING_)
		torc_enable_stealing();
#endif
		torc_waitall();
#if defined(_STEALING_)
		torc_disable_stealing();
#endif
		tt1 = torc_gettime();
		stt += (tt1-tt0);

		/* update the search distribution used for cmaes_SampleDistribution() */
		cmaes_UpdateDistribution(&evo, arFunvals);  

		/* read instructions for printing output or changing termination conditions */ 
		cmaes_ReadSignals(&evo, "signals.par");   
		fflush(stdout); /* useful in MinGW */

#if VERBOSE
		{
		const double *xbever = cmaes_GetPtr(&evo, "xbestever");
		double fbever = cmaes_Get(&evo, "fbestever");
	
		printf("BEST @ %5d: ", step);
		for (i = 0; i < dim; i++)
			printf("%25.16lf ", xbever[i]);
		printf("%25.16lf\n", fbever);
		}
		printf("Step %4d: time = %.3lf seconds\n", step, tt1-tt0);
#endif

#if defined(_RESTART_)
		cmaes_WriteToFile(&evo, "resume", "allresumes.dat");         /* write final results */
#endif
		step++;
	}

	gt2 = torc_gettime();
	printf("Stop:\n%s\n",  cmaes_TestForTermination(&evo)); /* print termination reason */
	cmaes_WriteToFile(&evo, "all", "allcmaes.dat");         /* write final results */
//	cmaes_WriteToFile(&evo, "resume", "allresumes.dat");         /* write final results */

	/* get best estimator for the optimum, xmean */
	xfinal = cmaes_GetNew(&evo, "xmean"); /* "xbestever" might be used as well */
	cmaes_exit(&evo); /* release memory */ 

	/* do something with final solution and finally release memory */
	free(xfinal); 

	gt3 = torc_gettime();
	printf("Total elapsed time = %.3lf  seconds\n", gt3-gt0);
	printf("Initialization time = %.3lf  seconds\n", gt1-gt0);
	printf("Processing time = %.3lf  seconds\n", gt2-gt1);
	printf("Funtion Evaluation time = %.3lf  seconds\n", stt);
	printf("Finalization time = %.3lf  seconds\n", gt3-gt2);

	torc_finalize();

	return 0;
}

