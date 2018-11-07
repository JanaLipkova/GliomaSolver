/*
 *  engine_ss2.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#define _XOPEN_SOURCE 500
#define _BSD_SOURCE

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <unistd.h>
#include <iostream>
using namespace std;
#include <mpi.h>
#include <omp.h>
#include <math.h>
extern "C"
{
#include <torc.h>
#include "engine_ss.h"
}

/**********************************************/
// PARAMS
/**********************************************/
data_t data;

void read_data()
{
	int i;

	data.Nth = 2;

	/* DEFAULT VALUES */
	data.N_init = 11000;
	data.N_seeds = 1000;
	data.N_steps = 10;

	data.lb = -6.0;	// Default LB, same for all
	data.ub = +6.0;	// Default UB, same for all
 
	data.lowerbound = (double *)malloc(data.Nth*sizeof(double));
	data.upperbound = (double *)malloc(data.Nth*sizeof(double));

	for (i = 0; i < data.Nth; i++) {
		data.lowerbound[i] = data.lb;
		data.upperbound[i] = data.ub;
	}

	data.MAXTHRESHOLDS = 10;
	data.threshold = (double *)malloc(data.MAXTHRESHOLDS*sizeof(double));
	for (i = 0; i < data.MAXTHRESHOLDS; i++) {
		data.threshold[i] = 0.0;
	}

	data.logval = 0;
	data.sigma = 0.1;
	data.seed = 280675;
	data.iplot = 0;	// gnuplot

	data.FACTOR = 0.50;

	/* USER-DEFINED VALUES */
	FILE *f = fopen("subset.par", "r");
	if (f == NULL) {
		return;
	}

	/*
	Nth		2
	N_init		11000
	N_seeds		1000
	N_steps		10
	Bdef		-4	4
	#B0		-6	6
	#B1		-6	6
	NTHRESHOLDS	10
	Thresholds	-64.0 -32.0 -16.0 -8.0 -4.0 -2.0 -1.0 -0.5 -0.25 -0.125
	logval		0
	sigma		0.1
	seed		280675
	iplot		0
	factor		0.50
	MAXTHRESHOLDS	20
	*/

	char line[256];

	int line_no = 0;
	while (fgets(line, 256, f)!= NULL) {
		line_no++;
		if ((line[0] == '#')||(strlen(line)==0)) {
			printf("ignoring line %d\n", line_no);
			continue;
		}

		if (strstr(line, "Nth")) {
			sscanf(line, "%*s %d", &data.Nth);
		}
		else if (strstr(line, "N_init")) {
			sscanf(line, "%*s %d", &data.N_init);
		}
		else if (strstr(line, "N_seeds")) {
			sscanf(line, "%*s %d", &data.N_seeds);
		}
		else if (strstr(line, "N_steps")) {
			sscanf(line, "%*s %d", &data.N_steps);
		}
		else if (strstr(line, "MAXTHRESHOLDS")) {
			sscanf(line, "%*s %d", &data.MAXTHRESHOLDS);
		}
		else if (strstr(line, "Bdef")) {
			sscanf(line, "%*s %lf %lf", &data.lb, &data.ub);
		}
		else if (strstr(line, "logval")) {
			sscanf(line, "%*s %d", &data.logval);
		}
		else if (strstr(line, "sigma")) {
			sscanf(line, "%*s %lf", &data.sigma);
		}
		else if (strstr(line, "seed")) {
			sscanf(line, "%*s %d", &data.seed);
		}
		else if (strstr(line, "iplot")) {
			sscanf(line, "%*s %d", &data.iplot);
		}
		else if (strstr(line, "factor")) {
			sscanf(line, "%*s %lf", &data.FACTOR);
		}
	}

	rewind(f);
	line_no = 0;

	free(data.lowerbound);
	free(data.upperbound);
	data.lowerbound = (double *)malloc(data.Nth*sizeof(double));
	data.upperbound = (double *)malloc(data.Nth*sizeof(double));

	for (i = 0; i < data.Nth; i++) {
		int found = 0;
		while (fgets(line, 256, f)!= NULL) {
			line_no++;

			if ((line[0] == '#')||(strlen(line)==0)) continue;

			char bound[8];
			sprintf(bound, "B%d", i);
			if (strstr(line, bound) != NULL) {
				sscanf(line, "%*s %lf %lf", &data.lowerbound[i], &data.upperbound[i]);
				found = 1;
				break;
			}	
		}
		if (!found) {
			data.lowerbound[i] = data.lb;	// Bdef value or Default LB
			data.upperbound[i] = data.ub;	// Bdef value of Default UB
		}
		rewind(f);
		line_no = 0;
	}

	free(data.threshold);
	data.threshold = (double *)calloc(1, data.MAXTHRESHOLDS*sizeof(double));	// goes to runinfo? 

	fclose(f);

}

void data_init()
{
	/* DATA: user's input parameters */
	read_data();
}

/**********************************************/
// RNG
/**********************************************/

void call_gsl_rand_init()
{
	printf("CALLING gsl_rand_init() on node %d\n", torc_node_id()); fflush(0);
	gsl_rand_init(data.seed);
}

void spmd_gsl_rand_init()
{
	int i;
	for (i = 0; i < torc_num_nodes(); i++) {
		torc_create_ex(i*torc_i_num_workers(), 1, (void *)call_gsl_rand_init, 0);
	}
	torc_waitall();
}

/**********************************************/
// FITFUN
/**********************************************/

#include "fitfun.c"

// 
double funeval(double *x, int n)
{
	inc_nfc();
	double f = fitfun(x, n, NULL, NULL);
	return f;
}

// # Demand
double demand(double sample[])
{
	return funeval(sample, data.Nth);
}

//# Parametrization
//# F_i = {demand > threshold_i}
int in_F(double fsample, int level)
{
	return fsample > data.threshold[level];
}

void taskf(double *point, double *fval)
{
	*fval = demand(point);// function evaluation
}

double fmodel(double *x, int n)
{
	inc_nfc();
	double f = fitfun(x, n, NULL, NULL);
	return f;
}

void taskf_model(double *x, int *pn, double *res)
{
	int n = *pn;
	*res = fmodel(x, n);
}

//# Modified Metropolis algorithm
void modified_metropolis(double seed[], double *pfseed, int *plevel, int *pN_steps)
{
	double fseed = *pfseed;
	int level = *plevel;
	int N_steps = *pN_steps;
	int n = data.Nth;
	
	double fleader = fseed;
	double leader[data.Nth];
	memcpy(leader, seed, data.Nth*sizeof(double));
	
//	samples = [seed]
	add_sample(leader, &fleader);

	double ftheta_leader = fmodel(leader, n);	// model evaluation, xxx: we must store it in the database (samples)!

	for (int i = 0; i < N_steps; i++) {
		//# Step 1: generate a candidate sample

		double xi[data.Nth];	// test point per direction
		for (int index = 0; index < data.Nth; index++) {
			//# generate xi from 1D uniform
			double spread = 2.0*data.sigma;
again:
			xi[index] = uniformrand(leader[index] - spread, leader[index] + spread);
			if ((xi[index] <= data.lowerbound[index]) || (xi[index] >= data.upperbound[index]))
				goto again;
		}

		double fthetaj[data.Nth];	// function value at each test point
		
		for (int index = 0; index < data.Nth; index++) {	// parallel
			double thetaj[data.Nth];	// test point

			memcpy(thetaj, leader, data.Nth*sizeof(double));
			thetaj[index] = xi[index];

                        int wid = torc_worker_id();     // -1
			torc_create(wid, (void *)taskf_model, 3,	// model evaluation
					data.Nth, MPI_DOUBLE, CALL_BY_COP,
					1, MPI_INT, CALL_BY_COP,
					1, MPI_DOUBLE, CALL_BY_RES,
					thetaj, &n, &fthetaj[index]);
			
		}
		torc_waitall();

		// build the candidate point and compute the function value there
		double candidate[data.Nth];
		for (int index = 0; index < data.Nth; index++) {
			//# compute the ratio                  
			double r, P;
			if (data.logval) {
				r = fthetaj[index] - ftheta_leader;
				P = log(uniformrand(0, 1));
			} else {
				r = (fthetaj[index] + 1e-6)/(ftheta_leader + 1e-6);
				if (r > 1) r = 1;
				P = (uniformrand(0, 1));
			}

			//# generate bernoulli rv to decide what to assign to candidate[index]
			if (P < r) {
				candidate[index] = xi[index];
			} else {
				candidate[index] = leader[index];
			}
		}
		double fcandidate = demand(candidate);

		//# Step 2: accept/reject the candidate
		if (in_F(fcandidate, level)) {
			memcpy(leader, candidate, data.Nth*sizeof(double));
			fleader = fcandidate;
			ftheta_leader = fmodel(leader, n);	// new leader, recompute ftheta_leader (not much to do, we can mix it with the parallel fevals that follow)
			add_sample(leader, &fleader);
		}
#if 0
		else {
			add_sample(leader, &fleader);
		}
#endif
	}
}


/**********************************************/
// EXTRA
/**********************************************/

void upd_threshold_task()
{
	MPI_Bcast(data.threshold, data.MAXTHRESHOLDS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void bcast_threshold()
{
	int i;
	
	if (torc_num_nodes() == 1) return;
	for (i = 0; i < torc_num_nodes(); i++) {
		torc_create_ex(i*torc_i_num_workers(), 1, (void *)upd_threshold_task, 0);
	}
	torc_waitall();
}

/**********************************************/
// MAIN
/**********************************************/

int main(int argc, char *argv[])
{
	double t0, t1;

	data_init();
	db_init();

	torc_register_task((void *)upd_threshold_task);
	torc_register_task((void *)reset_nfc_task);
	torc_register_task((void *)get_nfc_task);
	torc_register_task((void *)add_seed_task);
	torc_register_task((void *)add_sample_task);
	torc_register_task((void *)taskf_model);
	torc_register_task((void *)call_gsl_rand_init);
	torc_register_task((void *)taskf);
	torc_register_task((void *)modified_metropolis);

	torc_init(argc, argv, MODE_MW);

	spmd_gsl_rand_init();

	char *filename = (char *)"Pc.txt";
	unlink (filename);	// remove file 

	double P = 1.0; // # estimate of P(F), will be updated

//	# Step 1:
//	# Estimate P(F_1) using direct MCS
//	# P_1 := P(F_1) = #{samples in F_1} / #{all samples}
//	# Return samples q(*|F_1) which lie F_1

	int nthreads = torc_num_workers();
	
	reset_nfc();	// not required, already zero

	t0 = torc_gettime();
	set_nsamples(0);	// reset
	set_nseeds(0);		// reset

	system("shuf < curgen_db > init_samples.txt");  

	FILE *fps = fopen("init_samples.txt", "r");	// assumptions: the file should have been reshuffled, N_init <= #samples
	{
	for (int count = 0; count < data.N_init; count++) {
		double *a_sample = get_sample(count);
		if (fps == NULL) {
			// # generate candidate sample
			for (int d = 0; d < data.Nth; d++) {
				a_sample[d] = uniformrand(data.lowerbound[d], data.upperbound[d]);	// direct access to the database
			}
		}
		else {
			// # load candidate sample
			for (int d = 0; d < data.Nth; d++) {
				fscanf(fps, "%lf", &a_sample[d]);
			}
			double dummy;
			for (int d = 0; d < 1; d++) fscanf(fps, "%lf", &dummy);	// skip (the error and) the fitness value
		}

		torc_create(-1, (void *)taskf, 2,			// demand
				data.Nth, MPI_DOUBLE, CALL_BY_VAL,
				1, MPI_DOUBLE, CALL_BY_RES,
				a_sample, &a_sample[data.Nth]);

		set_nsamples(get_nsamples()+1);	// nsample++;

	}
	torc_waitall();
	}
	t1 = torc_gettime();

	if (fps != NULL) fclose(fps);
	int g_nfeval = get_nfc();
	printf("init. time = %lf ms for %d function calls (%lf ms/call)\n", 1000*(t1-t0), g_nfeval, nthreads*1000*(t1-t0)/g_nfeval);

#if 1   /* SPECIAL */
	dump_samples(0);
//	return 0;
#endif

	data.threshold[0] = 0;
	bcast_threshold(); // not necessary here

	// # check which samples are in F_0
	for (int count = 0; count < data.N_init; count++) {
		double *a_sample = get_sample(count);
		//if (in_F(samples[count][data.Nth], 0)) {	// all in
			add_seed(a_sample, &a_sample[data.Nth]);
		//}
	}

	// # update estimator
	P = P*get_nseeds()/data.N_init;
	P = pow(data.FACTOR, 1);
	
	FILE *fp;
	fp = fopen(filename, "a");
	fprintf(fp, "%.3e\n", P);
	
//	# Step 2: i=1 to len(threshold)-1
//	# Starting from each of q(*|F_i) simulate MC samples using modified Metropolis method
//	# Return samples q(*|F_{i+1}) which lie in F_{i+1}
//	# Estimate P_{i+1} := P(F_{i+1}|F_i) = #{samples in F_{i+1}} / #{all samples from q(*|F_i)}

//	# take the first N_seeds elements, ignore the rest

	sort_seeds(get_nseeds(), 1);
	set_nseeds((int) (data.FACTOR * get_nseeds()));
	printf("nseed = %d\n", get_nseeds());
	printf("L = %d, T = %lf, P = %.3e\n", 0, get_seed_f(get_nseeds()), P); //seeds[nseed*(data.Nth+1)+data.Nth], P);
	data.threshold[0] = get_seed_f(get_nseeds()); //seeds[nseed*(data.Nth+1)+data.Nth];	// set and broadcast (current position or even the full array)
	bcast_threshold();

	printf("Threshold at %lf\n", get_seed_f(get_nseeds())); //seeds[nseed*(data.Nth+1)+data.Nth]);
	if (get_nseeds() > data.N_seeds) {
		permute_seeds(get_nseeds());
		set_nseeds(data.N_seeds); 
	}

	int level;
	for (level = 0; level < data.MAXTHRESHOLDS-1; level++) {

		cout << "we have " << get_nseeds() << " seeds\n";
		dump_seeds(level);
		dump_samples(level);

		t0 = torc_gettime();
		set_nsamples(0);	// reset

		reset_nfc();

		for (int seed_i = 0; seed_i < get_nseeds(); seed_i++) {
			// # generate more samples with the same distribution
			double *a_seed = get_seed(seed_i);
			torc_create(-1, (void *)modified_metropolis, 4,
							data.Nth, MPI_DOUBLE, CALL_BY_VAL,
							1, MPI_DOUBLE, CALL_BY_VAL,
							1, MPI_INT, CALL_BY_VAL,
							1, MPI_INT, CALL_BY_VAL,
							a_seed, &a_seed[data.Nth], &level, &data.N_steps);
							//&seeds[seed_i*(data.Nth+1)], &seeds[seed_i*(data.Nth+1)+data.Nth], &level, &data.N_steps);

		}
		torc_enable_stealing();
		torc_waitall();
		torc_disable_stealing();

		t1 = torc_gettime();
		int g_nfeval = get_nfc();
		printf("level %d time = %lf ms for %d function calls (%lf ms/call)\n", level, 1000*(t1-t0), g_nfeval, 1000*(t1-t0)/g_nfeval);


		double N_total = get_nsamples();

		// # get seeds for the next MC
		set_nseeds(0);	// reset
		for (int sample_i = 0; sample_i < get_nsamples(); sample_i++) {
			// # check if the sample is in F_{i+1}
			double fsample_i = get_sample_f(sample_i);
			//if (in_F(fsample_i, level+1)) {	// all in
					add_seed(get_sample(sample_i), &fsample_i);
			//}
		}

		std::cout << "generated " << N_total << " samples; " << get_nseeds() << " of them can be new seeds\n";
		
		if (N_total == 0) break; // no samples, no seeds

		// # update estimator
		//P = P*nseed/N_total;
		P = pow(data.FACTOR, level+2);
		//printf("L = %d, T = %lf, P = %.3e\n", level+1, data.threshold[level+1], P);

		if (P == 0) break;

		// # take the first N_seeds elements, ignore the rest
		sort_seeds(get_nseeds(), 1);
		set_nseeds((int) (data.FACTOR * get_nseeds()));
		printf("nseed = %d\n", get_nseeds());
		printf("Threshold at %lf\n", get_seed_f(get_nseeds()));
		printf("L = %d, T = %lf, P = %.3e\n", level+1, get_seed_f(get_nseeds()), P); //seeds[nseed*(data.Nth+1)+data.Nth], P);
		data.threshold[level+1] = get_seed_f(get_nseeds());	//seeds[nseed*(data.Nth+1)+data.Nth];	// broadcast
		bcast_threshold();

		if (get_nseeds() > data.N_seeds) {
			permute_seeds(get_nseeds());
			set_nseeds(data.N_seeds); 
		}
		fprintf(fp, "%.3e\n", P);

		//# Now P(F) = P_1*...*P_m
	}
	
	dump_seeds(level);
	fclose(fp);

	printf("Total #nfeval = %d\n", get_tfc());

	torc_finalize();
	return 0;
}


