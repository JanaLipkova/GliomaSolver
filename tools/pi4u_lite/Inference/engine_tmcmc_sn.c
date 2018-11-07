/*
 *  engine_tmcmc_sn.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#define _XOPEN_SOURCE 500
#define _BSD_SOURCE
#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include "engine_tmcmc.h"
#include <pndlc.h>

//#define _RESTART_
#define _STEALING_
//#define _AFFINITY_

// extra parameters for sn-tmcmc
struct sn_data_t {
	double diffstep;
	int posdef_method;
} sn_data = { 1e-4, 4};

data_t data;
runinfo_t runinfo;
cgdb_t curgen_db;
db_t full_db;
resdb_t curres_db;

void read_data()
{
	int i;

	/* DEFAULT VALUES */
	data.Nth = 4; // Default PROBDIM
	data.MaxStages = 20; // Default MAXGENS;
	data.PopSize = 1024;	// Default DATANUM

	data.lb = -6;	// Default LB, same for all
	data.ub = +6;	// Default UB, same for all
 
	data.lowerbound = malloc(data.Nth*sizeof(double));
	data.upperbound = malloc(data.Nth*sizeof(double));

	for (i = 0; i < data.Nth; i++) {
		data.lowerbound[i] = data.lb;
		data.upperbound[i] = data.ub;
	}

	data.TolCOV = 1.0;	// 0.25, 0.5
	data.bbeta = 0.2;
	data.seed = 280675;

	data.options.MaxIter = 100;	// 
	data.options.Tol = 1e-6;
	data.options.Display = 0;

	data.iplot = 0;	// gnuplot

	data.Num = malloc(data.MaxStages*sizeof(double));
	for (i = 0; i < data.MaxStages; i++) {
		data.Num[i] = data.PopSize; // default DATANUM
	}
	data.LastNum = data.PopSize; // DATANUM;

//	sd_data.diffstep = 1e-4;
//	sd_data.posdef_method = 4;

	/* USER-DEFINED VALUES */
	FILE *f = fopen("tmcmc.par", "r");
	if (f == NULL) {
		return;
	}

	/*
	Nth             2
	MaxStages	200
	PopSize         1024
	Bdef		-4	4
	#B0              -6	6
	#B1              -6	6
	TolCOV          1
	bbeta           0.2
	seed            280675
	opt.MaxIter     100
	opt.Tol         1e-6
	opt.Display     0
	iplot           0
	diffstep	1e-4
	posdef		4
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
		else if (strstr(line, "MaxStages")) {
			sscanf(line, "%*s %d", &data.MaxStages);
		}
		else if (strstr(line, "PopSize")) {
			sscanf(line, "%*s %d", &data.PopSize);
		}
		else if (strstr(line, "TolCOV")) {
			sscanf(line, "%*s %lf", &data.TolCOV);
		}
		else if (strstr(line, "bbeta")) {
			sscanf(line, "%*s %lf", &data.bbeta);
		}
		else if (strstr(line, "seed")) {
			sscanf(line, "%*s %d", &data.seed);
		}
		else if (strstr(line, "opt.MaxIter")) {
			sscanf(line, "%*s %d", &data.options.MaxIter);
		}
		else if (strstr(line, "opt.Tol")) {
			sscanf(line, "%*s %lf", &data.options.Tol);
		}
		else if (strstr(line, "opt.Display")) {
			sscanf(line, "%*s %d", &data.options.Display);
		}
		else if (strstr(line, "iplot")) {
			sscanf(line, "%*s %d", &data.iplot);
		}
		else if (strstr(line, "Bdef")) {
			sscanf(line, "%*s %lf %lf", &data.lb, &data.ub);
		}
		else if (strstr(line, "diffstep")) {
			sscanf(line, "%*s %lf", &sn_data.diffstep);
		}
		else if (strstr(line, "posdef")) {
			sscanf(line, "%*s %d", &sn_data.posdef_method);
		}
	}

	rewind(f);
	line_no = 0;

	free(data.lowerbound);
	free(data.upperbound);
	data.lowerbound = malloc(data.Nth*sizeof(double));
	data.upperbound = malloc(data.Nth*sizeof(double));

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


	fclose(f);

	free(data.Num);
	data.Num = malloc(data.MaxStages*sizeof(double));
	for (i = 0; i < data.MaxStages; i++) {
		data.Num[i] = data.PopSize;
	}
	data.LastNum = data.PopSize;

}

void data_init()
{
	int i;

	/* DATA: user's input parameters */
	read_data();

#if 1
	init_curgen_db();
	init_curres_db();
	init_full_db();
#endif

	/* RUNINFO: running state */
	runinfo.CoefVar = calloc(1, data.MaxStages*sizeof(double));
	runinfo.p = calloc(1, data.MaxStages*sizeof(double));
	runinfo.currentuniques = calloc(1, data.MaxStages*sizeof(double));
	runinfo.logselection = calloc(1, data.MaxStages*sizeof(double));
	runinfo.acceptance = calloc(1, data.MaxStages*sizeof(double));

	double *SSmem = (double *)calloc(1, data.Nth*data.Nth*sizeof(double));
	runinfo.SS = (double **)malloc(data.Nth*sizeof(double *));
	for (i = 0; i < data.Nth; i++) {
		runinfo.SS[i] = SSmem + i*data.Nth; //&SSmem[i*data.Nth];
	}

	runinfo.meantheta = calloc(1, data.MaxStages*sizeof(double *));
	for (i = 0; i < data.MaxStages; i++) {
		runinfo.meantheta[i] = calloc(1, data.Nth*sizeof(double));
	}

	runinfo.Gen = 0;	// TODO, MPI: global shared memory (or per task)
	runinfo.CoefVar[0] = 10;

	/* already zero */
	runinfo.p[0] = 0;
	for (i = 0; i < data.MaxStages; i++) {
		runinfo.currentuniques[i] = 0;
		runinfo.logselection[i] = 0.0;
		runinfo.acceptance[i] = 0.0;
	}


	printf("runinfo = %p\n", &runinfo);
	printf("runinfo.p = %p\n", runinfo.p);
	printf("runinfo.SS = %p\n", runinfo.SS);
}


void save_runinfo()
{
	int i, j;

	/* allocate and initialize runinfo */
	FILE *f = fopen("runinfo.txt", "w");

	fprintf(f, "Gen=\n");
	fprintf(f, "%d\n", runinfo.Gen);
	
	fprintf(f, "CoefVar[%d]=\n", data.MaxStages);
	for (i = 0; i < data.MaxStages; i++) fprintf(f, "%.16lf\n", runinfo.CoefVar[i]);

	fprintf(f, "p[%d]=\n", data.MaxStages);
	for (i = 0; i < data.MaxStages; i++) fprintf(f, "%.16lf\n", runinfo.p[i]);

	fprintf(f, "currentuniques[%d]=\n", data.MaxStages);
	for (i = 0; i < data.MaxStages; i++) fprintf(f, "%d\n", runinfo.currentuniques[i]);

	fprintf(f, "logselection[%d]=\n", data.MaxStages);
	for (i = 0; i < data.MaxStages; i++) fprintf(f, "%.16lf\n", runinfo.logselection[i]);

	fprintf(f, "acceptance[%d]=\n", data.MaxStages);
	for (i = 0; i < data.MaxStages; i++) fprintf(f, "%.16lf\n", runinfo.acceptance[i]);

	fprintf(f, "SS[%d][%d]=\n", data.Nth, data.Nth);
	for (i = 0; i < data.Nth; i++)
		for (j = 0; j < data.Nth; j++)
			fprintf(f, "%.16lf\n", runinfo.SS[i][j]);

	fprintf(f, "meantheta[%d][%d]\n", data.MaxStages, data.Nth);
	for (i = 0; i < data.MaxStages; i++)
		for (j = 0; j < data.Nth; j++)
			fprintf(f, "%.16lf\n", runinfo.meantheta[i][j]);

	fclose(f);
}

int load_runinfo()
{
	int i, j;
	char header[256];

	/* allocate and initialize runinfo */
	FILE *f = fopen("runinfo.txt", "r");
	if (f == NULL) return 1;
	
	fscanf(f, "%s", header);
	fscanf(f, "%d", &runinfo.Gen);

	fscanf(f, "%s", header);
	for (i = 0; i < data.MaxStages; i++) fscanf(f, "%lf\n", &runinfo.CoefVar[i]);

	fscanf(f, "%s", header);
	for (i = 0; i < data.MaxStages; i++) fscanf(f, "%lf\n", &runinfo.p[i]);

	fscanf(f, "%s", header);
	for (i = 0; i < data.MaxStages; i++) fscanf(f, "%d\n", &runinfo.currentuniques[i]);

	fscanf(f, "%s", header);
	for (i = 0; i < data.MaxStages; i++) fscanf(f, "%lf\n", &runinfo.logselection[i]);

	fscanf(f, "%s", header);
	for (i = 0; i < data.MaxStages; i++) fscanf(f, "%lf\n", &runinfo.acceptance[i]);

	fscanf(f, "%s", header);
	for (i = 0; i < data.Nth; i++)
		for (j = 0; j < data.Nth; j++)
			fscanf(f, "%lf\n", &runinfo.SS[i][j]);

	fscanf(f, "%s", header);
	for (i = 0; i < data.MaxStages; i++)
		for (j = 0; j < data.Nth; j++)
			fscanf(f, "%lf\n", &runinfo.meantheta[i][j]);

	fclose(f);

	return 0;
}

#include <signal.h>
static int exit_signal_flag = 0;	// new

void check_for_exit()
{
	int val, exitgen = -1;
	char *s;

	s = (char *) getenv("EXIT_GEN");
	if (s != 0 && sscanf(s, "%d", &val) == 1 && val >= 0)
		exitgen = val;

	if (exitgen == runinfo.Gen) {
		printf("Read Exit Envrironment Variable!!!\n");
		torc_finalize();
		exit(1);
	}

	if (exit_signal_flag == 1) {
		printf("Received Exit Signal!!!\n");
		torc_finalize();
		exit(1);	// new 
	}

	FILE *fp;
	fp = fopen("exit.txt", "r");
	if (fp != NULL) {
		printf("Found Exit File!!!\n");
		unlink("exit.txt");
		torc_finalize();
		exit(1);
	}
}

// this is effective only on the master process
void handle_signal(int signal)
{
	switch (signal) {
	case SIGUSR1:
		printf("Caught SIGUSR1, exiting at next checkpoint!\n");
		exit_signal_flag = 1;
		break;
	default:
		fprintf(stderr, "Caught wrong signal: %d\n", signal);
		return;
	}
}

void setup_handler()
{
	struct sigaction sa;
	sa.sa_handler = &handle_signal;
	sa.sa_flags = SA_RESTART;
	sigfillset(&sa.sa_mask);

	if (sigaction(SIGUSR1, &sa, NULL) == -1) {
		perror("Error: cannot handle SIGUSR1"); // Should not happen
	}
}

#if 1	/* TORC-BASED DATA MANAGEMENT */
void torc_update_full_db_task(double point[], double *pF, double *G, int *pn, int *psurrogate)
{
	double F = *pF;
	int n = *pn;
	int surrogate = *psurrogate;

//	printf("torc_update_full_db_task: %f %f %f %f %f\n", point[0], point[1], point[2], F, G[0]); fflush(0);

	update_full_db(point, F, G, n, surrogate);
}

void torc_update_full_db(double point[], double F, double *G, int n, int surrogate)
{
	if (torc_node_id() ==0) {
		update_full_db(point, F, G, n, surrogate);
		return;
	}
	torc_create_direct(0, torc_update_full_db_task, 5,		/* message to the database manager (separate process?) or direct execution by server thread */
		data.Nth, MPI_DOUBLE, CALL_BY_VAL,
		1, MPI_DOUBLE, CALL_BY_COP,
		n, MPI_DOUBLE, CALL_BY_COP,	/* check this for CALL_BY_VAL: in the full-version of the library, with n=1 we had segv */
		1, MPI_INT, CALL_BY_COP,
		1, MPI_INT, CALL_BY_COP,
		point, &F, G, &n, &surrogate);
	torc_waitall3();
}

void torc_update_curgen_db_task(double point[], double *pF)
{
	double F = *pF;

//	printf("torc_update_curgen_db_task: %f %f %f %f\n", point[0], point[1], point[2], F); fflush(0);

	update_curgen_db(point, F);
}

void torc_update_curgen_db(double point[], double F)
{
	int me = torc_node_id();
//	printf("%d torc_update_curgen_db_task: %f %f %f %f\n", me, point[0], point[1], point[2], F); fflush(0);
	if (me == 0) {
		update_curgen_db(point, F);
		return;
	}
	torc_create_direct(0, torc_update_curgen_db_task, 2,		/* message to the database manager (separate process?) or direct execution by server thread */
		data.Nth, MPI_DOUBLE, CALL_BY_COP,
		1, MPI_DOUBLE, CALL_BY_COP,
		point, &F);
	torc_waitall3();	// wait without releasing the worker
}

void torc_update_curgen_db_der_task(double point[], double *pF, double grad[], double hes[])
{
	double F = *pF;

//	printf("torc_update_curgen_db_der_task: %f %f %f %f\n", point[0], point[1], point[2], F); fflush(0);

	update_curgen_db_der(point, F, grad, hes);
}

void torc_update_curgen_db_der(double point[], double F, double grad[], double hes[])
{
	int me = torc_node_id();
//	printf("%d torc_update_curgen_db_der_task: %f %f %f %f\n", me, point[0], point[1], point[2], F); fflush(0);
	if (me == 0) {
		update_curgen_db_der(point, F, grad, hes);
		return;
	}
	torc_create_direct(0, torc_update_curgen_db_task, 4,		/* message to the database manager (separate process?) or direct execution by server thread */
		data.Nth, MPI_DOUBLE, CALL_BY_COP,
		1, MPI_DOUBLE, CALL_BY_COP,
		data.Nth, MPI_DOUBLE, CALL_BY_COP,
		data.Nth*data.Nth, MPI_DOUBLE, CALL_BY_COP,
		point, &F, grad, hes);
	torc_waitall3();	// wait without releasing the worker
}

void torc_update_curres_db_task(double point[EXPERIMENTAL_RESULTS], double *pF)
{
	double F = *pF;

//	printf("torc_update_curres_db_task: %f %f %f\n", point[0], point[1], point[2]); fflush(0);

	update_curres_db(point, F);
}

void torc_update_curres_db(double point[EXPERIMENTAL_RESULTS], double F)
{
	int me = torc_node_id();
//	printf("%d torc_update_curres_db_task: %f %f %f\n", me, point[0], point[1], point[2]); fflush(0);
	if (me ==0) {
		update_curres_db(point, F);
		return;
	}
	torc_create_direct(0, torc_update_curres_db_task, 2,		/* message to the database manager (separate process?) or direct execution by server thread */
		EXPERIMENTAL_RESULTS, MPI_DOUBLE, CALL_BY_COP,
		1, MPI_DOUBLE, CALL_BY_COP,
		point, &F);
	torc_waitall3();
}

#endif


/*** TASK MANAGEMENT ***/
#include "fitfun.c" 

void taskfun(double /*const*/ *x, int *pN, double *res, int winfo[4])
{
	double f;
	int N = *pN;
//	printf("taskfun {%d,%d,%d,%d}\n", winfo[0], winfo[1], winfo[2], winfo[3]);

	inc_nfc();	// increment function call counter

	f = fitfun(x, N, (void *)NULL, winfo);
#if (EXPERIMENTAL_RESULTS > 0)	/* peh: decide about this (results should be passed as argument to fitfun) */
	int i;
	double results[EXPERIMENTAL_RESULTS];
	for (i = 0; i < EXPERIMENTAL_RESULTS; i++) {
		if (i < data.Nth)
			results[i] = x[i];
		else
			results[i] = 0.0;
	}
	torc_update_curres_db(results, f);
#endif

	*res = f;
	return;
}

double F(double *TP, int *pn)	/* for PNDL */
{
	double gres;

	taskfun(TP, pn, &gres, NULL);
//	gres = posterior(TP, PROBDIM, gres);

	return gres;
}

void evaluate_F(double point[], double *Fval, int worker_id, int gen_id, int chain_id, int step_id, int ntasks)
{
	int i;
	double G[64], F;	// maxtasks
	int winfo[4];
	int dim = data.Nth;

	for (i = 0; i < ntasks; i++) {
		winfo[0] = gen_id;
		winfo[1] = chain_id;
		winfo[2] = step_id;
		winfo[3] = i;

		if (ntasks == 1) {
			taskfun(point, &dim, &G[i], winfo);
		}
		else {
			torc_create(-1, taskfun, 4,
				data.Nth, MPI_DOUBLE, CALL_BY_VAL,
				1, MPI_INT, CALL_BY_COP,
				1, MPI_DOUBLE, CALL_BY_RES,
				4, MPI_INT, CALL_BY_COP,
				point, &dim, &G[i], winfo);
		}

	}
	if (ntasks > 1)
		torc_waitall();

	/* compute F from G's - how? : F fitness function = distance from ground truth */
	F = 0.0;
	for (i = 0; i < ntasks; i++) {
		F += G[i];
	}
	F = F/ntasks;

	//update_full_db(point, F, G, ntasks, 0); // not surrogate 
	//...torc_update_full_db(point, F, G, ntasks, 0); // not surrogate 

	*Fval = F;
}

void initchaintask(double in_tparam[], int *pdim, double *out_tparam, int winfo[4])
{
	int i;
	int gen_id, chain_id;
	gen_id = winfo[0];
	chain_id = winfo[1];

	long me = torc_worker_id();	//psthread_current_vp();
	double point[data.Nth], fpoint;

	for (i = 0; i < data.Nth; i++)
		point[i] = in_tparam[i];

	evaluate_F(point, &fpoint, me, gen_id, chain_id, 0, 1);	// 12 for namd

	/* update current db entry */
	torc_update_curgen_db(point, fpoint);
	*out_tparam = fpoint;	// currently not required, the result is already in the db

	return;
}


//~~~~~~~~~~~~~~~~~~~~~~~~MCMC with NEWTON UPDATE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//damianos
//==================================================//
//===========compute_moving_probab==================//
//==================================================//
/**
 * input(1): candidate point
 * input(2): current leader point
 * input(3): gradient on leader
 * input(4): hessian on leader
 * input(5): inverse of current hessian
 * input(3): random sigma (needed for stohastic Newton)
 * input(4): annealing coefficient
 * input(5): calculated probability q(xk,y) (equation 2.3)
 * input(6): calculated probability q(y,xk) (equation 2.4)
 * 
 * output(): none
 * 
 * remarks: function to calculate the moving probability q(xk -> y) and q(y -> xk) equation 2.3 and 2.4 of page 4 
 * from publication "A scaled Stohastic Newton algorithm for MCMC"
 */

void compute_grad(double point_to_compute[], double computed_gradient[]);	// forward declaration

void compute_moving_probab(double candidate[], double leader[], double current_gradient[],
	double *current_hessian/*2D*/, double *inv_current_hessian/*2D*/, double rand_sigma, double pj, double *q_xk_y, double *q_y_xk)
{
	int row;
	double expont_xk_y[data.Nth], expont_y_xk[data.Nth];
	double hessian_product_exp_xk_y[data.Nth], hessian_product_exp_y_xk[data.Nth];
	double weig_norm_exp_xk_y, weig_norm_exp_y_xk;
	double inv_hessian_product_cur_gradient[data.Nth]; //vector computed by inv_hessian * cur_gradient
	double inv_hessian_product_cand_gradient[data.Nth]; // vector computed by inv_hessian * cand_gradient
	double coef_expont = -0.5 * pow(rand_sigma,-2.0);
	double gradient_on_candidate[data.Nth];	// 1st derivatives on the candidate point

	memset(inv_hessian_product_cur_gradient, 0, data.Nth*sizeof(double));
	memset(inv_hessian_product_cand_gradient, 0, data.Nth*sizeof(double));
	memset(gradient_on_candidate, 0, data.Nth*sizeof(double));

//	Matlab
//	logProXtoY=-0.5*(data.bbeta^-2)*(thetac-thetaleader-((0.5*(data.bbeta^2)*A*dy')'))*H*(thetac-thetaleader-((0.5*(data.bbeta^2)*A*dy')'))';
//	logProbYtoX=-0.5*(data.bbeta^-2)*(thetaleader-thetac-((0.5*(data.bbeta^2)*A*dyc')'))*H*(thetaleader-thetac-((0.5*(data.bbeta^2)*A*dyc')'))';

	//compute inv_hessian * grad(log(target(leader) ))
	compute_mat_product_vect(inv_current_hessian, current_gradient, inv_hessian_product_cur_gradient, 1.0, data.Nth);

	//compute the gradient on candidate point and then inv_hessian * grad(log(target(cand_point)))
	compute_grad(candidate, gradient_on_candidate);
	compute_mat_product_vect(inv_current_hessian,  gradient_on_candidate, inv_hessian_product_cand_gradient, 1.0, data.Nth);

	// compute the two exponents in 2.3 and 2.4
	//expont_xk_y ->  candidate - leader - (sigma^2 /2) * inv(hessian(log(target_function))) * grad(log(target_function(xk))) 
	//expont_y_xk -> leader - candidate - (sigma^2 /2) * inv(hessian(log(target_function))) * grad(log(target_function(y)))

	for(row=0;row < data.Nth; ++row){
		expont_xk_y[row] = candidate[row] - leader[row] - (0.5*pow(rand_sigma,2) * inv_hessian_product_cur_gradient[row] );
		expont_y_xk[row] = leader[row] - candidate[row] - (0.5*pow(rand_sigma,2) * inv_hessian_product_cand_gradient[row]);
	}

	/*                                                    expont_q_xk_y
	//                                   /				      \  
	//compute the exp ( (-1/2*(sigma^2)) * ||y - xk - ((sigma^2)/2)* A * grad ||H )
	//                        |              \    /                 \       /
	//                     coef_expont        xk_y        inv_hessian_product_cur_gradient                   
	*/

	// || exponent_xk_y ||Hessian -> transpose(exponent_xk_y) * Hessian * exponent_xk_y

	//so first coef_expont compute Hessian * exponent_xk_y
	compute_mat_product_vect(current_hessian, expont_xk_y, hessian_product_exp_xk_y, coef_expont, data.Nth);
	compute_mat_product_vect(current_hessian, expont_y_xk, hessian_product_exp_y_xk, coef_expont, data.Nth);

	//then compute transpose(exponent) * hessian_product_exponent
	weig_norm_exp_xk_y = compute_dot_product(expont_xk_y, hessian_product_exp_xk_y, data.Nth);
	weig_norm_exp_y_xk = compute_dot_product(expont_y_xk, hessian_product_exp_y_xk, data.Nth);

	//apply exp to get practical q_xk_y and q_y_xk
	*q_xk_y = exp(weig_norm_exp_xk_y);
	*q_y_xk = exp(weig_norm_exp_y_xk);
}

//======================================================//
//=====================force_symm=======================//
//======================================================//
/**
 * 
 * input(1): matrix to be symmetric (input-output)
 * 
 * output(0): none
 * 
 * remarks: copies the upper triangular matrix to each lower forcing the matrix to be symmetric
 * the result is saved to the input matrix  [in-place]
 */
void force_symm(double *input_mat/*2D*/)
{
	int diag_coord, moving_coord;
	double value_to_copy;
  
	for(diag_coord=0; diag_coord<data.Nth; diag_coord++){ //loop through the non diagonal elements and copy each row to its' corresponding column
		for(moving_coord=diag_coord+1; moving_coord<data.Nth; moving_coord++){
			value_to_copy = input_mat[diag_coord*data.Nth+moving_coord];
			input_mat[moving_coord*data.Nth+diag_coord] = value_to_copy;
		}
	}
}
//======================================================//
//======================================================//

//int not_pos_def_times = 0;
extern int not_pos_def_times;

//==================================================//
//=============compute_grad=========================//
//==================================================//


/**
 * input(1): vector with point to compute gradient (either leader or candidate)
 *
 * input(2): vector to save the computed gradient 
 * output(): none
 * 
 * remarks: find the gradient in the given input point of log(target_function)
 */

void compute_grad(double point_to_compute[], double computed_gradient[])
{
	int i, pdim = data.Nth;
	double FEPS = 1e-3;
	int IPRINT = 0;
	int NOC, IERR;
	int IORD = 2;
	double XL[data.Nth], XU[data.Nth], UH[data.Nth];

	for (i = 0; i < data.Nth; i++) {
		XL[i] = data.lowerbound[i];
		XU[i] = data.upperbound[i];
		UH[i] = sn_data.diffstep;
	}
  
	c_pndlga(F, point_to_compute, &pdim, XL, XU, UH, &FEPS, &IORD, &IPRINT, computed_gradient, &NOC,&IERR);
//	print_matrix("pndl grad", computed_gradient, data.Nth);
}

//==================================================//
//=========compute_gradient_and_hessian=============//
//==================================================//

void compute_gradient_and_hessian(double point_to_compute[], double computed_gradient[], double *current_hessian/*2D*/)
{
	int i, pdim = data.Nth;
	double FEPS = 1e-3;
	int IPRINT = 0;
	int NOC, IERR;
	int IORD = 2;
	double XL[data.Nth], XU[data.Nth], UH[data.Nth];

	for (i = 0; i < data.Nth; i++) {
		XL[i] = data.lowerbound[i];
		XU[i] = data.upperbound[i];
		UH[i] = sn_data.diffstep;
	}
  
	c_pndlghfa (F, point_to_compute, &pdim, XL, XU, UH, &FEPS, &IORD, &IPRINT, computed_gradient, current_hessian, &pdim, &NOC,&IERR);
	//fill the other half of the computed Hessian matrixs
	force_symm(current_hessian);

}

//==================================================//
//=============compute_hessian======================//
//==================================================//

/**
 * input(1): vector with leaders value, leader[]
 * input(2): array to save the hessian matrix
 * output(): none
 * 
 * remarks: find the hessian in leader's point of log(target_function)
 */
void compute_hessian(double leader[], double *current_hessian/*2D*/)
{
	int i, pdim = data.Nth;

	double FEPS = 1e-3;
	int IPRINT = 0;
	int NOC, IERR;
	int IORD = 2;
	double XL[data.Nth], XU[data.Nth], UH[data.Nth];

	for (i = 0; i < data.Nth; i++) {
		XL[i] = data.lowerbound[i];
		XU[i] = data.upperbound[i];
		UH[i] = sn_data.diffstep;
	}

	c_pndlhfa(F, leader, &pdim, XL, XU, UH, &FEPS, &IORD, &IPRINT, current_hessian, &pdim, &NOC,&IERR);
	//fill the other half of the computed Hessian matrixs
	force_symm(current_hessian);

//	my_print_matrix_2d("pndl hess", current_hessian);
}

//==================================================//
//==================================================//

//damianos
//==================================================//
//========compute_deterministic_part================//
//==================================================//
/**
 * input(1): leader coordinates, leader[]
 * input(2): current gradient
 * input(3): inverse of current hessian
 * input(4): deterministic part of the proposed sample
 * input(5): random sigma
 * input(6): annealing coef
 * 
 * output: none
 * 
 * remarks: compute the deterministic part of the equation 2.1
 * from publication "A scaled Stohastic Newton algorithm for MCMC"
 */
void compute_deterministic_part(double leader[], double current_gradient[], double *inv_current_hessian/*2D*/,
				double deterministic_part[], double rand_sigma, double pj)
{
	int row;
	double inv_hessian_product_cur_gradient[data.Nth];

	memset(inv_hessian_product_cur_gradient, 0, data.Nth*sizeof(double));
 
	//compute inv(H) * gradient(log(target_function))
	compute_mat_product_vect(inv_current_hessian, current_gradient, inv_hessian_product_cur_gradient, 1.0, data.Nth); //the inverse of hessian SHOULD be multiplied by the annealing coefficient

	//now compute the deterministic part of equation 2.1
	for(row=0; row<data.Nth; row++){
		deterministic_part[row] = leader[row] + (0.5*pow(rand_sigma,2.0))*inv_hessian_product_cur_gradient[row];
//		printf("%d : %f = %f + %f x %f\n", row, deterministic_part[row], leader[row], (pow(rand_sigma,2) / 2.0) * (anneal_coef), inv_hessian_product_cur_gradient[row]);
	}
}

//==================================================//
//==================================================//

int newton_aborts = 0;
int not_good = 0;


//damianos
//==================================================//
//=============newton_compute_candidate=============//
//==================================================//
/**
 *
 * input(1): candidate vector
 * input(2): leader vector
 * input(3): gradient on the current point (leader)
 * input(4): hessian on current point (leader)
 * input(5): inverse hessian on current point
 * input(6): deterministic part of proposed sample
 * remarks: compute y = x_k + ((sigma^2)/2) * (-inv(hessian)) * gradient(log(target_function)) + sigma N(0,A)
 * publication "A scaled Stohastic Newton algorithm for MCMC" equation (2.1) page 4
 */
int newton_compute_candidate(double candidate[], double leader[], double grad[], double hes[], double current_gradient[], double *current_hessian/*2D*/,
				double *inv_current_hessian/*2D*/, double *deterministic_part, double rand_sigma, double pj)

{
	int i,j;
	double mean_nrm[data.Nth];
	double sigma_nrm[data.Nth*data.Nth];
	double stohastic_part[data.Nth];
	int inv_matrix_flag = 0;

//	double deterministic_part2[data.Nth];

	memset(mean_nrm, 0, data.Nth*sizeof(double));
	memset(sigma_nrm, 0, data.Nth*data.Nth*sizeof(double));
	memset(stohastic_part, 0, data.Nth*sizeof(double));

//	memset(deterministic_part2, 0, data.Nth*sizeof(double));

	//if we have a new leader then compute the gradient, the hessian and its inverse
	//also calculate the deterministic part of equation 2.1

	// adjust the Hessian matrix
	for (i = 0; i < data.Nth; i++)
		for (j = 0; j < data.Nth; j++)
			current_hessian[i*data.Nth+j] = -pj*hes[i*data.Nth+j];

	inv_matrix_flag = inv_matrix(1.0, current_hessian, inv_current_hessian, data.Nth);
	if (inv_matrix_flag == 0) {
		newton_aborts++;
		return 0; // early reject and fall back to the normal selection policy 
	}
		
	// adjust the gradient vector
	for (i = 0; i < data.Nth; i++)
		current_gradient[i] = pj*grad[i];

	compute_deterministic_part(leader, current_gradient, inv_current_hessian, deterministic_part, rand_sigma, pj);
//	compute_deterministic_part(leader, current_gradient, inv_current_hessian, deterministic_part2, rand_sigma, pj);

	//now compute the stohastic part of the equation
	for (i = 0; i < data.Nth; i++){
		for (j = 0; j < data.Nth; j++){
			sigma_nrm[i*data.Nth+j] = inv_current_hessian[i*data.Nth+j]; // the hessian and its inverse are NOT multiplied by
		}
	}

	//draw a sample from multivariate normal with given mean and sigma
	for (i = 0; i < data.Nth; i++) mean_nrm[i] = 0.0;
	mvnrnd(mean_nrm, (double *)sigma_nrm, stohastic_part, data.Nth);

	//for each problem dimension compute the candidate point 
	for( i=0; i < data.Nth; i++) {
		candidate[i] = deterministic_part[i] + rand_sigma * stohastic_part[i];
//		candidate[i] = deterministic_part2[i] + rand_sigma * stohastic_part[i];
	}

	int flag = 0;
	//for each problem dimension compute the candidate point 
	for(i=0; i < data.Nth; i++){
		if (candidate[i] < data.lowerbound[i]) { candidate[i] = data.lowerbound[i]; flag = 1; }
		if (candidate[i] > data.upperbound[i]) { candidate[i] = data.upperbound[i]; flag = 1; }
	}

	if (flag) return 0; // bad point, reject proposal
	return 1;
}
//==================================================//
//==================================================//


//~~~~~~~~~~~~~~~~~~~~~~~~MCMC with NEWTON UPDATE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



void compute_candidate(double candidate[], double leader[], double var)
{
	int i, j;
	double bSS[data.Nth*data.Nth];

	for (i = 0; i < data.Nth; i++)
		for (j = 0; j < data.Nth; j++)
			bSS[i*data.Nth+j]= data.bbeta*runinfo.SS[i][j];

retry:
	mvnrnd(leader, (double *)bSS, candidate, data.Nth);
	for (i = 0; i < data.Nth; i++) {
		if (isnan(candidate[i])) {
			printf("!!!!  isnan in candidate point!\n");
			exit(1);
			break;
		}
		if ((candidate[i] < data.lowerbound[i])||(candidate[i] > data.upperbound[i])) break;
	}
	if (i < data.Nth) goto retry;
}

void chaintask(double in_tparam[], int *pdim, int *pnsteps, double *out_tparam, double *grad, double *hes, int winfo[4])
{
	int i,step;
	int nsteps = *pnsteps;
	int gen_id = winfo[0];
	int chain_id = winfo[1];
	
	long me = torc_worker_id();

	double current_gradient[data.Nth];			// 1st derivatives on the leader point
	double current_hessian[data.Nth*data.Nth];		// hessian
	double inv_current_hessian[data.Nth*data.Nth];		// inverse of hessian
	double deterministic_part[data.Nth];			//deterministic part of equation (2.1)

	double bare_gradient[data.Nth];
	double bare_hessian[data.Nth*data.Nth];

	double leader[data.Nth], fleader, fpc_leader;			// fold
	double candidate[data.Nth], fcandidate, fpc_candidate;		// fnew
	double q_xk_y, q_y_xk;
	double rand_sigma = data.bbeta; // = 0.2;
	
	memset(current_gradient, 0, data.Nth*sizeof(double));
	memset(current_hessian, 0, data.Nth*data.Nth*sizeof(double));
	memset(inv_current_hessian, 0, data.Nth*data.Nth*sizeof(double));
	memset(deterministic_part, 0, data.Nth*sizeof(double));

	for (i = 0; i < data.Nth; i++) leader[i] = in_tparam[i]; //chainwork->in_tparam[i];	// get initial leader
	fleader = *out_tparam;									// and its value
	for (i = 0; i < data.Nth; i++) bare_gradient[i] = grad[i];					// and its gradient
	for (i = 0; i < data.Nth*data.Nth; i++) bare_hessian[i] = hes[i];				// and its hessian
	fpc_leader = posterior(leader, data.Nth, fleader);
	
	double pj = runinfo.p[runinfo.Gen];
	for (step = 0; step < nsteps; step++) {
		//damianos
		//compute the proposed sample by the stohastic Newton algorithm
		int good = newton_compute_candidate(candidate, leader, bare_gradient, bare_hessian, current_gradient, current_hessian, inv_current_hessian, deterministic_part, 
				rand_sigma, pj);
		if (!good) { not_good++;
			compute_candidate(candidate, leader, rand_sigma);

			evaluate_F(candidate, &fcandidate, me, gen_id, chain_id, step, 1);
			fpc_candidate = posterior(candidate, data.Nth, fcandidate);

			double prior_candidate = priorpdf(candidate, data.Nth); // from PanosA
			double prior_leader = priorpdf(leader, data.Nth);
			double L = exp((prior_candidate-prior_leader)+(fpc_candidate-fpc_leader)*pj);
			//double L = exp((fpc_candidate-fpc_leader)*pj);

			double P = uniformrand(0,1);
			if (L > 1) L = 1;
			if (P < L) {
				for (i = 0; i < data.Nth; i++) leader[i] = candidate[i]; // new leader!
				fleader = fcandidate;
				fpc_leader = fpc_candidate;
				compute_gradient_and_hessian(leader, bare_gradient, bare_hessian);
				torc_update_curgen_db_der(leader, fleader, bare_gradient, bare_hessian);
			}
			else {
				torc_update_curgen_db_der(leader, fleader, bare_gradient, bare_hessian);
			}
		}
		else {

			evaluate_F(candidate, &fcandidate, me, gen_id, chain_id, step, 1);
			fpc_candidate = posterior(candidate, data.Nth, fcandidate);

			/* Decide */
			{
			//damianos
			//compute the probability of moving from the current sample to the proposed and from the proposed to the current
			compute_moving_probab(candidate, leader, current_gradient, current_hessian, inv_current_hessian, rand_sigma, pj, &q_xk_y, &q_y_xk);

//			double log_target_product_moving_xk_y =  pj*fpc_leader + log(q_xk_y );
//			double log_target_product_moving_y_xk =  pj*fpc_candidate + log(q_y_xk );
//			double L = exp( (log_target_product_moving_y_xk - log_target_product_moving_xk_y) );
//			double L = exp((fpc_candidate-fpc_leader)*pj - log(q_xk_y) + log(q_y_xk)); 

//			double L = exp((fpc_candidate-fpc_leader)*pj - q_xk_y + q_y_xk); 
			double prior_candidate = priorpdf(candidate, data.Nth); // from PanosA
			double prior_leader = priorpdf(leader, data.Nth);
			double L = exp((prior_candidate-prior_leader)+(fpc_candidate-fpc_leader)*pj - q_xk_y + q_y_xk);

			double P = uniformrand(0,1);
			if (L > 1) L = 1;
			if (P < L) {
				for (i = 0; i < data.Nth; i++) leader[i] = candidate[i];	// new leader! 
				fleader = fcandidate;
				fpc_leader = fpc_candidate;

				compute_gradient_and_hessian(leader, bare_gradient, bare_hessian);
				torc_update_curgen_db_der(leader, fleader, bare_gradient, bare_hessian);
			}
			else {
				torc_update_curgen_db_der(leader, fleader, bare_gradient, bare_hessian);
			
			}
			} // decide

		}
	}
	return;
}

#if 1
typedef struct {
	int idx;
	int nsel;
	double F;
} sort_t;

int compar_desc(const void* p1, const void* p2)
{
	int dir = +1;   // -1: ascending order, +1: descending order
	sort_t *s1 = (sort_t *) p1;
	sort_t *s2 = (sort_t *) p2;

	if (s1->nsel < s2->nsel) return dir;
	if (s1->nsel > s2->nsel) return -dir;
//	if (s1->nsel == s2->nsel) return 0;
	return 0;
}
#endif

int prepare_newgen(int nchains, cgdbp_t *leaders)
{
	/* process curgen_db -> calculate statitics */
	/* compute probs based on F values */
	/* draw new samples (nchains or user-specified) */
	/* find unique samples: fill the (new) leaders table */
	/* count how many times they appear -> nsteps */
	/* return the new sample size (number of chains) */

	int i, p;
	int newchains; // = nchains;

	int n = curgen_db.entries;
	double fj[n];
	unsigned int sel[n];

	double **g_x;
	g_x = (double **)malloc(data.Nth*sizeof(double *));
	for (i = 0; i < data.Nth; i++)
		g_x[i] = (double *)malloc(n*sizeof(double));

	{//start block
	double **x = g_x;

	for (p = 0; p < data.Nth; p++) {
		for (i = 0; i < n; i++) {
			x[p][i] = curgen_db.entry[i].point[p];
		}
	}

	double meanx[data.Nth], stdx[data.Nth];
	for (p = 0; p < data.Nth; p++) {
		meanx[p] = compute_mean(x[p], n);
		stdx[p] = compute_std(x[p], n, meanx[p]);
	}

	printf("CURGEN DB (COMPLE) %d\n", runinfo.Gen);
	print_matrix("means", meanx, data.Nth);
	print_matrix("std", stdx, data.Nth);
	}//end block

	{//start block
	double **x = g_x;
	int un = 0, unflag, j;

	for (p = 0; p < data.Nth; p++) {
		x[p][un] = curgen_db.entry[0].point[p];	// un==0
	}
	un++;
	for (i = 1; i < n; i++) {
		double xi[data.Nth];
		for (p = 0; p < data.Nth; p++) {
			xi[p] = curgen_db.entry[i].point[p];
		}
		unflag = 1;	// is this point unique?
		for (j = 0; j < un; j++) {	// check all the previous unique points
			int compflag;
			compflag = 1;	// 
			for (p = 0; p < data.Nth; p++) {
				if (fabs(xi[p]-x[p][j]) > 1e-6) {
				//if (xi[p] != x[p][j]) {
					compflag = 0;	// they differ
					break;
				}
			}
			
			if (compflag == 1) { 
				unflag = 0;	// not unique, just found it in the unique points table
				break;
			}
		}
		if (unflag) {	// unique, put it in the table
			for (p = 0; p < data.Nth; p++) {
				x[p][un] = xi[p];
			}
			un++;
		}
	} // end block

	runinfo.currentuniques[runinfo.Gen] = un; //+ 1;
	runinfo.acceptance[runinfo.Gen] = (1.0*runinfo.currentuniques[runinfo.Gen])/data.Num[runinfo.Gen]; // check this

	double meanx[data.Nth], stdx[data.Nth];
	for (p = 0; p < data.Nth; p++) {
		meanx[p] = compute_mean(x[p], un);
		stdx[p] = compute_std(x[p], un, meanx[p]);
	}

	printf("CURGEN DB (UNIQUE) %d: [un = %d]\n", runinfo.Gen, un); // + 1);
	print_matrix("means", meanx, data.Nth);
	print_matrix("std", stdx, data.Nth);
	} // end block

	for (i = 0; i < n; i++) fj[i] = curgen_db.entry[i].F;	// separate point from F ?
	calculate_statistics(fj, n, data.Num[runinfo.Gen], runinfo.Gen, sel);

	newchains = 0;
	for (i = 0; i < n; i++) {
		if (sel[i] != 0) newchains++;
	}

	sort_t list[n];
	for (i = 0; i < n; i++) {
		list[i].idx = i;
		list[i].nsel = sel[i];
		list[i].F = curgen_db.entry[i].F;
	}

#if VERBOSE
	printf("Points before\n");
	for (i = 0; i < n; i++) {
		printf("%d: %d %d %f\n", i, list[i].idx, list[i].nsel, list[i].F);
	}
#endif

	qsort(list, n, sizeof(sort_t), compar_desc);

#if VERBOSE
	printf("Points after\n");
	for (i = 0; i < n; i++) {
		printf("%d: %d %d %f\n", i, list[i].idx, list[i].nsel, list[i].F);
	}
#endif
	int ldi;	// leader index
	ldi = 0;
	for (i = 0; i < n; i++) {	// newleader
		if (list[i].nsel != 0) {
			int idx = list[i].idx;
			for (p = 0; p < data.Nth; p++) {
				leaders[ldi].point[p] = curgen_db.entry[idx].point[p];
			}
			leaders[ldi].F = curgen_db.entry[idx].F;
			leaders[ldi].nsel = list[i].nsel;

			if (curgen_db.entry[idx].grad != NULL)
				for (p = 0; p < data.Nth; p++) {
					leaders[ldi].grad[p] = curgen_db.entry[idx].grad[p];
				}

			if (curgen_db.entry[idx].hes != NULL)
				for (p = 0; p < data.Nth*data.Nth; p++) {
					leaders[ldi].hes[p] = curgen_db.entry[idx].hes[p];
				}

			ldi++;
		}
	}

	for (i = 0; i < newchains; i++) leaders[i].queue = -1;	// rr

#if VERBOSE
	printf("Leaders before\n");
	for (i = 0; i < newchains; i++) {
		printf("%d %d %f %d\n", i, leaders[i].nsel, leaders[i].F, leaders[i].queue);
	}
#endif

	/* cool and greedy partitioning ala Panos-- ;-) */

	int nworkers = torc_num_workers();
	int *workload = calloc(1, nworkers*sizeof(int));	// workload[1..workers] = 0

	for (i = 0; i < newchains; i++) {
		int least_loader_worker = compute_min_idx_i(workload, nworkers);
		leaders[i].queue = least_loader_worker;
		workload[least_loader_worker] += leaders[i].nsel;
	}

	print_matrix_i("initial workload", workload, nworkers);
	free(workload);

#if VERBOSE
	printf("Leaders after\n");
	for (i = 0; i < newchains; i++) {
		printf("%d %d %f %d\n", i, leaders[i].nsel, leaders[i].F, leaders[i].queue);
	}
#endif


	{//start block
//	double x[data.Nth][n];
	double **x = g_x;
	for (i = 0; i < newchains; i++) {
		for (p = 0; p < data.Nth; p++) {
			x[p][i] = leaders[i].point[p];
		}
	}
	
	double meanx[data.Nth], stdx[data.Nth];
	for (p = 0; p < data.Nth; p++) {
		meanx[p] = compute_mean(x[p], newchains);
		stdx[p] = compute_std(x[p], newchains, meanx[p]);
	}

	printf("CURGEN DB (LEADER) %d: [nlead=%d]\n", runinfo.Gen, newchains);
	print_matrix("means", meanx, data.Nth);
	print_matrix("std", stdx, data.Nth);
	}//end block

	curgen_db.entries = 0;	// reset curgen db
	printf("calculate_statistics: newchains=%d\n", newchains);

	for (i = 0; i < data.Nth; i++) free(g_x[i]);
	free(g_x);

	return newchains;
}


void call_gsl_rand_init()
{
	printf("CALLING gsl_rand_init() on node %d\n", torc_node_id()); fflush(0);
	gsl_rand_init(data.seed);
}

void spmd_gsl_rand_init()
{
	int i;
	for (i = 0; i < torc_num_nodes(); i++) {
		torc_create_ex(i*torc_i_num_workers(), 1, call_gsl_rand_init, 0);
	}
	torc_waitall();
}

void call_print_matrix_2d()
{
	printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
	print_matrix_2d("runinfo.SS", runinfo.SS, data.Nth, data.Nth);
	printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
}

void spmd_print_matrix_2d()
{
	int i;
	for (i = 0; i < torc_num_nodes(); i++) {
		torc_create_ex(i*torc_i_num_workers(), 1, call_print_matrix_2d, 0);
	}
	torc_waitall();
}

void call_update_gdata()	// step for p[j]
{
	MPI_Bcast(runinfo.SS[0], data.Nth*data.Nth, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(runinfo.p, data.MaxStages, MPI_DOUBLE, 0, MPI_COMM_WORLD);	// just p[Gen]
	MPI_Bcast(&runinfo.Gen, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void spmd_update_gdata()	// step
{
	int i;
	if (torc_num_nodes() == 1) return;
	for (i = 0; i < torc_num_nodes(); i++) {
		torc_create_ex(i*torc_i_num_workers(), 1, call_update_gdata, 0);
	}
	torc_waitall();
}

#if defined(_AFFINITY_)	/* BRUTUS */
#include "affinity.c"

void call_setaffinity()
{
	int rank = torc_node_id();
        int numanodesize = 6;
        int start = rank*6 % 48;
        set_rankaff(start, numanodesize);

        get_rankaff(rank);
}

void spmd_setaffinity()
{
	int i;
	for (i = 0; i < torc_num_nodes(); i++) {
		torc_create_ex(i*torc_i_num_workers(), 1, call_setaffinity, 0);
	}
	torc_waitall();
}
#endif

int main(int argc, char *argv[])
{
	int i;
	double t0, gt0, gt1;
	int winfo[4];
	int nchains; // was below

	c_pndl_init();

	torc_register_task(initchaintask);
	torc_register_task(chaintask);
	torc_register_task(torc_update_full_db_task);
	torc_register_task(torc_update_curgen_db_task);
	torc_register_task(torc_update_curres_db_task);
	torc_register_task(reset_nfc_task);
	torc_register_task(get_nfc_task);
	torc_register_task(taskfun);
	torc_register_task(call_gsl_rand_init);
	torc_register_task(call_print_matrix_2d);
	torc_register_task(call_update_gdata);
#if defined(_AFFINITY_)
	torc_register_task(call_setaffinity);
#endif

	data_init();
//	gsl_rand_init();
	setup_handler();

	torc_init(argc, argv, MODE_MS);

#if defined(_AFFINITY_)
	spmd_setaffinity();
#endif

	spmd_gsl_rand_init();

	curgen_db.entries = 0; // peh+

#if defined(_RESTART_)
	int res;
	res = load_runinfo();
	if (res == 0)
	{
		load_curgen_db(runinfo.Gen);
		nchains = data.Num[0];
		printf("nchains = %d\n", nchains);
		gt0 = t0 = torc_gettime();
		goto next;
	}
#endif
	gt0 = t0 = torc_gettime();

	nchains = data.Num[0];
	double out_tparam[data.PopSize];	// nchains
	for (i = 0; i < nchains; i++) {
		winfo[0] = runinfo.Gen;
		winfo[1] = i;
		winfo[2] = -1;
		winfo[3] = -1;

		int d;
		double in_tparam[data.Nth];
		for (d = 0; d < data.Nth; d++) {
			in_tparam[d] = uniformrand(0,1);
			in_tparam[d] *= (data.upperbound[d]-data.lowerbound[d]);
			in_tparam[d] += data.lowerbound[d];
		}

		torc_create(-1, initchaintask, 4,
			data.Nth, MPI_DOUBLE, CALL_BY_COP,
			1, MPI_INT, CALL_BY_COP,
			1, MPI_DOUBLE, CALL_BY_RES,
			4, MPI_INT, CALL_BY_COP,
			in_tparam, &data.Nth, &out_tparam[i], winfo);
	}
#ifdef _STEALING_
	torc_enable_stealing();
#endif
	torc_waitall();
#ifdef _STEALING_
	torc_disable_stealing();
#endif

	gt1 = torc_gettime();
	int g_nfeval = get_nfc();
	printf("server: Generation %d: total elapsed time = %lf secs, generation elapsed time = %lf secs for function calls = %d\n", runinfo.Gen, gt1-t0, gt1-gt0, g_nfeval);
	reset_nfc();

//	print_full_db();
	print_curgen_db();
	dump_curgen_db(runinfo.Gen);
	display_curgen_db(runinfo.Gen);
//	dump_curres_db(runinfo.Gen);

	// save here
#if defined(_RESTART_)
	save_runinfo();
	check_for_exit();
#endif

#if defined(_RESTART_)
next:
	;
#endif
	static cgdbp_t *leaders; //[MAXCHAINS];
	leaders = calloc(1, data.PopSize*sizeof(cgdbp_t));
	for (i = 0; i < data.PopSize; i++) {
		leaders[i].point = calloc(1, data.Nth*sizeof(double));
		leaders[i].grad = calloc(1, data.Nth*sizeof(double));
		leaders[i].hes = calloc(1, data.Nth*data.Nth*sizeof(double));
	}

	curres_db.entries = 0;
	nchains = prepare_newgen(nchains, leaders);	// calculate statistics 

	spmd_update_gdata();
//	spmd_print_matrix_2d();
	call_print_matrix_2d();

	/* this can be moved above */
	if (runinfo.p[runinfo.Gen] == 1) {
		printf("p == 1 from previous run, nothing more to do\n");
		goto end;
	}

#if 0
	for (runinfo.Gen = 1; runinfo.Gen < data.MaxStages; runinfo.Gen++){
#else
	runinfo.Gen++;
	for (/*runinfo.Gen = 1*/; runinfo.Gen < data.MaxStages; runinfo.Gen++){
#endif
		/* process current generation, compute probs, find new chains */
		//leader[i]: { point[data.Nth], F, nsteps}

		int winfo[4];
		double in_tparam[data.Nth];
		int nsteps;
		gt0 = torc_gettime();
		
		/* compute gradient and hessian of each leader (todo: in parallel) */
		if (runinfo.Gen == 1) {
			for (i = 0; i < nchains; i++) {
				compute_gradient_and_hessian(leaders[i].point, leaders[i].grad, leaders[i].hes);
			}
		}

		for (i = 0; i < nchains; i++) {
			winfo[0] = runinfo.Gen;
			winfo[1] = i;
			winfo[2] = -1;	// not used
			winfo[3] = -1;	// not used

			int p;
			for (p = 0; p < data.Nth; p++)
				in_tparam[p] = leaders[i].point[p];
			nsteps = leaders[i].nsel;

			out_tparam[i] = leaders[i].F;	// fleader...

			torc_create(leaders[i].queue, chaintask, 7,
				data.Nth, MPI_DOUBLE, CALL_BY_COP,
				1, MPI_INT, CALL_BY_COP,
				1, MPI_INT, CALL_BY_COP,
				1, MPI_DOUBLE, CALL_BY_REF,
				data.Nth, MPI_DOUBLE, CALL_BY_VAL,
				data.Nth*data.Nth, MPI_DOUBLE, CALL_BY_VAL,
				4, MPI_INT, CALL_BY_COP,
				in_tparam, &data.Nth, &nsteps, &out_tparam[i], leaders[i].grad, leaders[i].hes, winfo);

		}
		/* wait for all chain tasks to finish */
#ifdef _STEALING_
		torc_enable_stealing();
#endif
		torc_waitall();
#ifdef _STEALING_
		torc_disable_stealing();
#endif

		gt1 = torc_gettime();
		int g_nfeval = get_nfc();
		printf("server: Generation %d: total elapsed time = %lf secs, generation elapsed time = %lf secs for function calls = %d\n", runinfo.Gen, gt1-t0, gt1-gt0, g_nfeval);
		reset_nfc();

//		print_full_db();
		print_curgen_db();
		dump_curgen_db(runinfo.Gen);
		display_curgen_db(runinfo.Gen);
//		dump_curres_db(runinfo.Gen);

		printf("newton_aborts = %d\n", newton_aborts);
		printf("not_good = %d\n", not_good);
//		printf("not_pos_def_times = %d\n", not_pos_def_times);

		// save here
#if defined(_RESTART_)
		save_runinfo();
		check_for_exit();
#endif

		curres_db.entries = 0;
		nchains = prepare_newgen(nchains, leaders);	// calculate statistics

		spmd_update_gdata();
		//spmd_print_matrix_2d();
		call_print_matrix_2d();

#if 0
		printf("=================\n");
		print_matrix("runinfo.p", runinfo.p, runinfo.Gen+1);
		print_matrix("runinfo.CoefVar", runinfo.CoefVar, runinfo.Gen+1);
		print_matrix_i("runinfo.currentuniques", runinfo.currentuniques, runinfo.Gen+1);
		print_matrix("runinfo.acceptance", runinfo.acceptance, runinfo.Gen+1);
		print_matrix("runinfo.logselection", runinfo.logselection, runinfo.Gen+1);
		printf("=================\n");
#endif

		if (runinfo.p[runinfo.Gen] == 1) {
			break;
		}
		if (runinfo.Gen+1 == data.MaxStages) {
			break;
		}
	}

	print_matrix("runinfo.p", runinfo.p, runinfo.Gen+1);
	print_matrix("runinfo.CoefVar", runinfo.CoefVar, runinfo.Gen+1);
	print_matrix_i("runinfo.currentuniques", runinfo.currentuniques, runinfo.Gen+1);
	print_matrix("runinfo.acceptance", runinfo.acceptance, runinfo.Gen+1);
	print_matrix("runinfo.logselection", runinfo.logselection, runinfo.Gen+1);

	double logEvidence[1];
	logEvidence[0] = compute_sum(runinfo.logselection, runinfo.Gen+1);
	print_matrix("logEvidence", logEvidence, 1);

	print_matrix_2d("runinfo.SS", runinfo.SS, data.Nth, data.Nth);

	for (i = 0; i < runinfo.Gen+1; i++) {
		char title[64];
		sprintf(title, "runinfo.meantheta(%d)", i);
		//print_matrix("runinfo.meantheta", runinfo.meantheta[i], data.Nth);
		print_matrix(title, runinfo.meantheta[i], data.Nth);
	}

	// last save here - do we need this? what happens if I restart the program with this saved data
#if defined(_RESTART_)
	save_runinfo();
#endif

end:
	/* terminate spawners */
	
	printf("newton_aborts = %d\n", newton_aborts);
	printf("not_good = %d\n", not_good);
	//printf("not_pos_def_times = %d\n", not_pos_def_times);

	printf("total function calls = %d\n", get_tfc());
	torc_finalize();

	return 0;
}
