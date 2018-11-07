/*
 *  fd_grad.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#define _XOPEN_SOURCE 600
#define _BSD_SOURCE
#include "fitfun.c"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <math.h>
#include <pndlc.h>
#include <torc.h>
void reset_nfc();
void inc_nfc();
int get_nfc();
int get_tfc();
void aux_init();

int romberg(double StepRatio, double der_init[], int N, double rombexpon[], int R, double *der_romb, double *errest);
int derivest(int id, double *der_init, int Nder, double H0, double StepRatio);
void pndlga_task(int *pid, double *X, int *pN, double *XL, double *XU, double *ph, int *pIORD, double *GRAD);

#ifndef MY_GETTIME
#define MY_GETTIME
static double my_gettime()
{
	struct timeval t;
	gettimeofday(&t, NULL);
	return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
}
#endif


struct data_s {
	/* read dimension, xl, xu, theta, iord, h0, stepratio, maxsteps  */
	int Nth;
	int iord;
	double h0;
	double stepratio;
	int maxsteps;
	double lb, ub;
	double *lowerbound;
	double *upperbound;
	double *theta;
} data;

void data_init()
{
	int i;

	data.Nth = 4;

	data.iord = 2;
	data.h0 = 0.1;
	data.stepratio = 1.5;
	data.maxsteps = 100;	

	data.lb = -1e38;
	data.ub = +1e38;

	data.lowerbound = malloc(data.Nth*sizeof(double));
	data.upperbound = malloc(data.Nth*sizeof(double));
	data.theta = malloc(data.Nth*sizeof(double));

	for (i = 0; i < data.Nth; i++) data.lowerbound[i] = data.lb;
	for (i = 0; i < data.Nth; i++) data.upperbound[i] = data.ub;
	for (i = 0; i < data.Nth; i++) data.theta[i] = 1.0;


/**********************************************/
// PARAMS
/**********************************************/

	/* USER-DEFINED VALUES */
	FILE *f = fopen("grad.par", "r");
	if (f == NULL) {
		return;
	}

	/*
	Nth		4
	order		2
	initstep	0.1
	stepratio	1.5
	maxsteps	100
	Bdef		-4	4

	lb		0 0 0 0 
	ub		2 2 2 2
	theta		1 1 1 1 
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
		else if (strstr(line, "order")) {
			sscanf(line, "%*s %d", &data.iord);
		}
		else if (strstr(line, "initstep")) {
			sscanf(line, "%*s %lf", &data.h0);
		}
		else if (strstr(line, "stepratio")) {
			sscanf(line, "%*s %lf", &data.stepratio);
		}
		else if (strstr(line, "maxsteps")) {
			sscanf(line, "%*s %d", &data.maxsteps);
		}
		else if (strstr(line, "Bdef")) {
			sscanf(line, "%*s %lf %lf", &data.lb, &data.ub);
		}
	}

	rewind(f);
	line_no = 0;

	free(data.lowerbound);
	free(data.upperbound);
	free(data.theta);
	data.lowerbound = (double *)malloc(data.Nth*sizeof(double));
	data.upperbound = (double *)malloc(data.Nth*sizeof(double));
	data.theta = (double *)malloc(data.Nth*sizeof(double));

	for (i = 0; i < data.Nth; i++) data.lowerbound[i] = data.lb;
	for (i = 0; i < data.Nth; i++) data.upperbound[i] = data.ub;
	for (i = 0; i < data.Nth; i++) data.theta[i] = 1.0;


	rewind(f);
	while (fgets(line, 256, f)!= NULL) {
		line_no++;

		if ((line[0] == '#')||(strlen(line)==0)) continue;

		if (strstr(line, "theta")) {
			char *p;
			p = strtok(line, " \t\n");
			for (i = 0; i < data.Nth; i++) {
				p = strtok(NULL, " \t\n");
				if (!p) break;
				data.theta[i] = atof(p);
                        }
                }
		else if (strstr(line, "lb")) {
			char *p;
			p = strtok(line, " \t\n");
			for (i = 0; i < data.Nth; i++) {
				p = strtok(NULL, " \t\n");
				if (!p) break;
				data.lowerbound[i] = atof(p);
                        }
                }
		else if (strstr(line, "ub")) {
			char *p;
			p = strtok(line, " \t\n");
			for (i = 0; i < data.Nth; i++) {
				p = strtok(NULL, " \t\n");
				if (!p) break;
				data.upperbound[i] = atof(p);
                        }
                }
        }

	for (i = 0; i < data.Nth; i++) {
                printf("param %d: %f < %f < %f\n", i, data.lowerbound[i], data.theta[i], data.upperbound[i]);
        }

	fclose(f);


}

#define NOISE	0

//Noise=@(x) 0.0*randn(size(x));
void Noise(double *x, double *n, int p, struct drand48_data *dr48_buffer)
{
	int i;
	for (i = 0; i < p; i++) {
		double r;
		drand48_r(dr48_buffer, &r);
#if NOISE
		n[i] = 1e-4*r;
#else
		n[i] = 0.0*r;
#endif
	}
}

static struct drand48_data *dr48_buffer;

void rand_init_task()
{
	int i;
	int node = torc_node_id();
	int workers = torc_i_num_workers(); 

	dr48_buffer = (struct drand48_data *) malloc(workers*sizeof(struct drand48_data));
	for (i = 0; i < workers; i++) {
		srand48_r((node+1)*time(0)+i, &dr48_buffer[i]);
	}
}

void rand_init()
{
	int i;

	for (i = 0; i < torc_num_nodes(); i++) {
		torc_create_ex(i*torc_i_num_workers(), 1, rand_init_task, 0);
	}
	torc_waitall();
}

double F(double *x, int *pn)
{
	int p = *pn;
	double sum = 0.0;
	double n[p];

	Noise(x, n, p, &dr48_buffer[torc_i_worker_id()]);

	inc_nfc();
      	sum = fitfun(x, p, NULL, NULL);
	sum += n[0];

//	usleep(10);

#if VERBOSE
	printf("F(%lf,%lf,%lf,%lf) = %lf\n", x[0], x[1], x[2], x[3], sum);
#endif
	return sum;

}

void run(double *TP, int *pn, double *res)
{
	*res = F(TP, pn);
	return;
}

void do_mytest();

int main(int argc, char *argv[])
{
	data_init();

	torc_register_task(rand_init_task);
	torc_register_task(pndlga_task);
	aux_init();
	c_pndl_init();
	torc_init(argc, argv, MODE_MW);

	rand_init();

	do_mytest();

	torc_finalize();

	return 0;
}


#define MAXGRADS	100

void pndlga_task(int *pid, double *X, int *pN, double *XL, double *XU, double *ph, int *pIORD, double *GRAD)
{
	/*int id = *pid;*/
	int N = *pN;
	double h = *ph;
	double FEPS = 1e-6;
	int IORD = *pIORD;
	int IPRINT = 0;
	int NOC, IERR;

	double UH[data.Nth];
	int i;
	
	for (i = 0; i < data.Nth; i++) UH[i] = h;
	c_pndlga(F,X,&N,XL,XU,UH,&FEPS,&IORD,&IPRINT,GRAD,&NOC,&IERR);
}

void do_mytest()
{
	int i, k;
	double t1, t2;

	double *XL = data.lowerbound;
	double *XU = data.upperbound;
	double *X = data.theta;
	int IORD = data.iord;
	int N = data.Nth;

	double GRADs[MAXGRADS][data.Nth];		// grouped by point  
	double allGRADs[data.Nth][MAXGRADS];	// grouped by direction
	double H[MAXGRADS];


	for (i = 0; i < data.Nth; i++) {
                printf("param %d: %e < %e < %e\n", i, data.lowerbound[i], data.theta[i], data.upperbound[i]);
        }

	reset_nfc();

//	double h0 = 0.1;
//	double step = 1.5;

	t1 = my_gettime();
	for (k = 0; k < MAXGRADS; k++) {
		double h = data.h0*pow(data.stepratio, -k);

		H[k] = h;
		if (h < 1e-6) break;

		// GRADIENT 
		torc_create(-1, pndlga_task, 8,
			1, MPI_INT, CALL_BY_COP,
			N, MPI_DOUBLE, CALL_BY_VAL,
			1, MPI_INT, CALL_BY_COP,
			N, MPI_DOUBLE, CALL_BY_VAL,
			N, MPI_DOUBLE, CALL_BY_VAL,
			1, MPI_DOUBLE, CALL_BY_COP,
			1, MPI_INT, CALL_BY_VAL,
			N, MPI_DOUBLE, CALL_BY_RES,
			//&k, X, &N, XL, XU, &h, &IORD, GRADs[k]);
			&k, X, &N, XL, XU, &H[k], &IORD, GRADs[k]);
	}
	torc_waitall();
	t2 = my_gettime();
	printf("t2-t1 = %lf seconds\n", t2-t1);

	int Nder = k;

#if VERBOSE
	for (k = 0; k < Nder; k++) {
		printf("GRADxS = %.6e\n", H[k]);
		printf("GRADIENT MATRIX :\n");
		for (i = 0; i < data.Nth; i++) {
			printf("GRADx%dV: %15.8lf\n", i, GRADs[k][i]);
		}
		printf("\n"); fflush(0);
	}
#endif

	for (i = 0; i < data.Nth; i++) {
		for (k = 0; k < Nder; k++) allGRADs[i][k] = GRADs[k][i];
		derivest(i, allGRADs[i], Nder, data.h0, data.stepratio);
	}

	get_nfc();
	printf("total function calls = %d\n", get_tfc());

}


#include <stdio.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_math.h>

//#define DBG 1

int	romberg(double StepRatio, double der_init[], int N, double rombexpon[], int R, double *der_romb, double *errest)
{
	/* do romberg extrapolation for each estimate
	 * StepRatio - Ratio decrease in step
	 * der_init - initial derivative estimates
	 * rombexpon - higher order terms to cancel using the romberg step
	 * der_romb - derivative estimate returned
	 * errest - error estimate
	 */

	double der_romb_, errest_;

	if (N-R != 2) {
		printf("Incorrect number of derivative estimates and romberg exponents\n");
		return 1;
	}

	double srinv = 1/StepRatio;
	int nexpon = R; //length(rombexpon);
	int i, j;
	
	//rmat = ones(nexpon+2,nexpon+1);
	gsl_matrix *rmat = gsl_matrix_alloc(nexpon+2, nexpon+1);	// ok, before return
	for (i = 0; i < nexpon+2; i++)
		for (j = 0; j < nexpon+1; j++)
			gsl_matrix_set(rmat, i, j, 1.0);

	/*
	 * switch nexpon
	 * case 0
	 * 	% rmat is simple: ones(2,1)
	 * case 1
	 * 	% only one romberg term
	 * 	rmat(2,2) = srinv^rombexpon;
	 * 	rmat(3,2) = srinv^(2*rombexpon);
	 * case 2
	 * 	% two romberg terms
	 * 	rmat(2,2:3) = srinv.^rombexpon;
	 * 	rmat(3,2:3) = srinv.^(2*rombexpon);
	 * 	rmat(4,2:3) = srinv.^(3*rombexpon);
	 * case 3
	 * 	% three romberg terms
	 * 	rmat(2,2:4) = srinv.^rombexpon;
	 * 	rmat(3,2:4) = srinv.^(2*rombexpon);
	 * 	rmat(4,2:4) = srinv.^(3*rombexpon);
	 * 	rmat(5,2:4) = srinv.^(4*rombexpon);
	 * end
	 */
	if (nexpon == 0) {
		/* nothing to do, rmat = ones(2,1) */
	}
	if (nexpon == 1) {
		gsl_matrix_set(rmat, 1, 1, pow(srinv,   rombexpon[0]));
		gsl_matrix_set(rmat, 2, 1, pow(srinv, 2*rombexpon[0]));
	}
	if (nexpon == 2) {
		gsl_matrix_set(rmat, 1, 1, pow(srinv,   rombexpon[0]));
		gsl_matrix_set(rmat, 1, 2, pow(srinv,   rombexpon[1]));
		gsl_matrix_set(rmat, 2, 1, pow(srinv, 2*rombexpon[0]));
		gsl_matrix_set(rmat, 2, 2, pow(srinv, 2*rombexpon[1]));
		gsl_matrix_set(rmat, 3, 1, pow(srinv, 3*rombexpon[0]));
		gsl_matrix_set(rmat, 3, 2, pow(srinv, 3*rombexpon[1]));
	}
	if (nexpon == 3) {
		gsl_matrix_set(rmat, 1, 1, pow(srinv,   rombexpon[0]));
		gsl_matrix_set(rmat, 1, 2, pow(srinv,   rombexpon[1]));
		gsl_matrix_set(rmat, 1, 3, pow(srinv,   rombexpon[2]));
		gsl_matrix_set(rmat, 2, 1, pow(srinv, 2*rombexpon[0]));
		gsl_matrix_set(rmat, 2, 2, pow(srinv, 2*rombexpon[1]));
		gsl_matrix_set(rmat, 2, 3, pow(srinv, 2*rombexpon[2]));
		gsl_matrix_set(rmat, 3, 1, pow(srinv, 3*rombexpon[0]));
		gsl_matrix_set(rmat, 3, 2, pow(srinv, 3*rombexpon[1]));
		gsl_matrix_set(rmat, 3, 3, pow(srinv, 3*rombexpon[2]));
		gsl_matrix_set(rmat, 4, 1, pow(srinv, 4*rombexpon[0]));
		gsl_matrix_set(rmat, 4, 2, pow(srinv, 4*rombexpon[1]));
		gsl_matrix_set(rmat, 4, 3, pow(srinv, 4*rombexpon[2]));
	}

#if DBG
	printf("rmat = \n"); gsl_matrix_fprintf (stdout, rmat, "%g");
#endif
	
	/*
	 * % qr factorization used for the extrapolation as well
	 * % as the uncertainty estimates
	 * [qromb,rromb] = qr(rmat,0)
	 */
	gsl_matrix *qromb = gsl_matrix_alloc(nexpon+2, nexpon+1);	// ok, before return
	gsl_matrix *rromb = gsl_matrix_alloc(nexpon+1, nexpon+1);	// ok, before return
	{
	gsl_vector *tau = gsl_vector_alloc(nexpon+1);	// ok
	gsl_matrix *rmat_qr = gsl_matrix_alloc(nexpon+2, nexpon+1);	// ok
	gsl_matrix_memcpy(rmat_qr, rmat);
	gsl_linalg_QR_decomp(rmat_qr, tau);
	
	gsl_matrix *Qtmp, *Rtmp;	// Q must be MxM
	Qtmp = gsl_matrix_alloc(nexpon+2, nexpon+2);	// ok
	Rtmp = gsl_matrix_alloc(nexpon+2, nexpon+1);	// ok

	gsl_linalg_QR_unpack (rmat_qr, tau, Qtmp, Rtmp);

	for(i = 0; i < nexpon+2; i++)
		for(j = 0; j < nexpon+1; j++)	// ignore last column 
			gsl_matrix_set(qromb, i, j, gsl_matrix_get(Qtmp, i, j));

	for(i = 0; i < nexpon+1; i++)
		for(j = 0; j < nexpon+1; j++)
			gsl_matrix_set(rromb, i, j, gsl_matrix_get(Rtmp, i, j));

	gsl_matrix_free(rmat_qr);
	gsl_vector_free(tau);
	gsl_matrix_free(Qtmp);
	gsl_matrix_free(Rtmp);
	}

	/*
	rombcoefs = rromb\(qromb'*der_init);
	der_romb = rombcoefs(1,:)';
	*/
	gsl_vector *rombcoefs = gsl_vector_alloc (nexpon+1); // ok, before return
	{
	gsl_vector *qxd = gsl_vector_alloc(nexpon+1); 	// ok
	gsl_vector_view der = gsl_vector_view_array(der_init, N);
	gsl_blas_dgemv (CblasTrans, 1.0, qromb, &der.vector, 0.0, qxd);	// qxd = qromb'*der_init

#if DBG
	printf ("qxd = \n"); gsl_vector_fprintf (stdout, qxd, "%g");
#endif

	gsl_permutation * perm = gsl_permutation_alloc (nexpon+1);	// ok
	int slu;
	gsl_linalg_LU_decomp (rromb, perm, &slu);
     
	gsl_linalg_LU_solve (rromb, perm, qxd, rombcoefs);
	
#if DBG
    printf ("rombcoefs = \n"); gsl_vector_fprintf (stdout, rombcoefs, "%g");
#endif

	der_romb_ = gsl_vector_get(rombcoefs, 0);	/* our result */

	gsl_permutation_free (perm);
	gsl_vector_free (qxd);
	}


	/*
	 * % uncertainty estimate of derivative prediction
	 * s = sqrt(sum((der_init - rmat*rombcoefs).^2,1));
	 */
	double s;
	{
	gsl_vector *rmxro = gsl_vector_alloc(nexpon+2); 	// ok
	gsl_blas_dgemv (CblasNoTrans, 1.0, rmat, rombcoefs, 0.0, rmxro);

#if DBG
	printf ("rmat*rombcoefs = \n"); gsl_vector_fprintf (stdout, rmxro, "%g");
#endif

	for (s = 0.0, i = 0; i < nexpon+2; i++)
		s += pow(der_init[i] - gsl_vector_get(rmxro, i), 2.0);
	s = sqrt(s);

	gsl_vector_free(rmxro);
	}

	/*
	 * % uncertainty estimate of derivative prediction (continued)
	 * rinv = rromb\eye(nexpon+1);
	 * cov1 = sum(rinv.^2,2); % 1 spare dof	
	 * errest = s'*12.7062047361747*sqrt(cov1(1));
	 */
	{
	// rinv: make LU decomposition and invert rromb 
	gsl_permutation * perm = gsl_permutation_alloc (nexpon+1);	// ok
	gsl_matrix * rinv = gsl_matrix_alloc (nexpon+1, nexpon+1);	// ok
	int slu;

	gsl_linalg_LU_decomp (rromb, perm, &slu);
	gsl_linalg_LU_invert (rromb, perm, rinv);

#if DBG
	printf ("rinv = \n"); gsl_matrix_fprintf (stdout, rinv, "%g");
#endif

	//cov1 = sum(rinv.^2,2); % 1 spare dof
	double covl[nexpon+1];
	for (i = 0; i < nexpon+1; i++) {
		covl[i] = 0;
		for (j = 0; j < nexpon+1; j++)
			covl[i] += pow(gsl_matrix_get(rinv, i, j), 2.0);		
	}
	
	// error estimation
	errest_ = s*12.7062047361747*sqrt(covl[0]);
	
	gsl_matrix_free(rinv);
	gsl_permutation_free(perm);
	}


	gsl_matrix_free(rmat);
	gsl_vector_free(rombcoefs);
	gsl_matrix_free(qromb);
	gsl_matrix_free(rromb);

	*der_romb = der_romb_;
	*errest = errest_;
	
	return 0;
}


#if 1
int derivest(int id, double *der_init, int Nder, double H0, double StepRatio)
{
	double der_romb, errest;

	int N = 5;
	int R = N-2;

	double *rombexpon;
	if (R > 0)
		rombexpon = malloc(R*sizeof(double));
	else
		rombexpon = NULL;
	
	int r;
	if (1)	// central
		for (r = 0; r < R; r++) rombexpon[r] = 2*(r+1); 
	else
		for (r = 0; r < R; r++) rombexpon[r] = 2*(r+1)-1; 

	int p ;
	
	double best_errest = 1e10;
	double best_stepsize = 0;
	double best_derest = 1e10;
	
	for (p = 0; p <= Nder-N; p++) 
	{
		double stepsize = H0*pow(StepRatio, -p);

		romberg(StepRatio, &der_init[p], N, rombexpon, R, &der_romb, &errest);
		if (p == 0) {
			best_errest = errest;
			best_stepsize = stepsize;
			best_derest = der_romb;
		}
		else {
			if (errest < best_errest) {
				best_errest = errest;
				best_stepsize = stepsize;
				best_derest = der_romb;
			}
		}

//		printf("der_romb = %.6lf\n", der_romb);
//		printf("errest = %.6lf\n", errest);
#if VERBOSE
		printf("[%3d] [der_romb errest stepsize] = [%20.15lf %20.15lf %20.15lf]\n", p, der_romb, errest, stepsize);
#endif
	}
	
	if (rombexpon != NULL) free(rombexpon);

	printf("%d: best:{errest,stepsize,derest} = {%lf,%lf,%lf}\n", id, best_errest, best_stepsize, best_derest);

//	printf("best_errest = %lf\n", best_errest);
//	printf("best_stepsize = %lf\n", best_stepsize);
//	printf("best_derest = %lf\n", best_derest);

	return 0;
}
#endif
