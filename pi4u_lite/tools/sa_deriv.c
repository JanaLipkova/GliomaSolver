/*
 *  sa_deriv.c
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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <torc.h>
#include "posdef.c"
#include <unistd.h>
int sqrtm(double **A, double **sA, int n);
double* allocate1D(int n);
double** allocate2D(int n);
void fprint_matrix_1d(FILE *fp, char *title, double *v, int n);
void fprint_matrix_2d(FILE *fp, char *title, double **v, int n1, int n2);
void reset_nfc();
void inc_nfc();
int get_nfc();
int get_tfc();
void aux_init();

struct data_s {
	/* read dimension, xl, xu, theta, iord, h0, stepratio, maxsteps  */
	int Nth;
	int iord;

//	double h0;
//	double stepratio;
//	int maxsteps;

	double lb, ub;
	double *lowerbound;
	double *upperbound;
	double *theta;

	double uh;
	int seed;
	int reps;	// total random vectors generated

	double *diffstep;
	int posdef;

} data;

void data_init()
{
	int i;

	data.Nth = 4;

	data.iord = 2;

	data.lb = -1e38;
	data.ub = +1e38;

	data.lowerbound = malloc(data.Nth*sizeof(double));
	data.upperbound = malloc(data.Nth*sizeof(double));
	data.theta = malloc(data.Nth*sizeof(double));

	for (i = 0; i < data.Nth; i++) data.lowerbound[i] = data.lb;
	for (i = 0; i < data.Nth; i++) data.upperbound[i] = data.ub;
	for (i = 0; i < data.Nth; i++) data.theta[i] = 1.0;

	data.uh = 1e-4;
	data.diffstep = malloc(data.Nth*sizeof(double));
	for (i = 0; i < data.Nth; i++) data.diffstep[i] = data.uh;
	data.posdef = -1;
	data.seed = 280675;

//	data.gH_avg = 5000;	

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
	Bdef		-4	4

	lb		0 0 0 0 
	ub		2 2 2 2
	theta		1 1 1 1

	Hdef		1e-4
	diffstep	1-e4 1e-4 1e-4 1e-4 
	posdef		-1
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
		else if (strstr(line, "Bdef")) {
			sscanf(line, "%*s %lf %lf", &data.lb, &data.ub);
		}
		else if (strstr(line, "Hdef")) {
			sscanf(line, "%*s %lf", &data.uh);
		}
		else if (strstr(line, "posdef")) {
			sscanf(line, "%*s %d", &data.posdef);
		}
		else if (strstr(line, "reps")) {
			sscanf(line, "%*s %d", &data.reps);
		}
		else if (strstr(line, "seed")) {
			sscanf(line, "%*s %d", &data.seed);
		}
	}

	rewind(f);
	line_no = 0;

	free(data.lowerbound);
	free(data.upperbound);
	free(data.theta);
	free(data.diffstep);
	data.lowerbound = (double *)malloc(data.Nth*sizeof(double));
	data.upperbound = (double *)malloc(data.Nth*sizeof(double));
	data.theta = (double *)malloc(data.Nth*sizeof(double));
	data.diffstep = (double *)malloc(data.Nth*sizeof(double));

	for (i = 0; i < data.Nth; i++) data.lowerbound[i] = data.lb;
	for (i = 0; i < data.Nth; i++) data.upperbound[i] = data.ub;
	for (i = 0; i < data.Nth; i++) data.theta[i] = 1.0;
	for (i = 0; i < data.Nth; i++) data.diffstep[i] = data.uh;

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
		else if (strstr(line, "diffstep")) {
			char *p;
			p = strtok(line, " \t\n");
			for (i = 0; i < data.Nth; i++) {
				p = strtok(NULL, " \t\n");
				if (!p) break;
				data.diffstep[i] = atof(p);
                        }
                }
        }

	for (i = 0; i < data.Nth; i++) {
                printf("param %d: %f < %f < %f: %f\n", i, data.lowerbound[i], data.theta[i], data.upperbound[i], data.diffstep[i]);
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

double Loss(double *x, int p, struct drand48_data *dr48_buffer)
{
	double n[p];
	double sum = 0.0;

	Noise(x, n, p, dr48_buffer);

	inc_nfc();
	sum = fitfun(x, p, NULL, NULL);
	sum += n[0];

#if VERBOSE
	printf("loss(%lf,%lf,%lf,%lf) = %lf\n", x[0], x[1], x[2], x[3], sum);
#endif
	return sum;
}


double *ghatinput;	// accumulated gradient vector
double **Hhatinput;	// accumulated hessian matrix
pthread_mutex_t hat_m = PTHREAD_MUTEX_INITIALIZER;	// lock for protecting multithreaded updates
//int gH_avg;		// random vectors generated (per task)

void do_reduce()
{
	if (torc_num_nodes() == 1) return;
#if VERBOSE
	int i;
	for (i = 0; i < data.Nth; i++) 
		printf("ghatinput[%d] = %f\n", i, ghatinput[i]);
#endif

	if (torc_node_id() == 0) {
		MPI_Reduce(MPI_IN_PLACE, ghatinput, data.Nth, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, Hhatinput[0], data.Nth*data.Nth, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#if VERBOSE
		printf("after\n");
		for (i = 0; i < data.Nth; i++) 
			printf("ghatinput[%d] = %f\n", i, ghatinput[i]);
#endif
	}
	else {
		MPI_Reduce(ghatinput, ghatinput, data.Nth, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(Hhatinput[0], Hhatinput[0], data.Nth*data.Nth, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
}

void do_work(int *pid)
{
	int p = data.Nth; 
	double c = data.uh; // same diffstep for every direction //1e-5;
//	int gH_avg = 20;	//%no. of random vectors generated 

	double *theta = data.theta;
	int gH_avg = data.reps / torc_num_workers();

//%%% GENERATION OF AVERAGED GRADIENT AND HESSIAN %%% %%% (NO AVERAGING IF gH_avg=1) %%%
	int id = *pid;

	struct drand48_data dr48_buffer;

	if (data.seed == 0) {
		data.seed = time(0);
	}
//	srand48_r(id, &dr48_buffer);
//	srand48_r((id+1)*time(0), &dr48_buffer);
	srand48_r((id+1)*data.seed, &dr48_buffer);

	int i, j;
	int m;
	for (m=1; m <=gH_avg; m++) {
		double delta[p];

		for (i = 0; i < p; i++) {
			double r;
			drand48_r(&dr48_buffer, &r);
			delta[i] = 2*round(r)-1;
		}

		double thetaplus[p], thetaminus[p];
		for (i = 0; i < p; i++) {
			thetaplus[i] = theta[i] + c*delta[i];		
			thetaminus[i] = theta[i] -c*delta[i];
		}

		double yplus  = Loss(thetaplus, p, &dr48_buffer);
		double yminus = Loss(thetaminus, p, &dr48_buffer);
            
		double ghat[p];
		for (i = 0; i < p; i++) 
			ghat[i]=(yplus-yminus)/(2*c*delta[i]);
            

		// update rank-level information
		pthread_mutex_lock(&hat_m);
		for (i = 0; i < p; i++) {
			ghatinput[i] += ghat[i];
		}
		pthread_mutex_unlock(&hat_m);

		//%%% GENERATE THE HESSIAN UPDATE %%% 
		double deltatilda[p];

		for (i = 0; i < p; i++) {
			double r;
			drand48_r(&dr48_buffer, &r);
			deltatilda[i]=2*round(r)-1;
		}

		double thetaplustilda[p], thetaminustilda[p];
		double thetaplusminustilda[p], thetaminusplustilda[p];

		for (i = 0; i < p; i++) {
			thetaplustilda[i] = thetaplus[i] + c*deltatilda[i];		
			thetaminustilda[i] = thetaminus[i] -c*deltatilda[i];

			thetaplusminustilda[i] = thetaplus[i] - c*deltatilda[i];		
			thetaminusplustilda[i] = thetaminus[i] + c*deltatilda[i];
		}

		double yplustilda=Loss(thetaplustilda, p, &dr48_buffer);
		double yminustilda=Loss(thetaminustilda, p, &dr48_buffer);
		double yplusminustilda=Loss(thetaplusminustilda, p, &dr48_buffer);
		double yminusplustilda=Loss(thetaminusplustilda, p, &dr48_buffer);
            
		double ghatplus[p], ghatminus[p];
		for (i = 0; i < p; i++) {
			ghatplus[i]  = (yplustilda-yplusminustilda)/(2*c*deltatilda[i]); 
			ghatminus[i] = (yminusplustilda-yminustilda)/(2*c*deltatilda[i]);
		}
            
		//%%% STATEMENT PROVIDING AN AVERAGE OF SP GRAD. APPROXS.PER ITERATION %%%
           
		double deltaghat[p];
		for (i = 0; i < p; i++) {
			deltaghat[i]=ghatplus[i]-ghatminus[i];
		}

		double Hhat[p][p];
		double HhatT[p][p];

		//for i=1:p
		//	 Hhat(:,i)=deltaghat(i)./(c*delta);
		//end
		for (j=0; j<p; j++) {
			for (i=0; i<p; i++) {
				Hhat[i][j]=deltaghat[j]/(2*c*delta[i]);
				HhatT[j][i] = Hhat[i][j];
			}
		}

		//Hhat=.5*(Hhat+Hhat');
		for (i=0; i<p; i++) {
			for (j=0; j<p; j++) {
				Hhat[i][j]=0.5*(Hhat[i][j] + HhatT[i][j]);
			}
		}

		// update rank-level information
		pthread_mutex_lock(&hat_m);
		for (i=0; i<p; i++) {
			for (j=0; j<p; j++) {
				Hhatinput[i][j]+=Hhat[i][j];
			}
		}
		pthread_mutex_unlock(&hat_m);


	}
}

int main(int argc, char *argv[])
{
	data_init();

	data.Nth = data.Nth;

	//theta = allocate1D(p);	//theta=zeros(p,1);

	ghatinput = allocate1D(data.Nth);
	Hhatinput = allocate2D(data.Nth);

	int i,j;

	torc_register_task(do_work);
	torc_register_task(do_reduce);
	aux_init();

	torc_init(argc, argv, 0);

	fprint_matrix_1d(stdout, "theta", data.theta, data.Nth);

//%%% Loop for number of hessians %%%

	reset_nfc();
	int ntasks = torc_num_workers();
	for (i = 0; i < ntasks; i++) {
		torc_create(-1, do_work, 1,
				1, MPI_INT, CALL_BY_COP, &i);
	}
	torc_waitall();


	for (i = 0; i < torc_num_nodes(); i++) {
		torc_create_ex(i*torc_i_num_workers(), 1, do_reduce, 0);
	}
	torc_waitall();

	//printf("after reduction\n");

	/* compute the average values */
	int gH_avg = data.reps / torc_num_workers();

	for (i=0; i<data.Nth; i++) {
		for (j=0; j<data.Nth; j++) {
			Hhatinput[i][j] /= (ntasks*gH_avg);
		}
		ghatinput[i] /= (ntasks*gH_avg);
	}

	/* store them */
	double *gbar = allocate1D(data.Nth);
	double **Hbar = allocate2D(data.Nth);

	for (i=0; i<data.Nth; i++) {
		for (j=0; j<data.Nth; j++) {
			Hbar[i][j] = Hhatinput[i][j];
		}
		gbar[i] = ghatinput[i];
	}


	//%%% THETA UPDATE BELOW USES GAUSSIAN ELIMINATION %%% %%% TO AVOID DIRECT COMPUTATION OF HESSIAN INVERSE %%%

	double **Hbarbar = allocate2D(data.Nth); //[PROBDIM][PROBDIM];
	double **sHbarbar = allocate2D(data.Nth);

//	Hbarbar=sqrtm(Hbar*Hbar+0.000001*eye(data.Nth)/1.0); 
	int l;
	for (i=0; i<data.Nth; i++) {
		for (j=0; j<data.Nth; j++) {
			double s=0;
			for (l=0; l<data.Nth; l++) {
				s+= Hbar[i][l]*Hbar[l][j];
			}
			Hbarbar[i][j] = s;
			if (i==j) Hbarbar[i][j] += 0.000001;
		}
	}
	if (sqrtm(Hbarbar, sHbarbar, data.Nth)) {
		return 1;
	}

//	fprint_matrix_2d(stdout, "HESSIAN", sHbarbar, data.Nth, data.Nth);

	FILE *fp = fopen("sHbarbar.txt", "w");
	fprint_matrix_2d(fp, "sHbarbar", sHbarbar, data.Nth, data.Nth);
	fclose(fp);


      	printf("GRADIENT VECTOR :\n");
        for (i = 0; i < data.Nth; i++) {
                printf("%15.8lf ", gbar[i]);
        }
	printf("\n"); fflush(0);
        printf("HESSIAN MATRIX :\n");
        for (i = 0; i < data.Nth; i++) {
                for (j = 0; j < data.Nth; j++)
                        printf("%15.8lf ", sHbarbar[i][j]);
                printf("\n");
        }
        printf("\n"); fflush(0);

	gsl_matrix *hessian_mat = gsl_matrix_alloc(data.Nth, data.Nth);
	for(i=0; i<data.Nth; i++){
		for(j=0; j<data.Nth; j++){
			gsl_matrix_set(hessian_mat, i, j, sHbarbar[i][j]);
		}
	}
	eigs(hessian_mat, data.Nth);

//	print(eig(Hbarbar));

	free(data.theta);
	free(Hbar);
	free(gbar);
	free(ghatinput);
	free(Hhatinput);


	get_nfc();
	printf("total function calls = %d\n", get_tfc());

	torc_finalize();
	return 0;
}


//http://yarchive.net/comp/sqrtm.html
int sqrtm(double **A, double **sA, int n)
{
	int i, j;


#if VERBOSE
	fprint_matrix_2d(stdout, "A", A, n, n);
#endif
	gsl_matrix *mA = gsl_matrix_alloc (n, n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			gsl_matrix_set(mA, i, j, A[i][j]);
		}
	}

#if 0
	printf("m =\n ");
	gsl_matrix_fprintf(stdout, mA, "%.4lf");
#endif


	gsl_vector *Eig = gsl_vector_alloc (n);
	gsl_matrix *V = gsl_matrix_alloc (n, n);

	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);
	gsl_eigen_symmv (mA, Eig, V, w);
	gsl_eigen_symmv_free (w);
//	gsl_eigen_symmv_sort (Ev, V, GSL_EIGEN_SORT_ABS_ASC);

	gsl_matrix *D = gsl_matrix_alloc (n, n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			gsl_matrix_set(D, i, j, 0.0);
		}
	}

	for (i = 0; i < n; i++) {
		double eval_i = gsl_vector_get (Eig, i);
		gsl_matrix_set(D, i, i, sqrt(eval_i));
	}
  
#if VERBOSE
	printf("m =\n ");
	gsl_matrix_fprintf(stdout, mA, "%.4lf");
	printf("V =\n ");
	gsl_matrix_fprintf(stdout, V, "%.4lf");
	printf("D =\n ");
	gsl_matrix_fprintf(stdout, D, "%.4lf");
#endif

	gsl_matrix *S0 = gsl_matrix_alloc (n, n);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
			1.0, V, D,
			0.0, S0);

	gsl_matrix *S = gsl_matrix_alloc (n, n);

	gsl_blas_dgemm (CblasNoTrans, CblasTrans,
			1.0, S0, V,
			0.0, S);

#if VERBOSE
	printf("S =\n ");
	gsl_matrix_fprintf(stdout, S, "%.4lf");
#endif

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			double v = gsl_matrix_get(S, i, j);
			sA[i][j] = v;
		}
	}

	gsl_vector_free (Eig);
	gsl_matrix_free (V);
	gsl_matrix_free (D);
	gsl_matrix_free(S0);
	gsl_matrix_free(S);

	return 0;	// everything ok
}

double** allocate2D(int n)
{
	int i;
#if 0
	double **a = (double **)calloc(1, n * sizeof(double *) + (n * n * sizeof(double)));
	double *mem = (double *)(a + n);
#else
	double *mem = (double *)calloc(1, n*n*sizeof(double));
	double **a =  (double **) malloc(n*sizeof(double *));
#endif
	for (i = 0; i < n; i++) {
		a[i] = mem + (i*n); 
	}

	return a;
}

double* allocate1D(int n)
{
	double *a = (double *)calloc(1, n * sizeof(double));
	return a;
}

