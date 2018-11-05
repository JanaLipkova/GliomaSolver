/*
 *  fd_deriv.c
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
#include "posdef.c" 
void reset_nfc();
void inc_nfc();
int get_nfc();
int get_tfc();
void aux_init();

void fprint_matrix_1d(FILE *fp, char *title, double *v, int n);

/*
static double my_gettime()
{
	struct timeval t;
	gettimeofday(&t, NULL);
	return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
}
*/

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
	int reps;
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

void do_mytest();

int main(int argc, char *argv[])
{
	/* read theta, xl, xu, uh, iord, reps */
	data_init();

	torc_register_task(rand_init_task);
	aux_init();
	c_pndl_init();

	torc_init(argc, argv, MODE_MW);

        fprint_matrix_1d(stdout, "theta", data.theta, data.Nth);

	rand_init();

	do_mytest();

	torc_finalize();

	return 0;
}


void do_mytest()  //double X[PROBDIM], int N, int IORD)
{
	double GRAD[data.Nth];	// 1st derivatives
	double sGRAD[data.Nth];	// 1st derivatives
	double HES[data.Nth][data.Nth];	// hessian
	double sHES[data.Nth][data.Nth];	// hessian
	double X[data.Nth];
	int i, j;
	double t1, t2;
	const int reps = 1; //200;

	double *theta = data.theta;
	double *XL = data.lowerbound;
	double *XU = data.upperbound;
	double *UH = data.diffstep;
	int N = data.Nth;
	int IORD = data.iord;

	double FEPS = 1e-6;
	int IPRINT = 0;
	int NOC;
	int IERR;

	for (i = 0; i < data.Nth; i++) X[i] = theta[i];

	reset_nfc();
#if 1	// GRADIENT
	printf("\n================================================\n");
	t1 = torc_gettime();
	{
	for (i = 0; i < data.Nth; i++) sGRAD[i] = 0.0;
	int t;
	for (t = 0; t < reps; t++) {
		if (t > 0) {
			double c = 1e-5;
			double delta;
			if (drand48() < 0.5) delta = -1; else delta = 1;
			for (i = 0; i < data.Nth; i++) X[i] = theta[i] + c*delta;
		}

		c_pndlga(F,X,&N,XL,XU,UH,&FEPS,&IORD,&IPRINT,GRAD,&NOC,&IERR);
		//printf("gradient -> %d\n", t);
		for (i = 0; i < data.Nth; i++) sGRAD[i] += GRAD[i];
	}
	for (i = 0; i < data.Nth; i++) GRAD[i] = sGRAD[i]/reps;
	}
	t2 = torc_gettime();

	printf("t2-t1 = %lf seconds\n", t2-t1);
#if 1 //VERBOSE
	printf("IORD = %d\n", IORD);
	printf("NOC = %d\n", NOC);
	printf("GRADIENT VECTOR :\n");
	for (i = 0; i < N; i++) {
		printf("%12.4lf ", GRAD[i]);
	}
	printf("\n"); fflush(0);
#endif

#endif

#if 1	// HESSIAN WITH FUNCTION CALLS
	if (IORD == 4) IORD = 2;
	printf("\n================================================\n");
	t1 = torc_gettime();
	{
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			sHES[i][j] = 0.0;
	int t;
	for (t = 0; t < reps; t++) {
		if (t > 0) {
			double c = 1e-5;
			double delta;
			if (drand48() < 0.5) delta = -1; else delta = 1;
			for (i = 0; i < data.Nth; i++) X[i] = theta[i] + c*delta;
		}
		c_pndlhfa(F,X,&N,XL,XU,UH,&FEPS,&IORD,&IPRINT,(double *)HES,&N,&NOC,&IERR); 
		//printf("hessian -> %d\n", t);

		for (i = 0; i < N; i++)
			for (j = i+1; j < N; j++)
				HES[j][i] = HES[i][j];

		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++)
				sHES[i][j] += HES[i][j];
	}

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			HES[j][i] = sHES[i][j]/reps;
	}
	t2 = torc_gettime();


	printf("t2-t1 = %lf seconds\n", t2-t1);
#if VERBOSE
	printf("IORD = %d\n", IORD);
	printf("NOC = %d\n", NOC);
	printf("HESSIAN MATRIX :\n");
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++)
			printf("%12.4lf ", HES[i][j]);
		printf("\n");
	}
	printf("\n"); fflush(0);
#endif

#endif

	get_nfc();
	printf("total function calls = %d\n", get_tfc());

	printf("GRADIENT VECTOR :\n");
	for (i = 0; i < N; i++) {
		printf("%15.8lf ", GRAD[i]);
	}
	printf("\n"); fflush(0);
	printf("HESSIAN MATRIX :\n");
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++)
			printf("%15.8lf ", HES[i][j]);
		printf("\n");
	}
	printf("\n"); fflush(0);

        gsl_matrix *hessian_mat = gsl_matrix_alloc(data.Nth, data.Nth);
        for(i=0; i<data.Nth; i++){
                for(j=0; j<data.Nth; j++){
                        gsl_matrix_set(hessian_mat, i, j, HES[i][j]);
                }
        }
	eigs(hessian_mat, data.Nth);

	if (data.posdef == -1) return;	// -1 or 0 to 3

	int res = check_mat_pos_def(hessian_mat, data.Nth);
	if (res == 0) {
		int m;
		//for (m = 0; m <=3; m++) {
		for (m = data.posdef; m <=data.posdef; m++) {
			printf(">>>  METHOD %d <<<<\n", m);
			for(i=0; i<data.Nth; i++){
        	        	for(j=0; j<data.Nth; j++){
                	 	       gsl_matrix_set(hessian_mat, i, j, HES[i][j]);
                		}
        		}
			gsl_matrix *hessian_mat2 = gsl_matrix_alloc(data.Nth, data.Nth);
			force_pos_def(hessian_mat, hessian_mat2, m, data.Nth); 

			//gsl_matrix_fprintf(stdout, hessian_mat2, "%lf"); 
			printf("mat2 = \n");
			for(i=0; i<data.Nth; i++){
        	        	for(j=0; j<data.Nth; j++){
                	 	       printf("%15.8lf ", gsl_matrix_get(hessian_mat2, i, j));
                		}
				printf("\n");
        		}

			eigs(hessian_mat2, data.Nth);
			//int res = check_mat_pos_def(hessian_mat2);
		}
	}
}

