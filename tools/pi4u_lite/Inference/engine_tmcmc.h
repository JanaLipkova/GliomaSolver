/*
 *  engine_tmcmc.h
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <torc.h>

/*
#ifndef MY_GETTIME
#define MY_GETTIME
#include <sys/time.h>
static double my_gettime()
{
	struct timeval t;
	gettimeofday(&t, NULL);
	return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
}
#endif
*/

#include "gsl_headers.h"

/*
#define DATANUM		(1024*1)
#define MAXCHAINS	(2*DATANUM)
extern int DATANUM;
extern int MAXCHAINS;  // = 2*1024; //DATANUM;

#define MAXGENS	20
#define PROBDIM	2
extern int PROBDIM;
*/

typedef struct data_s {
	int	Nth;		/* = PROBDIM*/
	int	MaxStages;	/* = MAXGENS*/
	int	PopSize;	/* = DATANUM*/

	double	*lowerbound;	/*[PROBDIM];*/
	double	*upperbound;	/*[PROBDIM];*/

#if 1
	double *prior_mu;
	double *prior_sigma;

	int auxil_size;
	double *auxil_data;
#endif

	int MinChainLength, MaxChainLength;

	double lb, ub;		/*generic lower and upper bound*/

	double	TolCOV;
	double	bbeta;
	int	seed;

	struct optim_options {
		int	MaxIter;
		double	Tol;
		int	Display;
	} options;

	int	iplot;

	int	*Num;		/*[MAXGENS];*/
	int	LastNum;

} data_t;


typedef struct runinfo_s {
	int 	Gen;
	double	*CoefVar;		/*[MAXGENS];*/
	double	*p;			/*[MAXGENS];		// cluster-wide*/
	int	*currentuniques;	/*[MAXGENS];*/
	double	*logselection;		/*[MAXGENS];*/
	double	*acceptance;		/*[MAXGENS];*/
	double	**SS;			/*[PROBDIM][PROBDIM];	// cluster-wide*/
	double	**meantheta; 		/*[MAXGENS][PROBDIM]*/
} runinfo_t;

extern data_t data;
extern runinfo_t runinfo;


/*** DATABASES ***/
/*
#define MAX_DB_ENTRIES	(2*DATANUM)
int MAX_DB_ENTRIES = 2*1024;
extern int MAX_DB_ENTRIES;
*/

typedef struct cgdbp_s {
	double *point; /*[PROBDIM];*/
	double F;
	int counter;	/* not used (?)*/
	int nsel;	/* for selection of leaders only*/
	int queue;	/* for submission of leaders only*/
#if defined(_TMCMC_SN_)
	int valid;
	double *grad;
	double *hes;
#endif
} cgdbp_t;

typedef struct cgdb_s {
	cgdbp_t *entry; /*[MAX_DB_ENTRIES];*/
	int entries;
	pthread_mutex_t m;
} cgdb_t;

typedef struct dbp_s {
	double *point; /*[PROBDIM];*/
	double F;
	int nG;
	double G[64];	/* maxG*/
	int surrogate;
} dbp_t;

typedef struct db_s {
	dbp_t *entry; /*[MAX_DB_ENTRIES];*/		/* */
	int entries;
	pthread_mutex_t m;
} db_t;

#define EXPERIMENTAL_RESULTS	0

typedef struct resdbp_s {
	double *point;	/*[EXPERIMENTAL_RESULTS+1]; // +1 for the result (F)*/
	double F;
	int counter;	/* not used (?)*/
	int nsel;	/* for selection of leaders only*/
} resdbp_t;

typedef struct resdb_s {
	resdbp_t *entry; /*[MAX_DB_ENTRIES];*/
	int entries;
	pthread_mutex_t m;
} resdb_t;

extern cgdb_t curgen_db;
extern db_t full_db;
extern resdb_t curres_db;

void update_full_db(double point[], double F, double *G, int n, int surrogate);
void init_full_db();

void update_curgen_db(double point[], double F);
#if defined(_TMCMC_SN_)
void update_curgen_db_der(double point[], double F, double grad[], double hes[]);
#endif
void init_curgen_db();

void update_curres_db(double point[], double F);
void init_curres_db();
void print_full_db();
void print_curgen_db();
void dump_curgen_db(int Gen);
void dump_curres_db(int Gen);
void display_curgen_db(int Gen);
int load_curgen_db(int Gen);

/*** UTILS ***/
double compute_sum(double *x, int n);
double compute_mean(double *x, int n);
double compute_std(double *x, int n, double mean);
double compute_min(double *x, int n);
int compute_min_idx_i(int *v, int n);
double compute_max(double *x, int n);
void print_matrix(char *name, double *x, int n);
void print_matrix_i(char *name, int *x, int n);
void print_matrix_2d(char *name, double **x, int n1, int n2);

/*** RNG ***/
void gsl_rand_init(int seed);
double normalrand(double mu, double var);
double uniformrand(double a, double b);
void multinomialrand(size_t K, unsigned int N, double q[], unsigned int nn[]);
void shuffle(int *perm, int N);
int mvnrnd(double *mean, double *var, double *res, int n);
double mvnpdf(int n, double *xv, double *mv, double *vm);
double logmvnpdf(int n, double *xv, double *mv, double *vm);

/*** STATISTICS ***/
void calculate_statistics(double flc[], int n, int nselections, int gen, unsigned int sel[]);

/*** PROBLEM FUNCTIONS ***/
double likelihood(double *x, int N);
double posterior(double *theta, int n, double LH);
double priorpdf(double *theta, int n);

/*** AUX ***/
void inc_nfc();
void get_nfc_task(int *);
int get_nfc();
void reset_nfc_task();
void reset_nfc();
int get_tfc();

/*** POSDEF ***/

void compute_mat_product_vect(double *mat/*2D*/, double vect[], double res_vect[], double coef, int PROBDIM);
double compute_dot_product(double row_vector[], double vector[], int PROBDIM);
int inv_matrix(double coef, double *current_hessian/*2D*/, double *inv_current_hessian/*2D*/, int PROBDIM);

