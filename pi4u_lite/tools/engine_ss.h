/*
 *  engine_ss.h
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

typedef struct data_s {
	int	N_init;		//  = 1100*10;
	int	N_seeds;	// = 100*10;
	int	N_steps;	// = 10;

	int	Nth;	// = 2;

	double	*lowerbound;	//[PROBDIM] = {-6.0, -6.0};
	double	*upperbound;	//[PROBDIM] = {+6.0, +6.0};
	double	lb, ub;		// generic lower and upper bound

	int	NTHRESHOLDS;	//	= 10;
	double	*threshold;	// [NTHRESHOLDS] = {-64.0, -32.0, -16.0, -8.0, -4.0, -2.0, -1.0, -1.0/2.0, -1.0/4.0, -1.0/8.0};

	int	logval;

	double	sigma;
	int	seed;
	int	iplot;


	double	FACTOR;		// for adaptive subset
	int	MAXTHRESHOLDS;	//

} data_t;


typedef struct runinfo_s {
} runinfo_t;

extern data_t data;
extern runinfo_t runinfo;


/*** DATABASES ***/

void db_init();
void add_seed_task(double s[], double *pfs);
void add_seed(double s[], double *pfs);
void add_sample_task(double s[], double *pfs);
void add_sample(double s[], double *pfs);
void permute_seeds(int N);
int compar_desc(const void* p1, const void* p2);
int compar_asc(const void* p1, const void* p2);
void sort_seeds(int N, int desc);
void dump_samples(int step);
void dump_seeds(int step);

void set_nsamples(int val);
void set_nseeds(int val);
int get_nsamples();
int get_nseeds();
double *get_sample(int i);
double *get_sample_f_addr(int i);
double get_sample_f(int i);
double *get_seed(int i);
double get_seed_f(int i);
double *get_seed_f_addr(int i);

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
void shuffle(int *data, int N);
int mvnrnd(double *mean, double *var, double *res, int n);

/*** AUX ***/
void inc_nfc();
void get_nfc_task(int *);
int get_nfc();
void reset_nfc_task();
void reset_nfc();
int get_tfc();
