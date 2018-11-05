/*
 *  ss_db.cpp
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
// DB
/**********************************************/

// Databases for samples and seeds
#define MAXSAMPLES	1000000	// peh:this should be max(data.N_init, data.N_seeds*data.N_steps)

typedef struct ss_db_s
{
	double *entry;
	int ncount;
	pthread_mutex_t db_mutex;	// = PTHREAD_MUTEX_INITIALIZER;
} ss_db_t;


static ss_db_t samples;
static ss_db_t seeds; 

void db_init()
{
	seeds.entry = (double *)calloc(1, MAXSAMPLES*(data.Nth+1)*sizeof(double));	// 2D in 1D due to qsort
	seeds.ncount = 0;
	pthread_mutex_init(&seeds.db_mutex, NULL);

	samples.entry = (double *)calloc(1, MAXSAMPLES*(data.Nth+1)*sizeof(double));	// similar format to seeds
	samples.ncount = 0;
	pthread_mutex_init(&samples.db_mutex, NULL);
}

void set_nsamples(int num)
{
	samples.ncount = num;
}

void set_nseeds(int num)
{
	seeds.ncount = num;
}

int get_nsamples()
{
	return samples.ncount;
}

int get_nseeds()
{
	return seeds.ncount;
}

double *get_sample(int i)
{
	return &samples.entry[i*(data.Nth+1)];
}

double get_sample_f(int i)
{
	return samples.entry[i*(data.Nth+1)+data.Nth];
}

double *get_seed(int i)
{
	return &seeds.entry[i*(data.Nth+1)];
}

double get_seed_f(int i)
{
	return seeds.entry[i*(data.Nth+1)+data.Nth];
}

double *get_sample_f_addr(int i)
{
	return &samples.entry[i*(data.Nth+1)+data.Nth];
}

double *get_seed_f_addr(int i)
{
	return &seeds.entry[i*(data.Nth+1)+data.Nth];
}


void add_seed_task(double s[], double *pfs)
{
	double fs = *pfs;
	int l_nseed;

	pthread_mutex_lock(&seeds.db_mutex);
	l_nseed = seeds.ncount++;
	pthread_mutex_unlock(&seeds.db_mutex);

	int offset = (data.Nth+1);
	memcpy(&seeds.entry[l_nseed*offset], s, data.Nth*sizeof(double));
	seeds.entry[l_nseed*offset+data.Nth] = fs;
}

void add_seed(double s[], double *pfs)
{
	if (torc_node_id() == 0)
		add_seed_task(s, pfs);
	else {
		torc_create_direct(0, (void *)add_seed_task, 2,
					data.Nth, MPI_DOUBLE, CALL_BY_VAL,
					1, MPI_DOUBLE, CALL_BY_VAL,
					s, pfs);
		torc_waitall3();
	}
}

void add_sample_task(double s[], double *pfs)
{
	double fs = *pfs;
	int l_nsample;

	pthread_mutex_lock(&samples.db_mutex);
	l_nsample = samples.ncount++;
	pthread_mutex_unlock(&samples.db_mutex);

	int offset = (data.Nth+1);
	memcpy(&samples.entry[l_nsample*offset], s, data.Nth*sizeof(double));
	samples.entry[l_nsample*offset+data.Nth] = fs;
}

void add_sample(double s[], double *pfs)
{
	if (torc_node_id() == 0)
		add_sample_task(s, pfs);
	else {
		torc_create_direct(0, (void *)add_sample_task, 2,
					data.Nth, MPI_DOUBLE, CALL_BY_VAL,
					1, MPI_DOUBLE, CALL_BY_VAL,
					s, pfs);
		torc_waitall3();
	}
}

// can be called earlier, to filter out the unnecessary seeds (if they are more than N_seeds)
void permute_seeds(int N)
{
	double **tmp_seeds;	// temporary storage for permutation-based rearrangement
	int i;

	if (N == 0) return;

	tmp_seeds = (double **)malloc(MAXSAMPLES*sizeof(double *));
	for (i = 0; i < MAXSAMPLES; i++)
		tmp_seeds[i] = (double *)calloc(1, (data.Nth+1)*sizeof(double));


	int *perm = (int *)malloc(N*sizeof(int));
	shuffle(perm, N);
	for (i = 0; i < N; i++) {
		memcpy(tmp_seeds[i], &seeds.entry[perm[i]*(data.Nth+1)], (data.Nth+1)*sizeof(double));
	}
	free(perm);
	
	for (i = 0; i < N; i++) {
		memcpy(&seeds.entry[i*(data.Nth+1)], tmp_seeds[i], (data.Nth+1)*sizeof(double));
	}

	for (i = 0; i < MAXSAMPLES; i++)
		free(tmp_seeds[i]);
	free(tmp_seeds);
}

int compar_desc(const void* p1, const void* p2)
{
	int dir = +1;	// -1: ascending order, +1: descending order
	double *s1 = (double *) p1;
	double *s2 = (double *) p2;
	
	if (s1[data.Nth] < s2[data.Nth]) return dir;
	if (s1[data.Nth] > s2[data.Nth]) return -dir;
//	if (s1[data.Nth] == s2[data.Nth]) return 0;
	return 0;
}

int compar_asc(const void* p1, const void* p2)
{
	int dir = -1;	// -1: ascending order, +1: descending order
	double *s1 = (double *) p1;
	double *s2 = (double *) p2;
	
	if (s1[data.Nth] < s2[data.Nth]) return dir;
	if (s1[data.Nth] > s2[data.Nth]) return -dir;
//	if (s1[data.Nth] == s2[data.Nth]) return 0;
	return 0;
}

void sort_seeds(int N, int desc)
{
	if (N == 0) return;

	if (desc)
		qsort(seeds.entry, N, (data.Nth+1)*sizeof(double), compar_desc);
	else
		qsort(seeds.entry, N, (data.Nth+1)*sizeof(double), compar_asc);
}

void dump_samples(int step)
{
	char fname[256];
	sprintf(fname, "samples_%03d.txt", step);
	FILE *fp = fopen(fname, "w");
	
	for (int i = 0; i < samples.ncount; i++) {
		for (int k = 0; k < data.Nth; k++) {
			fprintf(fp, "%lf ", samples.entry[i*(data.Nth+1)+k]);
		}
		fprintf(fp, "%lf\n", samples.entry[i*(data.Nth+1)+data.Nth]);
	}
	fclose(fp);
}

void dump_seeds(int step)
{
	char fname[256];
	sprintf(fname, "seeds_%03d.txt", step);
	FILE *fp = fopen(fname, "w");
	
	for (int i = 0; i < seeds.ncount; i++) {
		for (int k = 0; k < data.Nth; k++) {
			fprintf(fp, "%lf ", seeds.entry[i*(data.Nth+1)+k]);
		}
		fprintf(fp, "%lf\n", seeds.entry[i*(data.Nth+1)+data.Nth]);
	}
	fclose(fp);
}


