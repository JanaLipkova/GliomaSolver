/*
 *  tmcmc_db.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include "engine_tmcmc.h"

void update_full_db(double point[], double F, double *G, int n, int surrogate)
{
	int PROBDIM = data.Nth;
	int i, pos;

	pthread_mutex_lock(&full_db.m);
	pos = full_db.entries;
	full_db.entries++;
	pthread_mutex_unlock(&full_db.m);

	if (full_db.entry[pos].point == NULL) full_db.entry[pos].point = malloc(data.Nth*sizeof(double));

	for (i = 0; i < PROBDIM; i++) full_db.entry[pos].point[i] = point[i];
	full_db.entry[pos].F = F;
	full_db.entry[pos].nG = n;
	for (i = 0; i < n; i++) full_db.entry[pos].G[i] = G[i];
	full_db.entry[pos].surrogate = surrogate;
}

void init_full_db()
{
	pthread_mutex_init(&full_db.m, NULL);
	full_db.entries = 0;
	full_db.entry = calloc(1, data.PopSize*sizeof(dbp_t));	/* this should be something really large or expandable */
}


void update_curgen_db(double point[], double F)
{
	int PROBDIM = data.Nth;
	int i, pos;

	pthread_mutex_lock(&curgen_db.m);
	pos = curgen_db.entries;
	curgen_db.entries++;
	pthread_mutex_unlock(&curgen_db.m);

	if (curgen_db.entry[pos].point == NULL) curgen_db.entry[pos].point = malloc(data.Nth*sizeof(double));

	for (i = 0; i < PROBDIM; i++) curgen_db.entry[pos].point[i] = point[i];
	curgen_db.entry[pos].F = F;
}

#if defined(_TMCMC_SN_)
void update_curgen_db_der(double point[], double F, double *grad, double *hes)
{
	int PROBDIM = data.Nth;
	int i, pos;

	pthread_mutex_lock(&curgen_db.m);
	pos = curgen_db.entries;
	curgen_db.entries++;
	pthread_mutex_unlock(&curgen_db.m);

	if (curgen_db.entry[pos].point == NULL) curgen_db.entry[pos].point = malloc(data.Nth*sizeof(double));

	for (i = 0; i < PROBDIM; i++) curgen_db.entry[pos].point[i] = point[i];
	curgen_db.entry[pos].F = F;

	if (curgen_db.entry[pos].grad == NULL) curgen_db.entry[pos].grad = malloc(data.Nth*sizeof(double));
	for (i = 0; i < PROBDIM; i++) curgen_db.entry[pos].grad[i] = grad[i];

	if (curgen_db.entry[pos].hes == NULL) curgen_db.entry[pos].hes = malloc(data.Nth*data.Nth*sizeof(double));
	for (i = 0; i < PROBDIM*PROBDIM; i++) curgen_db.entry[pos].hes[i] = hes[i];
}
#endif

void init_curgen_db()
{
	pthread_mutex_init(&curgen_db.m, NULL);
	curgen_db.entries = 0;
        curgen_db.entry = calloc(1, (data.MinChainLength+1)*data.PopSize*sizeof(cgdbp_t));
}


void update_curres_db(double point[EXPERIMENTAL_RESULTS], double F)
{
	int i, pos;
	
#if (EXPERIMENTAL_RESULTS <=0)
	return; 
#endif
	pthread_mutex_lock(&curres_db.m);
	pos = curres_db.entries;
	curres_db.entries++;
	pthread_mutex_unlock(&curres_db.m);

	if (curres_db.entry[pos].point == NULL) curres_db.entry[pos].point = malloc((EXPERIMENTAL_RESULTS+1)*sizeof(double));

	for (i = 0; i < EXPERIMENTAL_RESULTS; i++) curres_db.entry[pos].point[i] = point[i];
	curres_db.entry[pos].F = F;	
}

void init_curres_db()
{
	pthread_mutex_init(&curres_db.m, NULL);
	curres_db.entries = 0;
	curres_db.entry = calloc(1, data.PopSize*sizeof(resdbp_t));
}


void print_full_db()
{
	int pos, i;

	printf("=======\n");
	printf("FULL_DB\n");
	for (pos = 0; pos < full_db.entries; pos++) {
		printf("ENTRY %d: POINT(%20.16lf,%20.16lf) F=%20.16lf SG=%d\n",
				pos, full_db.entry[pos].point[0], full_db.entry[pos].point[1],	/* extend it*/
				full_db.entry[pos].F, full_db.entry[pos].surrogate);
		printf("\tG=[");
		for (i = 0; i < full_db.entry[pos].nG-1; i++) printf("%20.16lf,", full_db.entry[pos].G[i]);
		printf("%20.16lf]\n", full_db.entry[pos].G[i]);
	}
	printf("=======\n");
}

void print_curgen_db()
{
/*	int pos;*/

	printf("=======\n");
	printf("CURGEN_DB [size=%d]\n", curgen_db.entries);
/*
	for (pos = 0; pos < curgen_db.entries; pos++) {
		printf("CHAIN %d: POINT(%20.16lf,%20.16lf) F=%20.16lf\n",
				pos, curgen_db.entry[pos].point[0], curgen_db.entry[pos].point[1], curgen_db.entry[pos].F);
	}
*/
	printf("=======\n");
}

void dump_curgen_db(int Gen)
{
	int PROBDIM = data.Nth;
	int pos;
	FILE *fp;
	char fname[256];

	sprintf(fname, "curgen_db_%03d.txt", Gen);
	fp = fopen(fname, "w");
	for (pos = 0; pos < curgen_db.entries; pos++) {
		if (PROBDIM == 2)
			fprintf(fp, "%20.16lf %20.16lf %20.16lf\n",
				curgen_db.entry[pos].point[0], curgen_db.entry[pos].point[1], curgen_db.entry[pos].F);
		else if (PROBDIM == 3) 
			fprintf(fp, "%20.16lf %20.16lf %20.16lf %20.16lf\n",
				curgen_db.entry[pos].point[0], curgen_db.entry[pos].point[1], curgen_db.entry[pos].point[2], curgen_db.entry[pos].F);
		else {
			int i;
			
			for (i = 0; i < PROBDIM; i++) {
				fprintf(fp, "%20.16lf ", curgen_db.entry[pos].point[i]);
			}
			fprintf(fp, "%20.16lf\n", curgen_db.entry[pos].F);
		}

	}
	fclose(fp);
}

int load_curgen_db(int Gen)
{
	int PROBDIM = data.Nth;	/* peh: be careful here */
	int pos;
	FILE *fp;
	char fname[256];

	sprintf(fname, "curgen_db_%03d.txt", Gen);
	fp = fopen(fname, "r");
	if (fp == NULL) {
		printf("DB file: %s not found!!!\n", fname);
		exit(1); 
		return 1;
	}

	curgen_db.entries = 0;
	char line[1024];
	while (fgets(line, 1024, fp) != NULL)
		curgen_db.entries++;

	fclose(fp);	/* fseek...*/
	fp = fopen(fname, "r");	

	for (pos = 0; pos < curgen_db.entries; pos++) {
		int i;
		for (i = 0; i < PROBDIM; i++) {
			if (curgen_db.entry[pos].point == NULL) curgen_db.entry[pos].point = malloc(PROBDIM*sizeof(double));
			fscanf(fp, "%lf", &curgen_db.entry[pos].point[i]);
		}
		fscanf(fp, "%lf", &curgen_db.entry[pos].F);
	}
	fclose(fp);

	return 0;
}

void dump_curres_db(int Gen)
{
	int pos;
	FILE *fp;
	char fname[256];

#if (EXPERIMENTAL_RESULTS <=0)
	return;
#endif
	sprintf(fname, "curres_db_%03d.txt", Gen);
	fp = fopen(fname, "w");
	for (pos = 0; pos < curres_db.entries; pos++) {
			int i;
			
			for (i = 0; i < EXPERIMENTAL_RESULTS; i++) {
				fprintf(fp, "%20.16lf ", curres_db.entry[pos].point[i]);
			}
			fprintf(fp, "%20.16lf\n", curres_db.entry[pos].F);
/*			fprintf(fp, "\n");*/
	}
	fclose(fp);
}


#include "gnuplot_i.h"
static gnuplot_ctrl * g = NULL;

void cleanup_gnuplot()
{
	if (g != NULL) {
		gnuplot_close(g);
	}
}

void display_curgen_db(int Gen)
{
	FILE *fp;
	char fname[256];

	if (!data.iplot) return;

	sprintf(fname, "curgen_db_%03d.txt", Gen);
	fp = fopen(fname, "r");
	if (fp == NULL) {
		printf("No file %s\n", fname);
	}
	fclose(fp);

/*	if (g != NULL) {
		gnuplot_close(g);
	}*/
	if (g == NULL) {
		g = gnuplot_init();
		atexit(cleanup_gnuplot);
	}

	gnuplot_cmd(g, "set view map");
	gnuplot_cmd(g, "set size ratio 1");
	gnuplot_cmd(g, "set palette rgbformulae 22,13,-31");
	gnuplot_cmd(g, "unset key");
/*	set terminal x11 [reset] <n> [[no]enhanced] [font <fontspec>] [title "<string>"] [[no]persist] [[no]raise] [close]*/
	gnuplot_cmd(g, "set term x11 %d persist", Gen);
	gnuplot_cmd(g, "set pointsize 1");
/*	gnuplot_cmd(g, "splot [-6:6][-6:6] \"%s\" with points pt 7 palette", fname);*/
	gnuplot_cmd(g, "splot \"%s\" with points pt 7 palette", fname);
/*	sleep(5);*/
}
