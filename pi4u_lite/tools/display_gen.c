/*
 *  display_gen.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "gnuplot_i.h"

int display = 1;

gnuplot_ctrl * g = NULL;

char *BASENAME;
// "curgen"
// "samples"
// "seeds"

void display_onegen_db_dim(int Gen, int Dim)
{
        FILE *fp;
        char fname[256];
        static int once = 0;

	sprintf(fname, "%s_%03d.txt", BASENAME, Gen);
        fp = fopen(fname, "r");
        if (fp == NULL) {
                printf("No file %s\n", fname);
		return;
        }
	fclose(fp);

//	if (g != NULL) {
//		gnuplot_close(g);
//	}
//	g = gnuplot_init();
	if (!once) {
		gnuplot_cmd(g, "set view map");
		gnuplot_cmd(g, "set size ratio 1");
		gnuplot_cmd(g, "set palette rgbformulae 22,13,-31");
		gnuplot_cmd(g, "unset key");
		gnuplot_cmd(g, "set pointsize 1");
		once = 1;
	}

//	set terminal x11 [reset] <n> [[no]enhanced] [font <fontspec>] [title "<string>"] [[no]persist] [[no]raise] [close]

	if (!display) {
		gnuplot_cmd(g, "set terminal png");
		gnuplot_cmd(g, "set output \"%s.png\"", fname);
	}

	int i, j;
//	gnuplot_cmd(g, "splot [-10:10][-10:10] \"%s\" using 2:3:4 with points pt 7 palette", fname);
//	gnuplot_cmd(g, "splot \"%s\" using 3:2:4 with points pt 7 palette", fname);

	for (i = 1; i <= Dim; i++) {
		for (j = i+1; j <=Dim; j++) {
			char using_str[64], title_str[64];
			sprintf(using_str, "using %d:%d:%d ", i,j,Dim+1);
			sprintf(title_str, "\"%d_%d_%d\"", Gen,i,j);
			gnuplot_cmd(g, "set term x11 %d persist title %s", i*Dim+j, title_str);
			gnuplot_cmd(g, "splot \"%s\" %s with points pt 7 palette", fname, using_str);
//			gnuplot_cmd(g, "splot [-6:6][-6:6] \"%s\" %s with points pt 7 palette", fname, using_str);
//			usleep(1000*500);
			sleep(2);
		}
	}

//	sleep(5);
}

void display_curgen_db(int Gen, int NGens)
{
        FILE *fp;
        char fname[256];
        static int once = 0;

	sprintf(fname, "%s_%03d.txt", BASENAME, Gen);
        fp = fopen(fname, "r");
        if (fp == NULL) {
                printf("No file %s\n", fname);
		return;
        }
	fclose(fp);

//	if (g != NULL) {
//		gnuplot_close(g);
//	}
//	g = gnuplot_init();
	if (!once) {
		gnuplot_cmd(g, "set view map");
		gnuplot_cmd(g, "set size ratio 1");
		gnuplot_cmd(g, "set palette rgbformulae 22,13,-31");
		gnuplot_cmd(g, "unset key");
		gnuplot_cmd(g, "set pointsize 1");
		once = 1;
	}

//	set terminal x11 [reset] <n> [[no]enhanced] [font <fontspec>] [title "<string>"] [[no]persist] [[no]raise] [close]
/*
	if ((Gen == 1) || (Gen == NGens))
		gnuplot_cmd(g, "set term x11 %d persist", Gen);
	else
		gnuplot_cmd(g, "set term x11 %d", Gen);

	if ((Gen > 2)
		gnuplot_cmd(g, "set term x11 %d close", Gen-1);

	gnuplot_cmd(g, "splot [-10:10][-10:10] \"%s\" with points pt 7 palette", fname);
	usleep(1000*500);
*/

	if (!display) {
		gnuplot_cmd(g, "set terminal png");
		gnuplot_cmd(g, "set output \"%s.png\"", fname);
	}

	if (Gen == 0) {
		if (display) gnuplot_cmd(g, "set term x11 %d persist", Gen);
		gnuplot_cmd(g, "splot [-10:10][-10:10] \"%s\" with points pt 7 palette", fname);
	}
	else if (Gen == 1) {
		//if (display) gnuplot_cmd(g, "set term x11 %d persist", Gen);
		gnuplot_cmd(g, "splot [-10:10][-10:10] \"%s\" with points pt 7 palette", fname);
	} else {
		gnuplot_cmd(g, "splot [-10:10][-10:10] \"%s\" with points pt 7 palette", fname);
	}
	usleep(1000*500);

//	sleep(5);
}

void display_curgen_db_single(int Gen, int p1, int p2, int Dim, int NGens)
{
        FILE *fp;
        char fname[256];
        static int once = 0;

	sprintf(fname, "%s_%03d.txt", BASENAME, Gen);
        fp = fopen(fname, "r");
        if (fp == NULL) {
                printf("No file %s\n", fname);
		return;
        }
	fclose(fp);

	if (!once) {
		gnuplot_cmd(g, "set view map");
		gnuplot_cmd(g, "set size ratio 1");
		gnuplot_cmd(g, "set palette rgbformulae 22,13,-31");
		gnuplot_cmd(g, "unset key");
		gnuplot_cmd(g, "set pointsize 1");
		once = 1;

//		gnuplot_cmd(g, "set logscale y");
	}

	if (!display) {
		gnuplot_cmd(g, "set terminal png");
		gnuplot_cmd(g, "set output \"%s.png\"", fname);
	}

	int i = p1;
	int j = p2;
	{
	char using_str[64], title_str[64];
	sprintf(using_str, "using %d:%d:%d ", i,j,Dim+1);
	sprintf(title_str, "\"%d_%d_%d\"", Gen,i,j);
#if 1
	gnuplot_cmd(g, "set term x11 %d persist title %s", i*Dim+j, title_str);
#else
	gnuplot_cmd(g, "set term x11 %d persist title %s", Gen, title_str);
#endif

//	gnuplot_cmd(g, "splot [1e2:1e6][1e-6:1e-3]\"%s\" %s with points pt 7 palette", fname, using_str);

	gnuplot_cmd(g, "splot \"%s\" %s with points pt 7 palette", fname, using_str);
//	gnuplot_cmd(g, "splot [-6:6][-6:6]\"%s\" %s with points pt 7 palette", fname, using_str);

	}


	sleep(1);
}


int main(int argc, char *argv[])
{
	int ngen = 0;
	int dim = 2;
	int p1 = 1, p2 = 2;

	if (argc == 1) {
		printf("usage: %s <basefilename> <ngen> <dim> [index i] [index j]\n", argv[0]);
		exit(1);
	}

	BASENAME = argv[1];
	if (argc >= 3) ngen = atoi(argv[2]);
	if (argc >= 4) dim = atoi(argv[3]);
	if (argc == 6) {
		p1 = atoi(argv[4]); 
		p2 = atoi(argv[5]); 
	}

	g = gnuplot_init();
#if 1
	int i;

	for (i = 0; i <= ngen; i++) {
//	for (i = ngen; i <= ngen; i++) {
//		display_curgen_db(i, ngen);
		display_curgen_db_single(i, p1, p2, dim, ngen);
		//sleep(5);
	}
#else
	display_onegen_db_dim(ngen, dim);
#endif

	//sleep(100);
	gnuplot_close(g);

	return 0;
}
