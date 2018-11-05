/*
 *  propagation_tool.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 29/6/15.
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */

#define _XOPEN_SOURCE 500
#define _BSD_SOURCE

//#include "fitfun.c" 
#include <stdio.h>
#include <stdlib.h> /* free() */
#include <sys/types.h>
#include <string.h>
#include <mpi.h>
#include <torc.h>
#include <unistd.h>
#include <math.h>
#include "fitfun_glioma.c"

#define PROBDIM	11

int main(int argc, char *argv[])
{
	int i;

	torc_register_task(taskfun);

	torc_init(argc, argv, MODE_MS);

	int n = PROBDIM;

	char ifilename[80];
	strcpy(ifilename, "points.txt");

	if (argc > 1) strcpy(ifilename, argv[1]);

	FILE *fp = fopen(ifilename, "r");
	if (fp == NULL) {
		printf("Input file %s does not exist!\n", ifilename);
		torc_finalize();
		return 1;
	}

/*
	int nlines=0;
	char line[256];	
	while (fgets(line, 256, fp)!= NULL) {
		nlines++;
	}
*/

	#define MAXPOINTS 6048

	char line[256];
	double TP[MAXPOINTS][PROBDIM], ref[MAXPOINTS], res[MAXPOINTS];
	int t = 0;
	while (fgets(line, 256, fp)!= NULL)
	{
		sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", &TP[t][0], &TP[t][1], &TP[t][2], &TP[t][3], &TP[t][4], &TP[t][5], &TP[t][6], &TP[t][7], &TP[t][8], &TP[t][9], &TP[t][10], &ref[t]);
		printf("line %d: %s", t, line);
		t++;
		if (t == MAXPOINTS) break;
	}
	fclose(fp);

	for (i = 0; i < t; i++) {
		int info[4];
		info[0] = 0; info[1] = 0; info[2] = 0; info[3] = i;

		torc_create(-1, taskfun, 4,
			PROBDIM, MPI_DOUBLE, CALL_BY_COP,
			1, MPI_INT, CALL_BY_COP,
			1, MPI_DOUBLE, CALL_BY_RES,
			4, MPI_INT, CALL_BY_COP,
			TP[i], &n, &res[i], info);
	}
	torc_waitall();

	for (i = 0; i < t; i++) {
		printf("RESULT %03d: %10.4f %10.4f %10.4e %10.4e %10.4lf\n", i, TP[i][0], TP[i][1], ref[i], res[i], fabs(ref[i]-res[i]));
	}

	torc_finalize();

	return 0;
}

