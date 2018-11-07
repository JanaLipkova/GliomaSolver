/*
 *  tmcmc_stats.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <math.h>
#include "engine_tmcmc.h"
#include <time.h>

int display = 0;
#define LARGE_SCALE_POPS
/** OBJLOGP FUNCTION **/
double Objlogp(double x, double *fj, int fn, double pj, double tol)
{
	int i;
	double fjmax = compute_max(fj, fn);

#ifdef LARGE_SCALE_POPS
	double *weight = (double *)malloc(fn*sizeof(double));
#else
	double weight[fn];
#endif
	for (i = 0; i < fn; i++)
		weight[i] = exp((fj[i]-fjmax)*(x-pj));

	double sum_weight = compute_sum(weight, fn);

#ifdef LARGE_SCALE_POPS
	double *q = (double *)malloc(fn*sizeof(double));
#else
	double q[fn];
#endif
	
	for (i = 0; i < fn; i++)
		q[i] = weight[i]/sum_weight;

	double mean_q = compute_mean(q, fn);
	double std_q = compute_std(q, fn, mean_q);

	double CoefVar = pow(std_q/mean_q-tol, 2);	/* result */

#ifdef LARGE_SCALE_POPS
	free(weight);
	free(q);
#endif
	return CoefVar;
}

typedef struct fparam_s {
	double *fj;
	int     fn;
	double  pj;
	double  tol;
} fparam_t;

fparam_t *sfp;

double Objlogp_s(double *x, int n)
{
	double *fj = sfp->fj;
	int fn = sfp->fn;
	double pj = sfp->pj;
	double tol = sfp->tol;

	return Objlogp(x[0], fj, fn, pj, tol);
}

double Objlogp_gsl(double x, void *param)
{
	fparam_t *fp = (fparam_t *) param;

	double *fj = fp->fj;
	int fn = fp->fn;
	double pj = fp->pj;
	double tol = fp->tol;

	double res = Objlogp(x, fj, fn, pj, tol);
/*	printf("Objlogp(%lf)=%lf\n", x, res); */
	return res;
}

double Objlogp_gsl2(const gsl_vector *v, void *param)
{
	double x;
	x = gsl_vector_get(v, 0);

	return Objlogp_gsl(x, param);
}

/*** OPTIMIZATION ***/
int fzerofind(double *fj, int fn, double pj, double tol, double *xmin, double *fmin)
{
	size_t iter = 0;
	/*size_t max_iter = data.options.MaxIter;*/	/* USER input - not used here */
	double Tol = data.options.Tol;
	int Display = data.options.Display;
	double Step = 1e-4; //data.options.Step;
	double x_lo = 0.0, x_hi = 4.0;
	int conv = 0;

	size_t niters;

	static int counter = -1;
	int first_try = 0;
	int dump = 0;
	FILE *fp = NULL;
	char fname[64];
	counter++;
retry:
	if (Display) printf("fminzero: x_lo = %e x_hi = %e Step = %e\n", x_lo, x_hi, Step);
	niters = (unsigned long) ((x_hi-x_lo) / Step);

	first_try++;
	if (first_try) dump=1;

	if (dump) {
		sprintf(fname, "fzero_%03d.txt", counter); 
		fp = fopen(fname, "w");
	}


	double m = 0;
	double fm = DBL_MAX;
	double t0 = torc_gettime();
	int found = 0;
#if !defined(_OPENMP)
	for (iter = 0; iter < niters; iter++)
	{
		double x = x_lo + iter*Step;
		double fx = Objlogp(x, fj, fn, pj, tol);
		if (dump) fprintf(fp, "%.16f %.16f\n", x, fx);

		if (fx < fm)
		{
			fm = fx;
			m = x;
		}
		if (fabs(fx) <= Tol) {
			found = 1;
			break;
		} 
	}
#else
	#pragma omp parallel
	{
	double lm = 0;
	double lfm = DBL_MAX; 
	
	#pragma omp for
	for (iter = 0; iter < niters; iter++)
	{
		double x, fx;

		if (found == 0)
		{
			x  = x_lo + iter*Step;
			fx = Objlogp(x, fj, fn, pj, tol);
			if (fx < lfm)
			{
				lfm = fx;
				lm = x;
			}
			if (fabs(fx) <= Tol) {
				found = 1;
				#pragma omp flush(found)
			}
		} /* task cancellation ? */
	}

	#pragma omp critical
	{
		if (lfm < fm)
		{
			fm = lfm;
			m = lm;
		}
	}

	}
#endif
	double t1 = torc_gettime();

	if (found) conv = 1;

	/* If fm is not within Tolerance, we can go back and retry with better refinement (more iterations) */
	if (!found) {
		x_lo = m - 10*Step;
		if (x_lo < 0) x_lo = 0;
		x_hi = m + 10*Step;
		if (x_hi > 4) x_hi = 4;
		Step = 0.1*Step; 
		if (Step < 1e-16) {
			return 1;
		} else {
			if (Display)
				printf("fzerofind (%e): m=%.16f fm=%.16f iter=%ld, time=%lf s\n", Step, m, fm, niters, t1-t0);
			goto retry;
		}
	}

	if (Display)
		printf("fzerofind: m=%.16f fm=%.16f iter=%ld, time=%lf s\n", m, fm, niters, t1-t0);

	*xmin = m;
	*fmin = fm;

	if (dump) fclose(fp);

	return (conv = 1);
}

int fzerofind1(double *fj, int fn, double pj, double tol, double *xmin, double *fmin)
{
	size_t iter = 0;
	/*size_t max_iter = data.options.MaxIter;*/	/* USER input - not used here */
	double Tol = data.options.Tol;
	int Display = data.options.Display;
	double Step = 1e-4; //data.options.Step;
	double x_lo = 0, x_hi = 4.0;
	int conv = 0;

	size_t niters;

retry:
	if (Display) printf("fminzero: Step = %e\n", Step);
	niters = (unsigned long) ((x_hi-x_lo) / Step);

	double m = 0;
	double fm = DBL_MAX;
	double t0 = torc_gettime();
	int found = 0;
#if !defined(_OPENMP)
	for (iter = 0; iter < niters; iter++)
	{
		double x = x_lo + iter*Step;
		double fx = Objlogp(x, fj, fn, pj, tol);
		if (fx < fm)
		{
			fm = fx;
			m = x;
		}
		if (fabs(fx) <= Tol) {
			found = 1;
			break;
		} 
	}
#else
	#pragma omp parallel
	{
	double lm = 0;
	double lfm = DBL_MAX; 
	
	#pragma omp for
	for (iter = 0; iter < niters; iter++)
	{
		double x, fx;

		if (found == 0)
		{
			x  = x_lo + iter*Step;
			fx = Objlogp(x, fj, fn, pj, tol);
			if (fx < lfm)
			{
				lfm = fx;
				lm = x;
			}
			if (fabs(fx) <= Tol) {
				found = 1;
				#pragma omp flush(found)
			}
		} /* task cancellation ? */
	}

	#pragma omp critical
	{
		if (lfm < fm)
		{
			fm = lfm;
			m = lm;
		}
	}

	}
#endif
	double t1 = torc_gettime();

	if (found) conv = 1;

	/* If fm is not within Tolerance, we can go back and retry with better refinement (more iterations) */
	if (!found) {
		Step = 0.1*Step; 
		if (Step < 1e-6) {
			return 1;
		} else {
			if (Display)
				printf("fzerofind (%e): m=%.16f fm=%.16f iter=%ld, time=%lf s\n", Step, m, fm, niters, t1-t0);
			goto retry;
		}
	}

	if (Display)
		printf("fzerofind: m=%.16f fm=%.16f iter=%ld, time=%lf s\n", m, fm, niters, t1-t0);

	*xmin = m;
	*fmin = fm;

	return (conv = 1);
}

int fminsearch(double *fj, int fn, double pj, double tol, double *xmin, double *fmin)
{
	const gsl_multimin_fminimizer_type *T; 
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;
	int conv = 0;

	size_t iter = 0, max_iter = data.options.MaxIter;	/* USER input*/
	double Tol = data.options.Tol;
	int Display = data.options.Display;
	double Step = 1e-4; //data.options.Step;
	int status;
	double size;

	fparam_t fp;
	fp.fj = fj; fp.fn = fn; fp.pj = pj; fp.tol = tol;

	/* Starting point */
	x = gsl_vector_alloc (1);
	gsl_vector_set (x, 0, pj);

	/* Set initial step sizes to Step */
	ss = gsl_vector_alloc (1);
	gsl_vector_set_all (ss, Step); /* input */

       /* Initialize method and iterate */
	minex_func.n = 1;
	minex_func.f = Objlogp_gsl2;
	minex_func.params = &fp;

/*	T = gsl_multimin_fminimizer_nmsimplex;*/
	T = gsl_multimin_fminimizer_nmsimplex2;
/*	T = gsl_multimin_fminimizer_nmsimplex2rand;*/
	s = gsl_multimin_fminimizer_alloc (T, 1);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	if (Display) {
		printf ("using %s method\n", gsl_multimin_fminimizer_name (s));
	}

	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, Tol);

		if (status == GSL_SUCCESS) {
			conv = 1;
			if (Display)
				printf ("converged to minimum at\n");
		}
#if 1
		else if (fabs(s->fval) <= Tol) {
			conv = 1;
			status = GSL_SUCCESS;
			if (Display)
				printf ("found minimum at\n");
		}
#endif
		if (Display)
			printf ("%3ld x =  %.16lf f() = %.16f size = %.16f\n",
				iter, gsl_vector_get (s->x, 0), s->fval, size);

	} while (status == GSL_CONTINUE && iter < max_iter);

#if 1
	/* double-check */
	if ((conv == 1) && (fabs(s->fval) > Tol)) {
		conv = 0;
		if (Display)
			printf ("fminsearch: converged but not found minimum.\n");
	}
#endif

	if (conv) {
		conv = 1;
		*fmin = s->fval;
		*xmin = gsl_vector_get(s->x, 0);
	} else {
		*fmin = 0;
		*xmin = 0.0;
	}

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return conv;
}

int fmincon(double *fj, int fn, double pj, double tol, double *xmin, double *fmin)
{
	int status;
	int iter = 0, max_iter = data.options.MaxIter;	/* USER input*/
	double Tol = data.options.Tol;
	int Display = data.options.Display;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	double x_lo = 0.0, x_hi = 4.0;    /* input */
	double m = 0.5, fm = 0.0;
	gsl_function F;
	int conv = 0;
	gsl_vector *x;
	int i;

	x = gsl_vector_alloc (1);

	fparam_t fp;

	fp.fj = fj; fp.fn = fn; fp.pj = pj; fp.tol = tol;

	F.function = Objlogp_gsl;
	F.params = &fp;

	T = gsl_min_fminimizer_brent; 
/*	T = gsl_min_fminimizer_goldensection;*/
/*	T = gsl_min_fminimizer_quad_golden;*/
	s = gsl_min_fminimizer_alloc (T);

	double f_lo = Objlogp_gsl(x_lo, &fp);
	double f_hi = Objlogp_gsl(x_hi, &fp);
	if (f_lo < f_hi) {
		m = x_lo;
		fm = f_lo;
	} else {
		m = x_hi;
		fm = f_hi;
	}

	for (i = 0; i < max_iter; i++) {
		double x = x_lo + i*(x_hi-x_lo)/max_iter;
		double fx = Objlogp_gsl(x, &fp);
		if (fx < fm) {
			m = x;
			fm = fx;
		}
	}

#if 1
	if (fabs(fm <= Tol)) {
		conv = 1;
		gsl_vector_free(x);
		gsl_min_fminimizer_free (s);
		if (Display)
			printf("fmincon: early return with m = %.16f fm = %.16f\n", m, fm);
		return conv;
	}
#endif

	if ((fm < f_lo) && (fm < f_hi)) {
		if (Display)
			printf("fmincon: initialized with %d tries and m = %f (fm = %f)\n", i, m, fm);
	} else {
		if (Display)
			printf("failed to initialize fmincon (%.16f, %.16f)!\n", f_lo, f_hi);
		return 0;
	}


	gsl_min_fminimizer_set (s, &F, m, x_lo, x_hi);

	if (Display) {
		printf ("using %s method\n", gsl_min_fminimizer_name (s));
		printf ("%5s [%18s, %18s] %18s %18s %18s\n", "iter", "lower", "upper", "min", "fmin", "err(est)");
		printf ("%5d [%.16f, %.16f] %.16f %.16f %.16f\n", iter, x_lo, x_hi, m, fm, x_hi - x_lo);
	}

	do {
		iter++;
		status = gsl_min_fminimizer_iterate (s);

		m  = gsl_min_fminimizer_x_minimum (s);
		x_lo = gsl_min_fminimizer_x_lower (s);
		x_hi = gsl_min_fminimizer_x_upper (s);

		status = gsl_min_test_interval (x_lo, x_hi, Tol, Tol);
		if (status == GSL_SUCCESS) {
			if (Display)
				printf ("Converged:\n");
			conv = 1;
		}
#if 1
		else if (fabs(gsl_min_fminimizer_f_minimum(s)) <= Tol) {
			conv = 1;
			status = GSL_SUCCESS;
			if (Display)
				printf ("found minimum at\n");
		}
#endif

		if (Display)
			printf ("%5d [%.16f, %.16f] %.16f f()=%.16f %.16f\n",
				iter, x_lo, x_hi, m, gsl_min_fminimizer_f_minimum(s), x_hi - x_lo);

	} while (status == GSL_CONTINUE && iter < max_iter);

#if 1
	/* double-check */
	if ((conv == 1) && (fabs(gsl_min_fminimizer_f_minimum(s)) > Tol)) {
		conv = 0;
		if (Display)
			printf ("converged but not found minimum.\n");
	}
#endif

	if (conv) {
		conv = 1;
		gsl_vector_set (x, 0, m);
		*fmin = Objlogp_gsl(m, &fp);
		*xmin = m;
	} else {
		*fmin = 0;
		*xmin = 0.0;
	}

	gsl_vector_free(x);
	gsl_min_fminimizer_free (s);

	return conv;
}


/*** STATISTICS ***/
#include "posdef.c" 	/* peh */

void calculate_statistics(double flc[], int n, int nselections, int gen, unsigned int sel[])
{
	int Display = data.options.Display;
	/*double pflag = 0;*/
	double tolCOV = data.TolCOV;
	double *CoefVar = runinfo.CoefVar;
	double *p = runinfo.p;
	int *Num = data.Num;
/*	int *currentuniques = runinfo.currentuniques;*/
	double *logselection = runinfo.logselection;
	double Step = 1e-4; //data.options.Step;
	double fmin = 0, xmin = 0;
	int conv = 0;
#if 1

#if 0
	conv = fmincon(flc, n, p[gen], tolCOV, &xmin, &fmin);
	if (Display)
		printf("fmincon: conv=%d xmin=%.16lf fmin=%.16lf\n", conv, xmin, fmin);
#endif
#if 1
	if (!conv) {
		conv = fminsearch(flc, n, p[gen], tolCOV, &xmin, &fmin);
		if (Display)
			printf("fminsearch: conv=%d xmin=%.16lf fmin=%.16lf\n", conv, xmin, fmin);
	}
#endif
#if 1
	if (!conv) {
		conv = fzerofind(flc, n, p[gen], tolCOV, &xmin, &fmin);
		if (Display)
			printf("fzerofind: conv=%d xmin=%.16lf fmin=%.16lf\n", conv, xmin, fmin);
	}
#endif

#else
	/* testing */
	if (Display) printf("\n");
	double t0, t1;
	t0 = torc_gettime();
	int conv_con = fmincon(flc, n, p[gen], tolCOV, &xmin, &fmin);
	t1 = torc_gettime();
	printf("statopt - fmincon: conv=%d xmin=%.16lf fmin=%.16lf in %f s\n", conv_con, xmin, fmin, t1-t0);
	t0 = torc_gettime();
	int conv_search = fminsearch(flc, n, p[gen], tolCOV, &xmin, &fmin);
	t1 = torc_gettime();
	printf("statopt - fminsearch: conv=%d xmin=%.16lf fmin=%.16lf in %f s\n", conv_search, xmin, fmin, t1-t0);
	t0 = torc_gettime();
	int conv_zero = fzerofind(flc, n, p[gen], tolCOV, &xmin, &fmin);
	t1 = torc_gettime();
	printf("statopt - fzerofind: conv=%d xmin=%.16lf fmin=%.16lf in %f s\n", conv_zero, xmin, fmin, t1-t0);
/*	conv = conv_search; */
	conv = conv_zero;
#endif
	/* gen: next generation number */
	int j = gen+1;

	if ((conv)&&(xmin > p[gen])) {
		p[j] = xmin;
		CoefVar[j] = fmin;
	} else {
		p[j] = p[gen] + 0.1*Step;
		CoefVar[j] = CoefVar[gen];
	}

	if (p[j] > 1) {
		/*pflag=p[j-1];*/
		p[j] = 1;
		Num[j]=data.LastNum;
	}

	/*print_matrix("p", p, j);*/

/*
	if (p[j] > 10) {
		 data.Rsq = data.Rsqstrict;
	}
*/
	/* Compute weights and normalize*/

	int i;
#if 1
	double *flcp = (double *)malloc(n*sizeof(double));
	for (i= 0; i<n; i++)
		flcp[i] = flc[i]*(p[j]-p[j-1]);


	double fjmax= compute_max (flcp,n );
	double *weight = (double *)malloc(n*sizeof(double));
	/*PA weight[i] = exp((flc[i]-fjmax)*(p[j]-p[j-1])); 23/06 */
	for (i = 0; i < n; i++)
		weight[i] = exp( flcp[i] - fjmax );

	if (display)
		print_matrix("weight", weight, n);

	double sum_weight = compute_sum(weight, n);

	double *q = (double *)malloc(n*sizeof(double));
	for (i = 0; i < n; i++)
		q[i] = weight[i]/sum_weight;

	if (display)
		print_matrix("runinfo_q", q, n);

	/*double sum_q = compute_sum(q, n);*/

	/*logselection[gen] = log(sum_weight/currentuniques[gen])+fjmax*(p[gen+1]-p[gen]); PA definition change for all types of resampling 23/06/15*/
	logselection[gen]= log(sum_weight) + fjmax -log(n);
#else
	double fjmax = compute_max(flc, n);

	double weight[n];
	for (i = 0; i < n; i++)
		weight[i] = exp((flc[i]-fjmax)*(p[j]-p[j-1]));

	if (display)
	print_matrix("weight", weight, n);

	double sum_weight = compute_sum(weight, n);

	double q[n];
	for (i = 0; i < n; i++)
		q[i] = weight[i]/sum_weight;

	if (display)
		print_matrix("runinfo_q", q, n);

	/*double sum_q = compute_sum(q, n);*/

	logselection[gen] = log(sum_weight/currentuniques[gen])+fjmax*(p[gen+1]-p[gen]);
#endif

	if (display)
		print_matrix("logselection", logselection, gen+1);

	double mean_q = compute_mean(q, n);
	double std_q = compute_std(q, n, mean_q);

	CoefVar[gen] = std_q/mean_q;

	if (display)
		print_matrix("CoefVar", CoefVar, gen+1);

	size_t K = n;
	unsigned int N = 1;

	unsigned int samples = n; /*1000;*/
	unsigned int *nn = (unsigned int *)malloc(samples*sizeof(unsigned int));

	for (i = 0; i < samples; i++) sel[i] = 0;

	if (nselections == 0) nselections = samples; /* n;*/
#if 0
	int k;
	for (k = 0; k < nselections; k++) {
		/*gsl_ran_multinomial (r, K, N, q, nn);*/
		multinomialrand (K, N, q, nn);
		for (i = 0; i < K; i++) sel[i]+=nn[i];
	}
#else
	N = nselections;
	/*gsl_ran_multinomial (r, K, N, q, nn);*/
	multinomialrand (K, N, q, nn);
	for (i = 0; i < K; i++) sel[i]+=nn[i];
#endif

	if (display) {
		printf("\n s = [");
		for (i = 0; i < K; i++) printf("%d ", sel[i]);
		printf("]\n");
	}

	/* compute SS */
	int PROBDIM = data.Nth;

	double mean_of_theta[PROBDIM];

	for (i = 0; i < PROBDIM; i++) {
		mean_of_theta[i] = 0;
		for (j = 0; j < n; j++) mean_of_theta[i]+=curgen_db.entry[j].point[i]*q[j];

		runinfo.meantheta[gen][i] = mean_of_theta[i];
	}

	if (display)
		print_matrix("mean_of_theta", mean_of_theta, PROBDIM);

	double meanv[PROBDIM];
	for (i = 0; i < PROBDIM; i++) {
		meanv[i] = mean_of_theta[i];
	}

	for (i = 0; i < PROBDIM; i++) {
		for (j = 0; j < PROBDIM; j++) {
			double s;
			int k;
			s = 0;
			for (k = 0; k < n; k++) {
				s += q[k]*(curgen_db.entry[k].point[i]-meanv[i])*(curgen_db.entry[k].point[j]-meanv[j]);
			}
			runinfo.SS[i][j] = runinfo.SS[j][i] = s;
		}
	}

#if 1	/* peh:check this */
	{
	int fixed = make_posdef(runinfo.SS[0], PROBDIM, 2);
	if (fixed) {
		printf("WARNING: runinfo.SS was forced to become positive definite\n");
	}
	}
#endif

	if (display)
		print_matrix_2d("runinfo.SS", runinfo.SS, PROBDIM, PROBDIM);

	free(flcp);
	free(weight);
	free(q);
	free(nn);
}

double priorpdf(double *theta, int n)
{
	/* peh:check this */
#if 1
	/* log-uniform */
	double res = 0;
	/*double res = 1;*/

	int i;
	for (i = 0; i < n; i++) {
		/*res += log(gsl_ran_flat_pdf(theta[i], data.lowerbound[i], data.upperbound[i])); PA EDIT
		res *= gsl_ran_flat_pdf(theta[i], data.lowerbound[i], data.upperbound[i]); */
		res += -log( data.upperbound[i]- data.lowerbound[i]);
	}
	if (res == 0) return 0;
	/*return log(res); PA EDIT*/
	return res;
#else
	/* gaussian */
	double res = logmvnpdf(n, theta, data.prior_mu, data.prior_sigma);
#endif
	return res;
}

double posterior(double *theta, int n, double LH)
{
	double res;
	/*
	double Prior = priorpdf(theta, n);
	printf("Prior = %lf\n", Prior);

	if (Prior > 0)
		res = LH + log(Prior);     // xxx
	else
	res = LH;
	*/

	res = LH;   /* Algorithm fix by PanosA */

	return res;
}
