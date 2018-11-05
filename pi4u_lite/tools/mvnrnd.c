#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <time.h>

/*
 *  Multivariate Normal density function and random number generator
 *  Multivariate Student t density function and random number generator
 *  Wishart random number generator
 *  Using GSL -> www.gnu.org/software/gsl
 *
 *  Copyright (C) 2006  Ralph dos Santos Silva
 */
int mvnrnd_silva(const gsl_rng *r, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result)
{
	size_t n = mean->size;
	/* multivariate normal distribution random number generator */
	/*
	*	n	dimension of the random vetor
	*	mean    vector of means of size n
	*	var     variance matrix of dimension n x n
	*	result  output variable with a sigle random vector normal distribution generation
	*/

	unsigned int k;
	gsl_matrix *work = gsl_matrix_alloc(n,n);

	gsl_matrix_memcpy(work,var);
	gsl_linalg_cholesky_decomp(work);

	for(k=0; k<n; k++)
        	gsl_vector_set( result, k, gsl_ran_ugaussian(r) );

	gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
	gsl_vector_add(result,mean);

	gsl_matrix_free(work);

	return 0;
}


/*
 *  @title multivariate normal random variables
 *  @author Carl Boettiger, <cboettig@gmail.com>
 *
 *  Based on the R function rmvnorm, from the mvtnorm package
 *  by Friedrich Leisch and Fabian Scheipl, implemented
 *  using the GSL libraries
 */

int mvnrnd_cboet(gsl_rng * rng, const gsl_vector * mean, gsl_matrix * covar, gsl_vector * ANS)
{
	unsigned int i;
	size_t n = mean->size;

	/* Calculate eigenvalues and eigenvectors of covar matrix */
	gsl_vector *eval = gsl_vector_alloc (n);
	gsl_matrix *evec = gsl_matrix_alloc (n, n);
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);
	gsl_eigen_symmv (covar, eval, evec, w);
	gsl_eigen_symmv_free (w);
//	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);



	/* Setup for: evec * matrix(diag(eval)) * transpose(evec)  */
	gsl_matrix *eval_mx = gsl_matrix_calloc (n, n);
	gsl_matrix * x_M = gsl_matrix_alloc (n,n);
	gsl_matrix * x_M_x = gsl_matrix_alloc (n,n);


	gsl_vector_view diagonal = gsl_matrix_diagonal(eval_mx);
	gsl_vector_memcpy(&diagonal.vector, eval);
	for(i=0;i<n;i++)
	{
		gsl_vector_set( &diagonal.vector, 
						i,  
						sqrt( gsl_vector_get(&diagonal.vector, i) )
					  );
	}



	/* evec * matrix(diag(eval)) * transpose(evec)  */
//	gsl_blas_dsymm (CblasLeft, CblasUpper, 
//					1.0, evec, eval_mx, 0.0, x_M);

	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 
					1.0, evec, eval_mx, 0.0, x_M);
	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 
					1.0, x_M, evec, 0.0, x_M_x);


	gsl_matrix_free(x_M);
	gsl_matrix_free(eval_mx);
	gsl_matrix_free(evec);
	gsl_vector_free(eval);

	gsl_vector * rnorms = gsl_vector_alloc(n);
	for(i=0;i<n;i++)
	{ 
		gsl_vector_set 
			( rnorms, i, 
			  gsl_ran_gaussian_ziggurat(rng, 1)
			);
	}

	gsl_blas_dgemv( CblasTrans, 1.0, x_M_x, rnorms, 0, ANS);
	gsl_vector_add(ANS, mean);
	gsl_matrix_free(x_M_x);
  gsl_vector_free(rnorms);

	return 0;
	/* answer provided through pass by reference */
}


int mvnrnd_gsl(gsl_rng * rng, const gsl_vector * mean, gsl_matrix * covar, gsl_vector * ANS)
{
#if 1
	return mvnrnd_silva(rng, mean, covar, ANS);
#else
	return mvnrnd_cboet(rng, mean, covar, ANS);
#endif
}



//gsl_rng *rng;

int mvnrnd(double *mean, double *sigma, double *out, int N)
{
	int res;

	gsl_vector_view mean_view = gsl_vector_view_array(mean, N);
	gsl_matrix_view sigma_view = gsl_matrix_view_array(sigma, N,N);
	gsl_vector_view out_view = gsl_vector_view_array(out, N);

//	res = mvnrnd_gsl(rng, &mean_view.vector, &sigma_view.matrix, &out_view.vector);
	res = mvnrnd_gsl(r, &mean_view.vector, &sigma_view.matrix, &out_view.vector);

	return res;
}
