#include <stdio.h>
#include <assert.h>
#include "gsl_headers.h"

#define POS_DEF_METHOD	4

//~~~~~~~~~~~~~~~~~~~~~~~~MCMC with NEWTON UPDATE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//damianos
//=======================================================//
//====================dot_product========================//
//=======================================================//
/**
 * 
 * input(1): transpose of first vector x^T
 * input(2): second vector (column vector) y
 * 
 * output(1): x^T * y
 * 
 * remarks: computes the dot product of input vectors
 */

double compute_dot_product(double row_vector[], double vector[], int PROBDIM)
{
	int row;
	double sum = 0.0;

	for(row=0; row<PROBDIM; row++){
		sum += row_vector[row] * vector[row];
	}
  
	return sum;
}

//=======================================================//
//=======================================================//

//damianos
//=======================================================//
//==============compute_mat_product_vect=================//
/**
 * 
 * input(1): matrix A
 * input(2): vector x
 * input(3): calculated A * x
 * 
 * output: none
 * 
 * remarks: compute the product of the given matrix with the given vector, which is multiplied by the given coefficient
 *		coef * (A * x)
 */
//=======================================================//
void compute_mat_product_vect(double *mat/*2D*/, double vect[], double res_vect[], double coef, int PROBDIM)
{
        int row, column;
        double current_dot_product;

	for(row=0; row<PROBDIM; row++){
                current_dot_product = 0.0;
                for(column=0; column<PROBDIM; column++) {
                        current_dot_product += mat[row*PROBDIM+column] * vect[column]; //row
                }
                res_vect[row] = coef * current_dot_product;
        }
}


//damianos
//==================================================//
//=================Force_Pos_Def====================//
//==================================================//
/**
 *
 * input(1): non positive definite matrix
 * input(2): resulted forced positive definite matix (output)
 * 
 * remarks: Force matrix to be positive definite following
 *  publication "Efficient stohastic generation of multi-site synthetic precipitation data"
 *  by F.P.Brissette et al.
 * "Modification of a negative eigenvalues to create a positive 
 * definite matrices and approximation for standard errors of correlation estimates by L.R. Schaeffer"
 */

void force_pos_def0(gsl_matrix *non_pos_def_mat, gsl_matrix *forced_pos_def_mat, int PROBDIM)
{
	int row, column;
//	int exist_pos_eig_values;
	int exist_zero_eig_value;
	double current_eig_value, sum_neg_eig_values, min_pos_eig_value;
	double normalization_factor, small_pos_eig_value;

	double eps = 1e-12;	//set a value to place instead of zero or negative eigenvalues
	double zero = 1e-15;	//set a value to look for zero

	gsl_matrix *working_mat, *temp_mat, *eig_vectors, *diag_mat, *intermediate_mat;//, *pos_def_mat;
	gsl_vector *eig_values;
	gsl_eigen_symmv_workspace *work_v;

	working_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);
	temp_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);
	diag_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);
	intermediate_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);

	eig_vectors = gsl_matrix_alloc(PROBDIM, PROBDIM);
	eig_values = gsl_vector_alloc(PROBDIM);
  
	//copy input non positive matrix twice
	gsl_matrix_memcpy(working_mat,non_pos_def_mat);
	gsl_matrix_memcpy(temp_mat, non_pos_def_mat);

	//Get the eigenystem of the input matrix
	//allocate the space for the eigenvalues and copy temp_mat to gsl matrix
	work_v = gsl_eigen_symmv_alloc(PROBDIM);
  
	gsl_eigen_symmv(temp_mat, eig_values, eig_vectors, work_v); //get eigenvalues and eigenvectors
  
	min_pos_eig_value = DBL_MAX;	//set the minimum positive eigenvalue to a big starting value
//	exist_pos_eig_values = 0; //set the flag for existence of positive eigenvalues to FALSE
	sum_neg_eig_values = 0.0;
	exist_zero_eig_value = 0; //set the flag for zero eigenvalue existence to FALSE
	//first_neg_value = 1;
   
	/**
	 * loop through the eigenvalues to find:
	 * the minimum positive value
	 * the sum of negative values
	 */

	if(gsl_vector_isneg(eig_values)){ //if all eigenvalues are negative
		min_pos_eig_value = -gsl_vector_max(eig_values); //if all values are negative assumed that the min positive value is the maximum of the negative values
		for(column=0; column<PROBDIM; column++) {	//sum of negative values equals the sum of all vector elements
			sum_neg_eig_values += gsl_vector_get(eig_values, column);
		}
		//exist_pos_eig_values = 1;
	}
	else {
		for(column=0; column < PROBDIM; column++){
			current_eig_value = gsl_vector_get(eig_values, column);
			if(current_eig_value > zero) { // if the current eigenvalue is positive, try to update the minimum positive eigenvalue
				if(current_eig_value <= min_pos_eig_value){
					min_pos_eig_value = current_eig_value;
				}
				//exist_pos_eig_values = 1; // positive eigenvalue was found, set the flag to TRUE
			}
			else if (current_eig_value < 0.0){
				sum_neg_eig_values += current_eig_value;
			}
			else {
				//handling zero eigenvalue
				printf("Warning: Zero eigenvalue\n");
				exist_zero_eig_value = 1; //as you find zero eigenvalue, set flag to true
			}
		}
	}

	if (sum_neg_eig_values == 0.0 && exist_zero_eig_value){ //let's check if only zero eigenvalues hinder matrix to be positive
		sum_neg_eig_values = 0.01; //set this value for correcting zero eigenvalues

		if (min_pos_eig_value == 1000000.0){
			printf("Warning: Correction of all zero eigenvalues \n");
			min_pos_eig_value = 101.0;
			sum_neg_eig_values = 1.0;
			//exist_pos_eig_values = 1; //modifying flag as we recover the situtation of all zero eigenvalues
		}
	}

	assert(sum_neg_eig_values != 0.0); //test that the sum of negative values is negative(!)
//	assert(exist_pos_eig_values && min_pos_eig_value != 1000000.0); //assert that the starting value for minimum is enough big

	sum_neg_eig_values = 2.0 * sum_neg_eig_values;
	//follow the publication on remarks to force the matrix to be positive definite
	normalization_factor = ( pow((sum_neg_eig_values),2) * 100.0) + 1.0; //compute the squared sum of negative values 2.0 * sum_neg
	//update the negative eigenvalues to create new small positive ones
	for(column=0; column< PROBDIM; column++){
		current_eig_value = gsl_vector_get(eig_values, column);

		if (current_eig_value <= zero){ //<= for zero or negative eigenvalue
			small_pos_eig_value = min_pos_eig_value * ( (pow(sum_neg_eig_values - current_eig_value,2)) / normalization_factor );
	
			while (small_pos_eig_value < eps){ //iterative find the closeste value to given eps
				small_pos_eig_value *= 10.0;
			}
			gsl_vector_set(eig_values, column, small_pos_eig_value); //update the eigenvalue
		}
	}

	//follows a Eigen_Vectors * diag(Eigen_Values) * (Eigen_Vectors)^T
	gsl_matrix_set_identity(diag_mat); //construct a diagonal matrix with eigenvalues in the first diagonal
	//hard copying the eigenvalues to the identity matrix to construct a diag[eig1 eig2 .. eign]
	for(row=0; row<PROBDIM; row++){
		gsl_matrix_set(diag_mat, row, row, gsl_vector_get(eig_values, row));
	}
  
	//now compute intermediate matrix = [v1 v2 .. vn] * diag[eig1 eig2 .. eign]
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eig_vectors, diag_mat, 0.0, intermediate_mat);
  
	//finally compute pos definite matrix = [v1 v2 .. vn]  * diag[eig1 eig2 .. eign] * [v1 v2 .. vn]^T
	gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, intermediate_mat, eig_vectors, 0.0, forced_pos_def_mat); //pos_def_mat

	//free the workspace
	gsl_eigen_symmv_free(work_v);
	//free the matrices
	gsl_matrix_free(working_mat);
	gsl_matrix_free(temp_mat);
	gsl_matrix_free(diag_mat);
	gsl_matrix_free(intermediate_mat);
	gsl_matrix_free(eig_vectors);
	//free the vector
	gsl_vector_free(eig_values);
}


void force_pos_def1(gsl_matrix *non_pos_def_mat, gsl_matrix *forced_pos_def_mat, int PROBDIM)
{
//	[v,d] = eig(a);
//	a_psd = v * diag(max(diag(d), eps))/v;

	int row;
//	double eps = 1e-3;  //set a value to place instead of zero or negative eigenvalues
	double eps = 1e-1;  //set a value to place instead of zero or negative eigenvalues

	gsl_matrix *working_mat, *temp_mat, *eig_vectors, *diag_mat, *intermediate_mat;
	gsl_vector *eig_values;
	gsl_eigen_symmv_workspace *work_v;

	working_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);
	temp_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);
	diag_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);
	intermediate_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);

	eig_vectors = gsl_matrix_alloc(PROBDIM, PROBDIM);
	eig_values = gsl_vector_alloc(PROBDIM);

	//copy input non positive matrix twice
	gsl_matrix_memcpy(working_mat,non_pos_def_mat);
	gsl_matrix_memcpy(temp_mat, non_pos_def_mat);

	//Get the eigenystem of the input matrix
	//allocate the space for the eigenvalues and copy temp_mat to gsl matrix
	work_v = gsl_eigen_symmv_alloc(PROBDIM);

	gsl_eigen_symmv(temp_mat, eig_values, eig_vectors, work_v); //get eigenvalues and eigenvectors

	//follows a Eigen_Vectors * diag(Eigen_Values) * (Eigen_Vectors)^T
	gsl_matrix_set_identity(diag_mat); //construct a diagonal matrix with eigenvalues in the first diagonal
	//hard copying the eigenvalues to the identity matrix to construct a diag[eig1 eig2 .. eign]

	for(row=0; row<PROBDIM; row++){
		double eig_val = gsl_vector_get(eig_values, row);
		if (eig_val < eps) eig_val = eps;
		gsl_matrix_set(diag_mat, row, row, eig_val); 
	}

	//now compute intermediate matrix = [v1 v2 .. vn] * diag[eig1 eig2 .. eign]
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eig_vectors, diag_mat, 0.0, intermediate_mat);

	//finally compute pos definite matrix = [v1 v2 .. vn]  * diag[eig1 eig2 .. eign] * [v1 v2 .. vn]^T
	gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, intermediate_mat, eig_vectors, 0.0, forced_pos_def_mat); //pos_def_mat

	//free the workspace
	gsl_eigen_symmv_free(work_v);
	//free the matrices and vector
	gsl_matrix_free(working_mat);
	gsl_matrix_free(temp_mat);
	gsl_matrix_free(diag_mat);
	gsl_matrix_free(intermediate_mat);
	gsl_matrix_free(eig_vectors);
	gsl_vector_free(eig_values);
}

int check_mat_pos_def(gsl_matrix *mat_to_check, int PROBDIM);	// forward declaration

//  Copyright 2003, Pekka Paalanen <pekka.paalanen@lut.fi>
void force_pos_def2(gsl_matrix *non_pos_def_mat, gsl_matrix *forced_pos_def_mat, int PROBDIM)
{
//	D = size(sigma, 1);
//	fixrate = 0.01;
//	covfixmat = ones(D) + fixrate*eye(D);
//	loops = 0;
//	min_limit = eps*10;

	int i; //, row;
	double fixrate = 0.01;
	int loops = 0;
	double min_limit = 1e-15;
	double convfix = 1 + fixrate;

	gsl_matrix *nsigma;

	nsigma = gsl_matrix_alloc(PROBDIM, PROBDIM);
	gsl_matrix_memcpy(nsigma, non_pos_def_mat);

	double d[PROBDIM];
//	gsl_matrix_fprintf (stdout, nsigma, "%10.3f");

	while (check_mat_pos_def(nsigma, PROBDIM) == 0) {
		loops++;
		printf("loops = %d\n", loops);

		int below_limit = 0;
		for (i = 0; i < PROBDIM; i++) {
			d[i] = gsl_matrix_get(nsigma, i, i);
			if (d[i] <= min_limit) below_limit = 1;
		}

		double maxE_abs = fabs(d[0]);
		double minE = d[0];
		for (i = 1; i < PROBDIM; i++) {
			if (fabs(d[i]) > maxE_abs) maxE_abs = fabs(d[i]);
			if (d[i] < minE) minE = d[i];
		}

		double m = maxE_abs * fixrate;
		double neg = minE;
		double addit;

//		printf("m = %f, neg = %f\n", m, neg);
		if (below_limit) {
			if (neg < 0) {
				addit = (m - neg);
			}
			else {
				if (m < min_limit) m = min_limit;
				addit = m;
			}
			for (i = 0; i < PROBDIM; i++) {
				double newv = gsl_matrix_get(nsigma, i, i);
				newv = newv + addit;
				gsl_matrix_set(nsigma, i, i, newv);
			}
		} else {
			for (i = 0; i < PROBDIM; i++) {
				double newv = gsl_matrix_get(nsigma, i, i);
				newv = newv * convfix;
				gsl_matrix_set(nsigma, i, i, newv);
			}
		}
	}

	gsl_matrix_memcpy(forced_pos_def_mat, nsigma);
	gsl_matrix_free(nsigma);
}


void force_pos_def3(gsl_matrix *non_pos_def_mat, gsl_matrix *forced_pos_def_mat, int PROBDIM)
{
	int row, column;
	for (row = 0; row < PROBDIM ; row++) {
		for (column = 0; column < PROBDIM; column++)
			gsl_matrix_set(forced_pos_def_mat, row, column, 0.0);
	}

	for (row = 0; row < PROBDIM ; row++)
		gsl_matrix_set(forced_pos_def_mat, row, row, 1.0);
}


void force_pos_def(gsl_matrix *non_pos_def_mat, gsl_matrix *forced_pos_def_mat, int method, int PROBDIM)
{
	switch (method) {
	case  0: force_pos_def0(non_pos_def_mat, forced_pos_def_mat, PROBDIM); break;
	case  1: force_pos_def1(non_pos_def_mat, forced_pos_def_mat, PROBDIM); break;
	case  2: force_pos_def2(non_pos_def_mat, forced_pos_def_mat, PROBDIM); break;
	case  3:
	default: force_pos_def3(non_pos_def_mat, forced_pos_def_mat, PROBDIM); break;
	}
}

//==================================================//
//==================================================//

//damianos
//======================================================//
//==============check_mat_pos_def===================//
//======================================================//
/**
 * 
 * input(1): matrix to check if positive definite
 * 
 * output(1): 1 or 0 if matrix is non negative or negative respectively
 * 
 * remarks: check if input matrix is non negative
 */
int check_mat_pos_def(gsl_matrix *mat_to_check, int PROBDIM)
{
	//printf("Check pos def START\n");
	gsl_set_error_handler_off();

	gsl_eigen_symm_workspace *work;
	gsl_vector *eig_values;
	int is_pos_def; //flag to show the matrix is positive definite
  
	//Get the eigenvalues of hessian and check if there are positive
	work = gsl_eigen_symm_alloc(PROBDIM);
	gsl_matrix *mat_copy = gsl_matrix_alloc(PROBDIM, PROBDIM);
	gsl_matrix_memcpy(mat_copy, mat_to_check);

	eig_values = gsl_vector_alloc(PROBDIM);
	gsl_eigen_symm(mat_copy, eig_values, work);//get the eigenvalues, the hessian_mat is going to be destroyed

//	int i;
//	for (i = 0; i < PROBDIM; i++) {
//		printf("eig(%d) = %lf\n", i, gsl_vector_get(eig_values, i));
//	}

	is_pos_def = gsl_vector_ispos(eig_values);
  
	gsl_vector_free(eig_values);
	gsl_eigen_symm_free(work);
	gsl_matrix_free(mat_copy);
  
	return is_pos_def;
}


int eigs(gsl_matrix *mat_to_check, int PROBDIM)
{
	//printf("Check pos def START\n");
	gsl_set_error_handler_off();

	gsl_eigen_symm_workspace *work;
	gsl_vector *eig_values;
	//int is_pos_def; //flag to show the matrix is positive definite
  
	//Get the eigenvalues of hessian and check if there are positive
	work = gsl_eigen_symm_alloc(PROBDIM);
	gsl_matrix *mat_copy = gsl_matrix_alloc(PROBDIM, PROBDIM);
	gsl_matrix_memcpy(mat_copy, mat_to_check);

	eig_values = gsl_vector_alloc(PROBDIM);
	gsl_eigen_symm(mat_copy, eig_values, work);//get the eigenvalues, the hessian_mat is going to be destroyed

	int i;
	for (i = 0; i < PROBDIM; i++) {
		printf("eig(%d) = %lf\n", i, gsl_vector_get(eig_values, i));
	}

	//is_pos_def = gsl_vector_ispos(eig_values);
  
	gsl_vector_free(eig_values);
	gsl_matrix_free(mat_copy);
	gsl_eigen_symm_free(work);
  
	return 0;
}

//==================================================//
//==================================================//


int not_pos_def_times = 0;

//damianos
//==================================================//
//==================inv_matrix======================//
//==================================================//
/**
 * input(1): coefficient to multiply the matrix
 * input(2): hessian matrix
 * input(3): matrix to store the inverse of current hessian
 * output(): none
 * 
 * remarks: compute the (inverse of the NEGATIVE hessian matrix)
 */
int inv_matrix(double coef, double *current_hessian/*2D*/, double *inv_current_hessian/*2D*/, int PROBDIM)
{
//	double current_hessian_copy[PROBDIM][PROBDIM]= {{0.0}};
	int row, column; 
	int is_pos_def; //flag to show if we need to force the hessian matrix to be semi positive definite
	int result = 1;

	gsl_matrix *hessian_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);
	gsl_matrix *pos_def_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);
	gsl_matrix *hessian_mat_cpy = gsl_matrix_alloc(PROBDIM,PROBDIM);

	for(row=0; row<PROBDIM; row++){
		for(column=0; column<PROBDIM; column++){
			gsl_matrix_set(hessian_mat, row, column, current_hessian[row*PROBDIM+column]);
		}
	}
	gsl_set_error_handler_off();

	gsl_matrix_memcpy(pos_def_mat, hessian_mat); //copy the current hessian to pos_def_mat
	gsl_matrix_memcpy(hessian_mat_cpy, hessian_mat);

	is_pos_def = check_mat_pos_def(hessian_mat_cpy, PROBDIM); //check if the hessian is non negative definite, hessian_mat_cpy will be destroyed

#if DBG
	printf("m = [\n");
	for (row = 0; row < PROBDIM ; row++) {
		printf("\t");
		for (column = 0; column < PROBDIM; column++)
			printf ("%g ", gsl_matrix_get(hessian_mat, row, column));
		printf("\n");
	}
	printf("];\n");
#endif

	if(!is_pos_def){//if the matrix is negative, force it to be positive
		not_pos_def_times++;
		//printf("not positive definitive\n");
#if (POS_DEF_METHOD == 4)
		result = 0;
		goto cleanup;
#endif

		force_pos_def(hessian_mat, pos_def_mat, POS_DEF_METHOD, PROBDIM);

#if DBG
		printf("p = [\n");
		for (row = 0; row < PROBDIM ; row++) {
			printf("\t");
			for (column = 0; column < PROBDIM; column++)
				printf ("%g ", gsl_matrix_get(pos_def_mat, row, column));
			printf("\n");
		}
		printf("];\n");
#endif
	}
	else {
		//printf("is positive definitive, p = m\n");
	}

	//as the matrix now is semipositive definite apply cholesky decomposition to invert it
	gsl_linalg_cholesky_decomp(pos_def_mat);
	gsl_linalg_cholesky_invert(pos_def_mat);

#if 0
	printf("h-1 = [\n");
	for (row = 0; row < PROBDIM ; row++) {
		printf("\t");
		for (column = 0; column < PROBDIM; column++)
			printf ("%g ", gsl_matrix_get(pos_def_mat, row, column));
		printf("\n");
	}
	printf("];\n");
#endif

	//pass the result to inv_current_hessian array in "hard way"
	for(row=0; row<PROBDIM; row++){
		for(column=0; column<PROBDIM; column++){
			inv_current_hessian[row*PROBDIM+column] = gsl_matrix_get(pos_def_mat, row, column);
		}
	}

#if (POS_DEF_METHOD == 4)
cleanup:
#endif
	//free used matrices
	gsl_matrix_free(hessian_mat);
	gsl_matrix_free(pos_def_mat);
	gsl_matrix_free(hessian_mat_cpy);

	return result;
}
//==================================================//
//==================================================//

void my_print_matrix_2d(char *msg, double **v, int PROBDIM)
{
	int n1 = PROBDIM;
	int n2 = PROBDIM;
	int i, j;

	printf("\n%s =\n\n", msg);
	for (i = 0; i < n1; i++) {
		for (j = 0; j < n2; j++) {
			printf("   %20.5lf", v[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

//==================================================//
//==================================================//

int make_posdef(double *mat, int dim, int method)
{
	int row, column; 
	int res = 1;

	gsl_matrix *hessian_mat = gsl_matrix_alloc(dim, dim);
	gsl_matrix *pos_def_mat = gsl_matrix_alloc(dim, dim);
	gsl_matrix *hessian_mat_cpy = gsl_matrix_alloc(dim, dim);

	for(row=0; row<dim; row++){
		for(column=0; column<dim; column++) {
			gsl_matrix_set(hessian_mat, row, column, mat[row*dim+column]);
		}
	}
	gsl_set_error_handler_off();

	gsl_matrix_memcpy(pos_def_mat, hessian_mat);    // create a copy
	gsl_matrix_memcpy(hessian_mat_cpy, hessian_mat);// create a copy

	int is_pos_def = check_mat_pos_def(hessian_mat_cpy, dim); //check if the matrix is positive definite

	if(!is_pos_def) {       // if not, apply the fix, store	the result in pos_def_mat
		force_pos_def(hessian_mat, pos_def_mat, method, dim);
		res = 1;
	} else {
		res = 0;
	}

	for(row=0; row<dim; row++){
		for(column=0; column<dim; column++){
			mat[row*dim+column] = gsl_matrix_get(pos_def_mat, row, column);
		}
	}

	gsl_matrix_free(hessian_mat);
	gsl_matrix_free(pos_def_mat);
	gsl_matrix_free(hessian_mat_cpy);

	return res;
}

