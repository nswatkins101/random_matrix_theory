#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_block.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>

#include<extern_var.h>
#include<matrix_utils.h>
#include<matrix_calc_utils.h>
#include<array_utils.h>

void gsl_matrix_inverse(gsl_matrix* z, gsl_matrix* m)
{
	int s = 0;
	/* gsl_vector* eval = gsl_vector_alloc(m->size1);
	gsl_matrix* temp_mat = gsl_matrix_alloc(m->size1, m->size2);
	printf("calculating eigenvalues...\n");
	gsl_matrix_memcpy(temp_mat, m);
	gsl_eigen_symm(temp_mat, eval, g_w_symm);
	array_fprintf_delim(stdout, eval->data, eval->size, ", ");
	printf("sizes: %i,%i and %i,%i with permutation of size %i\n", z->size1,z->size2,m->size1,m->size2,g_permutation->size);
	printf("calculating LU decomposition...\n");
	*/
	gsl_linalg_LU_decomp(m, g_permutation, &s);
	//printf("calculating inverse with LU...\n");
	gsl_linalg_LU_invert(m, g_permutation, z);
}

void gsl_matrix_complex_inverse(gsl_matrix_complex* z, gsl_matrix_complex* m)
{
	int s = 0;
	/*gsl_vector* eval = gsl_vector_alloc(m->size1);
	 * gsl_matrix* temp_mat = gsl_matrix_alloc(m->size1, m->size2);
	 * gsl_eigen_symm_workspace* w = gsl_eigen_symmv_alloc(m->size1);
	 * printf("calculating eigenvalues...\n");
	 * gsl_matrix_memcpy(temp_mat, m);
	 * gsl_eigen_symm(temp_mat, eval, w);
	 * print_array(eval->data, eval->size);
	 * printf("calculating LU decomposition...\n");
	*/
	gsl_linalg_complex_LU_decomp(m, g_permutation, &s);
	//printf("calculating inverse with LU...\n");
	gsl_linalg_complex_LU_invert(m, g_permutation, z);
}
	
// take the symmetric part of the matrix
void gsl_matrix_symmetric_part(gsl_matrix* mat)
{
	gsl_matrix* temp_mat = get_extern_temp_mat();
	gsl_matrix_memcpy(temp_mat, mat);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 0.5, temp_mat, identity_mat, 0.5, mat);	
	reset_extern_temp_mat(temp_mat);
}

void gsl_matrix_symmetrize(gsl_matrix* mat)
{
	// assumes a triangular matrix, doesn't matter if upper or lower
	gsl_matrix_symmetric_part(mat);
	// symm part sums the two transpose of the triangular parts and scales the matrix by 1/2, undo that
	gsl_matrix_scale(mat, 2);
}

double gsl_matrix_trace(const gsl_matrix* mat)
{
	if(mat->size1 != mat->size2)
			MY_GSL_ERROR(GSL_ENOTSQR);
	int dim = mat->size1;
	/* there is no gsl trace function */
	double tr = 0;
	for(int i=0; i<dim; i++)
		tr += gsl_matrix_get(mat, i, i);
	return tr;
}

double gsl_vector_l1_norm(const gsl_vector* x)
{
	double l1_norm = 0;
	for(int i=0; i<x->size; i++)
		l1_norm += gsl_vector_get(x, i);
	return l1_norm;
}