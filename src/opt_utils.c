#include<stdio.h>

#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>

#include<extern_var.h>
#include<opt_utils.h>

void gsl_matrix_grad_trace(gsl_matrix* grad, const gsl_matrix* mat)
{
	// gradient of the trace is the identity matrix
	gsl_matrix_set_identity(grad);
}

void gsl_vector_grad_l1_norm(gsl_vector* grad, const gsl_vector* x)
{
	gsl_vector_set_all(grad, 1);
}

// gets passed a symmetric matrix, calculates eigenvalues and K
double gsl_matrix_kemeny(const gsl_matrix* mat)
{
	double K =0;
	gsl_vector* temp_vec = get_extern_temp_vec();
	gsl_matrix* temp_mat = get_extern_temp_mat();
	gsl_matrix_memcpy(temp_mat, mat);
	gsl_eigen_symm(temp_mat, temp_vec, g_w_symm);	
	for(int i=0; i<mat->size1; i++)
	{
		if(gsl_vector_get(temp_vec, i) != 0)
			K += 1.0/gsl_vector_get(temp_vec, i);
	}
	reset_extern_temp_vec(temp_vec);
	reset_extern_temp_mat(temp_mat);
	return K;
}

// calculates the gradient of kemeny's constant at pos, returns a matrix
void grad_kemeny(gsl_matrix* pos, const gsl_matrix* grad)
{
	// result from Neumann series
}