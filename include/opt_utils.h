#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

void gsl_matrix_grad_trace(gsl_matrix* grad, const gsl_matrix* mat);

void gsl_vector_grad_l1_norm(gsl_vector* grad, const gsl_vector* x);

// gets passed a symmetric matrix, calculates eigenvalues and K
double gsl_matrix_kemeny(const gsl_matrix* mat);

// calculates the gradient of kemeny's constant at pos, returns a matrix
void gsl_matrix_grad_kemeny(gsl_matrix* pos, const gsl_matrix* grad);