#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>

void gsl_matrix_inverse(gsl_matrix* z, gsl_matrix* m);

void gsl_matrix_complex_inverse(gsl_matrix_complex* z, gsl_matrix_complex* m);

// take the symmetric part of the matrix
void gsl_matrix_symmetric_part(gsl_matrix* mat);

// completes a triangular matrix, assuming symmetry
void gsl_matrix_symmetrize(gsl_matrix* mat);
	
double gsl_matrix_trace(const gsl_matrix* mat);

double gsl_vector_l1_norm(const gsl_vector* x);
