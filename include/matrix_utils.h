#include<gsl/gsl_matrix.h>
#include<gsl/gsl_errno.h>

#define MY_GSL_ERROR(error) { \
				char* reason = malloc(100*sizeof(char)); \
				sprintf(reason, "In function '%s': %s", __func__, gsl_strerror(error)); \
				GSL_ERROR(reason, error); \
				free(reason); \
				}

int gsl_vector_fprintf_delim(FILE* stream, gsl_vector* vector, char* delimiter);

int gsl_matrix_fprintf_delim(FILE* stream, gsl_matrix* matrix, char* delimiter);

int gsl_matrix_complex_fprintf_delim(FILE* stream, gsl_matrix_complex* matrix, char* delimiter);

// C "overload" doesn't exist, style is to append _type to end of function
void gsl_matrix_set_ptr(gsl_matrix* m, const double* x);

// opposite of the above function, i.e. sets an array based on a matrix
void ptr_set_gsl_matrix(double* x, const gsl_matrix* m);

gsl_matrix* gsl_vector_diagonal_matrix(gsl_vector* x);

void gsl_matrix_complex_real_memcpy(gsl_matrix* mat, gsl_matrix_complex* z_mat);

void gsl_matrix_complex_imag_memcpy(gsl_matrix* mat, gsl_matrix_complex* z_mat);

gsl_matrix_view gsl_matrix_complex_real(gsl_matrix_complex* z);

gsl_matrix_view gsl_matrix_complex_imag(gsl_matrix_complex* z);
