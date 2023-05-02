/* header file for model_min.c */

#include<gsl/gsl_block.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

#ifndef MODEL_MIN_H
#define MODEL_MIN_H


/* functions that minimize the free energy wrt the internal parameters at the current position */
void model_minimize_global(MODEL* model);
void model_minimize_grad_descent(MODEL* model);
/* helper functions, static for encapsulation */
size_t multid_indices_to_index(gsl_vector_int* multid_indices, size_t size1);
void index_to_multid_indices(gsl_vector_int* multid_indices, size_t index, size_t size1);
size_t minimize_global(double* data, size_t size);
void gsl_vector_double_int(gsl_vector* dest, gsl_vector_int* src);
void gsl_vector_int_double(gsl_vector_int* dest, gsl_vector* src);
void gsl_vector_my_memcpy(gsl_vector* dest, gsl_vector* src);

#endif