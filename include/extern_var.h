/* header file for model.c */

#include<gsl/gsl_matrix.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_permutation.h>

/* external variables defined here, prefix g_ indicates a global variable */
extern int g_temp_num;
extern int g_dim;

extern double** g_temp_array;
extern int* g_temp_array_usage;

extern gsl_matrix** g_temp_mat;
extern int* g_temp_mat_usage;
extern gsl_matrix_complex** g_temp_cmat;
extern int* g_temp_cmat_usage;
extern gsl_vector** g_temp_vec;
extern int* g_temp_vec_usage;
extern gsl_vector_complex** g_temp_cvec;
extern int* g_temp_cvec_usage;

extern gsl_eigen_symm_workspace* g_w_symm;
extern gsl_eigen_nonsymm_workspace* g_w_nonsymm;
extern gsl_permutation* g_permutation;

extern gsl_matrix* identity_mat;

/* functions used to intiailize, allocate, and distribute temp memory */
int init_extern_var(int DIM);

int free_extern_var();

double* get_extern_temp_array();
int reset_extern_temp_array(double* temp_array);

gsl_matrix* get_extern_temp_mat();
int reset_extern_temp_mat(gsl_matrix* temp_mat);

gsl_matrix_complex* get_extern_temp_cmat();
int reset_extern_temp_cmat(gsl_matrix_complex* temp_cmat);

gsl_vector* get_extern_temp_vec();
int reset_extern_temp_vec(gsl_vector* temp_vec);

gsl_vector_complex* get_extern_temp_cvec();
int reset_extern_temp_cvec(gsl_vector_complex* temp_cvec);

