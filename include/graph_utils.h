#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>

gsl_matrix* compute_MFPT(gsl_matrix* g_mat, gsl_vector* pi_vec);

gsl_matrix* compute_matrix_laplacian(gsl_matrix* a_mat);

gsl_matrix* compute_matrix_laplacian_rw(gsl_matrix* a_mat);

gsl_matrix* compute_ever_passage(gsl_matrix* m_mat);

gsl_matrix* compute_first_passage(gsl_matrix* m_mat, double s);

// where g_mat is a generalized inverse with Ge=ge, see Hunter 2011
void compute_MFPT_inplace(gsl_matrix* mfpt_mat, gsl_matrix* g_mat, gsl_vector* pi_vec);

void compute_matrix_laplacian_inplace(gsl_matrix* m_mat, gsl_matrix* a_mat);

void compute_matrix_laplacian_rw_inplace(gsl_matrix* m_mat, gsl_matrix* a_mat);

void compute_ever_passage_inplace(gsl_matrix* ever_passage_mat, gsl_matrix* m_mat);

void compute_first_passage_inplace(gsl_matrix* mat, gsl_matrix* m_mat, double s);