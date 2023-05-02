#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

gsl_matrix* ER_RG(const int DIM, double p);

gsl_matrix* ER_RDG(const int DIM, double p);

gsl_matrix* ER_RDAG(const int DIM, double p);

void ER_RG_inplace(gsl_matrix* mat, double p);

void ER_RDG_inplace(gsl_matrix* mat, double p);

void ER_RDAG_inplace(gsl_matrix* mat, double p);

void bethe_lattice(gsl_matrix* a_mat, int coord);