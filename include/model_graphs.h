/* header file for model.c */

#include<gsl/gsl_block.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

#ifndef MODEL_GRAPH_H
#define MODEL__GRAPH_H

/* functions particular to the model I'm studying */
//static void compute_my_g_inverse(gsl_matrix* g_mat, gsl_matrix* q_mat);

void model_generate_bethe_erdos_renyi(MODEL* model);

double model_free_energy_kemeny(MODEL* model);


#endif