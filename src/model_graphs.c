/* Nat Watkins 02/01/2022 
 * This interfaces the model struct with graph operations that are necessary for the optimization project.
 */

#include<stdio.h>

#include<gsl/gsl_math.h>
#include<gsl/gsl_block.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

#include<matrix_utils.h>
#include<matrix_calc_utils.h>
#include<array_utils.h>
#include<graph_alloc_utils.h>
#include<graph_utils.h>
#include<random_matrix_utils.h>
#include<extern_var.h>
#include<opt_utils.h>

#include<model.h>
#include<model_graphs.h>

static void compute_my_g_inverse(gsl_matrix* g_mat, gsl_matrix* q_mat)
{
	FILE* q_fp = fopen("/root/my-documents/opt-ads/random-matrix-theory/data/q_mat", "w");
	FILE* temp_fp = fopen("/root/my-documents/opt-ads/random-matrix-theory/data/temp_mat", "w");
	FILE* g_fp = fopen("/root/my-documents/opt-ads/random-matrix-theory/data/g_mat", "w");
	gsl_matrix_fprintf_delim(q_fp, q_mat, ", ");
	fclose(q_fp);
	int n = q_mat->size1;
	// destroys q_mat
	gsl_matrix_add_constant(q_mat, 1.0/n);
	gsl_matrix_fprintf_delim(temp_fp, q_mat, ", ");
	fclose(temp_fp);
	gsl_matrix_inverse(g_mat, q_mat);
	gsl_matrix_fprintf_delim(g_fp, g_mat, ", ");
	fclose(g_fp);
}

void model_generate_bethe_erdos_renyi(MODEL* model)
{
	// get the parameters for the non local and local matrices from extern params, and mixing param alpha from intern param
	double p = gsl_vector_get(model->extern_params, 0);
	int z = gsl_vector_get(model->extern_params, 1);
	double alpha = gsl_vector_get(model->intern_params, 0);
	printf("model parameters: p=%f, z=%i, alpha=%f\n", p, z, alpha);
	
	gsl_matrix* nonlocal_a_mat = get_extern_temp_mat();
	gsl_matrix* local_a_mat = get_extern_temp_mat();
	
	// generate the non local and local matrices	
	ER_RG_inplace(nonlocal_a_mat, p);
	bethe_lattice(local_a_mat, z);
	
	// write
	FILE* nonlocal_fp = fopen("/root/my-documents/opt-ads/random-matrix-theory/data/nonlocal_mat", "w");
	FILE* local_fp = fopen("/root/my-documents/opt-ads/random-matrix-theory/data/local_mat", "w");
	gsl_matrix_fprintf_delim(nonlocal_fp , nonlocal_a_mat , ", ");
	fclose(nonlocal_fp );
	gsl_matrix_fprintf_delim(local_fp , local_a_mat, ", ");
	fclose(local_fp );
	
	// mix them according to alpha, internal param
	// Di Patti convention, alph=1 corresponds to a local graph, alpha=0 means a non-local lattice
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, alpha, identity_mat, local_a_mat,  (1-alpha), nonlocal_a_mat);
	gsl_matrix_set_zero(model->system);
	compute_matrix_laplacian_inplace(model->system, nonlocal_a_mat);
	FILE* sys_fp = fopen("/root/my-documents/opt-ads/random-matrix-theory/data/sys_mat", "w");
	gsl_matrix_fprintf_delim(sys_fp , model->system , ", ");
	fclose(sys_fp );
	
	reset_extern_temp_mat(nonlocal_a_mat);		
	reset_extern_temp_mat(local_a_mat);
}

double model_free_energy_kemeny(MODEL* model)
{
	// take the trace of the pseudo inverse of the system matrix, which will be the laplacian of the graph
	// actually, will be using my g-inverse because it is computationally less expensive
	
	gsl_matrix* g_mat = get_extern_temp_mat();	
	// compute_my_g_inverse destroys model->system, but that doesn't matter because we only needed it for one calculation
	gsl_matrix_set_zero(g_mat);
	compute_my_g_inverse(g_mat, model->system);
	double kemeny = gsl_matrix_trace(g_mat);
	
	// generate the model once more so that model->system isn't garbage
	model->generate(model);
	reset_extern_temp_mat(g_mat);

	return kemeny;
}