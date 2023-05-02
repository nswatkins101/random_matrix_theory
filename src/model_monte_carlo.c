/* Nat Watkins 03/03/2021. 
 * This creates the model struct, which encompasses a number of gsl data blocks and several functions to allocate them and calculate their properties.
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

void model_minimize_monte_carlo(MODEL* model)
{
	// minimizes the free energy wrt the internal parameters at the current position
	// phase diagram does minimize over all external parameters, i.e. generates the phase diagram
	
	/* the energy landscape is determined by the particular form of the generated matrices. We want
	* average over the ensemble when calculating the critical internal parameters.
	*/
	const int MAX_NUM_STEPS = 50;
	const double FTOL = 0.0001;
	const double dx = 0.1;

	gsl_vector* point = get_extern_temp_vec();
	gsl_vector* pros_point = get_extern_temp_vec();
	
	for(int i=0; i<MAX_NUM_STEPS; i++)
	{
		model_minimize_monte_carlo_proposal(pros_point, point);
		if(model_minimize_monte_carlo_accept(model, pros_point, point))
			gsl_vector_memcpy(point, pros_point);
		// calculate the derivative at point model->pos
		model->set_grad_free_energy(model);
		if(gsl_blas_dnrm2(model->grad) < FTOL*FTOL)
		{
			//printf("reached FTOL\n");
			break;
		}
		// updates the intern_params at the point j in phase space, moves down the gradient
		//printf("dp is %f\n", dx*model->grad->data[0]);
		gsl_blas_daxpy(-dx, model->grad, model->intern_params);
		// check that intern_params is still in the boundary
		if(!model->params_in_boundary(model))
		{
			//printf("ERROR: model_minimize: out of bounds\n");
			// if not in boundary, reverse update step and break
			gsl_blas_daxpy(dx, model->grad, model->intern_params);
			break;
		}
	}
}

void model_phase_diagram(MODEL* model)
{
	// runs minimize over all positions, i.e. all external parameters, hence creating a phase diagram
	for(int i=0; i<model->len; i++)
	{
		model->set_pos(model, i);
		model->minimize(model);
	}		
}

void model_generate_bethe_erdos_renyi(MODEL* model)
{
	// get the parameters for the non local and local matrices from extern params, and mixing param alpha from intern param
	double p = gsl_vector_get(model->extern_params, 0);
	int z = gsl_vector_get(model->extern_params, 1);
	double alpha = gsl_vector_get(model->intern_params, 0);
	
	gsl_matrix* nonlocal_a_mat = get_extern_temp_mat();
	gsl_matrix* local_a_mat = get_extern_temp_mat();
	
	// generate the non local and local matrices	
	ER_RG_inplace(nonlocal_a_mat, p);
	bethe_lattice(local_a_mat, z);
	
	// mix them according to alpha, internal param
	// Di Patti convention, alph=1 corresponds to a local graph, alpha=0 means a non-local lattice
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, alpha, identity_mat, local_a_mat,  (1-alpha), nonlocal_a_mat);
	gsl_matrix_set_zero(model->system);
	compute_matrix_laplacian_inplace(model->system, nonlocal_a_mat);
	
	reset_extern_temp_mat(nonlocal_a_mat);		
	reset_extern_temp_mat(local_a_mat);
}

double model_free_energy_kemeny(MODEL* model)
{
	// take the trace of the pseudo inverse of the system matrix, which will be the laplacian of the graph
	// actually, will be using my g-inverse because it is computationally less expensive
	// returns the ensemble averaged kemeny
	const int NUM_SAMPLES = 100;
	
	double kemeny = 0;
	gsl_matrix* g_mat = get_extern_temp_mat();
	
	for(int i=0; i<NUM_SAMPLES; i++)
	{
		model->generate(model);
		// compute_my_g_inverse destroys model->system, but that doesn't matter because we only needed it for one calculation
		gsl_matrix_set_zero(g_mat);
		compute_my_g_inverse(g_mat, model->system);
		kemeny += gsl_matrix_trace(g_mat)/NUM_SAMPLES;
		//kemeny += gsl_matrix_kemeny(model->system)/NUM_SAMPLES;
	}
	// generate the model once more so that model->system isn't garbage
	model->generate(model);
	reset_extern_temp_mat(g_mat);
	//printf("kemeny %f\n", kemeny);

	return kemeny;
}

void model_set_grad_free_energy(MODEL* model)
{
	const double dx = 0.1;
	double forward_fe = 0;
	double backward_fe = 0;
	double derivative = 0;
	// calculating the gradient to <F>, the averaging is done in free_energy(MODEL*)
	for(int i=0; i<model->intern_params->size; i++)
	{
		model->intern_params->data[i] += dx;
		forward_fe = model->free_energy(model);
		model->intern_params->data[i] -= dx;
		backward_fe = model->free_energy(model);
		derivative = (forward_fe-backward_fe)/dx;
		//printf("derivative %f\n", derivative);
		gsl_vector_set(model->grad, i, derivative);
	}
}