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

void compute_my_g_inverse(gsl_matrix* g_mat, gsl_matrix* q_mat)
{
	int n = q_mat->size1;
	// destroys q_mat
	gsl_matrix_add_constant(q_mat, 1.0/n);
	gsl_matrix_inverse(g_mat, q_mat);
}

void model_init(MODEL* model, int len, int num_intern_params, int num_extern_params, int num_states, double* input_data)
{
	// initialize function pointers
	model->init = &model_init;
	model->set_pos = &model_set_pos;
	model->params_in_boundary = &model_params_in_boundary;
	model->minimize = &model_minimize;
	model->expected_spectrum = &model_expected_spectrum;
	model->phase_diagram = &model_phase_diagram;
	
	model->generate = NULL;
	model->free_energy = NULL;
	model->set_grad_free_energy = &model_default_set_grad_free_energy;
	
	model->len = len;
	model->num_intern_params = num_intern_params;
	model->num_extern_params = num_extern_params;
	model->num_states = num_states;
	
	// new initialization using blocks
	model->block = malloc(sizeof(gsl_block));
	model->block->data = input_data;
	
	// alloc params and have them reference slices of block
	list_extern_params = malloc(sizeof(gsl_matrix));
	list_intern_params = malloc(sizeof(gsl_matrix));
	
	list_extern_params.size1 = num_extern_params;
	list_extern_params.size2 = len;
	list_extern_params.tda = num_extern_params+num_intern_params;
	list_extern_params.data = model->raw_data->data;
	
	list_intern_params.size1 = num_intern_params;
	list_intern_params.size2 = len;
	list_intern_params.tda = num_extern_params+num_intern_params;
	list_intern_params.data = model->raw_data->data+num_extern_params;							// offset is the number of extern params before first intern param
	
	model->list_intern_params = (gsl_vector**)malloc(len*sizeof(gsl_vector*));						// free parameters that can change, i.e. alpha, a value of alpha for every external parameter tuple
	for(int i=0; i<len; i++)
		model->list_intern_params[i] = gsl_vector_calloc(num_intern_params);
	model->list_grad = (gsl_vector**)malloc(len*sizeof(gsl_vector*));																			// list of gradients at every point in phase space
	for(int i=0; i<len; i++)
		model->list_grad[i] = gsl_vector_calloc(num_intern_params);
	// doesn't allocate memory for extern params or boundaries, those should be defined in main
	model->list_extern_params = NULL;
	model->intern_params_lower_boundary = NULL;
	model->intern_params_upper_boundary = NULL;
	
	model->intern_params = gsl_vector_calloc(num_intern_params);																	// these are pointers that point to a particular element in list_intern_params, don't have to modify list_intern_params directly
	model->grad = gsl_vector_calloc(num_intern_params);																					// gradient of a scalar field on intern_params
	model->extern_params = gsl_vector_calloc(num_extern_params);

	model->state = gsl_vector_calloc(num_states);
	model->system = gsl_matrix_calloc(num_states, num_states);										
}

void model_set_pos(MODEL* model, int i)
{
	if(i < model->len)
	{
		model->pos = i;
		model->intern_params = model->list_intern_params[i];
		model->grad = model->list_grad[i];
		model->extern_params = model->list_extern_params[i];
	}
	else
		printf("ERROR: index in model_set_pos(MODEL*, int) exceeds length.\n");
}

int model_params_in_boundary(MODEL* model)
{
	double element;
	double upper;
	double lower;
	for(int i=0; i<model->intern_params->size; i++)
	{
		element = gsl_vector_get(model->intern_params, i);
		upper = gsl_vector_get(model->intern_params_upper_boundary, i);
		lower = gsl_vector_get(model->intern_params_lower_boundary, i);
		if(lower <= element && element <= upper)
		{
			// pass
		}
		// if an element lies outside the interval in any dimension, return false
		else
			return 0;														// here return 0 means failure, i.e. returns false
	}
	return 1;
}

void model_expected_spectrum(MODEL* model)
{
	// not implemented
}

void model_minimize(MODEL* model)
{
	// minimizes the free energy wrt the internal parameters at the current position
	// phase diagram does minimize over all external parameters, i.e. generates the phase diagram
	
	/* the energy landscape is determined by the particular form of the generated matrices. We want
	* average over the ensemble when calculating the critical internal parameters.
	*/
	const int MAX_NUM_STEPS = 100;
	const double TOL = 0.01;
	const double dx = 0.01;

	for(int j=0; j<MAX_NUM_STEPS; j++)
	{
		// calculate the derivative at point model->pos
		model->set_grad_free_energy(model);
		if(gsl_blas_dnrm2(model->grad) < TOL)
			break;
		// updates the intern_params at the point j in phase space
		gsl_blas_daxpy(dx, model->grad, model->intern_params);
		// check that intern_params is still in the boundary
		if(!model->params_in_boundary(model))
		{
			// if not in boundary, reverse update step and break
			gsl_blas_daxpy(-dx, model->grad, model->intern_params);
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

void generate_bethe_erdos_renyi(MODEL* model)
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
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, (1.0-alpha), identity_mat, local_a_mat,  alpha, nonlocal_a_mat);
	compute_matrix_laplacian_inplace(model->system, nonlocal_a_mat);
	
	reset_extern_temp_mat(nonlocal_a_mat);
	reset_extern_temp_mat(local_a_mat);
}

double free_energy_kemeny(MODEL* model)
{
	// take the trace of the pseudo inverse of the system matrix, which will be the laplacian of the graph
	// actually, will be using my g-inverse because it is computationally less expensive
	// returns the ensemble averaged kemeny
	const int NUM_SAMPLES = 100;
	
	double kemeny;
	gsl_matrix* g_mat = get_extern_temp_mat();
	
	for(int i=0; i<NUM_SAMPLES; i++)
	{
		model->generate(model);
		// compute_my_g_inverse destroys model->system, but that doesn't matter because we only needed it for one calculation
		compute_my_g_inverse(g_mat, model->system);
		kemeny += gsl_matrix_kemeny(g_mat)/NUM_SAMPLES;
	}
	// generate the model once more so that model->system isn't garbage
	model->generate(model);
	reset_extern_temp_mat(g_mat);
	
	return kemeny;
}

void model_default_set_grad_free_energy(MODEL* model)
{
	const double dx = 0.01;
	double forward_fe;
	double backward_fe;
	double derivative = 0;
	// calculating the gradient to <F>, the averaging is done in free_energy(MODEL*)
	for(int i=0; i<model->intern_params->size; i++)
	{
		model->intern_params->data[i] += dx;
		forward_fe = model->free_energy(model);
		model->intern_params->data[i] -= dx;
		backward_fe = model->free_energy(model);
		derivative = (forward_fe-backward_fe)/dx;
		gsl_vector_set(model->grad, i, derivative);
	}
}