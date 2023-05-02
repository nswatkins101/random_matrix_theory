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

void model_init(MODEL* model, size_t len, size_t num_intern_params, size_t num_extern_params, size_t num_states, gsl_block* input_data, gsl_vector* lower_boundary, gsl_vector* upper_boundary)
{
	// initialize function pointers
	model->init = &model_init;
	model->free = &model_free;
	model->set_pos = &model_set_pos;
	model->params_in_boundary = &model_params_in_boundary;
	model->expected_spectrum = &model_expected_spectrum;
	model->phase_diagram = &model_phase_diagram;
	model->set_grad_free_energy = &model_set_grad_free_energy;		// default function, calculates finite difference grad
	model->mean_free_energy = &model_mean_free_energy;
	
	model->generate = NULL;
	model->free_energy = NULL;
	
	model->len = len;
	model->num_intern_params = num_intern_params;
	model->num_extern_params = num_extern_params;
	
	// set block using passed block*
	model->block = input_data;
	
	// grad block is allocated using gsl function
	model->grad_block = gsl_block_calloc(len*num_intern_params);
	
	// malloc block and set each member manually using function args
	model->intern_params_lower_boundary = lower_boundary;
	model->intern_params_upper_boundary = upper_boundary;
	
	// params at a particular position will be slices of the data block
	model->pos =0;
	model->extern_params = malloc(sizeof(gsl_vector));
	model->extern_params->size = num_extern_params;
	model->extern_params->stride = 1;											// elements of the vector are contigous in the data block
	model->extern_params->data = model->block->data;
	model->extern_params->block = model->block;
	model->extern_params->owner = 0;										// don't deallocate block when etern_params is deallocated
	
	model->intern_params = malloc(sizeof(gsl_vector));
	model->intern_params->size = num_intern_params;
	model->intern_params->stride = 1;											// elements of the vector are contigous in the data block
	model->intern_params->data = model->block->data + num_extern_params;
	model->intern_params->block = model->block;
	model->intern_params->owner = 0;
	
	model->grad = malloc(sizeof(gsl_vector));
	model->grad->size = num_intern_params;
	model->grad->stride = 1;															// elements of the vector are contigous in the data block
	model->grad->data = model->grad_block->data;
	model->grad->block = model->grad_block;
	model->grad->owner = 0;															// don't deallocate block when intern_params is deallocated

	model->system = gsl_matrix_calloc(num_states, num_states);										
}

void model_free(MODEL* model)
{
	gsl_matrix_free(model->system);
	
	// free in reverse order in which I allocated
	// only frees the data it allocated, not the data passed to it
	gsl_vector_free(model->extern_params);
	gsl_vector_free(model->intern_params);
	gsl_vector_free(model->grad);
}

void model_set_pos(MODEL* model, size_t i)
{
	if(i < model->len)
	{
		model->pos = i;
		model->extern_params->data = model->block->data + model->pos*(model->num_extern_params+model->num_intern_params);
		model->intern_params->data = model->block->data + model->num_extern_params + model->pos*(model->num_extern_params+model->num_intern_params);
		model->grad->data = model->grad_block->data + model->pos*model->num_intern_params;
	}
	else
		printf("ERROR: index in model_set_pos(MODEL*, int) exceeds length.\n");
}

int model_params_in_boundary(MODEL* model)
{
	double element;
	double upper;
	double lower;
	for(size_t i=0; i<model->intern_params->size; i++)
	{
		element = gsl_vector_get(model->intern_params, i);
		upper = model->intern_params_upper_boundary->data[i];
		lower = model->intern_params_lower_boundary->data[i];
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
	// not implemented, would calculate the spectral density of the system (i.e. the matri systme) averaging over the ensemble
	// would probably set or return a histogram
}

void model_phase_diagram(MODEL* model)
{
	// runs minimize over phase space, i.e. all external parameters, setting the values of internal parameters and hence creating a phase diagram
	for(size_t i=0; i<model->len; i++)
	{
		model->set_pos(model, i);
		model->minimize(model);
	}		
}

void model_set_grad_free_energy(MODEL* model)
{
	const double dx = 0.1;
	double forward_fe = 0;
	double backward_fe = 0;
	double derivative = 0;
	// calculating the gradient to <F>, the averaging is done in free_energy(MODEL*)
	for(size_t i=0; i<model->intern_params->size; i++)
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

double model_mean_free_energy(MODEL* model)
{
	// returns the ensemble averaged free energy according to the ensemble specified in generate
	const unsigned NUM_SAMPLES = 10;
	double FE = 0;	
	for(unsigned i=0; i<NUM_SAMPLES; i++)
	{
		model->generate(model);
		FE += model->free_energy(model);
	}
	// generate the model once more so that model->system isn't garbage, free energy calculation operates on system data in place
	model->generate(model);

	return FE;
}