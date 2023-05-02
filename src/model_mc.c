/* Nat Watkins 10/08/2021. 
 * This creates the monte carlo function and related functions for minimization of params in model.
 */

#include<stdio.h>
#include<stdlib.h>
#ifndef min
	#define min(a,b) ((a) < (b) ? (a) : (b))
#endif

#include<gsl/gsl_math.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_block.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

#include<matrix_utils.h>
#include<model.h>
#include<model_mc.h>

/* returns a pointer to the value which is the mode of the histogram */
void* histogram_mode(HISTOGRAM* hist)
{
	int mode_index = 0;
	for(int i=0; i<hist->num_bins; i++)
	{
		if(hist->bin[mode_index] < hist->bin[i])
			mode_index = i;
	}
	return ((gsl_vector**)(hist->value))[mode_index];
}

/* sets the internal parameter using the mode of the bin returned from mc */
void model_minimize_mc(MODEL* model)
{
	HISTOGRAM* hist = model_mc(model);
	// calculate the mode of the histogram
	gsl_vector* mode = histogram_mode(hist);
	gsl_vector_memcpy(model->intern_params, mode);
	// model_mc allocates a hist, have to free the memory
	for(int i=0; i<hist->num_bins; i++)
	{
		gsl_vector_free(((gsl_vector**)(hist->value))[i]);
	}
	free(hist);
}

HISTOGRAM* model_mc(MODEL* model)
{
	// minimizes the free energy wrt the internal parameters at the current position
	// phase diagram does minimize over all external parameters, i.e. generates the phase diagram
	
	/* the energy landscape is determined by the particular form of the generated matrices. We want
	* average over the ensemble when calculating the critical internal parameters.
	*/
	const int MAX_NUM_STEPS = 50;
	const double dx = 0.1;

	// not using extern var because they are tensors of size num_states, need vectors in intern param space
	gsl_vector* point = gsl_vector_calloc(model->num_intern_params);
	gsl_vector* pros_point = gsl_vector_calloc(model->num_intern_params);
	
	// allocate and struct hist and allocate members, rare instance where memory is allocated in a function	
	HISTOGRAM* hist = malloc(sizeof(HISTOGRAM));
	hist->num_bins = 10;
	hist->dx = 1.0/hist->num_bins;
	hist->value = malloc(hist->num_bins*sizeof(gsl_vector*));
	for(int i=0; i<hist->num_bins; i++)
	{
		((gsl_vector**)(hist->value))[i] = gsl_vector_calloc(model->num_intern_params);
	}
	
	// define the value of each bin
	for(int i=0; i<hist->num_bins; i++)
		((gsl_vector**)(hist->value))[i]->data[0] = i*hist->dx + dx/2;
	
	// main loop for monte carlo algorithm
	for(int i=0; i<MAX_NUM_STEPS; i++)
	{
		model_mc_proposal(model, pros_point, point);
		if(model_mc_accept(model, pros_point, point))
		{
			gsl_vector_memcpy(point, pros_point);
			model_mc_bin(hist, point);
		}
	}
	
	gsl_vector_free(point);
	gsl_vector_free(pros_point);
	return hist;
}

void model_mc_proposal(MODEL* model, gsl_vector* pros_point, gsl_vector* point)
{
	double range = 0;
	double offset = 0;
	for(int i=0; i<point->size; i++)
	{
		range = model->intern_params_upper_boundary->data[i] - model->intern_params_lower_boundary->data[i];
		offset = model->intern_params_lower_boundary->data[i];
		// the simplest proposal is to uniformly sample the parameter space
		pros_point->data[i] = rand()/RAND_MAX*range + offset;
		printf("pros alpha %f\n", pros_point->data[i]);
	}
}

int model_mc_accept(MODEL* model, gsl_vector* pros_point, gsl_vector* point)
{
	const int kT = 1.0;
	int success = 0;
	gsl_vector_memcpy(model->intern_params, pros_point);
	gsl_vector_fprintf_delim(stdout, model->intern_params, ", ");
	double fe_new = model->free_energy(model);
	gsl_vector_memcpy(model->intern_params, point);
	double fe_old = model->free_energy(model);
	double acceptance_prob = min(exp(-(fe_new-fe_old)/kT), 1);
	success = (int)(acceptance_prob+rand()/RAND_MAX);
	return success;
}

void model_mc_bin(HISTOGRAM* hist, gsl_vector* point)
{
	gsl_vector* delta = gsl_vector_calloc(point->size);
	int in_bin = 0;
	gsl_vector* value = NULL;
	for(int i=0; i<hist->num_bins; i++)
	{
		// doesn't change point
		gsl_vector_memcpy(delta, point);
		value = ((gsl_vector**)(hist->value))[i];
		gsl_blas_daxpy(-1.0, value, delta);
		// assume it is in the bin
		in_bin = 1;
		for(int j=0; j<point->size; j++)
		{
			if(abs(delta->data[j]) > hist->dx)
			{
				// if it finds the delta has a component greater the dx, it is not in bin, break
				in_bin = 0;
				break;
			}
		}
		if(in_bin)
		{
			hist->bin[i]++;
			break;
		}
	}
	gsl_vector_free(delta);
}