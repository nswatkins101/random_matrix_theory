/* Nat Watkins 10/08/2021. 
 * This creates the  monte carlo function and related functions for minimization of params in model.
 */

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

#include<model.h>

#ifndef MODEL_MC_H
#define MODEL_MC_H

typedef struct tagHISTOGRAM
{
	int num_bins;
	double dx;
	
	int* bin;
	void* value;
} HISTOGRAM;

/* returns a pointer to the value which is the mode of the histogram */
void* histogram_mode(HISTOGRAM* hist);

/* sets the internal parameter using the mode of the bin returned from mc */
void model_minimize_mc(MODEL* model);

HISTOGRAM* model_mc(MODEL* model);

void model_mc_proposal(MODEL* model, gsl_vector* pros_point, gsl_vector* point);

int model_mc_accept(MODEL* model, gsl_vector* pros_point, gsl_vector* point);

void model_mc_bin(HISTOGRAM* hist, gsl_vector* point);

#endif