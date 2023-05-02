/* header file for model.c */

#include<gsl/gsl_block.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

#ifndef MODEL_H
#define MODEL_H

// forward declaraction so I can create pointers to said structures in following function definitions
typedef struct tagMODEL MODEL;

/* helper functions for manipulating data in model struct, independent of particulars of the model */
void model_init(MODEL* model, size_t len, size_t num_intern_params, size_t num_extern_params, size_t num_state, gsl_block* data, gsl_vector* lower_boundary, gsl_vector* upper_boundary);

void model_free(MODEL* model);

void model_set_pos(MODEL* model, size_t i);

int model_params_in_boundary(MODEL* model);

void model_expected_spectrum(MODEL* model);

/* runs minimize over all positions, i.e. all external parameters, hence creating a phase diagram */
void model_phase_diagram(MODEL* model);

void model_set_grad_free_energy(MODEL* model);

double model_mean_free_energy(MODEL* model);
		
struct tagMODEL
{
	void (*init)(MODEL*, size_t, size_t, size_t, size_t, gsl_block*, gsl_vector*, gsl_vector*);
	void (*free)(MODEL*);
	void (*set_pos)(MODEL*, size_t);
	int (*params_in_boundary)(MODEL*);
	void (*expected_spectrum)(MODEL*) ;
	void (*phase_diagram)(MODEL*);
	void (*set_grad_free_energy)(MODEL*);			// sets the gradient of the mean field at the current point in parameter space, set by set_pos
	double (*mean_free_energy)(MODEL*);			// returns the value of the mean field at the current point in parameter space, set by set_pos
	
	// function ptrs that are set to NULL in init, defines the behavior of the model and necessary for most of the above functions, necessary to define in main
	void (*generate)(MODEL*);
	double (*free_energy)(MODEL*);
	void (*minimize)(MODEL*);
	
	size_t len;															// the number of points to consider in the phase space
	size_t num_intern_params;								// number of internal state parameters or degrees of freedom
	size_t num_extern_params;								//number of eternal state parameters, i.e. 2 for PV models
	gsl_block* block;													// stores all the params together in a block, each state, i.e. (extern_params, state/intern_params, FE) is concentated, slice data into managable chunks, i.e. matrices
	gsl_block* grad_block;										// list of gradients wrt DoF at every point in phase space, c style array
	
	size_t pos;															// sets the position of the model in external parameter space, so that functions don't have to do calculations at every point in parameter space
	gsl_vector* extern_params;								// point to a slice of block, don't have to modify list_intern_params directly
	gsl_vector* intern_params;									// point to a slice of block, don't have to modify list_intern_params directly
	gsl_vector* grad;													// gradient of FE wrt intern_params at current position in phase space

	gsl_vector* intern_params_upper_boundary;		// defines domain for the internal parameters
	gsl_vector* intern_params_lower_boundary;
	
	gsl_matrix* system;											// typically the laplacian or hamiltonian of the system, generator of time evolution of the chain
};

#endif