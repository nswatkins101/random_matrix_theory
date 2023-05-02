/* header file for model.c */

#include<gsl/gsl_math.h>
#include<gsl/gsl_block.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

// forward declaraction so I can create pointers to said structures in following function definitions
typedef struct tagMODEL MODEL;

void compute_my_g_inverse(gsl_matrix* g_mat, gsl_matrix* q_mat);

void model_init(MODEL* model, int len, int num_intern_params, int num_extern_params, int num_state);

void model_free(MODEL* model);

void model_set_pos(MODEL* model, int i);

int model_params_in_boundary(MODEL* model);

void model_expected_spectrum(MODEL* model);

/* minimizes the free energy wrt the internal parameters at the current position */
void model_minimize(MODEL* model);

/* runs minimize over all positions, i.e. all external parameters, hence creating a phase diagram */
void model_phase_diagram(MODEL* model);

void generate_bethe_erdos_renyi(MODEL* model);

double free_energy_kemeny(MODEL* model);

void model_default_set_grad_free_energy(MODEL* model);
		
struct tagMODEL
{
	void (*init)(MODEL*, int, int, int, int);
	void (*free)(MODEL*);
	void (*set_pos)(MODEL*, int);
	int (*params_in_boundary)(MODEL*);
	void (*expected_spectrum)(MODEL*) ;
	void (*minimize)(MODEL*);
	void (*phase_diagram)(MODEL*);
	
	void (*generate)(MODEL*);
	double (*free_energy)(MODEL*);
	void (*set_grad_free_energy)(MODEL*);			// sets the gradient of the scalar field at a particular point
	
	int len;																// the number of points to consider in the phase space, i.e. length of list_extern_params
	int num_intern_params;
	int num_extern_params;
	gsl_block* block;													// stores all the params together in a block, each state, i.e. (extern_params, intern_params) is concentated, slice data into managable chunks, i.e. matrices
	gsl_vector** list_intern_params;							// free parameters that can change, i.e. alpha. This is a vector on phase space, i.e. extern param space
	gsl_vector** grad_block;										// list of gradients at every point in phase space
	const gsl_vector** list_extern_params;				// points to a const block of memory, but ptr is not const, fixed parameters set at the beginning of the program, i.e. p. The mesh for phase space.
	
	int pos;																// sets the position of the model in external parameter space, so that functions don't have to do calculations at every point in parameter space
	gsl_vector* intern_params;									// these are pointers that point to a particular element in list_intern_params, don't have to modify list_intern_params directly
	gsl_vector* grad;													// gradient of a scalar field on intern_params
	const gsl_vector* extern_params;

	const gsl_block* intern_params_upper_boundary;	
	const gsl_block* intern_params_lower_boundary;
	
	int num_states;
	gsl_vector* state;
	gsl_matrix* system;												// typically the laplacian or hamiltonian of the system
};