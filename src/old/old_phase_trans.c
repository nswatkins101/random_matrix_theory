/* Nat Watkins 03/03/2021. 
 * This project seeks to generate random graphs, calculate, and plot the empirical limiting distribution of their spectra.
 */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#include<gsl/gsl_math.h>
#include<gsl/gsl_block.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>

#include<matrix_utils.h>
#include<matrix_calc_utils.h>
#include<array_utils.h>
#include<opt_utils.h>
#include<graph_alloc_utils.h>
#include<graph_utils.h>
#include<random_matrix_utils.h>
#include<extern_var.h>

void compute_my_g_inverse(gsl_matrix* g_mat, gsl_matrix* q_mat)
{
	int n = q_mat->size1;
	// destroys q_mat
	gsl_matrix_add_constant(invertible_mat, 1.0/n);
	gsl_cblas_dgemm(q_mat, invertible_mat);
	gsl_matrix_inverse(g_mat, invertible_mat);
	reset_extern_temp_mat(invertible_mat);
}

/* the derivative requires (in my special case) two matrices, the local
 * and nonlocal adjacency matrices. 
 */
double deriviative_kemeny(gsl_matrix* local_a_mat, gsl_matrix* nonlocal_a_mat)
{
	double derivative;
	q_mat = get_extern_temp_mat();
	g_mat = get_exter
	temp_mat = get_extern_temp_mat();
	temp_product = get_extern_temp_mat();
	temp_array = get_extern_temp_array();
	int n = temp_mat->size;
	double* temp_array = calloc(sizeof(double), n);
	
	compute_matrix_laplacian_inplace(q_mat, nonlocal_a_mat, temp_array);		// laplacian on the Erdos-Renyi graph (nonlocal graph)
	compute_my_g_inverse(g_mat,q_mat);																// g-inverse of laplacian on Erdos-Renyi graph
	compute_matrix_laplacian_inplace(q_mat, local_a_mat, temp_arrayr);			// laplacian on the Bethe lattice (local graph)
	gsl_matrix_memcpy(temp_mat, g_mat);
	gsl_cblas_dgemm(temp_mat, q_mat, temp_product);
	gsl_matrix_memcpy(temp_mat, g_mat);
	gsl_matrix_dgemm(temp_mat, temp_product, product);
	
	reset_extern_temp_mat(temp_mat);
	reset_extern_temp_mat(temp_product);
	derivative = gsl_matrix_trace(g_mat) - gsl_matrix_trace(product)
	return derivative;
}

double op_derivative_kemeny(void** op_args)
{
	gsl_matrix* local_a_mat = op_args[0];
	gsl_matrix* nonlocal_a_mat = op_args[1];
	return deriviative_kemeny(local_a_mat, nonlocal_a_mat);
}

void gen_bethe_erdos_renyi(void** gen_args)
{
	gsl_matrix* local_a_mat = gsl_matrix*(gen_args[0]);
	int z = *(int*(gen_args[1]));
	gsl_matrix* nonlocal_a_mat = gsl_matrix*(gen_args[2]);
	double p = *(double*(gen_args[3]));
	
	bethe_lattice(local_a_mat, z);
	Erdos_Renyi_inplace(nonlocal_a_mat, p);
}

void minimize(MODEL*)
{
	const int MAX_NUM_STEPS = 100;
	const double TOL = 0.01;
	const double EPS = 0.01;
	for(int i=0; i<model->len; i++)
	{
		for(int j=0; j<MAX_NUM_STEPS; j++)
		{
			// calculate the derivative at point j
			model->deriviative_free_energy(model, j);
			if(gsl_blas_dnrm2(model->gradient) << TOL)
				break;
			// updates the intern_params at the point j in phase space
			gsl_blas_daxpy(EPS, &gradient[j], &intern_params[j])
		}
	}
}
		
typedef struct tagMODEL
{
	void calloc(MODEL*) = &model_calloc;
	void minimize(MODEL*) = &model_minimize;
	void expected_spectrum(MODEL*) = &model_expected_spectrum;
	
	void generate(MODEL*);
	void free_energy(MODEL*);
	void derivative_free_energy(MODEL*, int i);	// sets the gradient of the scalar field at a particular point
	
	int len;																// the number of points to consider in the phase space
	gsl_vector* intern_params;									// free parameters that can change, i.e. alpha
	gsl_vector gradient;											// gradient of a scalar field on intern_params
	const gsl_vector* extern_param;						// points to a const block of memory, but ptr is not const, fixed parameters set at the beginning of the program, i.e. p
	
	gsl_vector state;
	gsl_matrix system;												// typically the laplacian or hamiltonian of the system
} MODEL;


void calc_phase_diagram(MODEL* model)
{
	for(int i=0; i<model->len; i++)
	{
		
}

int main()
{
	void* cmd_args;
	srand(time(NULL));											// sets seed using time
	const int VERBOSE = 0;
	/* disable buffering for testing */
	setbuf(stdout, NULL);
	
	const int DIM = 1000;
	const int MAXTIME = 1800;									// in seconds, I should add a function that keeps track of how long the program took
	int num_pnts = 10;
	int num_samples = 5;
	
	/* open the data files for writing matrix elements and esd */
	char data_dir[] = "/root/my-documents/opt-ads/random-matrix-theory/data";
	char curve_file[100];
	sprintf(curve_file, "%s/curve.csv", data_dir);
	FILE* curve_fp = fopen(esd_file, "w");
	
	/* allocate memory for gsl data types */
	printf("creating gsl structs...");
	gsl_matrix* local_a_mat = gsl_matrix_alloc(DIM, DIM);
	gsl_matrix* nonlocal_a_mat = gsl_matrix_alloc(DIM, DIM);
	printf("finished.\n");
	
	/* allocating memory for external variables */
	printf("allocating memory for external variables...");
	init_extern_var(DIM);
	printf("finished.\n");
	
	/* storage for kemeny's constant and critical alpha */
	double* curve = calloc(2*num_pnts, sizeof(double));		
		
	/* calculate the curve for which a phase transition occurs along */
	/* raster across delta,z space, writing coordinates for delta,z such that <f(delta, z)>=0 (or within a tolerance) (<> denotes an ensemble avg) */
	void** op_args = calloc(sizeof(void*), 2);
	void** gen_args = calloc(sizeof(void*), 4);
	int z;
	double delta;
	op_args[0] = local_a_mat;
	op_args[1] = nonlocal_a_mat;
	
	gen_args[0] = local_a_mat;
	gen_args[1] = &z;
	gen_args[2] = nonlocal_a_mat;
	gen_args[3] = &delta;
	int k=0;
	for(z=0; z<num_pnts; z++)
	{
		op_args[1] = z;
		for(delta=0; delta<1; delta+=1.0/num_steps)
		{
			f = random_matrix_expectation(num_samples, op_derivative_kemeny, op_args, generator, gen_args);
			if(abs(f) < tol)
			{
				curve[2*k] = delta;
				curve[2*k+1] = z;
				k++;
			}
		}
	}
	
	/* write curve data */
	printf("writing data...");
	array2d_fprintf_delim(curve_fp, curve, k, 2, 2, ", ");
	printf("finished.\n");
	
	/* close the files */
	fclose(curve_fp);

	/* free gsl memory */
	printf("freeing gsl memory...\n");
	gsl_matrix_free(local_a_mat);
	gsl_matrix_free(nonlocal_a_mat);
	
	free_extern_var();

	printf("end of program.\n");
	return 0;
}


	
	

