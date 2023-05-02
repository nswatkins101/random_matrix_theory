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
#include<model.h>
#include<model_graphs.h>
#include<model_mc.h>
#include<model_min.h>
#include<color_macros.h>

int main()
{
	void* cmd_args;
	srand(time(NULL));											// sets seed using time
	const int VERBOSE = 0;
	/* disable buffering for testing */
	setbuf(stdout, NULL);
	
	const int DIM = 100;
	const double MAXTIME = 1800;									// in seconds, I should add a function that keeps track of how long the program took, not implemented
	
	const size_t num_extern_params = 2;
	const size_t num_intern_params = 1;
	const size_t tda = num_extern_params+num_intern_params;

	const size_t num_pnts = 10;
	const size_t len = pow(num_pnts, num_extern_params);
	const double dx = 1.0/num_pnts;
	
	/* open the data files for writing data */
	char data_dir[] = "/root/my-documents/opt-ads/random-matrix-theory/data";
	char block_file[100];
	sprintf(block_file, "%s/phase_diagram.csv", data_dir);
	FILE* block_fp = fopen(block_file, "w");
		
	/* allocating memory for external variables */
	printf("allocating extern memory...");
	init_extern_var(DIM);
	printf(GRN "done.\n" RESET);
	
	/* create model, define the boundaries and external params */
	printf("allocating model memory...");
	MODEL* model = calloc(1, sizeof(MODEL));
	printf(GRN "done.\n" RESET);
	printf("creating interal parameters boundary vectors...");
	gsl_vector* lower_boundary = gsl_vector_calloc(num_intern_params);
	gsl_vector* upper_boundary = gsl_vector_calloc(num_intern_params);
	// boundary for internal parameters, only free parameter is alpha, which is in (0,1)
	lower_boundary->data[0] = 0.01;
	upper_boundary->data[0] = 0.99;
	printf(GRN "done.\n" RESET);
	printf("allocating data block..");
	gsl_block* block = gsl_block_calloc(len*tda);
	printf(GRN "done.\n" RESET);
	printf("initializing elements of data block...");
	// define the points of the mesh, (p, z)
	double p = 1;
	int z = 2;
	// the values of alpha are init to 0.5, but minimization will change there values
	double init_alpha = 0.5;
	size_t row = 0;
	for(size_t i=0; i<num_pnts; i++)
	{
		p = 5*(double)(i)/DIM;
		for(size_t j=0; j<num_pnts; j++)
		{
			z = 5*j+2;
			row = i*num_pnts+j;
			block->data[row*tda+0] = p;
			block->data[row*tda+1] = z;
			block->data[row*tda+2] = init_alpha;
		}
	}
	printf(GRN "done.\n" RESET);
	// initialize the members of the model struct
	printf("initializing model...");
	// strangely enough, model->init is undefined, so I am defining it then running it, to define the other function ptrs
	model->init = &model_init;
	model->init(model, len, num_intern_params, num_extern_params, DIM, block, lower_boundary, upper_boundary);
	printf(GRN "done.\n" RESET);
	// have to still define the particular functions that define how the model behaves
	model->generate = model_generate_bethe_erdos_renyi;
	model->free_energy = model_free_energy_kemeny;
	// redefining the default minimization function to be the global minimization scheme
	model->minimize = &model_minimize_global;
	// generates the minimized structures for each condition, stores them in block, aka minimizes FE wrt to internal parameters at each set of external parameters
	printf("generating phase diagram...");
	model->phase_diagram(model);
	printf(GRN "done.\n" RESET);
	
	/* write data block */
	printf("writing data...");
	array2d_fprintf_delim(block_fp, model->block->data, len, tda, tda, ", ");
	printf(GRN "done.\n" RESET);
	
	/* close the files */
	fclose(block_fp);

	/* free gsl memory */
	printf("freeing model's internal memory...");
	model->free(model);
	printf(GRN "done.\n" RESET);
	printf("freeing boundary vectors...");
	gsl_vector_free(lower_boundary);
	gsl_vector_free(upper_boundary);
	printf(GRN "done.\n" RESET);
	printf("freeing data block...");
	gsl_block_free(block);
	printf(GRN "done.\n" RESET);
	printf("freeing model struct...");
	free(model);
	printf(GRN "done.\n" RESET);
	printf("freeing extern memory...");
	free_extern_var();
	printf(GRN "done.\n" RESET);

	printf("end of program.\n");
	return 0;
}


	
	

