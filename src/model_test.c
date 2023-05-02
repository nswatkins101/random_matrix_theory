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

int main()
{
	void* cmd_args;
	srand(time(NULL));											// sets seed using time
	const int VERBOSE = 0;
	/* disable buffering for testing */
	setbuf(stdout, NULL);
	
	const int DIM = 1000;
	const int MAXTIME = 1800;									// in seconds, I should add a function that keeps track of how long the program took
	const double PERCOLATION = (double)(1.0/DIM);
	const int num_pnts = 10;
	const double dp = 1.0/num_pnts;
	const int num_extern_params = 2;
	const int num_intern_params = 1;
	const double tda = num_extern_params+num_intern_params;
	
	/* open the data files for writing matrix elements and esd */
	char data_dir[] = "/root/my-documents/opt-ads/random-matrix-theory/data";
	char block_file[100];
	sprintf(block_file, "%s/phase_diagram.csv", data_dir);
	FILE* block_fp = fopen(block_file, "w");
	
	/* allocate memory for gsl data types */
	
	/* allocating memory for external variables */
	printf("creating extern memory...\n");
	init_extern_var(DIM);
	printf("finished.\n");
	
	/* create model, define the external params and boundaries, minimize FE over interal params */
	MODEL* model = malloc(sizeof(MODEL));
	double* data = calloc(sizeof(double));
	double lower_boundary[] = {1.0};
	double upper_boundary[] = {0.0};
	// define the points of the mesh, (p, z)
	for(int i=0; i<len; i++)
	{
		for(int j=0; j<num_pnts; j++)
		{
			p = i*dp;
			z = j;
			data[(i+j)*tda+0] = p;
			data[(i+j)*tda+1] = z;
		}
	}
	model->init(model, num_intern_params, num_extern_params, num_states, data, lower_boundary, upper_boundary);
	// generates the minimized structures for each condition, stores them in list_intern_params aka block
	model->phase_diagram(model);
	
	/* write data block */
	printf("writing data...");
	array2d_fprintf_delim(block_fp, model->block->data, len, tda, tda, ", ");
	printf("finished.\n");
	
	/* close the files */
	fclose(block_fp);

	/* free gsl memory */
	printf("freeing model memory...\n");
	model->free(model);
	printf("freeing gsl memory...\n");
	gsl_matrix_free(a_mat);
	printf("freeing extern memory...\n");
	free_extern_var();

	printf("end of program.\n");
	return 0;
}


	
	

