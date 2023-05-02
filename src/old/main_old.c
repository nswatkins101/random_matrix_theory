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
	int num_pnts = 10;
	double log_max = -4;
	int num_samples = 5;
	double span[] = {-5, 5};
	double dx = 5e-2;
	double num_evals = (int)((span[1]-span[0])/dx);
	
	/* open the data files for writing matrix elements and esd */
	char data_dir[] = "/root/my-documents/opt-ads/random-matrix-theory/data";
	char mat_file[100];
	sprintf(mat_file, "%s/matrix.csv", data_dir);
	FILE* mat_fp = fopen(mat_file, "w");
	char esd_file[100];
	sprintf(esd_file, "%s/esd.csv", data_dir);
	FILE* esd_fp = fopen(esd_file, "w");
	
	/* allocate memory for gsl data types */
	printf("creating gsl structs...");
	gsl_matrix* a_mat = gsl_matrix_alloc(DIM, DIM);
	
	/* allocating memory for external variables */
	init_extern_var(DIM);
	printf("finished.\n");
	
	/* spectrum will be stored as a 2d array, i.e. {(eval1, density1),...}*/
	double* esd = calloc(num_evals*2, sizeof(double));
	for(int i=0; i < num_evals; i++)
		esd[2*i+0] = span[0] + i*dx;
	
	/* storage for kemeny's constant and critical alpha */
	double* KvsA = calloc(2*num_pnts, sizeof(double));
	double* AvsLogDelta = calloc(2*num_pnts, sizeof(double));
		
	double coord = 3;
	double delta = 0;
	double log_delta = 0;
	double a = 0;
	double expectation = 0;
	double (*op_fptr)(void**) = &op_kemeny;
	void** op_args[] = {a_mat};
	void (*gen_fptr)(void**) = &gen_L_NL;
	void** gen_args[] = {a_mat, &coord, &delta, &a};
	double min[2];
	for(int i=1; i<num_pnts; i++)
	{
		log_delta = i/num_pnts*log_max;
		delta = pow(10.0, log_delta);
		for(int j=1; j<num_pnts; j++)
		{
			a = j/num_pnts;
			expectation = random_matrix_expectation(num_samples, op_fptr, op_args, gen_fptr, gen_args);
			KvsA[2*j] = a;
			KvsA[2*j+1] = expectation;
		}
		array_min(min, num_pnts, KvsA);
		AvsLogDelta[2*i] = log_delta;
		AvsLogDelta[2*i+1] = min[0];
	}
	
	
	/* write esd */
	printf("writing data...");
	array2d_fprintf_delim(esd_fp, esd, num_evals, 2, 2, ", ");
	printf("finished.\n");
	
	/* close the files */
	fclose(mat_fp);
	fclose(esd_fp);

	/* free gsl memory */
	printf("freeing gsl memory...\n");
	gsl_matrix_free(a_mat);
	
	free_extern_var();

	printf("end of program.\n");
	return 0;
}


	
	

