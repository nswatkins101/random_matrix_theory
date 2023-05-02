#include<stdio.h>
#include<math.h>

#include<matrix_calc_utils.h>
#include<graph_alloc_utils.h>
#include<matrix_utils.h>

/* all of the graph_alloc functions return a gsl_matrix* pointing to the adjacency matrix
 * of such a randomly generated graph
 */

gsl_matrix* ER_RG(const int DIM, double p)
{
	gsl_matrix* a_mat = gsl_matrix_alloc(DIM, DIM);
	ER_RG_inplace(a_mat, p);
	return a_mat;
}

gsl_matrix* ER_RDG(const int DIM, double p)
{
	gsl_matrix* a_mat = gsl_matrix_alloc(DIM, DIM);
	ER_RDG_inplace(a_mat, p);
	return a_mat;
}

gsl_matrix* ER_RDAG(const int DIM, double p)
{
	gsl_matrix* mat = gsl_matrix_alloc(DIM, DIM);
	for(int i=0; i<DIM; i++)
	{
		for(int j=0; j<i; j++)
		{
			mat->data[i*mat->tda+j] = (int)(p+(double)rand()/RAND_MAX);							// generates pseudo-random numbers between 0 and 1, adds p and floors
		}
	}
	return mat;
}

void ER_RG_inplace(gsl_matrix* mat, double p)
{
	gsl_matrix_set_zero(mat);
	for(int i=0; i<mat->size1; i++)
	{
		for(int j=0; j<i; j++)
		{
			mat->data[i*mat->tda+j] = (int)(p+(double)rand()/RAND_MAX);							// generates pseudo-random numbers between 0 and 1, adds p and floors
			mat->data[j*mat->tda+i] = mat->data[i*mat->tda+j];
		}
	}
}

void ER_RDG_inplace(gsl_matrix* mat, double p)
{
	gsl_matrix_set_zero(mat);
	for(int i=0; i<mat->size1; i++)
	{
		for(int j=0; j<mat->size2; j++)
		{
			mat->data[i*mat->tda+j] = (int)(p+(double)rand()/RAND_MAX);							// generates pseudo-random numbers between 0 and 1, adds p and floors
		}
	}
}

void ER_RDAG_inplace(gsl_matrix* mat, double p)
{
	gsl_matrix_set_zero(mat);
	for(int i=0; i<mat->size1; i++)
	{
		for(int j=0; j<i; j++)
		{
			mat->data[i*mat->tda+j] = (int)(p+(double)rand()/RAND_MAX);							// generates pseudo-random numbers between 0 and 1, adds p and floors
		}
	}
}

/* Creates branches until the size of the matrix hits the fixed size of a_mat.
 * Therefore, it could happen that the outermost shell of the tree is incomplete.
 */
void bethe_lattice(gsl_matrix* a_mat, int coord)
{
	if(coord >= a_mat->size1)
		printf("ERROR: bethe_lattice: coordination is equal or greater than number of vertices.\n");
	gsl_matrix_set_zero(a_mat);
	if(coord >= a_mat->size1)
		return;
	// let the root by the first labelled node
	for(int i=0; i<coord; i++)
		gsl_matrix_set(a_mat, 0, i+1, 1);

	for(int i=1; i<a_mat->size1; i++)
	{
		for(int k=0; k<coord; k++)
		{
			int j = 2+i*(coord-1)+k;
			if(j<a_mat->size2)
				gsl_matrix_set(a_mat, i, j, 1);
		}
	}
	// the above only defines the upper triangular part
	gsl_matrix_symmetrize(a_mat);
}