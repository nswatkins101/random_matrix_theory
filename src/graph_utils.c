#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include<graph_utils.h>
#include<matrix_utils.h>
#include<matrix_calc_utils.h>
#include<extern_var.h>

// where g_mat is a generalized inverse with Ge=ge, see Hunter 2011
gsl_matrix* compute_MFPT(gsl_matrix* g_mat, gsl_vector* pi_vec)
{
	if(g_mat->size1 != g_mat->size2)
		printf("ERROR: matrix input to compute_MFPT(gsl_matrix*, gsl_vector*) is not square.\n");
	const int DIM = pi_vec->size;
	gsl_matrix* mfpt_mat  = gsl_matrix_alloc(DIM, DIM);
	compute_MFPT_inplace(mfpt_mat, g_mat, pi_vec);
	return mfpt_mat;
}

gsl_matrix* compute_matrix_laplacian(gsl_matrix* a_mat)
{
	if(a_mat->size1 != a_mat->size2)
		printf("ERROR: matrix input to compute_matrix_laplacian(gsl_matrix*) is not square.\n");
	const int DIM = a_mat->size1;
	gsl_matrix* m_mat = gsl_matrix_alloc(DIM, DIM);
	compute_matrix_laplacian_inplace(m_mat, a_mat);
	return m_mat;
}

gsl_matrix* compute_matrix_laplacian_rw(gsl_matrix* a_mat)
{
	if(a_mat->size1 != a_mat->size2)
		printf("ERROR: matrix input to compute_matrix_laplacian_rw(gsl_matrix*) is not square.\n");
	const int DIM = a_mat->size1;
	gsl_matrix* m_mat = gsl_matrix_alloc(DIM, DIM);	
	compute_matrix_laplacian_rw_inplace(m_mat, a_mat);
	return m_mat;
}

gsl_matrix* compute_ever_passage(gsl_matrix* m_mat)
{
	// ever-passage matrix is the limit of F(s) as s->0+
	double s = 1e-8;
	gsl_matrix* ever_passage_mat = compute_first_passage(m_mat, s);
	return ever_passage_mat;
}

/* calculates the first passage matrix in fourier space */
gsl_matrix* compute_first_passage(gsl_matrix* m_mat, double s)
{
	if(m_mat->size1 != m_mat->size2)
		printf("ERROR: matrix input to compute_first_passage(gsl_matrix*, double) is not square.\n");
	const int DIM = m_mat->size1;
	gsl_matrix* first_passage_mat = gsl_matrix_alloc(DIM, DIM);
	compute_first_passage_inplace(first_passage_mat, m_mat, s);	
	return first_passage_mat;
}
	
// where g_mat is a generalized inverse with Ge=ge, see Hunter 2011
void compute_MFPT_inplace(gsl_matrix* mfpt_mat, gsl_matrix* g_mat, gsl_vector* pi_vec)
{
	if(g_mat->size1 != g_mat->size2)
		printf("ERROR: matrix input to compute_MFPT(gsl_matrix*, gsl_vector*) is not square.\n");
	const int DIM = pi_vec->size;
	for(int i=0; i<DIM; i++)
	{
		for(int j=0; j<DIM; j++)
		{
			// inline gsl_matrix set and get
			mfpt_mat->data[i*mfpt_mat->tda+j] = (g_mat->data[j*g_mat->tda+j] - g_mat->data[i*g_mat->tda+j] + !(i-j))/(pi_vec->data[j]);
		}
	}
}

void compute_matrix_laplacian_inplace(gsl_matrix* m_mat, gsl_matrix* a_mat)
{
	// using global temp_array for temp structures
	double* temp_array = get_extern_temp_array();
	gsl_matrix_set_zero(m_mat);
	if(a_mat->size1 != a_mat->size2)
		printf("ERROR: matrix input to compute_matrix_laplacian(gsl_matrix*) is not square.\n");
	const int DIM = a_mat->size1;
	double* degree = temp_array;
	for(int i=0; i<DIM; i++)
	{
		degree[i] = 0;
		for(int j=0; j<DIM; j++)
		{
			degree[i] += a_mat->data[i*a_mat->tda+j];
		}
	}
	for(int i=0; i<DIM; i++)
	{
		for(int j=0; j<DIM; j++)
		{
			m_mat->data[i*m_mat->tda+j] = degree[i]*(!(i-j)) - a_mat->data[i*a_mat->tda+j];
		}
	}
	reset_extern_temp_array(temp_array);
}

void compute_matrix_laplacian_rw_inplace(gsl_matrix* m_mat, gsl_matrix* a_mat)
{
	double* temp_array = get_extern_temp_array();
	gsl_matrix_set_zero(m_mat);
	if(a_mat->size1 != a_mat->size2)
		printf("ERROR: matrix input to compute_matrix_laplacian(gsl_matrix*) is not square.\n");
	const int DIM = a_mat->size1;
	double* degree = temp_array;
	for(int i=0; i<DIM; i++)
	{
		degree[i] = 0;
		for(int j=0; j<DIM; j++)
		{
			degree[i] += a_mat->data[i*a_mat->tda+j];
		}
	}
	for(int i=0; i<DIM; i++)
	{
		for(int j=0; j<DIM; j++)
		{
			m_mat->data[i*m_mat->tda+j] = a_mat->data[i*a_mat->tda+j]/degree[i];
		}
	}
	reset_extern_temp_array(temp_array);
}

void compute_ever_passage_inplace(gsl_matrix* ever_passage_mat, gsl_matrix* m_mat)
{
	// ever-passage matrix is the limit of F(s) as s->0+
	double s = 1e-10;
	compute_first_passage_inplace(ever_passage_mat, m_mat, s);
}

/* calculates the first passage matrix in fourier space */
void compute_first_passage_inplace(gsl_matrix* mat, gsl_matrix* m_mat, double s)
{
	// using global variable for temp structures
	gsl_matrix* temp_mat = get_extern_temp_mat();
	
	gsl_matrix_set_zero(mat);
	if(m_mat->size1 != m_mat->size2)
		printf("ERROR: matrix input to compute_first_passage(gsl_matrix*, double) is not square.\n");
	gsl_matrix_memcpy(mat, m_mat);

	gsl_matrix* identity = temp_mat;
	gsl_matrix_set_identity(identity);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, s, identity, identity, -1.0, mat);

	gsl_matrix_inverse(temp_mat, mat);
	gsl_vector diag = gsl_matrix_diagonal(temp_mat).vector;
	gsl_matrix* diag_mat = gsl_vector_diagonal_matrix(&diag);
	gsl_blas_dtrsm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, diag_mat, temp_mat);
	gsl_matrix_set_zero(mat);
	gsl_matrix_memcpy(mat, temp_mat);
	
	// clean global extern memory
	reset_extern_temp_mat(temp_mat);
}

/* writes the elements of a gsl_matrix to a file with csv formatting */
int gsl_matrix_write(gsl_matrix* mat, const char* filename)
{
	FILE* fp = fopen(filename, "w");
	for(int i=0; i<mat->size1; i++)
	{
		for(int j=0; j<mat->size2-1; j++)
		{
			if(fprintf(fp, "%f, ", mat->data[i*mat->tda+j]) < 0)
				return -1;
		}
		if(fprintf(fp, "%f\n", mat->data[i*mat->tda+mat->size2-1]) < 0)
			return -1;
	}
	fclose(fp);
	return 1;
}