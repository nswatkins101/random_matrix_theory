#include<stdio.h>

#include<gsl/gsl_eigen.h>
#include<gsl/gsl_permutation.h>

#include<extern_var.h>

/* external variables defined here, prefix g_ indicates a global variable */
int g_temp_num;
int g_dim;

double** g_temp_array;
int* g_temp_array_usage;

gsl_matrix** g_temp_mat;
int* g_temp_mat_usage;
gsl_matrix_complex** g_temp_cmat;
int* g_temp_cmat_usage;
gsl_vector** g_temp_vec;
int* g_temp_vec_usage;
gsl_vector_complex** g_temp_cvec;
int* g_temp_cvec_usage;

gsl_eigen_symm_workspace* g_w_symm;
gsl_eigen_nonsymm_workspace* g_w_nonsymm;
gsl_permutation* g_permutation;

gsl_matrix* identity_mat;

int init_extern_var(int DIM)
{
	g_temp_num = 4;
	g_temp_array_usage = calloc(g_temp_num, sizeof(int));
	g_temp_mat_usage = calloc(g_temp_num, sizeof(int));
	g_temp_cmat_usage = calloc(g_temp_num, sizeof(int));
	g_temp_vec_usage = calloc(g_temp_num, sizeof(int));
	g_temp_cvec_usage = calloc(g_temp_num, sizeof(int));
	g_temp_array = calloc(g_temp_num, sizeof(void*));
	for(int i=0; i<g_temp_num; i++)
		g_temp_array[i] = calloc(sizeof(double), DIM);
	g_temp_mat = calloc(g_temp_num, sizeof(void*));
	for(int i=0; i<g_temp_num; i++)
		g_temp_mat[i] = gsl_matrix_alloc(DIM, DIM);
	g_temp_cmat = calloc(g_temp_num, sizeof(void*));
	for(int i=0; i<g_temp_num; i++)
		g_temp_cmat[i] = gsl_matrix_complex_alloc(DIM, DIM);
	g_temp_vec = calloc(g_temp_num, sizeof(void*));
	for(int i=0; i<g_temp_num; i++)
		g_temp_vec[i] = gsl_vector_alloc(DIM);
	g_temp_cvec = calloc(g_temp_num, sizeof(void*));
	for(int i=0; i<g_temp_num; i++)
		g_temp_cvec[i] = gsl_vector_complex_alloc(DIM);
		
	g_w_symm = gsl_eigen_symm_alloc(DIM);
	g_w_nonsymm = gsl_eigen_nonsymm_alloc(DIM);
	g_permutation = gsl_permutation_alloc(DIM);
	
	identity_mat = gsl_matrix_alloc(DIM, DIM);
	gsl_matrix_set_identity(identity_mat);
	return 0;
}

int free_extern_var()
{
	free(g_temp_mat_usage);
	for(int i=0; i<g_temp_num; i++)
		gsl_matrix_free(g_temp_mat[i]);
	for(int i=0; i<g_temp_num; i++)
		gsl_matrix_complex_free(g_temp_cmat[i]);
	for(int i=0; i<g_temp_num; i++)
		gsl_vector_free(g_temp_vec[i]);
	for(int i=0; i<g_temp_num; i++)
		gsl_vector_complex_free(g_temp_cvec[i]);
	gsl_eigen_symm_free(g_w_symm);
	gsl_eigen_nonsymm_free(g_w_nonsymm);
	gsl_permutation_free(g_permutation);
	
	gsl_matrix_free(identity_mat);
	
	return 0;
}

double* get_extern_temp_array()
{
	double* temp_array = NULL;
	for(int i=0; i<g_temp_num; i++)
	{
		if(!g_temp_array_usage[i])
		{
			temp_array = g_temp_array[i];
			g_temp_array_usage[i] = 1;
			break;
		}
	}
	if(!temp_array)
		printf("ERROR: not enough temporary matrices.\n");
	return temp_array;
}

int reset_extern_temp_array(double* temp_array)
{
	for(int i=0; i<g_temp_num; i++)
	{
		if(temp_array == g_temp_array[i])
		{
			for(int j=0; j<g_temp_num; j++)
				temp_array[j] = 0;
			g_temp_array_usage[i] = 0;
			temp_array = NULL;
			break;
		}
	}
	if(temp_array)
	{
		printf("ERROR: cannot reset temp_mat, could not find a match for the address.\n");
		return -1;
	}
	return 0;
}

gsl_matrix* get_extern_temp_mat()
{
	gsl_matrix* temp_mat = NULL;
	for(int i=0; i<g_temp_num; i++)
	{
		if(!g_temp_mat_usage[i])
		{
			temp_mat = g_temp_mat[i];
			g_temp_mat_usage[i] = 1;
			break;
		}
	}
	if(!temp_mat)
		printf("ERROR: not enough temporary matrices.\n");
	return temp_mat;
}

int reset_extern_temp_mat(gsl_matrix* temp_mat)
{
	for(int i=0; i<g_temp_num; i++)
	{
		if(temp_mat == g_temp_mat[i])
		{
			gsl_matrix_set_zero(temp_mat);
			g_temp_mat_usage[i] = 0;
			temp_mat = NULL;
			break;
		}
	}
	if(temp_mat)
	{
		printf("ERROR: cannot reset temp_mat, could not find a match for the address.\n");
		return -1;
	}
	return 0;
}

gsl_matrix_complex* get_extern_temp_cmat()
{
	gsl_matrix_complex* temp_cmat = NULL;
	for(int i=0; i<g_temp_num; i++)
	{
		if(!g_temp_cmat_usage[i])
		{
			temp_cmat = g_temp_cmat[i];
			g_temp_cmat_usage[i] = 1;
			break;
		}
	}
	if(!temp_cmat)
		printf("ERROR: not enough temporary complex matrices.\n");
	return temp_cmat;
}

int reset_extern_temp_cmat(gsl_matrix_complex* temp_cmat)
{
	for(int i=0; i<g_temp_num; i++)
	{
		if(temp_cmat == g_temp_cmat[i])
		{
			gsl_matrix_complex_set_zero(temp_cmat);
			g_temp_cmat_usage[i] = 0;
			temp_cmat = NULL;
			break;
		}
	}
	if(temp_cmat)
	{
		printf("ERROR: cannot reset temp_cmat, could not find a match for the address.\n");
		return -1;
	}
	return 0;
}

gsl_vector* get_extern_temp_vec()
{
	gsl_vector* temp_vec = NULL;
	for(int i=0; i<g_temp_num; i++)
	{
		if(!g_temp_vec_usage[i])
		{
			temp_vec = g_temp_vec[i];
			g_temp_vec_usage[i] = 1;
			break;
		}
	}
	if(!temp_vec)
		printf("ERROR: not enough temporary vectors.\n");
	return temp_vec;
}

int reset_extern_temp_vec(gsl_vector* temp_vec)
{
	for(int i=0; i<g_temp_num; i++)
	{
		if(temp_vec == g_temp_vec[i])
		{
			gsl_vector_set_zero(temp_vec);
			g_temp_vec_usage[i] = 0;
			temp_vec = NULL;
			break;
		}
	}
	if(temp_vec)
	{
		printf("ERROR: cannot reset temp_vec, could not find a match for the address.\n");
		return -1;
	}
	return 0;
}

gsl_vector_complex* get_extern_temp_cvec()
{
	gsl_vector_complex* temp_cvec = NULL;
	for(int i=0; i<g_temp_num; i++)
	{
		if(!g_temp_cvec_usage[i])
		{
			temp_cvec = g_temp_cvec[i];
			g_temp_cvec_usage[i] = 1;
			break;
		}
	}
	if(!temp_cvec)
		printf("ERROR: not enough temporary complex vectors.\n");
	return temp_cvec;
}

int reset_extern_temp_cvec(gsl_vector_complex* temp_cvec)
{
	for(int i=0; i<g_temp_num; i++)
	{
		if(temp_cvec == g_temp_cvec[i])
		{
			gsl_vector_complex_set_zero(temp_cvec);
			g_temp_cvec_usage[i] = 0;
			temp_cvec = NULL;
			break;
		}
	}
	if(temp_cvec)
	{
		printf("ERROR: cannot reset temp_vec, could not find a match for the address.\n");
		return -1;
	}
	return 0;
}