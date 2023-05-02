#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include"array_utils.h"
#include"matrix_utils.h"

int gsl_vector_fprintf_delim(FILE* stream, gsl_vector* vector, char* delimiter)
{
	return array_fprintf_delim(stream, vector->data, vector->size, delimiter);
}

int gsl_matrix_fprintf_delim(FILE* stream, gsl_matrix* matrix, char* delimiter)
{
	return array2d_fprintf_delim(stream, matrix->data, matrix->size1, matrix->size2, matrix->tda, delimiter);
}

int gsl_matrix_complex_fprintf_delim(FILE* stream, gsl_matrix_complex* matrix, char* delimiter)
{
	int size1 = matrix->size1;
	int size2 = matrix->size2;
	gsl_complex z;
	for(int i=0; i<size1; i++)
	{
		for(int j=0; j<size2-1; j++)
		{
			z = gsl_matrix_complex_get(matrix, i, j);
			if(fprintf(stream, "%g+%gi%s", GSL_REAL(z), GSL_IMAG(z), delimiter) != 0)
				return -1;
		}
		z = gsl_matrix_complex_get(matrix, i, size2-1);
		if(fprintf(stream, "%g+%gi\n", GSL_REAL(z), GSL_IMAG(z)) != 0)
			return -1;
	}
	return 0;
}

// C "overload" doesn't exist, style is to append _type to end of function
void gsl_matrix_set_ptr(gsl_matrix* m, const double* x)
{
	gsl_matrix_view view_array = gsl_matrix_view_array(x, m->size1, m->size2);
	//print_matrix(&(view_array.matrix));
	//printf("\n");
	gsl_matrix_memcpy(m, &(view_array.matrix));
	//print_matrix(m);
	//printf("\n");	
}

// opposite of the above function, i.e. sets an array based on a matrix
void ptr_set_gsl_matrix(double* x, const gsl_matrix* m)
{
	for(int i=0; i<m->size1; i++)
	{
		for(int j=0; j<m->size2; j++)
		{
			x[i*m->size1+j] = gsl_matrix_get(m, i, j);
		}
	}
}

// input is a vector, allocs and returns a matrix with those vector elements on the diagonal
gsl_matrix* gsl_vector_diagonal_matrix(gsl_vector* x)
{
	gsl_matrix* mat = gsl_matrix_alloc(x->size, x->size);
	gsl_vector_view diag = gsl_matrix_diagonal(mat);
	gsl_matrix_set_all(mat, 0.0);
	gsl_vector_memcpy(&diag.vector, x);
	return mat;
}

void gsl_matrix_complex_real_memcpy(gsl_matrix* mat, gsl_matrix_complex* z_mat)
{
	gsl_complex z;
	for(int i=0; i< mat->size1; i++)
	{
		for(int j=0; j<(mat->size2); j++)
		{
			z = gsl_matrix_complex_get(z_mat, i, j);
			gsl_matrix_set(mat, i, j, GSL_REAL(z));
		}
	}
}

void gsl_matrix_complex_imag_memcpy(gsl_matrix* mat, gsl_matrix_complex* z_mat)
{
	gsl_complex z;
	for(int i=0; i< mat->size1; i++)
	{
		for(int j=0; j<(mat->size2); j++)
		{
			z = gsl_matrix_complex_get(z_mat, i, j);
			gsl_matrix_set(mat, i, j, GSL_IMAG(z));
		}
	}
}

gsl_matrix_view gsl_matrix_complex_real(gsl_matrix_complex* z)
{
	double array[(z->size1)*(z->size2)];
	for(int i=0; i<z->size1; i++)
	{
		for(int j=0; j<(z->size2); j++)
			array[i*(z->size2)+j] = z->data[2*(i*(z->tda)+j)];
	}
	return gsl_matrix_view_array(array, z->size1, z->size2);
}

gsl_matrix_view gsl_matrix_complex_imag(gsl_matrix_complex* z)
{
	double array[(z->size1)*(z->size2)];
	for(int i=0; i<z->size1; i++)
	{
		for(int j=0; j<(z->size2); j++)
			array[i*(z->size2)+j] = z->data[2*(i*(z->tda)+j)+1];
	}
	return gsl_matrix_view_array(array, z->size1, z->size2);
}
