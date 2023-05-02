/* Nat Watkins 02/01/2022. 
 * This defines several different minimization functions for use by the model struct.
 */

#include<stdio.h>

#include<gsl/gsl_math.h>
#include<gsl/gsl_block.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

#include<matrix_utils.h>
#include<matrix_calc_utils.h>
#include<array_utils.h>
#include<graph_alloc_utils.h>
#include<graph_utils.h>
#include<random_matrix_utils.h>
#include<extern_var.h>
#include<opt_utils.h>

#include<model.h>
#include<model_min.h>

void model_minimize_global(MODEL* model)
{
	/* minimizes the free energy wrt the internal parameters at the current point in phase space
	 * (i.e. for a set of eternal parameters) by calculating and then comparing the value of the FE at all 
	 * points in the internal space. This is a simplest global minimization technique. 
	 */
	 
	 /* a clearer method might be to create a multidimensional array of coordinates (a mesh), 
	  * then flatten that array, then iterate through that array
	  */
	  // gsl_vector* span is the span of the internal parameter space in its units
	 gsl_vector* span = gsl_vector_alloc(model->num_intern_params);
	 gsl_vector_memcpy(span, model->intern_params_upper_boundary);
	 gsl_vector_sub(span, model->intern_params_lower_boundary);
	 
	 // gsl_vector* coord are coordinates in the internal paramter space in its units
	 gsl_vector* coord = gsl_vector_alloc(model->num_intern_params);
	 
	 // gsl_vector_int* multid_ind is a vector of integers containing the indices for a point in a multidimensional array
	 gsl_vector_int* multid_indices = gsl_vector_int_calloc(model->num_intern_params);
	 
	 // allocate memory for free energy values
	 double* flat_data = malloc(sizeof(double)*size);
	 
	 for(int i=0; 
	 
	 
	  
	 // first define an internal parameter space mesh indeed by i
	 const int size1 = 10;															// number of points along each axis in the mesh, in my style I call size1,size2,... the number of rows, columns, etc. Size is the total number of elements
	 int size = pow(size1, model->num_intern_params);				// total number of points
	 // scale and shift mesh
	 // gsl_vector* span is the span of the internal parameter space in its units
	 gsl_vector* span = gsl_vector_alloc(model->num_intern_params);
	 gsl_vector_memcpy(span, model->intern_params_upper_boundary);
	 gsl_vector_sub(span, model->intern_params_lower_boundary);
	 
	// gsl_vector* coord are coordinates in the internal paramter space in its units
	 gsl_vector* coord = gsl_vector_alloc(model->num_intern_params);
	 
	 // gsl_vector_int* multid_ind is a vector of integers containing the indices for a point in a multidimensional array
	 gsl_vector_int* multid_indices = gsl_vector_int_calloc(model->num_intern_params);
	 
	 // allocate memory for free energy values
	 double* flat_data = malloc(sizeof(double)*size);
	 
	 // rasters across set of internal parameters while keeping eternal parameters fied
	 // data contains the FEs evaluated on a mesh as a 1d array
	 for(int index=0; index<size; index++)
	 {
		 // casts the index a c-array (index) to the inde of a multidimensional array (multid_indices)
		 index_to_multid_indices(multid_indices, index, size1);
		 // casts the index of the multidimensional array (multid_indices) to a coordinate of the internal parameter space in its units (coord)
		 gsl_vector_double_int(coord, multid_indices);
		 gsl_vector_scale(coord, 1.0/size1);
		 gsl_vector_mul(coord, span);
		 gsl_vector_add(coord, model->intern_params_lower_boundary);
		 gsl_vector_memcpy(model->intern_params, coord);
		 // calculates and returns the FE at the current point in parameter space set above
		 flat_data[index] = model->mean_free_energy(model);
	}
	 int index_min = minimize_global(flat_data, size);
	 // casts the index (of a c-array) of the minimum (index_min) to the indices of a multidimensional array (multid_indices)
	 index_to_multid_indices(multid_indices, index_min, size1);
	 // casts the index of the multidimensional array (multid_indices) to a coordinate of the internal parameter space in its units (coord)
	 gsl_vector_double_int(coord, multid_indices);
	gsl_vector_scale(coord, 1.0/size1);
	gsl_vector_mul(coord, span);
	gsl_vector_add(coord, model->intern_params_lower_boundary);
	// the result of the function is to set the internal params to the value that minimizes FE	 
	gsl_vector_memcpy(model->intern_params, coord);	 
}

/* Calculates the inde of a multidimensional array cast as a c-array
 * gsl_vector* coord points a struct containing the multidimensional indices
 * size1 is the number of rows, columns, etc. Assumes a hypercubic mulitdimensional array
 */
int multid_indices_to_index(gsl_vector_int* multid_indices, int size1)
{
	int index = 0;
	int dim = multid_indices->size;
	for(int i=0; i<dim; i++)
	{
		index += gsl_vector_int_get(multid_indices, i)*pow(size1, dim-1-i);
	}
	return index;
}

/* Calculates the coords of a c-array cast as a multidimensional array
 * gsl_vector* coord points a struct containing the multidimensional indices, assumes coord is empty and fills with the indices for a multidim array
 * int ind is the inde of the c-array
 * size1 is the number of elements along each ais, i.e. assumes mapping an element of an interval to a hypercube, i.e. [0, len^dim-1] -> [0,len-1]^dim
 */
void index_to_multid_indices(gsl_vector_int* multid_indices, int index, int size1)
{ 
	int k = 0;
	int dim = multid_indices->size;
	for(int i=dim-1; i>=0; i--)
	{
		k = index/pow(size1,i);
		index  -= k*pow(size1, i);		
		gsl_vector_int_set(multid_indices, i, k);
	}
}

/*
 * double* data is a c-array containing the data, finds the smallest element of data
 * int size is the number of elements in double* data
 */
int minimize_global(double* data, int size)
{
	// initially assumes first element is the minimum
	int index_min = 0;
	for(int i=0; i<size; i++)
	{
		if(data[i] < data[index_min])
			index_min = i;
	}
	return index_min;
}

void model_minimize_grad_descent(MODEL* model)
{
	// minimizes the free energy wrt the internal parameters at the current position
	// phase diagram does minimize over all external parameters, i.e. generates the phase diagram
	
	/* the energy landscape is determined by the particular form of the generated matrices. We want
	* average over the ensemble when calculating the critical internal parameters.
	*/
	const int MAX_NUM_STEPS = 50;
	const double FTOL = 0.0001;
	const double dx = 0.1;

	for(int i=0; i<MAX_NUM_STEPS; i++)
	{
		// calculate the derivative at point model->pos
		model->set_grad_free_energy(model);
		if(gsl_blas_dnrm2(model->grad) < FTOL*FTOL)
		{
			//printf("reached FTOL\n");
			break;
		}
		// updates the intern_params at the point j in phase space, moves down the gradient
		//printf("dp is %f\n", dx*model->grad->data[0]);
		gsl_blas_daxpy(-dx, model->grad, model->intern_params);
		// check that intern_params is still in the boundary
		if(!model->params_in_boundary(model))
		{
			//printf("ERROR: model_minimize: out of bounds\n");
			// if not in boundary, reverse update step and break
			gsl_blas_daxpy(dx, model->grad, model->intern_params);
			break;
		}
	}
}

/* Casts/Copies a integer valued vector into a double valued vector
 */
void gsl_vector_double_int(gsl_vector* dest, gsl_vector_int* src)
{
	for(int i=0; i<dest->size; i++)
	{
		dest->data[i] = (double)(src->data[i]);
	}	
}

/* Casts/Copies a double valued vector into a integer valued vector
 */
void gsl_vector_int_double(gsl_vector_int* dest, gsl_vector* src)
{
	for(int i=0; i<dest->size; i++)
	{
		dest->data[i] = (int)(src->data[i]);
	}	
}