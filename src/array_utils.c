#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include<array_utils.h>

/* allocates a nested array, i.e. a 2d array */
void** array2d_calloc(int m, int n, size_t size)
 {
	void** rtn = calloc(m, sizeof(void*));
	for(int i=0; i<m; i++)	
		rtn[i] = calloc(n,size);
	return rtn;
}

 /* frees a nested array, i.e. a 2d array */
void array2d_free(void* arr, int m, int n)
{
	for(int i=0; i<m; i++)
		free(((void**)arr)[i]);
	free(arr);
	arr = NULL;
}

/* allocates an array and copies a 2d array onto it, concentates rows */
double* array2d_flatten(double** arr, int m, int n)
{
	double* rtn = calloc(m*n, sizeof(double*));
	for(int i=0; i<m; i++)
		memcpy(rtn+i, arr[i], n*sizeof(double));
	// destroys the input array
	array2d_free(arr, m, n);
	arr = NULL;
	return rtn;
}

/* finds the tuple in a 2d array for which y is minimized */
void array_min(double* min, int size, double* arr)
{
	min[0] = 0;
	min[1] = 0;
	for(int i=0; i<size; i++)
	{
		if(arr[2*i+1] <min[1])
		{
			min[0] = arr[2*i];
			min[1] = arr[2*i+1];
		}
	}
}

/* writes an array into a file with formatting according to delimiter
 * the first element of the record is the size, followed by the elements of the list
 */
 int array_fprintf_delim(FILE* stream, double* data, const int size, char* delimiter)
{
	// write the parameters of the data set
	//fprintf(stream, "%i%s", size, delimiter);
	// write the elements of the list
	for(int i=0; i<size-1; i++)
	{
		if(fprintf(stream, "%f%s", data[i], delimiter) < 0)
			return -1;
	}
	// write the final element of the list and the newline char
	if(fprintf(stream, "%f\n", data[size-1]))
		return -1;
	return 0;
}

/* writes an array formatted as a 2d array into a file with csv formatting */
int array2d_fprintf_delim(FILE* stream, double* data, const int size1, int size2, int tda, char* delimiter)
{
	// first write the parameters of the data set
	//fprintf(stream, "2D ARRAY. size1: %i%ssize2: %i%stda: %i\n", size1, delimiter, size2, delimiter, tda);
	for(int i=0; i<size1; i++)
	{
		for(int j=0; j<size2-1; j++)
		{
			if(fprintf(stream, "%f%s", data[i*tda+j], delimiter) < 0)
				return -1;
		}
		if(fprintf(stream, "%f\n", data[i*tda+size2-1]) < 0)
			return -1;
	}
	return 0;
}