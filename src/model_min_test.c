/* Nat Watkins 02/01/2022. 
 * This program tests the program model_min independently of the rest of the project.
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include<gsl/gsl_math.h>
#include<gsl/gsl_block.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

#include<array_utils.h>
#include<model.h>
#include<model_min.h>

int main()
{
	// first try using the methods from model_min without the model part
	double data[] = {1, 2, 3, 4, 3, 2, 1};
	printf("Minimum is at %i", minimize_global(data, 7));
	
	// test the index packing methods
	double multi[] = {{0, 1, 2}, {3, 4, 5}};
	
	return 0;
}