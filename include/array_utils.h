/* allocates a nested array, i.e. a 2d array */
void** array2d_calloc(int m, int n, size_t size);

/* frees an nested array, i.e. a 2d array */
void array2d_free(void* arr, int m, int n);

/* allocates an array and copies a 2d array onto it, concentates rows */
double* array2d_flatten(double** arr, int m, int n);

/* finds the tuple in a 2d array for which y is minimized */
void array_min(double* min, int size, double* arr);

/* writes an array into a file with formatting according to delimiter
 * the first element of the record is the size, followed by the elements of the list
 */
 int array_fprintf_delim(FILE* stream, double* data, const int size, char* delimiter);
 
/* writes an array formatted as a 2d array into a file with csv formatting 
 * tda == size2 typically
 */
int array2d_fprintf_delim(FILE* stream, double* data, const int size1, int size2, int tda, char* delimiter);
