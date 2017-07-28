/*


*/


#ifndef _ARRAY_H_
#define _ARRAY_H_ 

#include <stdlib.h>
// Structures

typedef struct {
	int ni;
	int nj;
	int* data;
} Array2d;

typedef struct {
	int ni;
	int nj;
	int nk;
	double* data;
} Array3d;

// Function prototypes

// Allocate contiguous memory for 2d array
extern int malloc2d(double ***a, int n, int m);

// Allocate contiguous memory for 3d array
extern int malloc3d(double ****a, int p, int q, int r);

// Allocate contigous memory for 2d array with variable row length
extern int malloc2dY(double ***a, int n, int *m);

// Free the allocated memory of 2d/3d arrays
extern void free2dArray(double ***a);
extern void free3dArray(double ****a);

#endif
