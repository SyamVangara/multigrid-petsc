#include "array.h"

#define SUM_ARRAY(result, a, size ) {int i; result=0; for (i=0;i<size;i++)\
	{result += a[i];}}

void CreateArrayInt2d(int ni, int nj, ArrayInt2d *a) {
	//Allocates memory to a logical integer 2d array		
	a->ni = ni;
	a->nj = nj;
	a->data = malloc(a->ni*a->nj*sizeof(int));
}

void DeleteArrayInt2d(ArrayInt2d *a) {
	// Safely delete a logical 2d integer array
	
	if (a->data != NULL) free(a->data);
}

int malloc2d(double ***a, int n, int m) {
	// Allocate contiguous memory for 2d array "A" of (nxm) 
	//
	// Pointer to the pointer of array - "&A" (= a)
	// Number of rows (n), columns - m; 
	// Location of function call - fle_name, line_num
	//
	// Array "A" can be accessed using A[i][j]                           

	*a = (double **)malloc(n*sizeof(double *));
	**a = (double *)malloc(n*m*sizeof(double));
	if (*a==NULL||**a==NULL) return 1;
	
	// Assign a pointer to each row
	for(int i=1;i<n;i++){
		*(*a+i) = *(*a+i-1) + m;
	}

	return 0;
}

int malloc3d(double ****a, int p, int q, int r) {
	// Allocate contiguous memory for 3d array "A" of (pxqxr) 
	//
	// Pointer to the pointer of array - "&A" (= a) 
	// Location of function call - fle_name, line_num
	//
	// Array "A" can be accessed using A[i][j]                           

	*a = (double ***)malloc(p*sizeof(double **));
	**a = (double **)malloc(p*q*sizeof(double *));
	***a = (double *)malloc(p*q*r*sizeof(double));
	if (*a==NULL||**a==NULL || ***a==NULL) return 1;
	
	// Assign a pointers
//	for(int i=1;i<p;i++){
//		*(*a+i) = *(*a+i-1) + q;
//	}
//	for (int i=0;i<p;i++) {
//		for(int j=1;j<q;j++) {
//			*(*(*a+i)+j) = (*a)[i][j-1] + r;
//		}
//	}

	for (int i=1;i<p*q;i++) {
		*(**a+i) = *(**a+i-1) + r;
	}
	for (int i=1;i<p;i++) {
		*(*a+i) = *(*a+i-1) + q;
	}
	return 0;
}

int malloc2dY(double ***a, int n, int *m) {
	// Allocate contiguous memory for 2d array "A" with variable row length (nxm(i))
	//
	// Pointer to the pointer of array - "&A" (= a)
	// Number of rows (n), columns - m(i) in i^th row 
	// Location of function call - fle_name, line_num
	//
	// Array "A" can be accessed using A[i][j]                           
	
	int aTotal;

	// aTotal - total number of elements in A
	aTotal=0; 
	for (int i=0;i<n;i++) aTotal += m[i];

	//SUM_ARRAY(aTotal, m, n); // aTotal - total number of elements in A
	
	*a = (double **)malloc(n*sizeof(double *));
	**a = (double *)malloc(aTotal*sizeof(double));
	if (*a==NULL||**a==NULL) return 1;
	
	// Assign a pointer to each row
	for(int i=1;i<n;i++){
		*(*a+i) = *(*a+i-1) + *(m+i-1);
	}

	return 0;
}

void free2dArray(double ***a) {
	//free the allocated memory
	free(**a);
	free(*a);
}

void free3dArray(double ****a) {
	//free the allocated memory
	free(***a);
	free(**a);
	free(*a);
}

