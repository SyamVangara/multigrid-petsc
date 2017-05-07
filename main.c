#include <math.h>
#include "petscksp.h"

#define PI 3.14159265358979323846

double func(double x, double y);
Mat matrixA(int n);
Vec vecb(int n, double (*f)(double, double));

int main(int argc, char *argv[]) {
	
	KSP	solver;
	Mat	A;
	Vec	b, x;
	int	iters, n;

	PetscInitialize(&argc, &argv, 0, 0);
	n = 10;
	//PetscOptionsGetInt(PETSC_NULL, "-n", &n, 0);
	PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);

	A = matrixA(n);
	b = vecb(n, func);
	VecDuplicate(b, &x);
	KSPCreate(PETSC_COMM_WORLD, &solver);
	KSPSetOperators(solver, A, A);
	KSPSetTolerances(solver, 1.e-15, 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT);
	KSPSetFromOptions(solver);
	KSPSolve(solver, b, x);
	KSPGetIterationNumber(solver, &iters);
	PetscPrintf(PETSC_COMM_WORLD, "Solution in %d iterations is: \n", iters);
	//VecView(x, PETSC_VIEWER_STDOUT_WORLD);

	MatDestroy(&A); VecDestroy(&b); VecDestroy(&x);
	KSPDestroy(&solver);
	PetscFinalize();
	return 0;
}

double func(double x, double y) {
	return -2.0*PI*PI*sin(PI*x)*sin(PI*y);
}

Mat matrixA(int n) {
	Mat	A;
	int	r, rowStart, rowEnd, i, j;
	double	h, invh2;

	h	= 1.0/(n+1);
	invh2	= 1.0/(h*h);

	MatCreate(PETSC_COMM_WORLD, &A);
	MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n*n, n*n);
	MatSetFromOptions(A);
	MatSetUp(A);
	//MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);
	MatGetOwnershipRange(A, &rowStart, &rowEnd);

	for (r=rowStart; r<rowEnd; r++) {
		i = r%n; j = r/n;
		if (j-1>0) {
			MatSetValue(A, r, r-n, invh2, INSERT_VALUES);
		}
		if (i-1>0) {
			MatSetValue(A, r, r-1, invh2, INSERT_VALUES); 
		}
		MatSetValue(A, r, r, -4*invh2, INSERT_VALUES);
		if (i+1<n-1) {
			MatSetValue(A, r, r+1, invh2, INSERT_VALUES);
		}
		if (j+1<n-1) {
			MatSetValue(A, r, r+n, invh2, INSERT_VALUES);
		}
	}
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	return A;
}

Vec vecb(int n, double (*f)(double, double)) {
	Vec	V;
	int	r, rowStart, rowEnd, i, j;
	double	h;

	h = 1.0/(n+1);
	VecCreate(PETSC_COMM_WORLD, &V);
	VecSetSizes(V, PETSC_DECIDE, n*n);
	VecSetFromOptions(V);
	VecGetOwnershipRange(V, &rowStart, &rowEnd);
	for (r=rowStart; r<rowEnd; r++) {
		i = (r%n)+1;
		j = (r/n)+1;
		VecSetValue(V, r, (*f)(i*h, j*h), INSERT_VALUES);
	}
	VecAssemblyBegin(V);
	VecAssemblyEnd(V);

	return V;
}
