#include <math.h>
#include "petscksp.h"

#define PI 3.14159265358979323846

double func(double x, double y);

int main(int argc, char *argv[]) {
	
	KSP	solver;
	Mat	A;
	Vec	b, x;
	int	iters, n;

	PetscInitialize(&argc, &argv, 0, 0);
	n = 10;
	PetscOptionsGetInt(PETSC_NULL, "-n", &n, 0);

	A = matrixA();
	b = vecb();
	VecDuplicate(b, &x);
	KSPCreate(PETSC_COMM_WORLD, &solver);
	KSPSetOperators(solver, A, A, DIFFERENT_NONZERO_PATTERN);
	KSPSetFromOptions(solver);
	KSPSolve(solver, b, x);
	KSPGetIterationNumber(solver, &iters);
	PetscPrintf(PETSC_COMM_WORLD, "Solution in %d iterations is: \n");
}

double func(double x, double y) {
	return -2.0*PI*PI*sin(PI*x)*sin(PI*y);
}


