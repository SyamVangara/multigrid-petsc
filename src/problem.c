#include "problem.h"

void OpA(double *A, double *metrics, double *h) {
	//Computes the coefficients
	//
	//i - row    - y coord
	//j - column - x coord
	//A[0]*u(i,j-1) + A[1]*u(i-1,j) + A[2]*u(i,j) + A[3]*u(i+1,j) + A[4]*u(i,j+1) = f(i,j)
	//
	//metrics[5]	- metrics at a point
	//h[2]		- mesh width in computational domain in each direction
	
	double hy2, hx2;
	
	hx2 = h[0]*h[0];
	hy2 = h[1]*h[1];
	A[0] = (metrics[1]/hy2) - (metrics[3]/(2*h[1]));
	A[1] = (metrics[0]/hx2) - (metrics[2]/(2*h[0]));
	A[2] = -2.0*((metrics[0]/hx2) + (metrics[1]/hy2));
	A[3] = (metrics[0]/hx2) + (metrics[2]/(2*h[0]));
	A[4] = (metrics[1]/hy2) + (metrics[3]/(2*h[1]));
}

double Ffunc(double x, double y) {
	// Gives the f(x,y)
	
	return -2*PI*PI*sin(PI*x)*sin(PI*y);
}

double SOLfunc(double x, double y) {
	// Gives the f(x,y)
	
	return sin(PI*x)*sin(PI*y);
}

void SetUpProblem(Problem *prob) {
	// Sets up the function, solution, domain boundaries and discrete operator of the given problem type
	
	prob->Ffunc   = &Ffunc;
	prob->SOLfunc = &SOLfunc;
	prob->OpA     = &OpA;
//	for (int i=0;i<DIMENSION;i++) {
//		prob.bounds[i*2] = 0.0;
//		prob.bounds[i*2] = 0.0;
//	}
}
