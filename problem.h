



#ifndef _PROBLEM_H_
#define _PROBLEM_H_

#include <stdio.h>
#include <math.h>

#include "array.h"

#define PI 3.14159265358979323846
#define DIMENSION 2

//#define FUNC(i,j) (-2*PI*PI*sin(PI*coord[0][(j)])*sin(PI*coord[1][(i)]))
//#define SOL(i,j) (sin(PI*coord[0][(j)])*sin(PI*coord[1][(i)]))

typedef struct {
	double (*Ffunc)(double x, double y);
	double (*SOLfunc)(double x, double y);
	void   (*OpA)(double *A, double *metrics, double *h);
//	double bounds[DIMENSION*2];
} Problem;

void SetUpProblem(Problem *prob);
#endif
