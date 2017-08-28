



#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <stdio.h>
#include <math.h>
#include "matbuild.h"
#include "array.h"
#include <time.h>

typedef enum {VCYCLE, ICYCLE} Cycle;

typedef struct {
	Cycle		cycle;
	int		numIter;
	int		v[2];
	double		*rnorm;
} Solver;

typedef struct {
	double	error[3];
	FILE	*solData;
	FILE	*errData;
	FILE	*resData;
} PostProcess;

void SetUpSolver(Indices *indices, Solver *solver, Cycle c);
void DestroySolver(Solver *solver);
void Solve(Assembly *assem, Solver *solver);

void SetUpPostProcess(PostProcess *pp);
void DestroyPostProcess(PostProcess *pp);
void Postprocessing(Problem *prob, Mesh *mesh, Indices *indices, Assembly *assem, Solver *solver, PostProcess *pp);

#endif

