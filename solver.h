



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
	double		*rnorm;
	int		(*range)[2];
} Solver;

typedef struct {
	double	error[3];
	FILE	*solData;
	FILE	*errData;
	FILE	*resData;
} PostProcess;

void SetUpSolver(Indices *indices, Solver *solver, Cycle c);
void DestroySolver(Solver *solver);
void SetPostProcess(PostProcess *pp);
void DestroyPostProcess(PostProcess *pp);
void Postprocessing(Problem *prob, Mesh *mesh, Indices *indices, Assembly *assem, Solver *solver, PostProcess *pp);

void UpdateRHS(double *A, double **u, double **r, int *n);
double Residual(double **u, double **f, double **r, double *As, int *n);
void JacobiStep(double **u, double **f, double *As, double w, int *n);
void Jacobi(double **u, double **f, double **r, double *As, double w, double *rnorm, int v,int *n);
void ResidualRestriction(double **f, double **r, int *n);
void ErrorCorrection(double **u, int *n, int flag);
void Vcycle(double **u, double **f, double **r, double *As, double w, int *v,int levels,int *n);
void Multigrid(double **u, double **f, double **r, double *As, double w, double *rnorm, int levels, int *n,int m);
void MultigridPetsc(Array2d u, Array2d metrics, double *f, double **opIH2h, double **opIh2H, ArrayInt2d *IsStencil, ArrayInt2d *IsResStencil, ArrayInt2d *IsProStencil, IsRange *range, double *rnorm, int levels, int *fulln, int *m);
//void MultigridPetsc(Array2d u, Array2d metrics, double *f, double **opIH2h, double **opIh2H, ArrayInt2d *IsStencil, double *rnorm, int levels, int *fulln, int *m);
void PMultigrid(double **u, double **f, double **r, double *As, double w, double *rnorm, int levels, int*n, int m);
//void AsyncMultigrid(double **u, double **f, double **r, double *As, double w, double *rnorm, int*n, int m);
double L2norm(double *a, int n);
double L1Norm(double *a, int n);
double LiNorm(double *a, int n);
double norm(double *a, int n);
void Initialization(double **u, int *n);
int JacobiMalloc(double ***f, double ***u, double ***r, int *n);
int MultigridMalloc(double ***f, double ***u, double ***r, int *n, int levels);
int AsyncMultigridMalloc(double ***f, double ***u, double ***r,int *n, int levels);

#endif

