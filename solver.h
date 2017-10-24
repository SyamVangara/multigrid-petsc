



#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <stdio.h>
#include <math.h>
//#include "matbuild.h"
#include "array.h"
#include <time.h>

#include "petscksp.h"
#include "mesh.h"

typedef struct {
	int		grids;   // num of grids in this level
	int		*gridId; // Grid Id of each grid in a given level
	double		(*h)[2]; // Delta h in reference domain
	int		*ranges; // ranges of global indices processes
	ArrayInt2d	global;  // global to grid map
	ArrayInt2d	*grid;   // grid to global map
} Level;

typedef struct {
	int	levels;
	int	totalGrids;
	int	coarseningFactor;
	Level	*level;
} Indices;

typedef struct {
	int	totalGrids;
	Array2d	*res;
	Array2d	*pro;
} Operator;

typedef struct {
	int	levels;
	Mat 	*res;
	Mat 	*pro;
	Mat	*A;
	Mat	*A2;
	Vec	*b;
	Vec	*u;
	IS	*bottomIS;
	IS	*topIS;
	IS	*subFineIS;
} Assembly;

typedef enum {VCYCLE, ICYCLE, ECYCLE, D1CYCLE, D2CYCLE, D3CYCLE, D4CYCLE} Cycle;
//typedef enum {False, True} CustomBool;

typedef struct {
	Cycle		cycle;
//	CustomBool	moreInfo;
	int		moreInfo; // 0: False; 1: True
	int		numIter;
	int		v[2];
	// For more info; move them to a different struct?
	int		grids;
	double		**rNormGrid;
	double		*rNormGlobal;
	// more info ends
	double		*rnorm;
	Assembly	*assem;
} Solver;

typedef struct {
	double	error[3];
	FILE	*solData;
	FILE	*errData;
	FILE	*resData;
} PostProcess;

void SetUpIndices(Mesh *mesh, Indices *indices);
void DestroyIndices(Indices *indices);
void mapping(Indices *indices, int mappingStyleflag);

void SetUpOperator(Indices *indices, Operator *op);
void DestroyOperator(Operator *op);
void GridTransferOperators(Operator op, Indices indices);

void SetUpSolver(Indices *indices, Solver *solver, Cycle c);
//void SetUpSolver(Indices *indices, Solver *solver, Cycle cyc, CustomBool moreInfo);
void DestroySolver(Solver *solver);
void Solve(Solver *solver);

void SetUpPostProcess(PostProcess *pp);
void DestroyPostProcess(PostProcess *pp);
void Postprocessing(Problem *prob, Mesh *mesh, Indices *indices, Solver *solver, PostProcess *pp);

void Assemble(Problem *prob, Mesh *mesh, Indices *indices, Operator *op, Solver *solver);

#endif

