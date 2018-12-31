



#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <stdio.h>
#include <math.h>
#include <omp.h>
//#include "matbuild.h"
#include "array.h"
#include <time.h>

#include "petscksp.h"
#include "mesh.h"

typedef struct {
	int		ngrids;   // num of grids in this level
	int		*gridId; // Grid Id of each grid in a given level
	double		(*h)[2]; // Delta h in reference domain // ! Remove
	long int	(*ranges)[2]; // ranges of global indices for each grid
	ArrayInt2d	global;  // global to grid map// ! Remove
	ArrayInt2d	*grid;   // grid to global map// ! Remove
} Level;

typedef struct {
	int	nlevels;
	int	totalGrids; // ! Remove
	int	coarseningFactor; // ! Remove
	Mat 	*res;
	Mat 	*pro;
	Mat	*A;
	Vec	*b;
	Vec	*u;
	Level	*level;
} Levels;

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
	int	moreInfo; // 0:False; 1:True
//	IS	*subFineIS;
	IS	**gridIS; // moreNorm flag related info begins
} Assembly;

//typedef enum {VCYCLE, ICYCLE, ECYCLE, D1CYCLE, D2CYCLE, D3CYCLE, D4CYCLE, D1PSCYCLE, PetscPCMG, FILTER, VFILTER, ADDITIVE, ADDITIVEScaled} Cycle;
//typedef enum {False, True} CustomBool;

typedef struct {
	int		cycle;
//	Cycle		cycle;
//	CustomBool	moreInfo;
	int		moreInfo; // 0: False; 1: True
	int		numIter;
	int		v[2];
	// For more info; move them to a different struct?
	int		grids;
	double		**rNormGrid; // !Remove
	double		*rNormGlobal; // !Remove
	// more info ends
	double		*rnorm;
	Levels		*levels;
//	Assembly	*assem;
} Solver;

typedef struct {
	double	error[3];
	FILE	*solData;
	FILE	*errData;
	FILE	*resData;
	FILE	*XgridData;
	FILE	*YgridData;
} PostProcess;

int CreateLevels(Grids *grids, Levels *levels);
void DestroyLevels(Levels *levels);

//int CreateIndices(Grids *grids, Indices *indices);
////void SetUpIndices(Mesh *mesh, Indices *indices);
//void DestroyIndices(Indices *indices);
//void mapping(Indices *indices, int mappingStyleflag);

int CreateOperator(Grids *grids, Operator *op);
//void SetUpOperator(Indices *indices, Operator *op);
void DestroyOperator(Operator *op);
void GridTransferOperators(Operator op, Levels levels);

int CreateSolver(Grids *grids, Solver *solver); 
//void SetUpSolver(Indices *indices, Solver *solver, Cycle c);
//void SetUpSolver(Indices *indices, Solver *solver, Cycle cyc, CustomBool moreInfo);
void DestroySolver(Solver *solver);
int Solve(Solver *solver);

void SetUpPostProcess(PostProcess *pp);
void DestroyPostProcess(PostProcess *pp);
void Postprocessing(Problem *prob, Mesh *mesh, Levels *levels, Solver *solver, PostProcess *pp);

void Assemble(Problem *prob, Mesh *mesh, Levels *levels, Operator *op, Solver *solver);

#endif

