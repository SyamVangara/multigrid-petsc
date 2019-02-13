
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
	int		rank; // -ve if BC block 
	int		blockID[MAX_DIMENSION]; // -ve in a dir => it is BC block in that dir
	long int	bcStartIndex; // Global start index on BC cells on neighboring block
	long int	bcGStartIndex; // Global grid start index on BC cells on neighboring block
	long int	bcInc[MAX_DIMENSION]; // Increments on BC cells on neigbhoring block
} BCindices;

typedef struct {
	int		ngrids;   // num of grids in this level
	int		*gridId; // Grid Id of each grid in a given level
	double		(*h)[2]; // Delta h in reference domain // ! Remove
	long int	*ranges; // ranges of global level indices for each grid
	long int	(*granges)[2]; // ranges of global grid indices for each grid
	long int	(*inc)[MAX_DIMENSION]; // increments of global index 
						// in each direction for each grid
	BCindices	(*bcindices)[MAX_DIMENSION][2]; // BC indices for all grids in each dir
	BCindices	(*ebcindices)[MAX_DIMENSION][2][2]; // Edge BC indices for all grids in each dir; [dim][2][2] -> Cyclic
	BCindices	(*cbcindices)[2][2][2]; // Corner BC indices for all grids in each dir
	ArrayInt2d	global;  // global to grid map// ! Remove
	ArrayInt2d	*grid;   // grid to global map// ! Remove
} Level;

typedef struct {
	int	dimension;
	int	nlevels;
	int	totalGrids; // ! Remove
	int	coarseningFactor; // ! Remove
	Mat 	*res; // only between successive grids (not levels)
	Mat 	*pro; // only between successive grids (not levels)
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
	double		**rNormGrid; // !Remove
	double		*rNormGlobal; // !Remove
	// more info ends
	double		error[3];
	double		*rnorm;
	Levels		*levels;
//	OutFiles	*outfiles;
//	Assembly	*assem;
} Solver;

typedef struct {
	double	error[3];
	FILE	*solData;
	FILE	*errData;
	FILE	*resData;
	FILE	*XgridData;
	FILE	*YgridData;
} OutFiles;

int CreateLevels(Grids *grids, Levels *levels);
void DestroyLevels(Levels *levels);

int CreateOperator(Grids *grids, Operator *op);
void DestroyOperator(Operator *op);
void GridTransferOperators(Operator op, Levels levels);

int CreateSolver(Grids *grids, Solver *solver); 
void DestroySolver(Solver *solver);
int Solve(Solver *solver);

void PostProcessing(Grids *grids, Solver *solver);

#endif

