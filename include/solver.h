
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

typedef struct BCblocks{
	int		rank; // -ve if block is beyond domain
	int		blockID[MAX_DIMENSION]; // -ve in a dir => it is beyond domain in that dir
	long int	bcStartIndex; // Global start index on BC cells on neighboring block
	long int	bcGStartIndex; // Global grid start index on BC cells on neighboring block
	long int	bcInc[MAX_DIMENSION]; // Increments on BC cells on neigbhoring block
	struct BCblocks	*sbcindices; // Second neighoring block
} BCindices;

typedef struct {
	int		prob;
	int		ngrids;   // num of grids in this level
	int		*gridId; // Grid Id of each grid in a given level
	double		eps;
	double		(*h)[2]; // Delta h in reference domain // ! Remove
	long int	*ranges; // ranges of global level indices for each grid
	long int	(*granges)[2]; // ranges of global grid indices for each grid
	long int	(*inc)[MAX_DIMENSION]; // increments of global index 
						// in each direction for each grid
	BCindices	(*bcindices)[MAX_DIMENSION][2]; // BC indices for all grids in each dir
	BCindices	(*ebcindices)[MAX_DIMENSION][2][2]; // Edge BC indices for all grids in each dir; [dim][2][2] -> Cyclic
	BCindices	(*cbcindices)[2][2][2]; // Corner BC indices for all grids in each dir
	IS		*is; // index sets to extract vectors of each grid
	ArrayInt2d	global;  // global to grid map// ! Remove
	ArrayInt2d	*grid;   // grid to global map// ! Remove
} Level;

typedef struct {
	int	dimension;
	int	prob;
	int	nlevels;
	double	eps;
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
	int	cycle;
	int	prob;
	int	numIter;
	int	v[2];
	double	eps;
	double	rtol;
	double	error[3];
	double	*rnorm;
	Levels	*levels;
} Solver;

int CreateLevels(Grids *grids, Levels *levels);
void DestroyLevels(Levels *levels);

void GetSubIS(int lg, Grid *grid, Level *level, IS *indexSet);

int CreateOperator(Grids *grids, Operator *op);
void DestroyOperator(Operator *op);
void GridTransferOperators(Operator op, Levels levels);

int CreateSolver(Grids *grids, Solver *solver); 
void DestroySolver(Solver *solver);
int Solve(Solver *solver);

int PostProcessing(Grids *grids, Solver *solver);

#endif

