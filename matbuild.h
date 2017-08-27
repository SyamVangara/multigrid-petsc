




#ifndef _MATBUILD_H_
#define _MATBUILD_H_

#include <stdio.h>
#include <math.h>
#include "petscksp.h"
#include "array.h"
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
	Vec	*b;
	Vec	*u;
} Assembly;

void SetUpIndices(Mesh *mesh, Indices *indices);
void DestroyIndices(Indices *indices);
void mapping(Indices *indices, int mappingStyleflag);

void SetUpOperator(Indices *indices, Operator *op);
void DestroyOperator(Operator *op);
void GridTransferOperators(Operator op, Indices indices);

void SetUpAssembly(Indices *indices, Assembly *assem);
void DestroyAssembly(Assembly *assem);

void Assemble(Problem *prob, Mesh *mesh, Indices *indices, Operator *op, Assembly *assem);

extern Mat matrixA(double *metrics, double **opIH2h, double **opIh2H, int n0, int levels);
extern Mat restrictionMatrix(double **Is, int m, int nh, int nH);
extern Mat prolongationMatrix(double **Is, int m, int nh, int nH);
Mat restrictionMatrixMPI(double **Is, int m, IsRange rangeh, IsRange rangeH, ArrayInt2d IsResStencil);
Mat prolongationMatrixMPI(double **Is, int m, IsRange rangeh, IsRange rangeH, ArrayInt2d IsProStencil);
extern Mat GridTransferMatrix(double **Is, int m, int nh, int nH, char *type);
extern void vecb(Vec *b, Array2d f, double **opIh2H, int n0, int levels);

#endif


