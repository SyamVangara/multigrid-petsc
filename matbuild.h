




#ifndef _MATBUILD_H_
#define _MATBUILD_H_

#include <stdio.h>
#include <math.h>
#include "petscksp.h"
#include "array.h"

typedef struct {

} Operations;

typedef struct {
	int		grids;   // num of grids per level
	int		*gridId; // Grid Id of each grid in a given level
	ArrayInt2d	global;  // global to grid map
	ArrayInt2d	*grid;   // grid to global map
} Indices;

typedef struct {
	Indices		indices;
	Operations	ops;
} Level;

typedef struct {
	int	levels;
	Level	*level;
} Assembly;

extern Mat matrixA(double *metrics, double **opIH2h, double **opIh2H, int n0, int levels);
Mat levelMatrixA(Array2d metrics, ArrayInt2d IsStencil, int n, int l); 
extern Mat restrictionMatrix(double **Is, int m, int nh, int nH);
extern Mat prolongationMatrix(double **Is, int m, int nh, int nH);
Mat restrictionMatrixMPI(double **Is, int m, IsRange rangeh, IsRange rangeH, ArrayInt2d IsResStencil);
//extern Mat restrictionMatrixMPI(double **Is, int m, int nh, int nH);
Mat prolongationMatrixMPI(double **Is, int m, IsRange rangeh, IsRange rangeH, ArrayInt2d IsProStencil);
//extern Mat prolongationMatrixMPI(double **Is, int m, int nh, int nH);
extern Mat GridTransferMatrix(double **Is, int m, int nh, int nH, char *type);
void levelvecb(Vec *b, double *f);
extern void vecb(Vec *b, Array2d f, double **opIh2H, int n0, int levels);
//extern int ipow(int base, int exp);

#endif


