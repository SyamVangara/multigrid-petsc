




#ifndef _MATBUILD_H_
#define _MATBUILD_H_

#include <stdio.h>
#include <math.h>
#include "petscksp.h"
#include "array.h"

extern Mat matrixA(double *metrics, double **opIH2h, double **opIh2H, int n0, int levels);
Mat levelMatrixA(Array3d metrics, ArrayInt2d IsGlobalToGrid, ArrayInt2d IsGridToGlobal, int n, int l);
//extern Mat levelMatrixA(Array3d metrics, int n, int level);
extern Mat restrictionMatrix(double **Is, int m, int nh, int nH);
extern Mat prolongationMatrix(double **Is, int m, int nh, int nH);
extern Mat restrictionMatrixMPI(double **Is, int m, int nh, int nH);
extern Mat prolongationMatrixMPI(double **Is, int m, int nh, int nH);
extern Mat GridTransferMatrix(double **Is, int m, int nh, int nH, char *type);
void levelvecb(Vec *b, double *f);
extern void vecb(Vec *b, Array2d f, double **opIh2H, int n0, int levels);
//extern int ipow(int base, int exp);

#endif


