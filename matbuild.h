




#ifndef _MATBUILD_H_
#define _MATBUILD_H_

#include <stdio.h>
#include <math.h>
#include "petscksp.h"
extern Mat matrixA(double ***metrics, double **opIH2h, double **opIh2H, int n0, int levels);
extern Mat levelMatrixA(double ***metrics, int n, int level);
extern Mat restrictionMatrix(double **Is, int m, int nh, int nH);
extern Mat prolongationMatrix(double **Is, int m, int nh, int nH);
extern Mat restrictionMatrixMPI(double **Is, int m, int nh, int nH);
extern Mat prolongationMatrixMPI(double **Is, int m, int nh, int nH);
extern Mat GridTransferMatrix(double **Is, int m, int nh, int nH, char *type);
extern void vecb(Vec *b, double **f, double **opIh2H, int n0, int levels);
//extern int ipow(int base, int exp);

#endif


