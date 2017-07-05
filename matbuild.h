




#ifndef _MATBUILD_H_
#define _MATBUILD_H_

#include <stdio.h>
#include <math.h>
#include "petscksp.h"
Mat matrixA(double ***metrics, double **opIH2h, double **opIh2H, int n0, int levels);
Mat restrictionMatrix(double **Is, int m, int nh, int nH);
Mat prolongationMatrix(double **Is, int m, int nh, int nH);
Mat GridTransferMatrix(double **Is, int m, int nh, int nH, char *type);
void vecb(Vec *b, double **f, double **opIh2H, int n0, int levels);

#endif


