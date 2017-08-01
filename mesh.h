



#ifndef _MESH_H_
#define _MESH_H_

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <petscksp.h>

#include "array.h"

extern int UniformMesh(double ***pcoord, int *n, double *bounds, double *h, int dimension);

extern int NonUniformMeshY(double ***pcoord, int *n, double *bounds, double *h, int dimension, double (*Transform)(double *bounds, double range, double s) );

int MetricCoefficients2D(Array2d *metrics, double **coord, ArrayInt2d *IsGlobalToGrid, IsRange *range, double *bounds, int dimension, void (*MetricCoefficientsFunc)(double *metricsAtPoint, double *bounds, double *lengths, double x, double y));

#endif
