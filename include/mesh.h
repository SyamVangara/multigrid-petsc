



#ifndef _MESH_H_
#define _MESH_H_

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <petscksp.h>

#include "array.h"
#include "problem.h"

#define PI 3.14159265358979323846
#define MAX_DIMENSION 3
#define MIN_POINTS 3

//typedef enum {UNIFORM, NONUNIFORM1, NONUNIFORM2} MeshType;
//
//typedef struct {
//	int	dimension;
//	MeshType type;
//	int	n[MAX_DIMENSION];
//	double	bounds[MAX_DIMENSION][2];
////	double	**coord;
////	double	h;
////	void	(*MetricCoefficients)(void *mesh, double x, double y, double *metrics);
//} Topo;

typedef struct {
	int	dimension;
	int	blockID[MAX_DIMENSION];
	int	type[MAX_DIMENSION];
	int	n[MAX_DIMENSION];
	int	dimProcs[MAX_DIMENSION];
	double	bounds[MAX_DIMENSION][2];
	double	range[MAX_DIMENSION][2];
	double	**coord;
	double	h;
	void	(*MetricCoefficients)(void *mesh, double x, double y, double *metrics);
} Mesh;

int CreateMesh(Mesh *mesh);
void DestroyMesh(Mesh *mesh);

#endif
