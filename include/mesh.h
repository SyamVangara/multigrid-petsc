



#ifndef _MESH_H_
#define _MESH_H_

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <petscksp.h>

#include "array.h"
#include "problem.h"

#define PI 3.14159265358979323846
#define DIMENSION 2

typedef enum {UNIFORM, NONUNIFORM} MeshType;

typedef struct {
	int	n[DIMENSION];
	double	bounds[DIMENSION*2];
	double	**coord;
	double	h;
	void	(*MetricCoefficients)(void *mesh, double x, double y, double *metrics);
} Mesh;

void SetUpMesh(Mesh *mesh, MeshType type);
void DestroyMesh(Mesh *mesh);

#endif
